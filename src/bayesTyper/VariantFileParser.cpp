
/*
VariantFileParser.cpp - This file is part of BayesTyper (v1.1)


The MIT License (MIT)

Copyright (c) 2016 Jonas Andreas Sibbesen and Lasse Maretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <string>
#include <iostream>
#include <assert.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <list>
#include <fstream>
#include <stdlib.h> 
#include <math.h> 
#include <mutex>
#include <thread>
#include <algorithm>

#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "VariantFileParser.hpp"
#include "Utils.hpp"
#include "VariantClusterGroup.hpp"
#include "VariantCluster.hpp"
#include "VariantClusterGraph.hpp"
#include "Sequence.hpp"

using namespace std;

static const uint max_chromosome_substring_size = 10000;
static const uint variant_cluster_group_batch_size = 1000;


VariantFileParser::VariantFileParser(const string genome_filename, const uint max_allele_length_in, const double copy_number_variant_threshold_in, const ushort num_threads_in, const uint prng_seed) : max_allele_length(max_allele_length_in), copy_number_variant_threshold(copy_number_variant_threshold_in), num_threads(num_threads_in) {

    number_of_variants = 0;
    max_alternative_alleles = 0;

    prng = mt19937(prng_seed);

    cout << "[" << Utils::getLocalTime() << "] Parsing reference genome ..." << endl;

    ulong genome_size = parseFasta(genome_filename, &genome_sequences);

    cout << "[" << Utils::getLocalTime() << "] Parsed " << genome_sequences.size() << " sequence(s) from the reference genome (" << genome_size << " nucleotides)\n" << endl;
}


VariantFileParser::~VariantFileParser() {

    for (auto &cit: genome_sequences) {

        delete cit.second;
    }

    for (auto &cit: decoy_sequences) {

        delete cit.second;
    }
}


ulong VariantFileParser::parseFasta(const string filename, unordered_map<string, string *> * fasta_sequences) {

    assert(fasta_sequences->empty());

    ifstream fasta_file(filename.c_str());
    assert(fasta_file.is_open());

    string fasta_line;
    string sequence_name = "";
    string * sequence = new string();

    ulong fasta_size = 0;

    while (fasta_file.good()) {

        getline(fasta_file, fasta_line);

        if (fasta_line.substr(0,1) == ">") {

            if (sequence_name != "") {
            
                fasta_size += sequence->size();    

                transform(sequence->begin(), sequence->end(), sequence->begin(), ::toupper);

                assert(!(sequence->empty()));
                assert(sequence->size() <= Utils::uint_overflow);   
                assert(fasta_sequences->insert({sequence_name, sequence}).second);
                
                sequence = new string();
            }

            vector<string> sequence_name_split;

            boost::split(sequence_name_split, fasta_line, boost::is_any_of("\t "));
            assert(sequence_name_split.front().at(0) == '>');

            sequence_name = simplifyChromosomeId(sequence_name_split.front().substr(1));
            assert(!(sequence_name.empty()));

        } else {

            assert(!(sequence_name.empty()));
            sequence->append(fasta_line);
        }
    }

    fasta_size += sequence->size();    

    transform(sequence->begin(), sequence->end(), sequence->begin(), ::toupper);

    assert(!(sequence_name.empty()));
    assert(!(sequence->empty()));
    assert(sequence->size() <= Utils::uint_overflow);   
    assert(fasta_sequences->insert({sequence_name, sequence}).second);
    
    fasta_file.close();

    return fasta_size;
}


template <int kmer_size>
void VariantFileParser::addDecoys(const string decoy_filename) {

    cout << "[" << Utils::getLocalTime() << "] Parsing decoy sequence(s) ..." << endl;

    ulong decoy_size = 0;

    if (!(decoy_filename.empty())) {

        decoy_size = parseFasta(decoy_filename, &decoy_sequences);

        for (auto &cit: decoy_sequences) {

            assert(genome_sequences.find(cit.first) == genome_sequences.end());
        }

        unordered_set<string> empty_set; 
        assert(empty_set.empty());

        addSequencesToInterclusterRegions<kmer_size>(decoy_sequences, &empty_set);
        assert(empty_set.size() == decoy_sequences.size());
    }

    cout << "[" << Utils::getLocalTime() << "] Parsed " << decoy_sequences.size() << " decoy sequence(s) (" << decoy_size << " nucleotides)\n" << endl;
}


template <int kmer_size>
void VariantFileParser::readVariantFile(const string vcf_filename, vector<VariantClusterGraph*> * variant_cluster_graphs, vector<VariantClusterGroup*> * variant_cluster_groups) { 

    ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > variant_cluster_group_queue(Utils::queue_size_thread_scaling * num_threads);
    mutex counting_mutex;

    pair<vector<uint>, vector<uint> > variant_stats;

    vector<thread> processing_threads(num_threads);

    for (int i=0; i < num_threads; i++) {

        processing_threads.at(i) = thread(&VariantFileParser::processVariantClusterGroupsCallback<kmer_size>, this, &variant_cluster_group_queue, &counting_mutex, variant_cluster_graphs, variant_cluster_groups);
    }  

    if ((vcf_filename.substr(vcf_filename.size()-3,3) == ".gz") or (vcf_filename.substr(vcf_filename.size()-5,5) == ".gzip")) {

        ifstream gz_vcf_file(vcf_filename.c_str(), ifstream::binary);
        assert(gz_vcf_file.is_open());            
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(gz_vcf_file);

        variant_stats = parseVariants<boost::iostreams::filtering_istream, kmer_size>(&in, &variant_cluster_group_queue);

        gz_vcf_file.close();

    } else {

        assert((vcf_filename.substr(vcf_filename.size()-4,4) == ".vcf"));
        ifstream vcf_file(vcf_filename.c_str());
        assert(vcf_file.is_open());

        variant_stats = parseVariants<ifstream, kmer_size>(&vcf_file, &variant_cluster_group_queue);
    
        vcf_file.close();
    }

    variant_cluster_group_queue.pushedLast();

    for (auto& th: processing_threads) {

        th.join();
    } 

    cout << "\n[" << Utils::getLocalTime() << "] Parsed " << variant_stats.first.at(uchar(AlleleCount::Total)) << " alternative alleles and excluded:\n" << endl; 
    cout << "\t- Alleles on chromosome(s) not in genome: " << variant_stats.first.at(uchar(AlleleCount::Excluded_genome)) << endl;  
    cout << "\t- Alleles with reference not equal to genome sequence: " << variant_stats.first.at(uchar(AlleleCount::Excluded_match)) << endl;
    cout << "\t- Alleles with non-canonical bases (not ACGTN): " << variant_stats.first.at(uchar(AlleleCount::Excluded_canon)) << endl;
    cout << "\t- Alleles within " << kmer_size << " bases of chromosome end: " << variant_stats.first.at(uchar(AlleleCount::Excluded_end)) << endl;
    cout << "\t- Alleles longer than " << max_allele_length << " bases: " << variant_stats.first.at(uchar(AlleleCount::Excluded_length)) << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Out of " << number_of_variants << " variants:\n" << endl; 
    cout << "\t- Single nucleotides polymorphism: " << variant_stats.second.at(uchar(Utils::VariantType::SNP)) << endl;
    cout << "\t- Insertion: " << variant_stats.second.at(uchar(Utils::VariantType::Insertion)) << endl;
    cout << "\t- Deletion: " << variant_stats.second.at(uchar(Utils::VariantType::Deletion)) << endl;
    cout << "\t- Complex: " << variant_stats.second.at(uchar(Utils::VariantType::Complex)) << endl;
    cout << "\t- Mixture: " << variant_stats.second.at(uchar(Utils::VariantType::Mixture)) << endl;
    cout << "\t- Unsupported (excluded): " << variant_stats.second.at(uchar(Utils::VariantType::Unsupported)) << endl;

    uint sum_variant_types = 0;

    for (auto i = 0; i < static_cast<uchar>(Utils::VariantType::VARIANT_TYPE_SIZE); i++) {

        sum_variant_types += variant_stats.second.at(i);
    }

    assert(number_of_variants == sum_variant_types);

    cout << "\n[" << Utils::getLocalTime() << "] Merged variants into " << variant_cluster_graphs->size() << " clusters and further into " << variant_cluster_groups->size() << " groups\n" << endl;


    cout << "[" << Utils::getLocalTime() << "] Shuffling intercluster regions ..." << endl;

    shuffle(intercluster_regions.begin(), intercluster_regions.end(), prng);

    cout << "[" << Utils::getLocalTime() << "] Finished shuffling\n" << endl;
}


template <typename FileType, int kmer_size>
pair<vector<uint>, vector<uint> > VariantFileParser::parseVariants(FileType * vcf_file, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue) {

    cout << "[" << Utils::getLocalTime() << "] Parsing variant file ..." << endl; 

    vector<uint> allele_counter = vector<uint>(static_cast<uchar>(AlleleCount::ALLELE_COUNT_SIZE), 0);
    vector<uint> variant_type_counter = vector<uint>(static_cast<uchar>(Utils::VariantType::VARIANT_TYPE_SIZE), 0);

    string cur_chromosome_id = "";
    string prev_chromosome_id = "";
    string cur_chromosome_name = "";
    string prev_chromosome_name = "";
    
    unordered_map<string, string *>::iterator chromosome_seq_it = genome_sequences.end();

    long cur_position = 0;
    long prev_position = -1;

    long prev_var_end_position = -1;
    long cur_group_end_position = -1;     

    vector<unordered_map<uint, VariantCluster*> * > * variant_cluster_group_batch = new vector<unordered_map<uint, VariantCluster*> * >();  
    variant_cluster_group_batch->reserve(variant_cluster_group_batch_size);
    
    uint variant_cluster_group_batch_complexity = 0;

    unordered_map<uint, VariantCluster*> * variant_cluster_group = new unordered_map<uint, VariantCluster*>();
    
    uint variant_cluster_group_counter = 0;

    map<uint, VariantCluster*> variant_cluster_group_flanks;
    list<unordered_set<uint> > variant_cluster_group_merge_sets;

    vector<string> variant_line_split(5);

    ulong chromosome_variant_counter = 0;
    unordered_set<string> counted_chromosomes;

    set<uint> variant_depedencies;

    while (vcf_file->good()) {

        getline(*vcf_file, variant_line_split.at(0), '\t');

        if (variant_line_split.at(0).size() == 0) {

            vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');

            continue;              
        }

        if (variant_line_split.at(0).substr(0,1) == "#") {

            assert(variant_line_split.at(0).substr(variant_line_split.at(0).size() - 6, 6) == "#CHROM");

            vector<string> header_string_split;
            boost::split(header_string_split, variant_line_split.at(0), boost::is_any_of("\n"));

            for (auto &header_line: header_string_split) {

                if (header_line.substr(0,8) == "##contig") {

                    header_line.pop_back();
                    vector<string> contig_line_split;
                    boost::split(contig_line_split, header_line, boost::is_any_of(","));

                    assert(contig_line_split.size() >= 2);

                    string contig_id = "";
                    uint contig_length = 0;  

                    for (auto &contig_info: contig_line_split) {

                        if (contig_info.substr(0,13) == "##contig=<ID=") {
        
                            assert(contig_info.size() > 13);
                            contig_id = contig_info.substr(13);

                        } else if (contig_info.substr(0,7) == "length=") {

                            assert(contig_info.size() > 7);
                            contig_length = stoi(contig_info.substr(7));
                        }
                    }                  

                    assert(!(contig_id.empty()));
                    assert(contig_length > 0);

                    auto genome_sequences_it = genome_sequences.find(simplifyChromosomeId(contig_id)); 

                    if (genome_sequences_it != genome_sequences.end()) {

                        assert(contig_length == genome_sequences_it->second->size());
                    }
                }
            }

            string last_header_string;
            getline(*vcf_file, last_header_string);
            assert(count(last_header_string.begin(), last_header_string.end(), '\t') >= 6);

            continue;   
        } 

        for (uchar i = 1; i < 5; i++) {

            getline(*vcf_file, variant_line_split.at(i), '\t');
        }

        vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');

        cur_chromosome_name = variant_line_split.at(0);
        cur_chromosome_id = simplifyChromosomeId(variant_line_split.at(0));

        cur_position = stoi(variant_line_split.at(1)) - 1;

        if ((prev_chromosome_id != "") and (cur_chromosome_id != prev_chromosome_id)) {

            auto pgit = genome_sequences.find(prev_chromosome_id);

            if (pgit != genome_sequences.end()) {

                processVariantClusters(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);                   
                intercluster_regions.emplace_back(pgit->second->cbegin() + prev_var_end_position + 1, pgit->second->cend(), classifyChromosome(prev_chromosome_id));

                prev_var_end_position = -1;
                cur_group_end_position = -1;
                variant_cluster_group_counter = 0;

                assert(counted_chromosomes.insert(prev_chromosome_id).second);

            } else {

                assert(prev_var_end_position == -1);
                assert(cur_group_end_position == -1);
                assert(variant_cluster_group_counter == 0);
            }

            variant_depedencies.clear();
   
            cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosome_variant_counter << " variants on chromosome " << prev_chromosome_name << endl; 
            chromosome_variant_counter = 0;

        } else {

            assert(prev_position < cur_position);
        }

        VariantCluster::Variant cur_variant;
        cur_variant.id = number_of_variants;
        cur_variant.type = Utils::VariantType::Unsupported;
        cur_variant.max_reference_length = 0;

        number_of_variants += 1;
        chromosome_variant_counter += 1;

        prev_chromosome_id = cur_chromosome_id;
        prev_chromosome_name = cur_chromosome_name;
        prev_position = cur_position;

        assert(variant_line_split.at(3).find(",") == string::npos);        

        string var_ref_seq = variant_line_split.at(3);

        transform(var_ref_seq.begin(), var_ref_seq.end(), var_ref_seq.begin(), ::toupper);

        vector<string> alt_alleles;
        boost::split(alt_alleles, variant_line_split.at(4), boost::is_any_of(","));
        assert(alt_alleles.size() < (Utils::ushort_overflow - 2));

        auto vit = variant_depedencies.begin();
        
        while (vit != variant_depedencies.end()) {
            
            if (*vit >= cur_position) {

                break;
            }
                
            variant_depedencies.erase(vit);
            vit = variant_depedencies.begin();
        }

        cur_variant.has_dependency = !(variant_depedencies.empty());

        if (alt_alleles.back() == "*") {

            assert(cur_variant.has_dependency);
            alt_alleles.pop_back();
        }

        assert(!(alt_alleles.empty()));
        allele_counter.at(uchar(AlleleCount::Total)) += alt_alleles.size();

        chromosome_seq_it = genome_sequences.find(cur_chromosome_id);

        if (chromosome_seq_it == genome_sequences.end()) {

            assert(variant_cluster_group->empty());

            allele_counter.at(uchar(AlleleCount::Excluded_genome)) += alt_alleles.size();
            variant_type_counter.at(uchar(Utils::VariantType::Unsupported))++;

            continue;
        }

        vector<string> ref_alleles(alt_alleles.size(), var_ref_seq);

        for (ushort i = 0; i < alt_alleles.size(); i++) {

            transform(alt_alleles.at(i).begin(), alt_alleles.at(i).end(), alt_alleles.at(i).begin(), ::toupper);
            assert(count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0);

            rightTrimAllele(&ref_alleles.at(i), &alt_alleles.at(i));
        }

        assert(ref_alleles.size() == alt_alleles.size());

        bool is_excluded = false;

        string gen_ref_seq = chromosome_seq_it->second->substr(cur_position, variant_line_split.at(3).size()); 
        assert(gen_ref_seq.size() == var_ref_seq.size());

        if (var_ref_seq.compare(gen_ref_seq) != 0) {

            allele_counter.at(uchar(AlleleCount::Excluded_match)) += alt_alleles.size();
            is_excluded = true;
        }

        if (var_ref_seq.find_first_not_of("ACGTN") != string::npos) {

            allele_counter.at(uchar(AlleleCount::Excluded_canon)) += alt_alleles.size();
            is_excluded = true;
        }

        if (cur_position < (kmer_size - 1)) {

            allele_counter.at(uchar(AlleleCount::Excluded_end)) += alt_alleles.size();
            is_excluded = true;
        }

        unordered_set<ushort> excluded_alleles;

        if (!is_excluded) {

            for (ushort i = 0; i < alt_alleles.size(); i++) {

                if (alt_alleles.at(i).find_first_not_of("ACGTN") != string::npos) {

                    allele_counter.at(uchar(AlleleCount::Excluded_canon))++;
                    assert(excluded_alleles.insert(i).second);

                } else if ((cur_position + static_cast<long>(ref_alleles.at(i).size()) - 1 + kmer_size) > static_cast<long>(chromosome_seq_it->second->size())) {

                    allele_counter.at(uchar(AlleleCount::Excluded_end))++;
                    assert(excluded_alleles.insert(i).second);

                } else if ((ref_alleles.at(i).size() > max_allele_length) or (alt_alleles.at(i).size() > max_allele_length)) {

                    allele_counter.at(uchar(AlleleCount::Excluded_length))++;
                    assert(excluded_alleles.insert(i).second);

                } else {

                    variant_depedencies.insert(cur_position + ref_alleles.at(i).size() - 1);
                }
            }
        }

        assert(excluded_alleles.size() <= alt_alleles.size());

        if (is_excluded or excluded_alleles.size() == alt_alleles.size()) {

            variant_type_counter.at(uchar(Utils::VariantType::Unsupported))++;            
            continue;     
        }

        assert(!(variant_depedencies.empty()));
        assert(prev_var_end_position <= cur_group_end_position);

        if ((cur_position - cur_group_end_position) >= kmer_size) {

            processVariantClusters(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);
            variant_cluster_group_counter = 0;
        } 

        if ((cur_position - (prev_var_end_position + 1)) >= kmer_size) {

            intercluster_regions.emplace_back(chromosome_seq_it->second->cbegin() + prev_var_end_position + 1, chromosome_seq_it->second->cbegin() + cur_position, classifyChromosome(cur_chromosome_id));
        }

        set<uint> cur_end_positions;

        for (ushort i = 0; i < alt_alleles.size(); i++) {

            assert(!(ref_alleles.at(i).empty()));            
            assert(!(alt_alleles.at(i).empty()));

            if (excluded_alleles.count(i) > 0) {

                cur_variant.excluded_alt_alleles.push_back(i + 1);
            
            } else {

                parseAllele(ref_alleles.at(i), alt_alleles.at(i), cur_position, i + 1, &cur_variant, *(chromosome_seq_it->second));
                
                uint allele_start_position = 0;
                uint chromosome_start_position = cur_position + ref_alleles.at(i).size();

                if (ref_alleles.at(i).front() == alt_alleles.at(i).front()) {

                    allele_start_position++;
                }

                uint ref_copy_number_variant_length = copyNumberVariantLength<kmer_size>(ref_alleles.at(i), allele_start_position, *chromosome_seq_it->second, chromosome_start_position);
                uint alt_copy_number_variant_length = copyNumberVariantLength<kmer_size>(alt_alleles.at(i), allele_start_position, *chromosome_seq_it->second, chromosome_start_position);

                cur_end_positions.insert(cur_position + ref_alleles.at(i).size() - 1);
                cur_group_end_position = max(cur_group_end_position, static_cast<long>(cur_position + ref_alleles.at(i).size() - 1 + max(ref_copy_number_variant_length, alt_copy_number_variant_length)));
            }
        }

        assert(!(cur_end_positions.empty()));
        assert(!cur_variant.alternative_alleles.empty());
        assert((cur_variant.alternative_alleles.size() + cur_variant.excluded_alt_alleles.size()) == alt_alleles.size());
        assert(cur_variant.excluded_alt_alleles.size() == excluded_alleles.size());

        prev_var_end_position = max(prev_var_end_position, static_cast<long>(*cur_end_positions.rbegin()));
        assert(cur_group_end_position >= prev_var_end_position);

        max_alternative_alleles = max(max_alternative_alleles, ushort(alt_alleles.size()));

        clusterVariants<kmer_size>(cur_variant, cur_position, cur_end_positions, cur_chromosome_id, cur_chromosome_name, &variant_cluster_group_counter, &variant_cluster_group_flanks, variant_cluster_group, &variant_cluster_group_merge_sets);

        variant_cluster_group_batch_complexity++;   
        variant_type_counter.at(uchar(cur_variant.type))++;
    }

    variant_cluster_group_batch_complexity = variant_cluster_group_batch_size;

    processVariantClusters(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);

    assert(variant_cluster_group_batch->empty());
    delete variant_cluster_group_batch;

    assert(variant_cluster_group->empty());
    delete variant_cluster_group;
    
    if (chromosome_seq_it != genome_sequences.end()) {
        
        intercluster_regions.emplace_back(chromosome_seq_it->second->cbegin() + prev_var_end_position + 1, chromosome_seq_it->second->cend(), classifyChromosome(cur_chromosome_id));                                 
        assert(counted_chromosomes.insert(cur_chromosome_id).second);
        assert(cur_chromosome_id == prev_chromosome_id);
    }

    cout << "[" << Utils::getLocalTime() << "] Parsed " << chromosome_variant_counter << " variants on chromosome " << cur_chromosome_name << endl; 

    addSequencesToInterclusterRegions<kmer_size>(genome_sequences, &counted_chromosomes);
    assert(counted_chromosomes.size() == genome_sequences.size());

    uint variant_counter_sanity_check = 0;

    for (auto &type_count: variant_type_counter) {

        variant_counter_sanity_check += type_count;
    }

    assert(number_of_variants == variant_counter_sanity_check);

    return make_pair(allele_counter, variant_type_counter);
}


void VariantFileParser::rightTrimAllele(string * ref_allele, string * alt_allele) {

    while ((ref_allele->size() > 1) and (alt_allele->size() > 1)) {

        if (ref_allele->back() == alt_allele->back()) {

            ref_allele->pop_back();
            alt_allele->pop_back();

        } else {

            break;
        }
    }

    assert(!(ref_allele->empty()));
    assert(!(alt_allele->empty()));
}


void VariantFileParser::parseAllele(const string & ref_allele, const string & alt_allele, const uint cur_position, const ushort alternative_allele_counter, VariantCluster::Variant * cur_variant, const string & chromosome_sequence) {

    auto ref_allele_it = ref_allele.begin();
    auto alt_allele_it = alt_allele.begin();

    uint identical_left_nucleotides = 0;

    while ((ref_allele_it != ref_allele.end()) and (alt_allele_it != alt_allele.end())) {

        if (*alt_allele_it == *ref_allele_it) {

            identical_left_nucleotides++;

        } else {

            break;
        }

        ref_allele_it++;
        alt_allele_it++;
    }

    assert((ref_allele.size() > identical_left_nucleotides) or (alt_allele.size() > identical_left_nucleotides));

    cur_variant->alternative_alleles.emplace_back(alternative_allele_counter, pair<uint, const string>(ref_allele.size(), alt_allele));
    cur_variant->max_reference_length = max(cur_variant->max_reference_length, static_cast<uint>(ref_allele.size()));

    Utils::VariantType variant_type = classifyAllele(ref_allele.size() - identical_left_nucleotides, alt_allele.size() - identical_left_nucleotides);

    if (cur_variant->type == Utils::VariantType::Unsupported) {

        cur_variant->type = variant_type;

    } else {

        if (cur_variant->type != variant_type) {

            cur_variant->type = Utils::VariantType::Mixture;
        }
    }
}


Utils::VariantType VariantFileParser::classifyAllele(const int reference_size, const int allele_size) {

    assert((reference_size > 0) or (allele_size > 0));

    if ((reference_size == 1) and (allele_size == 1)) {

        return Utils::VariantType::SNP;
    
    } else if ((reference_size == 0) or (allele_size == 0)) {

        if ((allele_size - reference_size) > 0) {

            return Utils::VariantType::Insertion;

        } else {

            assert((allele_size - reference_size) < 0);
            
            return Utils::VariantType::Deletion;
        }
    }

    return Utils::VariantType::Complex;
}


template <int kmer_size>
uint VariantFileParser::copyNumberVariantLength(const string & allele_sequence, const uint allele_start_position, const string & chromosome_sequence, const uint chromosome_start_position) {

    assert(allele_start_position <= allele_sequence.size());
    assert(chromosome_start_position <= chromosome_sequence.size());

    uint copy_number_variant_length = 0;

    if (allele_sequence.size() == allele_start_position) {

        return copy_number_variant_length;
    }

    const uint allele_length = allele_sequence.size() - allele_start_position;

    if (allele_length < static_cast<uint>(kmer_size)) {

        while (chromosome_sequence.compare(chromosome_start_position + copy_number_variant_length, allele_length, allele_sequence, allele_start_position, allele_length) == 0) {

            copy_number_variant_length += allele_length;
        }

    } else {

        unordered_set<bitset<kmer_size * 2> > allele_kmers;
        KmerPair<kmer_size> kmer_pair = KmerPair<kmer_size>();

        auto allele_sequence_it = allele_sequence.begin();
        advance(allele_sequence_it, allele_start_position);

        assert(allele_sequence_it <= allele_sequence.end());

        while (allele_sequence_it != allele_sequence.end()) {

            auto nt_bits = Sequence::ntToBit<1>(*allele_sequence_it);

            if (nt_bits.second) {

                if (kmer_pair.move(nt_bits.first)) {

                    allele_kmers.emplace(kmer_pair.getLexicographicalLowestKmer());
                }

            } else {

                kmer_pair.reset();
            }

            allele_sequence_it++;
        } 

        uint chromosome_window_end_position = min(chromosome_start_position + copy_number_variant_length + allele_length, static_cast<uint>(chromosome_sequence.size()));

        while (true) {

            kmer_pair.reset();

            uint num_bases = 0;
            uint num_identical_kmers = 0;

            pair<double, uint> highest_scoring_window(0, 0);

            for (uint chromosome_position = chromosome_start_position + copy_number_variant_length; chromosome_position < chromosome_window_end_position; chromosome_position++) {

                auto nt_bits = Sequence::ntToBit<1>(chromosome_sequence.at(chromosome_position));

                if (nt_bits.second) {

                    if (kmer_pair.move(nt_bits.first)) {

                        if (allele_kmers.count(kmer_pair.getLexicographicalLowestKmer()) > 0) {

                            num_identical_kmers++;
                        }
                    }

                } else {

                    kmer_pair.reset();
                }

                num_bases++;

                if (num_identical_kmers > 0) {

                    assert(num_bases >= kmer_size);
                    double identical_kmer_fraction = num_identical_kmers/static_cast<double>(num_bases - kmer_size + 1);

                    if (Utils::doubleCompare(identical_kmer_fraction, highest_scoring_window.first) or (identical_kmer_fraction > highest_scoring_window.first)) {

                        highest_scoring_window.first = identical_kmer_fraction;
                        highest_scoring_window.second = num_bases;
                    }
                }
            }

            if (highest_scoring_window.first < copy_number_variant_threshold) {

                break;
            }

            copy_number_variant_length += highest_scoring_window.second;

            if (chromosome_window_end_position == chromosome_sequence.size()) {

                break;
            }

            chromosome_window_end_position = min(chromosome_start_position + copy_number_variant_length + allele_length, static_cast<uint>(chromosome_sequence.size()));
        }
    }

    return copy_number_variant_length;
}


template <int kmer_size>
void VariantFileParser::clusterVariants(VariantCluster::Variant & cur_variant, const uint cur_position, const set<uint> cur_end_positions, const string & cur_chromosome_id, const string & cur_chromosome_name, uint * variant_cluster_group_counter, map<uint, VariantCluster*> * variant_cluster_group_flanks, unordered_map<uint, VariantCluster*> * variant_cluster_group, list<unordered_set<uint> > * variant_cluster_group_merge_sets) {

    if (!variant_cluster_group_flanks->empty()) {

        auto lit = variant_cluster_group_flanks->begin();            
        assert(!variant_cluster_group->empty());

        while (static_cast<int>(cur_position - lit->first) >= kmer_size) {

            lit = variant_cluster_group_flanks->erase(lit);

            if (lit == variant_cluster_group_flanks->end()) {

                assert(variant_cluster_group_flanks->empty());
                break;
            }
        }

    } else {

        assert(variant_cluster_group->empty());
    }

    set<VariantCluster*> second_overlaps;
    VariantCluster * variant_cluster_group_first = nullptr;

    auto vit = variant_cluster_group_flanks->begin();

    while (vit != variant_cluster_group_flanks->end()) {

        if ((abs(static_cast<int>(cur_position - vit->first)) + 1) <= kmer_size) {

            if (variant_cluster_group_first == nullptr) {

                variant_cluster_group_first = vit->second;

                if (cur_position >= vit->first) {

                    vit = variant_cluster_group_flanks->erase(vit);
                    continue;                    
                }
            
            } else if (variant_cluster_group_first != vit->second) {   

                second_overlaps.insert(vit->second);
            }  
        } 

        for (auto &cit: cur_end_positions) {

            if ((abs(static_cast<int>(cit - vit->first)) + 1) <= kmer_size) {

                if (variant_cluster_group_first == nullptr) {

                    variant_cluster_group_first = vit->second;

                    if (cur_position >= vit->first) {

                        vit = variant_cluster_group_flanks->erase(vit);
                        continue;                    
                    }
                
                } else if (variant_cluster_group_first != vit->second) {   

                    second_overlaps.insert(vit->second);
                }   
        
            } else if ((cur_position < vit->first) and (cit > vit->first)) {

                if (variant_cluster_group_first == nullptr) {

                    variant_cluster_group_first = vit->second;
                
                } else if (variant_cluster_group_first != vit->second) {   

                    second_overlaps.insert(vit->second);
                }            
            }
        }

        vit++;
    }

    if (variant_cluster_group_first == nullptr) {

        assert(second_overlaps.empty());

        VariantCluster * variant_cluster = new VariantCluster();
        variant_cluster->cluster_idx = *variant_cluster_group_counter;
        variant_cluster->variants = map<uint, VariantCluster::Variant>();
        variant_cluster->left_flank = cur_position;
        variant_cluster->right_flank = *cur_end_positions.rbegin();

        variant_cluster->chromosome_id = cur_chromosome_id;
        variant_cluster->chromosome_name = cur_chromosome_name;
        variant_cluster->chromosome_class = classifyChromosome(cur_chromosome_id);

        assert(variant_cluster->variants.insert(pair<uint, VariantCluster::Variant>(cur_position, cur_variant)).second);

        for (auto &cit: cur_end_positions) {

            assert(variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cit, variant_cluster)).second);
        }

        if (*cur_end_positions.rbegin() - cur_position >= kmer_size) {

            if (!(variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cur_position, variant_cluster)).second)) {

                assert(cur_end_positions.count(cur_position) > 0);   
            }          
        }

        assert(variant_cluster_group->insert({*variant_cluster_group_counter, variant_cluster}).second);

        (*variant_cluster_group_counter)++;
        assert(*variant_cluster_group_counter < Utils::uint_overflow);

    } else {

        assert(variant_cluster_group_first->chromosome_id == cur_chromosome_id);
        assert(variant_cluster_group_first->chromosome_name == cur_chromosome_name);
        assert(variant_cluster_group_first->chromosome_class == classifyChromosome(cur_chromosome_id));

        assert(variant_cluster_group_first->variants.insert(pair<uint, VariantCluster::Variant>(cur_position, cur_variant)).second);
        variant_cluster_group_first->right_flank = max(*cur_end_positions.rbegin(), variant_cluster_group_first->right_flank);

        for (auto &cit: cur_end_positions) {

            if (!(variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cit, variant_cluster_group_first)).second)) {

                if (variant_cluster_group_flanks->at(cit) != variant_cluster_group_first) {

                    bool overlap_check = false;
                    assert(!second_overlaps.empty());

                    for (auto &lit: second_overlaps) {
                        
                        if (lit == variant_cluster_group_flanks->at(cit)) {

                           overlap_check = true; 
                        }
                    }

                    assert(overlap_check);
                }
            }
        }
      
        if (*cur_end_positions.rbegin() - cur_position >= kmer_size) {

            if (!(variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cur_position, variant_cluster_group_first)).second)) {

                if (variant_cluster_group_flanks->at(cur_position) != variant_cluster_group_first) {

                    bool overlap_check = false;
                    assert(!second_overlaps.empty());

                    for (auto &lit: second_overlaps) {
                        
                        if (lit == variant_cluster_group_flanks->at(cur_position)) {

                           overlap_check = true; 
                        }
                    }

                    assert(overlap_check);
                }
            }        
        }
    }


    if (!second_overlaps.empty()) {

        assert(variant_cluster_group_first != nullptr);

        auto found_set = variant_cluster_group_merge_sets->end();
        auto sit = variant_cluster_group_merge_sets->begin();

        while (sit != variant_cluster_group_merge_sets->end()) {

            assert(!(sit->empty()));

            if (sit->count(variant_cluster_group_first->cluster_idx) > 0) {

                if (found_set == variant_cluster_group_merge_sets->end())

                    found_set = sit;

                else {

                    if (sit != found_set) {

                        found_set->insert(sit->begin(), sit->end());
                        sit = variant_cluster_group_merge_sets->erase(sit);  
                        continue;
                    }
                }
            }

            bool merged_cluster_merge_sets = false;

            for (auto &lit: second_overlaps) {

                if (sit->count(lit->cluster_idx) > 0) {

                    if (found_set == variant_cluster_group_merge_sets->end())

                        found_set = sit;

                    else {

                        if (sit != found_set) {

                            found_set->insert(sit->begin(), sit->end());
                            sit = variant_cluster_group_merge_sets->erase(sit);
                            merged_cluster_merge_sets = true; 
                            break;
                        }
                    }
                } 
            } 

            if (!merged_cluster_merge_sets) {

                sit++;
            }
        }

        if (found_set != variant_cluster_group_merge_sets->end()) {

            found_set->insert(variant_cluster_group_first->cluster_idx);

            for (auto &lit: second_overlaps) {
                
                found_set->insert(lit->cluster_idx);
            }
        
        } else {

            variant_cluster_group_merge_sets->push_back(unordered_set<uint>());
            assert(variant_cluster_group_merge_sets->back().insert(variant_cluster_group_first->cluster_idx).second);

            for (auto &lit: second_overlaps) {
                
                variant_cluster_group_merge_sets->back().insert(lit->cluster_idx);
            }
        }
    }
}


void VariantFileParser::processVariantClusters(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue, vector<unordered_map<uint, VariantCluster*> * > ** variant_cluster_group_batch, uint * variant_cluster_group_batch_complexity, unordered_map<uint, VariantCluster*> ** variant_cluster_group, list<unordered_set<uint> > * variant_cluster_group_merge_sets, map<uint, VariantCluster*> * variant_cluster_group_flanks) {

    if (!(*variant_cluster_group)->empty()) {

        mergeVariantClusters(*variant_cluster_group, *variant_cluster_group_merge_sets);
        
        (*variant_cluster_group_batch)->push_back(*variant_cluster_group);            
        *variant_cluster_group = new unordered_map<uint, VariantCluster*>();
    }

    variant_cluster_group_merge_sets->clear();
    variant_cluster_group_flanks->clear();

    if (*variant_cluster_group_batch_complexity >= variant_cluster_group_batch_size) {

        *variant_cluster_group_batch_complexity = 0;

        variant_cluster_group_queue->push(*variant_cluster_group_batch);
        *variant_cluster_group_batch = new vector<unordered_map<uint, VariantCluster*> * >();
        (*variant_cluster_group_batch)->reserve(variant_cluster_group_batch_size);
    }
}


void VariantFileParser::mergeVariantClusters(unordered_map<uint, VariantCluster*> * variant_cluster_group, list<unordered_set<uint> > & variant_cluster_group_merge_sets) {

    unordered_set<uint> already_clustered;

    for (auto &cit: variant_cluster_group_merge_sets) {

        assert(cit.size() > 1);

        auto cur_cluster = cit.begin();
        uint first_cluster = *cit.begin();

        assert(already_clustered.emplace(*cur_cluster).second);
        cur_cluster++;

        while (cur_cluster != cit.end()) {

            assert(already_clustered.emplace(*cur_cluster).second);
            variant_cluster_group->at(first_cluster)->mergeVariantClusters(*variant_cluster_group->at(*cur_cluster));

            delete variant_cluster_group->at(*cur_cluster);
            assert(variant_cluster_group->erase(*cur_cluster));

            cur_cluster++;
        }
    }
}


template <int kmer_size>
void VariantFileParser::processVariantClusterGroupsCallback(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue, mutex * process_mutex, vector<VariantClusterGraph*> * variant_cluster_graphs, vector<VariantClusterGroup*> * variant_cluster_groups) {

    vector<VariantClusterGraph*> variant_cluster_graphs_local;
    vector<VariantClusterGroup*> variant_cluster_groups_local;

    vector<unordered_map<uint, VariantCluster *> * > * variant_cluster_group_batch = nullptr;

    while (variant_cluster_group_queue->pop(&variant_cluster_group_batch)) {

        for (auto & cur_variant_cluster_group: *variant_cluster_group_batch) {

            assert(!(cur_variant_cluster_group->empty()));
            assert(cur_variant_cluster_group->size() < Utils::uint_overflow);

            auto cur_variant_cluster_group_dependencies = getVariantClusterGroupDependencies(cur_variant_cluster_group);

            vector<VariantCluster *> cur_variant_clusters;
            vector<VariantClusterGraph *> cur_variant_cluster_graphs;

            cur_variant_clusters.reserve(cur_variant_cluster_group->size());
            cur_variant_cluster_graphs.reserve(cur_variant_cluster_group->size());

            for (auto & variant_cluster: *cur_variant_cluster_group) {

                assert(variant_cluster.first == variant_cluster.second->cluster_idx);

                assert(cur_variant_cluster_group->begin()->second->chromosome_id == variant_cluster.second->chromosome_id);
                assert(cur_variant_cluster_group->begin()->second->chromosome_name == variant_cluster.second->chromosome_name);
                assert(cur_variant_cluster_group->begin()->second->chromosome_class == variant_cluster.second->chromosome_class);

                assert(!(variant_cluster.second->variants.empty()));
                assert(variant_cluster.second->variants.size() < Utils::ushort_overflow);
                assert(variant_cluster.second->left_flank == variant_cluster.second->variants.begin()->first);

                cur_variant_clusters.push_back(variant_cluster.second);
                cur_variant_cluster_graphs.emplace_back(new VariantClusterKmerGraph<kmer_size>(variant_cluster.second, *genome_sequences.at(variant_cluster.second->chromosome_id)));       

                variant_cluster_graphs_local.push_back(cur_variant_cluster_graphs.back());
            }

            variant_cluster_groups_local.emplace_back(new VariantClusterGroup(cur_variant_clusters, cur_variant_cluster_graphs, cur_variant_cluster_group_dependencies, cur_variant_cluster_group->begin()->second->chromosome_class, cur_variant_cluster_group->begin()->second->chromosome_name));

            for (auto &vit: *cur_variant_cluster_group) {
                
                delete vit.second;
            }
 
            delete cur_variant_cluster_group;
        }

        delete variant_cluster_group_batch;
    }

    lock_guard<mutex> process_lock(*process_mutex);

    for (auto & variant_cluster_graph: variant_cluster_graphs_local) {

        assert(variant_cluster_graph->getInfo().size() > 0);
        
        variant_cluster_graphs->push_back(variant_cluster_graph);
    }

    for (auto & variant_cluster_group: variant_cluster_groups_local) {

        assert(variant_cluster_group->numberOfVariants() > 0);
        assert(variant_cluster_group->numberOfVariantClusters() > 0);
        assert(variant_cluster_group->numberOfVariantClusterGroupTrees() > 0);

        variant_cluster_groups->push_back(variant_cluster_group);
    }
}


unordered_map<uint, uint> VariantFileParser::getVariantClusterGroupDependencies(unordered_map<uint, VariantCluster*> * variant_cluster_group) {

    unordered_map<uint, uint> variant_cluster_depedencies;

    auto vit_first = variant_cluster_group->begin();
    
    while (vit_first != variant_cluster_group->end()) {

        assert(vit_first->first == vit_first->second->cluster_idx);
        assert(vit_first->second->left_flank <= vit_first->second->right_flank);
        
        auto vit_second = variant_cluster_group->begin();
        auto nested_variant_cluster = variant_cluster_group->end();

        while (vit_second != variant_cluster_group->end()) {

            if (vit_first != vit_second) {

                if ((vit_first->second->left_flank > vit_second->second->left_flank) and (vit_first->second->right_flank < vit_second->second->right_flank)) {

                    if (nested_variant_cluster == variant_cluster_group->end()) {

                        nested_variant_cluster = vit_second;
                    
                    } else if ((vit_second->second->left_flank > nested_variant_cluster->second->left_flank) and (vit_second->second->right_flank < nested_variant_cluster->second->right_flank)) {

                        nested_variant_cluster = vit_second;
                    } 
                
                } else if (!((vit_first->second->left_flank < vit_second->second->left_flank) and (vit_first->second->right_flank > vit_second->second->right_flank))) {

                    assert((vit_first->second->right_flank < vit_second->second->left_flank) or (vit_second->second->right_flank < vit_first->second->left_flank));
                }
            }
    
            vit_second++;
        }

        if (nested_variant_cluster != variant_cluster_group->end()) {  
            
            assert(variant_cluster_depedencies.emplace(vit_first->first, nested_variant_cluster->second->cluster_idx).second);  
        }          

        vit_first++;
    }

    for (auto & vit: variant_cluster_depedencies) {  

        assert(variant_cluster_group->at(vit.second)->contained_clusters.emplace(variant_cluster_group->at(vit.first)->cluster_idx, variant_cluster_group->at(vit.first)->left_flank, variant_cluster_group->at(vit.first)->right_flank).second);
    }

    return variant_cluster_depedencies;
}


template <int kmer_size>
void VariantFileParser::addSequencesToInterclusterRegions(unordered_map<string, string *> & sequences, unordered_set<string> * counted_sequences) {

    for (auto &cit: sequences) {

        auto counted_sequences_insert = counted_sequences->insert(cit.first);

        if (counted_sequences_insert.second) {

            auto sit = cit.second->begin();

            while (sit <= (cit.second->end() - max_chromosome_substring_size)) {

                intercluster_regions.emplace_back(sit, sit + max_chromosome_substring_size, classifyChromosome(cit.first));              
                
                sit += max_chromosome_substring_size - (kmer_size - 1);
            }

            intercluster_regions.emplace_back(sit, cit.second->end(), classifyChromosome(cit.first));              
        }
    }
}


const vector<VariantFileParser::InterClusterRegion> & VariantFileParser::getInterclusterRegions() {

    return intercluster_regions;
}


ulong VariantFileParser::getNumberOfVariants() {

    return number_of_variants;
}


ushort VariantFileParser::getMaxAlternativeAlleles() {

    return max_alternative_alleles;
}


string VariantFileParser::simplifyChromosomeId(const string & chromosome_id) {

    string simple_chromosome_id = chromosome_id;

    transform(simple_chromosome_id.begin(), simple_chromosome_id.end(), simple_chromosome_id.begin(), ::tolower);

    if (simple_chromosome_id.substr(0,3) == "chr") {

        return simple_chromosome_id.substr(3);
    }

    return simple_chromosome_id;
}


Utils::ChromosomeClass VariantFileParser::classifyChromosome(const string & chromosome_id) {

    if (decoy_sequences.find(chromosome_id) != decoy_sequences.end()) {

        return Utils::ChromosomeClass::Decoy;

    } else {

        return classifyGenomeChromosome(chromosome_id);
    }
}


Utils::ChromosomeClass VariantFileParser::classifyGenomeChromosome(const string & chromosome_id) {

    if (chromosome_id == "x") {

        return Utils::ChromosomeClass::X;

    } else if (chromosome_id == "y") {

        return Utils::ChromosomeClass::Y;
    }

    return Utils::ChromosomeClass::Autosomal;
}


template void VariantFileParser::addDecoys<31>(const string);
template void VariantFileParser::addDecoys<39>(const string);
template void VariantFileParser::addDecoys<47>(const string);
template void VariantFileParser::addDecoys<55>(const string);
template void VariantFileParser::addDecoys<63>(const string);

template void VariantFileParser::readVariantFile<31>(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *);
template void VariantFileParser::readVariantFile<39>(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *);
template void VariantFileParser::readVariantFile<47>(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *);
template void VariantFileParser::readVariantFile<55>(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *);
template void VariantFileParser::readVariantFile<63>(const string, vector<VariantClusterGraph*> *, vector<VariantClusterGroup*> *);

