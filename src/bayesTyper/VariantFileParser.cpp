
/*
VariantFileParser.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "VariantFileParser.hpp"
#include "Utils.hpp"
#include "VariantClusterGroup.hpp"
#include "VariantCluster.hpp"
#include "VariantClusterGraph.hpp"
#include "Nucleotide.hpp"


using namespace std;

static const uint variant_cluster_group_batch_size = 1000;

bool InterclusterRegionCompare(const VariantFileParser::InterClusterRegion & first_intercluster_region, const VariantFileParser::InterClusterRegion & second_intercluster_region) {
    
    assert(first_intercluster_region.start_position <= first_intercluster_region.end_position);
    assert(second_intercluster_region.start_position <= second_intercluster_region.end_position);

    return ((first_intercluster_region.end_position - first_intercluster_region.start_position) > (second_intercluster_region.end_position - second_intercluster_region.start_position));
}

VariantFileParser::VariantFileParser(const OptionsContainer & options_container) : num_threads(options_container.getValue<ushort>("threads")), max_allele_length(options_container.getValue<uint>("max-allele-length")), copy_number_variant_threshold(options_container.getValue<float>("copy-number-variant-threshold")) {

    num_variants = 0;
    num_variant_clusters = 0;
    num_variant_cluster_groups = 0;

    allele_type_counter = vector<uint>(static_cast<uchar>(AlleleCount::ALLELE_COUNT_SIZE), 0);
    variant_type_counter = vector<uint>(static_cast<uchar>(VariantCluster::VariantType::VARIANT_TYPE_SIZE), 0);

    intercluster_regions_length = 0;

    prev_chrom_name = "";

    prev_position = -1;
    prev_var_end_position = -1;

    total_num_variants = 0;


    openVariantFile(options_container.getValue<string>("variant-file"));

    for (string variants_line; getline(variants_infile_fstream, variants_line, '\t');) {

        assert(!variants_line.empty());

        if (variants_line.substr(0,1) != "#") {

            total_num_variants++;
        }

        variants_infile_fstream.ignore(numeric_limits<streamsize>::max(), '\n'); 
    }

    variants_infile_fstream.reset();


    openVariantFile(options_container.getValue<string>("variant-file"));

    string header_string;
    assert(getline(variants_infile_fstream, header_string, '\t'));

    assert(!header_string.empty());
    assert(header_string.substr(header_string.size() - 6, 6) == "#CHROM");

    assert(getline(variants_infile_fstream, header_string));
    
    const uint number_of_columns = count(header_string.begin(), header_string.end(), '\t') + 2;

    assert(number_of_columns >= 8);
    has_format = number_of_columns > 8;

    variant_line = vector<string>(6, "");
    updateVariantLine();
}

void VariantFileParser::openVariantFile(const string & variant_filename) {

    assert(!variants_infile.is_open());
    assert(variants_infile_fstream.empty());
    
    if (variant_filename.substr(variant_filename.size() - 7, 7) == ".vcf.gz") {

        variants_infile_fstream.push(boost::iostreams::gzip_decompressor());        
        variants_infile.open(variant_filename, ios::binary);
    
    } else {

        assert((variant_filename.substr(variant_filename.size() - 4, 4) == ".vcf"));
        variants_infile.open(variant_filename);
    }

    if (!variants_infile.is_open()) {

        cerr << "\nERROR: Unable to open file " << variant_filename << "\n" << endl;
        exit(1);
    }

    variants_infile_fstream.push(boost::ref(variants_infile));
    assert(variants_infile_fstream.is_complete());    
}

bool VariantFileParser::updateVariantLine() {

    getline(variants_infile_fstream, variant_line.at(0), '\t');

    for (ushort i = 1; i < 5; i++) {

        getline(variants_infile_fstream, variant_line.at(i), '\t');
    }

    variants_infile_fstream.ignore(numeric_limits<streamsize>::max(), '\t');  
    variants_infile_fstream.ignore(numeric_limits<streamsize>::max(), '\t');

    if (has_format) {

        getline(variants_infile_fstream, variant_line.at(5), '\t');
        variants_infile_fstream.ignore(numeric_limits<streamsize>::max(), '\n');

    } else {

        getline(variants_infile_fstream, variant_line.at(5));
    }

    return variants_infile_fstream.good();
}

void VariantFileParser::addSequenceToInterclusterRegions(const string & chrom_name, const bool is_decoy, const uint start_position, const uint end_position) {

    assert(start_position <= end_position);

    intercluster_regions_length += end_position - start_position + 1;

    if ((end_position - start_position + 1) >= Utils::kmer_size) {

        intercluster_regions.emplace_back(chrom_name, is_decoy, start_position, end_position); 
    }             
}

bool VariantFileParser::constructVariantClusterGroups(InferenceUnit * inference_unit, const uint min_unit_variants, const Chromosomes & chromosomes) { 

    cout << "[" << Utils::getLocalTime() << "] Parsing variants in unit " << inference_unit->index << " ..." << endl; 

    assert(inference_unit->num_variants == 0);
    assert(inference_unit->variant_cluster_groups.empty());

    const uint num_variants_pre = num_variants;
    const uint num_variant_clusters_pre = num_variant_clusters;
    const uint num_variant_cluster_groups_pre = num_variant_cluster_groups;

    ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > variant_cluster_group_queue(Utils::queue_size_thread_scaling * num_threads);

    mutex process_mutex;

    vector<thread> processing_threads;
    processing_threads.reserve(num_threads);

    for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        processing_threads.push_back(thread(&VariantFileParser::processVariantClusterGroupsCallback, this, &(inference_unit->variant_cluster_groups), &variant_cluster_group_queue, &process_mutex, ref(chromosomes)));
    }  

    parseVariants(&variant_cluster_group_queue, min_unit_variants, chromosomes);

    variant_cluster_group_queue.pushedLast();

    for (auto & processing_thread: processing_threads) {

        processing_thread.join();
    } 

    assert(num_variants_pre < num_variants);
    assert(num_variant_clusters_pre < num_variant_clusters);
    assert(num_variant_cluster_groups_pre < num_variant_cluster_groups);

    inference_unit->num_variants = num_variants - num_variants_pre;
    inference_unit->num_variant_clusters = num_variant_clusters - num_variant_clusters_pre;

    assert(inference_unit->variant_cluster_groups.size() == (num_variant_cluster_groups - num_variant_cluster_groups_pre));

    assert(num_variant_cluster_groups <= num_variant_clusters);
    assert(num_variant_clusters <= num_variants);
    
    assert(num_variants <= total_num_variants);
    assert((total_num_variants != num_variants) == variants_infile_fstream.good());

    cout << "[" << Utils::getLocalTime() << "] Parsed " << inference_unit->num_variants << " variants as " << inference_unit->num_variant_clusters << " clusters" << endl; 

    return (num_variants == total_num_variants);
}

void VariantFileParser::parseVariants(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue, const uint min_unit_variants, const Chromosomes & chromosomes) {

    bool is_first_unit_variant = true;
    uint unit_variant_counter = 0;
    
    string cur_chrom_name = "";

    auto chromosomes_it = chromosomes.cend();

    if (prev_chrom_name != "") {

        chromosomes_it = chromosomes.find(prev_chrom_name);
    }

    int cur_position = 0;
    int cur_group_end_position = prev_var_end_position;     

    uint variant_cluster_group_batch_complexity = 0;

    vector<unordered_map<uint, VariantCluster*> * > * variant_cluster_group_batch = new vector<unordered_map<uint, VariantCluster*> * >();  
    variant_cluster_group_batch->reserve(variant_cluster_group_batch_size);
    
    unordered_map<uint, VariantCluster*> * variant_cluster_group = new unordered_map<uint, VariantCluster*>();

    map<uint, VariantCluster*> variant_cluster_group_flanks;
    list<unordered_set<uint> > variant_cluster_group_merge_sets;

    set<uint> variant_depedencies;

    while (is_first_unit_variant or updateVariantLine()) {

        is_first_unit_variant = false;
        assert(!variant_line.at(0).empty());

        cur_chrom_name = variant_line.at(0);
        cur_position = stoi(variant_line.at(1)) - 1;

        if (cur_chrom_name != prev_chrom_name) {

            if (prev_chrom_name != "") {

                auto prev_chromosomes_it = chromosomes.find(prev_chrom_name);

                if (prev_chromosomes_it != chromosomes.cend()) {

                    processVariantClusterGroups(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);

                    assert((prev_var_end_position + 1) <= static_cast<int>(prev_chromosomes_it->second.size()));
                    addSequenceToInterclusterRegions(prev_chrom_name, chromosomes.isDecoy(prev_chrom_name), prev_var_end_position + 1, prev_chromosomes_it->second.size() - 1);

                    assert(intercluster_chromosomes.insert(prev_chrom_name).second);

                } else {

                    assert(prev_var_end_position == -1);
                    assert(cur_group_end_position == -1);
                }
            }

            assert(!cur_chrom_name.empty());

            chromosomes_it = chromosomes.find(cur_chrom_name);

            prev_var_end_position = -1;
            cur_group_end_position = -1;

            variant_depedencies.clear();

        } else {

            assert(prev_position < cur_position);
        }

        prev_chrom_name = cur_chrom_name;

        auto vit = variant_depedencies.begin();
        
        while (vit != variant_depedencies.end()) {
            
            if (static_cast<int>(*vit) >= cur_position) {

                break;
            }
                
            variant_depedencies.erase(vit);
            vit = variant_depedencies.begin();
        }

        if ((unit_variant_counter >= min_unit_variants) and ((cur_position - cur_group_end_position) >= static_cast<int>(Utils::kmer_size))) {

            assert(variant_depedencies.empty());
            break;
        }

        prev_position = cur_position;

        assert(variant_line.at(3).find(",") == string::npos);        

        string var_ref_seq = variant_line.at(3);

        transform(var_ref_seq.begin(), var_ref_seq.end(), var_ref_seq.begin(), ::toupper);

        vector<string> alt_alleles;
        boost::split(alt_alleles, variant_line.at(4), boost::is_any_of(","));
        assert(alt_alleles.size() < (Utils::ushort_overflow - 2));

        vector<string> origin_allele_att;
        origin_allele_att.reserve(alt_alleles.size());

        auto origin_att_str = getInfoAttributeString(variant_line.at(5), "ACO");

        if (origin_att_str.second) {

            boost::split(origin_allele_att, origin_att_str.first, boost::is_any_of(","));        
        
        } else {

            origin_allele_att = vector<string>(alt_alleles.size(), "");          
        }

        assert(origin_allele_att.size() == alt_alleles.size());

        VariantCluster::Variant cur_variant(variant_line.at(2), !variant_depedencies.empty());

        num_variants += 1;
        unit_variant_counter += 1;

        if (alt_alleles.back() == "*") {

            assert(cur_variant.has_dependency);
            alt_alleles.pop_back();
        }

        assert(!alt_alleles.empty());
        allele_type_counter.at(uchar(AlleleCount::Total)) += alt_alleles.size();

        if (chromosomes.isDecoy(cur_chrom_name)) {

            assert(variant_cluster_group->empty());

            allele_type_counter.at(uchar(AlleleCount::Excluded_decoy)) += alt_alleles.size();
            variant_type_counter.at(uchar(VariantCluster::VariantType::Unsupported))++;

            continue;
        }

        if (chromosomes_it == chromosomes.cend()) {

            assert(variant_cluster_group->empty());

            allele_type_counter.at(uchar(AlleleCount::Excluded_genome)) += alt_alleles.size();
            variant_type_counter.at(uchar(VariantCluster::VariantType::Unsupported))++;

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

        string gen_ref_seq = chromosomes_it->second.substr(cur_position, variant_line.at(3).size()); 
        transform(gen_ref_seq.begin(), gen_ref_seq.end(), gen_ref_seq.begin(), ::toupper);

        assert(gen_ref_seq.size() == var_ref_seq.size());

        if (var_ref_seq.compare(gen_ref_seq) != 0) {

            allele_type_counter.at(uchar(AlleleCount::Excluded_match)) += alt_alleles.size();
            is_excluded = true;
        }

        if (cur_position < (static_cast<int>(Utils::kmer_size) - 1)) {

            allele_type_counter.at(uchar(AlleleCount::Excluded_end)) += alt_alleles.size();
            is_excluded = true;
        }

        unordered_set<ushort> excluded_alleles;

        if (!is_excluded) {

            for (ushort i = 0; i < alt_alleles.size(); i++) {

                if ((cur_position + ref_alleles.at(i).size() - 1 + Utils::kmer_size) > chromosomes_it->second.size()) {

                    allele_type_counter.at(uchar(AlleleCount::Excluded_end))++;
                    assert(excluded_alleles.insert(i).second);

                } else if ((ref_alleles.at(i).size() > max_allele_length) or (alt_alleles.at(i).size() > max_allele_length)) {

                    allele_type_counter.at(uchar(AlleleCount::Excluded_length))++;
                    assert(excluded_alleles.insert(i).second);

                } else {

                    variant_depedencies.insert(cur_position + ref_alleles.at(i).size() - 1);
                }
            }
        }

        assert(excluded_alleles.size() <= alt_alleles.size());

        if (is_excluded or excluded_alleles.size() == alt_alleles.size()) {

            variant_type_counter.at(uchar(VariantCluster::VariantType::Unsupported))++;
            
            continue;     
        }

        assert(!variant_depedencies.empty());
        assert(prev_var_end_position <= cur_group_end_position);

        if ((cur_position - cur_group_end_position) >= static_cast<int>(Utils::kmer_size)) {

            processVariantClusterGroups(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);
        } 

        if (cur_position > (prev_var_end_position + 1)) {

            addSequenceToInterclusterRegions(cur_chrom_name, chromosomes.isDecoy(cur_chrom_name), prev_var_end_position + 1, cur_position - 1);
        }

        set<uint> cur_end_positions;

        for (ushort alt_allele_idx = 0; alt_allele_idx < alt_alleles.size(); alt_allele_idx++) {

            assert(!ref_alleles.at(alt_allele_idx).empty());            
            assert(!alt_alleles.at(alt_allele_idx).empty());

            if (excluded_alleles.count(alt_allele_idx) == 0) {

                addAlternativeAllele(&cur_variant, ref_alleles.at(alt_allele_idx), alt_alleles.at(alt_allele_idx), origin_allele_att.at(alt_allele_idx));
                
                uint ref_copy_number_variant_length = copyNumberVariantLength(ref_alleles.at(alt_allele_idx), chromosomes_it->second, cur_position + ref_alleles.at(alt_allele_idx).size());
                uint alt_copy_number_variant_length = copyNumberVariantLength(alt_alleles.at(alt_allele_idx), chromosomes_it->second, cur_position + ref_alleles.at(alt_allele_idx).size());

                cur_end_positions.insert(cur_position + ref_alleles.at(alt_allele_idx).size() - 1);
                cur_group_end_position = max(cur_group_end_position, static_cast<int>(cur_position + ref_alleles.at(alt_allele_idx).size() - 1 + max(ref_copy_number_variant_length, alt_copy_number_variant_length)));
            }
        }

        assert(!cur_end_positions.empty());
        assert(!cur_variant.alt_alleles.empty());
        assert((cur_variant.alt_alleles.size() + excluded_alleles.size()) == alt_alleles.size());

        prev_var_end_position = max(prev_var_end_position, static_cast<int>(*cur_end_positions.rbegin()));
        assert(cur_group_end_position >= prev_var_end_position);

        clusterVariants(cur_variant, cur_position, cur_end_positions, cur_chrom_name, &variant_cluster_group_flanks, variant_cluster_group, &variant_cluster_group_merge_sets);

        variant_cluster_group_batch_complexity++;   
        variant_type_counter.at(uchar(cur_variant.type))++; 
    }

    variant_cluster_group_batch_complexity = variant_cluster_group_batch_size;

    processVariantClusterGroups(variant_cluster_group_queue, &variant_cluster_group_batch, &variant_cluster_group_batch_complexity, &variant_cluster_group, &variant_cluster_group_merge_sets, &variant_cluster_group_flanks);

    assert(variant_cluster_group_batch->empty());
    delete variant_cluster_group_batch;

    assert(variant_cluster_group->empty());
    delete variant_cluster_group;

    if (total_num_variants == num_variants) {

        if (chromosomes_it != chromosomes.cend()) {
            
            assert(cur_chrom_name == chromosomes_it->first);
            addSequenceToInterclusterRegions(cur_chrom_name, chromosomes.isDecoy(cur_chrom_name), prev_var_end_position + 1, chromosomes_it->second.size() - 1);
     
            assert(intercluster_chromosomes.insert(cur_chrom_name).second);
            assert(cur_chrom_name == prev_chrom_name);
        }

        chromosomes_it = chromosomes.cbegin();

        while (chromosomes_it != chromosomes.cend()) {

            if (intercluster_chromosomes.insert(chromosomes_it->first).second) {

                addSequenceToInterclusterRegions(chromosomes_it->first, chromosomes.isDecoy(chromosomes_it->first), 0, chromosomes_it->second.size() - 1);
            }   

            chromosomes_it++;                     
        }

        assert(intercluster_chromosomes.size() == chromosomes.getTotalCount());
    }
}

pair<string, bool> VariantFileParser::getInfoAttributeString(const string & info_str, const string & att_name) {

    stringstream info_ss(info_str);
    string att_str = "";

    while (getline(info_ss, att_str, ';')) {

        if ((att_str.substr(0, att_name.size()) == att_name) and (att_str.substr(att_name.size(), 1) == "=")) {

            return make_pair(att_str.substr(att_name.size() + 1), true);
        }
    }

    return make_pair("", false);
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

    assert(!ref_allele->empty());
    assert(!alt_allele->empty());
}

void VariantFileParser::addAlternativeAllele(VariantCluster::Variant * cur_variant, const string & ref_allele, const string & alt_allele, const string & origin_att) {

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
    cur_variant->num_redundant_nucleotides = min(cur_variant->num_redundant_nucleotides, identical_left_nucleotides);

    cur_variant->alt_alleles.emplace_back(ref_allele.size(), alt_allele, origin_att);

    VariantCluster::VariantType variant_type = classifyAllele(ref_allele.size() - identical_left_nucleotides, alt_allele.size() - identical_left_nucleotides);

    if (cur_variant->type == VariantCluster::VariantType::Unsupported) {

        cur_variant->type = variant_type;

    } else {

        if (cur_variant->type != variant_type) {

            cur_variant->type = VariantCluster::VariantType::Mixture;
        }
    }
}

VariantCluster::VariantType VariantFileParser::classifyAllele(const int reference_size, const int allele_size) {

    assert((reference_size > 0) or (allele_size > 0));

    if ((reference_size == 1) and (allele_size == 1)) {

        return VariantCluster::VariantType::SNV;
    
    } else if ((reference_size == 0) or (allele_size == 0)) {

        if ((allele_size - reference_size) > 0) {

            return VariantCluster::VariantType::Insertion;

        } else {

            assert((allele_size - reference_size) < 0);
            
            return VariantCluster::VariantType::Deletion;
        }
    }

    return VariantCluster::VariantType::Complex;
}

uint VariantFileParser::copyNumberVariantLength(const string & allele_sequence, const string & chrom_sequence, const uint chrom_start_position) {

    assert(chrom_start_position <= chrom_sequence.size());

    uint copy_number_variant_length = 0;

    if (allele_sequence.size() < Utils::kmer_size) {

        return copy_number_variant_length;
    }

    unordered_set<bitset<Utils::kmer_size * 2> > allele_kmers;
    KmerPair<Utils::kmer_size> kmer_pair = KmerPair<Utils::kmer_size>();

    auto allele_sequence_it = allele_sequence.begin();
    assert(allele_sequence_it != allele_sequence.end());

    while (allele_sequence_it != allele_sequence.end()) {

        if (kmer_pair.move(Nucleotide::ntToBit<1>(*allele_sequence_it))) {

            allele_kmers.emplace(kmer_pair.getLexicographicalLowestKmer());
        }

        allele_sequence_it++;
    } 

    if (allele_kmers.empty()) {

        return copy_number_variant_length;        
    }

    uint chrom_window_end_position = min(chrom_start_position + copy_number_variant_length + static_cast<uint>(allele_sequence.size()), static_cast<uint>(chrom_sequence.size()));

    while (true) {

        kmer_pair.reset();

        uint num_bases = 0;
        uint num_identical_kmers = 0;

        pair<double, uint> highest_scoring_window(0, 0);

        for (uint chrom_position = chrom_start_position + copy_number_variant_length; chrom_position < chrom_window_end_position; chrom_position++) {

            if (kmer_pair.move(Nucleotide::ntToBit<1>(chrom_sequence.at(chrom_position)))) {

                if (allele_kmers.count(kmer_pair.getLexicographicalLowestKmer()) > 0) {

                    num_identical_kmers++;
                }
            }

            num_bases++;

            if (num_identical_kmers > 0) {

                assert(num_bases >= static_cast<uint>(Utils::kmer_size));
                double identical_kmer_fraction = num_identical_kmers/static_cast<double>(num_bases - Utils::kmer_size + 1);

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

        if (chrom_window_end_position == chrom_sequence.size()) {

            break;
        }

        chrom_window_end_position = min(chrom_start_position + copy_number_variant_length + static_cast<uint>(allele_sequence.size()), static_cast<uint>(chrom_sequence.size()));
    }

    return copy_number_variant_length;
}

void VariantFileParser::clusterVariants(VariantCluster::Variant & cur_variant, const uint cur_position, const set<uint> cur_end_positions, const string & cur_chrom_name, map<uint, VariantCluster*> * variant_cluster_group_flanks, unordered_map<uint, VariantCluster*> * variant_cluster_group, list<unordered_set<uint> > * variant_cluster_group_merge_sets) {

    if (!variant_cluster_group_flanks->empty()) {

        auto lit = variant_cluster_group_flanks->begin();            
        assert(!variant_cluster_group->empty());

        while (static_cast<int>(cur_position - lit->first) >= static_cast<int>(Utils::kmer_size)) {

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

        if ((abs(static_cast<int>(cur_position - vit->first)) + 1) <= static_cast<int>(Utils::kmer_size)) {

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

            if ((abs(static_cast<int>(cit - vit->first)) + 1) <= static_cast<int>(Utils::kmer_size)) {

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
        variant_cluster->cluster_idx = variant_cluster_group->size();
        variant_cluster->variants = map<uint, VariantCluster::Variant>();
        variant_cluster->left_flank = cur_position;
        variant_cluster->right_flank = *cur_end_positions.rbegin();

        variant_cluster->chrom_name = cur_chrom_name;

        assert(variant_cluster->variants.insert(pair<uint, VariantCluster::Variant>(cur_position, cur_variant)).second);

        for (auto &cit: cur_end_positions) {

            assert(variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cit, variant_cluster)).second);
        }

        if (*cur_end_positions.rbegin() - cur_position >= Utils::kmer_size) {

            if (!variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cur_position, variant_cluster)).second) {

                assert(cur_end_positions.count(cur_position) > 0);   
            }          
        }

        assert(variant_cluster_group->insert({variant_cluster_group->size(), variant_cluster}).second);
        assert(variant_cluster_group->size() < Utils::uint_overflow);

    } else {

        assert(variant_cluster_group_first->chrom_name == cur_chrom_name);

        assert(variant_cluster_group_first->variants.insert(pair<uint, VariantCluster::Variant>(cur_position, cur_variant)).second);
        variant_cluster_group_first->right_flank = max(*cur_end_positions.rbegin(), variant_cluster_group_first->right_flank);

        for (auto &cit: cur_end_positions) {

            if (!variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cit, variant_cluster_group_first)).second) {

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
      
        if (*cur_end_positions.rbegin() - cur_position >= Utils::kmer_size) {

            if (!variant_cluster_group_flanks->insert(pair<uint, VariantCluster *>(cur_position, variant_cluster_group_first)).second) {

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

            assert(!sit->empty());

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

void VariantFileParser::processVariantClusterGroups(ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue, vector<unordered_map<uint, VariantCluster*> * > ** variant_cluster_group_batch, uint * variant_cluster_group_batch_complexity, unordered_map<uint, VariantCluster*> ** variant_cluster_group, list<unordered_set<uint> > * variant_cluster_group_merge_sets, map<uint, VariantCluster*> * variant_cluster_group_flanks) {

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

            assert(variant_cluster_group->at(first_cluster)->cluster_idx != variant_cluster_group->at(*cur_cluster)->cluster_idx);
            assert(variant_cluster_group->at(first_cluster)->chrom_name == variant_cluster_group->at(*cur_cluster)->chrom_name);
            
            assert(variant_cluster_group->at(first_cluster)->contained_clusters.empty());
            assert(variant_cluster_group->at(*cur_cluster)->contained_clusters.empty());

            variant_cluster_group->at(first_cluster)->left_flank = min(variant_cluster_group->at(first_cluster)->left_flank, variant_cluster_group->at(*cur_cluster)->left_flank);
            variant_cluster_group->at(first_cluster)->right_flank = max(variant_cluster_group->at(first_cluster)->right_flank, variant_cluster_group->at(*cur_cluster)->right_flank);

            for (auto & variant: variant_cluster_group->at(*cur_cluster)->variants) {

                assert(variant_cluster_group->at(first_cluster)->variants.insert(variant).second);
            }


            delete variant_cluster_group->at(*cur_cluster);
            assert(variant_cluster_group->erase(*cur_cluster));

            cur_cluster++;
        }
    }
}

void VariantFileParser::processVariantClusterGroupsCallback(vector<VariantClusterGroup*> * variant_cluster_groups, ProducerConsumerQueue<vector<unordered_map<uint, VariantCluster*> * > * > * variant_cluster_group_queue, mutex * process_mutex, const Chromosomes & chromosomes) {

    vector<VariantClusterGroup*> variant_cluster_groups_local;
    vector<unordered_map<uint, VariantCluster *> * > * variant_cluster_group_batch = nullptr;

    while (variant_cluster_group_queue->pop(&variant_cluster_group_batch)) {

        for (auto & cur_variant_cluster_group: *variant_cluster_group_batch) {

            assert(!cur_variant_cluster_group->empty());
            assert(cur_variant_cluster_group->size() < Utils::uint_overflow);

            auto cur_variant_cluster_group_dependencies = getVariantClusterGroupDependencies(cur_variant_cluster_group);

            vector<VariantCluster *> cur_variant_clusters;
            vector<VariantClusterGraph *> cur_variant_cluster_graphs;

            cur_variant_clusters.reserve(cur_variant_cluster_group->size());
            cur_variant_cluster_graphs.reserve(cur_variant_cluster_group->size());

            for (auto & variant_cluster: *cur_variant_cluster_group) {

                assert(variant_cluster.first == variant_cluster.second->cluster_idx);

                assert(cur_variant_cluster_group->begin()->second->chrom_name == variant_cluster.second->chrom_name);

                assert(!variant_cluster.second->variants.empty());
                assert(variant_cluster.second->variants.size() < Utils::ushort_overflow);

                assert(variant_cluster.second->left_flank == variant_cluster.second->variants.begin()->first);
                assert(variant_cluster.second->right_flank >= variant_cluster.second->variants.rbegin()->first);

                cur_variant_clusters.push_back(variant_cluster.second);
                cur_variant_cluster_graphs.emplace_back(new VariantClusterGraph(variant_cluster.second, chromosomes.find(variant_cluster.second->chrom_name)->second));
            }

            variant_cluster_groups_local.emplace_back(new VariantClusterGroup(cur_variant_clusters, cur_variant_cluster_graphs, cur_variant_cluster_group_dependencies));

            for (auto &vit: *cur_variant_cluster_group) {
                
                delete vit.second;
            }
 
            delete cur_variant_cluster_group;
        }

        delete variant_cluster_group_batch;
    }

    lock_guard<mutex> process_lock(*process_mutex);

    for (auto & variant_cluster_group: variant_cluster_groups_local) {

        assert(variant_cluster_group->numberOfVariants() > 0);
        assert(variant_cluster_group->numberOfVariantClusters() > 0);
        assert(variant_cluster_group->numberOfVariantClusterGroupTrees() > 0);

        num_variant_clusters += variant_cluster_group->numberOfVariantClusters();
        num_variant_cluster_groups++;

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

string VariantFileParser::getVariantStatsString(const uint num_units) {

    stringstream variant_stats_ss;

    variant_stats_ss << "[" << Utils::getLocalTime() << "] Parsed " << allele_type_counter.at(uchar(AlleleCount::Total)) << " alternative alleles across " << num_units << " units and excluded:\n" << endl; 
    variant_stats_ss << "\t- Alleles on a decoy sequence: " << allele_type_counter.at(uchar(AlleleCount::Excluded_decoy)) << endl;  
    variant_stats_ss << "\t- Alleles on a chromosome not in the reference genome: " << allele_type_counter.at(uchar(AlleleCount::Excluded_genome)) << endl;  
    variant_stats_ss << "\t- Alleles with a reference not equal to the reference genome sequence: " << allele_type_counter.at(uchar(AlleleCount::Excluded_match)) << endl;
    variant_stats_ss << "\t- Alleles within " << Utils::kmer_size << " bases of a chromosome end: " << allele_type_counter.at(uchar(AlleleCount::Excluded_end)) << endl;
    variant_stats_ss << "\t- Alleles longer than " << max_allele_length << " bases: " << allele_type_counter.at(uchar(AlleleCount::Excluded_length)) << endl;

    variant_stats_ss << "\n[" << Utils::getLocalTime() << "] Out of " << num_variants << " variants:\n" << endl; 
    variant_stats_ss << "\t- Single nucleotide variant: " << variant_type_counter.at(uchar(VariantCluster::VariantType::SNV)) << endl;
    variant_stats_ss << "\t- Insertion: " << variant_type_counter.at(uchar(VariantCluster::VariantType::Insertion)) << endl;
    variant_stats_ss << "\t- Deletion: " << variant_type_counter.at(uchar(VariantCluster::VariantType::Deletion)) << endl;
    variant_stats_ss << "\t- Complex: " << variant_type_counter.at(uchar(VariantCluster::VariantType::Complex)) << endl;
    variant_stats_ss << "\t- Mixture: " << variant_type_counter.at(uchar(VariantCluster::VariantType::Mixture)) << endl;
    variant_stats_ss << "\t- Unsupported (excluded): " << variant_type_counter.at(uchar(VariantCluster::VariantType::Unsupported)) << endl;

    variant_stats_ss << "\n[" << Utils::getLocalTime() << "] Variants merged into " << num_variant_clusters << " clusters and further into " << num_variant_cluster_groups << " groups" << endl;

    return variant_stats_ss.str();
}

void VariantFileParser::sortInterclusterRegions() {

    sort(intercluster_regions.begin(), intercluster_regions.end(), InterclusterRegionCompare);
}

void VariantFileParser::writeInterclusterRegions(const string & output_prefix) {

    ofstream regions_outfile(output_prefix + ".txt.gz", ios::binary);

    if (!regions_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << output_prefix + ".txt.gz" << "\n" << endl;
        exit(1);
    }

    boost::iostreams::filtering_ostream regions_outfile_fstream;

    regions_outfile_fstream.push(boost::iostreams::gzip_compressor());
    regions_outfile_fstream.push(boost::ref(regions_outfile));

    assert(regions_outfile_fstream.is_complete());    
    
    for (auto & intercluster_region: intercluster_regions) {

        regions_outfile_fstream << intercluster_region.chrom_name << "\t" << intercluster_region.is_decoy << "\t" << intercluster_region.start_position << "\t" << intercluster_region.end_position << endl;
    }
}

const vector<VariantFileParser::InterClusterRegion> & VariantFileParser::getInterclusterRegions() {

    return intercluster_regions;
}

ulong VariantFileParser::getInterclusterRegionLength() {

    return intercluster_regions_length;
}

ulong VariantFileParser::getNumberOfInterclusterRegionKmers() {

    assert(intercluster_regions_length > (intercluster_regions.size() * (Utils::kmer_size - 1)));
    return (intercluster_regions_length - (intercluster_regions.size() * (Utils::kmer_size - 1)));
}

uint VariantFileParser::getNumberOfVariants() {

    return total_num_variants;
}


