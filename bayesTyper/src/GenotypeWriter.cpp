
/*
GenotypeWriter.cpp - This file is part of BayesTyper (v0.9)


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
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <unordered_map> 
#include <unordered_set> 
#include <math.h>
#include <algorithm>
#include <thread>
#include <random>

#include "boost/algorithm/string.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filter/gzip.hpp"
#include "boost/functional/hash.hpp"
#include "ProducerConsumerQueue.hpp"

#include "GenotypeWriter.hpp"
#include "Genotypes.hpp"
#include "Sample.hpp"
#include "Utils.hpp"
#include "VariantKmerStats.hpp"
#include "VariantInfo.hpp"
#include "VariantFileParser.hpp"
#include "OptionsContainer.hpp"
#include "ChromosomePloidy.hpp"

using namespace std;

static const string format_column = "\tGT:GPP:MAP:NK:NOK:NUK:AFF";
static const string empty_variant_sample = ".:.:.:.:.:.";

GenotypeWriter::GenotypeWriter(const vector<Sample> & samples, const OptionsContainer & options_container) : num_samples(samples.size()), output_prefix(options_container.getValue<string>("output-prefix")), num_gibbs_samples(options_container.getValue<ushort>("gibbs-samples") * options_container.getValue<ushort>("number-of-gibbs-chains")), min_observed_kmers(options_container.getValue<ushort>("minimum-number-of-observed-informative-kmers")), chromosome_ploidy(samples) {

    genotype_queue = new ProducerConsumerQueue<vector<Genotypes *> *>(Utils::queue_size_scaling * options_container.getValue<ushort>("threads"));
    writer_threads.push_back(thread(&GenotypeWriter::writeTemporaryFile, this));
}


void GenotypeWriter::writeTemporaryFile() {

    string output_file_name = output_prefix;
    output_file_name.append("_tmp.txt");

    ofstream output_file(output_file_name);
    assert(output_file.is_open());

    vector<Genotypes *> * genotyped_variants = nullptr;

    while (genotype_queue->pop(&genotyped_variants)) {

        for (auto &genotyped_variant: *genotyped_variants) {

            assert(genotyped_variant);
            assert(genotyped_variant->variant_info.num_alleles > 1);
            assert(genotyped_variant->variant_info.num_alleles < Utils::ushort_overflow);

            vector<vector<string> > estimated_genotypes(num_samples);
            
            vector<vector<float> > genotype_posteriors;
            genotype_posteriors.reserve(num_samples);

            vector<vector<float> > allele_posteriors;
            allele_posteriors.reserve(num_samples);

            AlleleCounts allele_counts(genotyped_variant->variant_info.num_alleles + genotyped_variant->variant_info.has_dependency - 1);
            auto allele_filter_flags = getAlleleFilterFlags(genotyped_variant->variant_kmer_stats, genotyped_variant->variant_info);
            
            estimateGenotypesAndAlleleCounts(genotyped_variant->posteriors, genotyped_variant->variant_info, &estimated_genotypes, &genotype_posteriors, &allele_posteriors, &allele_counts, &allele_filter_flags);
            vector<float> allele_call_probabilities = calculateAlleleCallProbabilities(allele_posteriors, genotyped_variant->variant_info);
            
            output_file << genotyped_variant->variant_info.variant_id << "\t" << genotyped_variant->variant_info.has_dependency << "\t";

            writeQualityAndFilter(&output_file, genotyped_variant->filter_status, allele_call_probabilities, genotyped_variant->variant_info);
            writeAlleleCounts(&output_file, allele_counts);
            writeAlleleCallProbabilities(&output_file, allele_call_probabilities);

            output_file << ";VT=" << genotyped_variant->variant_info.variant_type << ";VCS=" << genotyped_variant->variant_cluster_size << ";VCI=" << genotyped_variant->variant_cluster_id;

            output_file << ";VCGS=" << genotyped_variant->variant_cluster_group_size << ";VCGI=" << genotyped_variant->variant_cluster_group_id;

            output_file << ";HC=" << genotyped_variant->num_candidates;

            if (genotyped_variant->has_complex_region) {

                output_file << ";HCR";
            }

            if (genotyped_variant->has_redundant_sequence) {

                output_file << ";HRS";
            }

            assert(genotyped_variant->variant_info.num_alleles > genotyped_variant->variant_info.excluded_alternative_alleles.size());

            if (!(genotyped_variant->variant_info.excluded_alternative_alleles.empty())) {

                output_file << ";AE=";
                
                sort(genotyped_variant->variant_info.excluded_alternative_alleles.begin(), genotyped_variant->variant_info.excluded_alternative_alleles.end());
                auto eait = genotyped_variant->variant_info.excluded_alternative_alleles.begin(); 
                
                assert(*eait > 0);    
                assert(genotyped_variant->covered_alleles.count(*eait) < 1);
                    
                output_file << *eait;
                eait++;

                while (eait != genotyped_variant->variant_info.excluded_alternative_alleles.end()) {

                    assert(genotyped_variant->covered_alleles.count(*eait) < 1);

                    output_file << "," << *eait;
                    eait++;
                }
            }
            
            bool is_first_not_covered = true;
            assert(genotyped_variant->covered_alleles.size() <= genotyped_variant->variant_info.num_alleles);

            for (ushort allele_idx = 0; allele_idx < genotyped_variant->variant_info.num_alleles; allele_idx++) {

                if (find(genotyped_variant->variant_info.excluded_alternative_alleles.begin(), genotyped_variant->variant_info.excluded_alternative_alleles.end(), allele_idx) == genotyped_variant->variant_info.excluded_alternative_alleles.end()) {

                    if (genotyped_variant->covered_alleles.count(allele_idx) < 1) {

                        if (is_first_not_covered) {

                            output_file << ";ANC=" << allele_idx;
                            is_first_not_covered = false;
                        
                        } else {

                            output_file << "," << allele_idx;
                        }
                    }
                }
            }

            writeSamples(&output_file, estimated_genotypes, genotype_posteriors, allele_posteriors, genotyped_variant->variant_kmer_stats, genotyped_variant->variant_info, &allele_filter_flags);

            output_file << "\n";

            delete genotyped_variant;
        }

        delete genotyped_variants;
    }

    output_file.close();
}


uint GenotypeWriter::genotypeToOneDimensionalIndex(const pair<ushort, ushort> genotype, const ushort num_alleles) {

    assert(genotype.first <= genotype.second);
    return (genotype.second + num_alleles * genotype.first - genotype.first * (genotype.first + 1) / 2);
}


vector<vector<string> > GenotypeWriter::getAlleleFilterFlags(const VariantKmerStats & variant_kmer_stats, const VariantInfo & variant_info) {

    vector<vector<string> > allele_filter_flags(num_samples, vector<string>(variant_info.num_alleles + variant_info.has_dependency, "P"));

    assert(variant_kmer_stats.mean_num_kmers.size() == num_samples);
    assert(variant_kmer_stats.mean_num_observed_kmers.size() == num_samples);
    assert(variant_kmer_stats.mean_num_unique_kmers.size() == num_samples);

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

        assert(variant_kmer_stats.mean_num_kmers.at(sample_idx).size() == variant_info.num_alleles);
        assert(variant_kmer_stats.mean_num_observed_kmers.at(sample_idx).size() == variant_info.num_alleles);
        assert(variant_kmer_stats.mean_num_unique_kmers.at(sample_idx).size() == variant_info.num_alleles);

        for (ushort allele_idx = 0; allele_idx < variant_info.num_alleles; allele_idx++) {

            if (!(Utils::floatCompare(variant_kmer_stats.mean_num_observed_kmers.at(sample_idx).at(allele_idx), -1))) {

                assert(variant_kmer_stats.mean_num_observed_kmers.at(sample_idx).at(allele_idx) >= 0);

                if (variant_kmer_stats.mean_num_observed_kmers.at(sample_idx).at(allele_idx) < min_observed_kmers) {

                    allele_filter_flags.at(sample_idx).at(allele_idx) = "F";
                }             
            }
        }
    }

    return allele_filter_flags;
}


void GenotypeWriter::estimateGenotypesAndAlleleCounts(const PosteriorContainer & posteriors, const VariantInfo & variant_info, vector<vector<string> > * estimated_genotypes, vector<vector<float> > * genotype_posteriors, vector<vector<float> > * allele_posteriors, AlleleCounts * allele_counts, vector<vector<string> > * allele_filter_flags) {

    const uint num_all_alleles = variant_info.num_alleles + variant_info.has_dependency;
    const uint num_potential_genotypes = (num_all_alleles * (num_all_alleles - 1)) / 2 + num_all_alleles;

    assert(num_samples == posteriors.size());
    assert(num_samples == estimated_genotypes->size());

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

        if (posteriors.at(sample_idx).first == Utils::Ploidy::Diploid) {

            genotype_posteriors->emplace_back(vector<float>(num_potential_genotypes, 0));
            allele_posteriors->emplace_back(vector<float>(num_all_alleles, 0));

            assert(posteriors.at(sample_idx).second.size() <= num_potential_genotypes);
        
        } else if (posteriors.at(sample_idx).first == Utils::Ploidy::Haploid) {

            genotype_posteriors->emplace_back(vector<float>(num_all_alleles, 0));
            allele_posteriors->emplace_back(vector<float>(num_all_alleles, 0));

            assert(posteriors.at(sample_idx).second.size() <= num_all_alleles);
        
        } else {

            assert(posteriors.at(sample_idx).first == Utils::Ploidy::Null);

            genotype_posteriors->emplace_back(vector<float>());
            allele_posteriors->emplace_back(vector<float>());
            allele_filter_flags->at(sample_idx).clear();

            assert(posteriors.at(sample_idx).second.size() == 1);
        }

        vector<pair<ushort, ushort> > max_genotypes;
        max_genotypes.reserve(posteriors.at(sample_idx).second.size());

        double max_num_gibbs_samples = 0;
        double sum_num_gibbs_samples = 0;

        for (auto &pit: posteriors.at(sample_idx).second) {

            pair<ushort, ushort> cur_genotype = pit.first;

            if (cur_genotype.first == Utils::ushort_overflow) {

                cur_genotype.first = num_all_alleles - 1;
                assert((posteriors.at(sample_idx).first == Utils::Ploidy::Null) or variant_info.has_dependency);
            } 

            if (cur_genotype.second == Utils::ushort_overflow) {

                cur_genotype.second = num_all_alleles - 1;
                assert((posteriors.at(sample_idx).first != Utils::Ploidy::Diploid) or variant_info.has_dependency);
            } 

            if (posteriors.at(sample_idx).first == Utils::Ploidy::Diploid) {

                uint genotype_number = genotypeToOneDimensionalIndex(cur_genotype, num_all_alleles);
                
                assert(genotype_number < Utils::uint_overflow);
                assert(Utils::floatCompare(genotype_posteriors->at(sample_idx).at(genotype_number), 0));

                genotype_posteriors->at(sample_idx).at(genotype_number) = pit.second;
                allele_posteriors->at(sample_idx).at(cur_genotype.first) += pit.second;

                if (cur_genotype.first != cur_genotype.second) {
                    
                    allele_posteriors->at(sample_idx).at(cur_genotype.second) += pit.second;
                }

            } else if (posteriors.at(sample_idx).first == Utils::Ploidy::Haploid) {

                assert(cur_genotype.second == (num_all_alleles - 1));
                assert(Utils::floatCompare(genotype_posteriors->at(sample_idx).at(cur_genotype.first), 0));
                
                genotype_posteriors->at(sample_idx).at(cur_genotype.first) = pit.second;
                allele_posteriors->at(sample_idx).at(cur_genotype.first) += pit.second;
            }

            assert(cur_genotype.first <= cur_genotype.second);

            sum_num_gibbs_samples += pit.second;

            if (Utils::doubleCompare(pit.second, max_num_gibbs_samples)) {

                max_genotypes.push_back(cur_genotype);

            } else if (pit.second > max_num_gibbs_samples) {

                max_genotypes.clear();
                max_genotypes.push_back(cur_genotype);

                max_num_gibbs_samples = pit.second;
            }
        }

        assert(static_cast<uint>(round(sum_num_gibbs_samples)) == static_cast<uint>(round(num_gibbs_samples)));

        assert(!(max_genotypes.empty()));
        assert(max_num_gibbs_samples > 0);

        assert(max_genotypes.front().first < num_all_alleles);
        assert(max_genotypes.front().second < num_all_alleles);

        assert(estimated_genotypes->at(sample_idx).empty());

        if (posteriors.at(sample_idx).first == Utils::Ploidy::Diploid) {

            if (max_genotypes.size() == 1) {

                if (allele_filter_flags->at(sample_idx).at(max_genotypes.front().first) == "P") {

                    estimated_genotypes->at(sample_idx).push_back(to_string(max_genotypes.front().first));
                    allele_counts->addAllele(max_genotypes.front().first);
                
                } else {

                    assert(allele_filter_flags->at(sample_idx).at(max_genotypes.front().first) == "F");
                    estimated_genotypes->at(sample_idx).push_back(".");
                }

                if (allele_filter_flags->at(sample_idx).at(max_genotypes.front().second) == "P") {

                    estimated_genotypes->at(sample_idx).push_back(to_string(max_genotypes.front().second));
                    allele_counts->addAllele(max_genotypes.front().second);
                
                } else {

                    assert(allele_filter_flags->at(sample_idx).at(max_genotypes.front().second) == "F");
                    estimated_genotypes->at(sample_idx).push_back(".");
                }

            } else {

                estimated_genotypes->at(sample_idx).push_back(".");
                estimated_genotypes->at(sample_idx).push_back(".");
            }
        
        } else if (posteriors.at(sample_idx).first == Utils::Ploidy::Haploid) {

            if (max_genotypes.size() == 1) {

                if (allele_filter_flags->at(sample_idx).at(max_genotypes.front().first) == "P") {

                    estimated_genotypes->at(sample_idx).push_back(to_string(max_genotypes.front().first));
                    allele_counts->addAllele(max_genotypes.front().first);
                
                } else {

                    assert(allele_filter_flags->at(sample_idx).at(max_genotypes.front().first) == "F");
                    estimated_genotypes->at(sample_idx).push_back(".");
                }

            } else {

                estimated_genotypes->at(sample_idx).push_back(".");             
            }    
        } 

        for (auto &posterior: genotype_posteriors->at(sample_idx)) {

            posterior /= num_gibbs_samples;
        }

        for (auto &posterior: allele_posteriors->at(sample_idx)) {

            posterior /= num_gibbs_samples;
        }
    }
}


vector<float> GenotypeWriter::calculateAlleleCallProbabilities(const vector<vector<float> > & allele_posteriors, const VariantInfo & variant_info) {

    const uint num_all_alleles = variant_info.num_alleles + variant_info.has_dependency;

    vector<float> allele_call_prob(num_all_alleles, 0);

    for (ushort j = 0; j < num_all_alleles; j++) {

        double not_called_prob = 0;
        bool has_perfect = false;

        for (ushort i = 0; i < allele_posteriors.size(); i++) {

            if (!(allele_posteriors.at(i).empty())) {

                assert(allele_posteriors.at(i).size() == num_all_alleles);

                if (!(Utils::floatCompare(allele_posteriors.at(i).at(j), 1))) {

                    assert(allele_posteriors.at(i).at(j) < 1);
                    not_called_prob += log(1 - allele_posteriors.at(i).at(j));
                
                } else {

                    has_perfect = true;
                    break;
                }
            }
        }

        if (has_perfect) {

            allele_call_prob.at(j) = 1;

        } else {

            allele_call_prob.at(j) = 1 - exp(not_called_prob);
        }
    }

    return allele_call_prob;
}


void GenotypeWriter::writeQualityAndFilter(ofstream * output_file, Utils::FilterStatus filter_status, const vector<float> & allele_call_probabilities, const VariantInfo & variant_info) {

    assert(allele_call_probabilities.size() >= variant_info.num_alleles);

    float max_probability = 0;

    for (ushort i = 1; i < variant_info.num_alleles; i++) {

        max_probability = max(max_probability, allele_call_probabilities.at(i));  
    }

    if (Utils::floatCompare(max_probability, 1)) {

        *output_file << "999";
    
    } else if (Utils::floatCompare(max_probability, 0)) {

        *output_file << "0";

    } else {

        assert(max_probability > 0);
        assert(max_probability < 1);

        *output_file << -10*log10(1 - max_probability);
    }

    *output_file << "\t" << filter_status << "\t";   
}


void GenotypeWriter::writeAlleleCounts(ofstream * output_file, const AlleleCounts & allele_counts) {

    vector<float> allele_frequencies;
    allele_frequencies.reserve(allele_counts.counts.size());

    auto cit = allele_counts.counts.begin();
    assert(cit != allele_counts.counts.end());
        
    *output_file << "AC=" << *cit;

    if (allele_counts.num_alleles == 0) {

        allele_frequencies.push_back(0);

    } else {

        allele_frequencies.push_back(*cit/static_cast<float>(allele_counts.num_alleles));
    }
    
    cit++;

    while (cit != allele_counts.counts.end()) {

        *output_file << "," << *cit;
    
        if (allele_counts.num_alleles == 0) {

            allele_frequencies.push_back(0);

        } else {

            allele_frequencies.push_back(*cit/static_cast<float>(allele_counts.num_alleles));
        }

        cit++; 
    }

    auto fit = allele_frequencies.begin();
    assert(fit != allele_frequencies.end());

    *output_file << ";AF=" << *fit;    
    fit++;

    while (fit != allele_frequencies.end()) {

        *output_file << "," << *fit;
        fit++; 
    }    

    *output_file << ";AN=" << allele_counts.num_alleles; 
}


void GenotypeWriter::writeAlleleCallProbabilities(ofstream * output_file, const vector<float> & allele_call_probabilities) {

    auto fit = allele_call_probabilities.begin();
    assert(fit != allele_call_probabilities.end());

    *output_file << ";ACP=" << *fit;
    fit++;

    while (fit != allele_call_probabilities.end()) {

        *output_file << "," << *fit;
        fit++; 
    }
}


void GenotypeWriter::writeSamples(ofstream * output_file, const vector<vector<string> > & estimated_genotypes, const vector<vector<float> > & genotype_posteriors, const vector<vector<float> > & allele_posteriors, const VariantKmerStats & variant_kmer_stats, const VariantInfo & variant_info, vector<vector<string> > * allele_filter_flags) {

    assert(num_samples == estimated_genotypes.size());
    assert(num_samples == genotype_posteriors.size());
    assert(num_samples == allele_posteriors.size());

    assert(variant_kmer_stats.mean_num_kmers.size() == num_samples);
    assert(variant_kmer_stats.mean_num_observed_kmers.size() == num_samples);
    assert(variant_kmer_stats.mean_num_unique_kmers.size() == num_samples);

    for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

        assert(estimated_genotypes.at(sample_idx).empty() == genotype_posteriors.at(sample_idx).empty());
        assert(estimated_genotypes.at(sample_idx).empty() == allele_posteriors.at(sample_idx).empty());
        assert(estimated_genotypes.at(sample_idx).empty() == allele_filter_flags->at(sample_idx).empty());    

        *output_file << "\t";

        if (!(estimated_genotypes.at(sample_idx).empty())) {

            auto egit = estimated_genotypes.at(sample_idx).begin();         
            assert(egit != estimated_genotypes.at(sample_idx).end());         
            
            *output_file << *egit;
            egit++;

            while (egit != estimated_genotypes.at(sample_idx).end()) {

                *output_file << "/" << *egit;
                egit++;
            }

            const uint num_all_alleles = variant_info.num_alleles + variant_info.has_dependency;
            const uint num_potential_genotypes = (num_all_alleles * (num_all_alleles - 1)) / 2 + num_all_alleles;

            assert(genotype_posteriors.at(sample_idx).size() == num_potential_genotypes);    
            assert(allele_posteriors.at(sample_idx).size() == num_all_alleles);    

            writeAlleleFormatField<float>(output_file, genotype_posteriors.at(sample_idx), false);
            writeAlleleFormatField<float>(output_file, allele_posteriors.at(sample_idx), false);

            assert(variant_kmer_stats.mean_num_kmers.at(sample_idx).size() == variant_info.num_alleles);
            assert(variant_kmer_stats.mean_num_observed_kmers.at(sample_idx).size() == variant_info.num_alleles);
            assert(variant_kmer_stats.mean_num_unique_kmers.at(sample_idx).size() == variant_info.num_alleles);

            writeAlleleFormatField<float>(output_file, variant_kmer_stats.mean_num_kmers.at(sample_idx), variant_info.has_dependency);
            writeAlleleFormatField<float>(output_file, variant_kmer_stats.mean_num_observed_kmers.at(sample_idx), variant_info.has_dependency);
            writeAlleleFormatField<float>(output_file, variant_kmer_stats.mean_num_unique_kmers.at(sample_idx), variant_info.has_dependency);

            assert(allele_filter_flags->at(sample_idx).size() == num_all_alleles);    

            for (ushort allele_idx = 0; allele_idx < allele_posteriors.at(sample_idx).size(); allele_idx++) {

                if (Utils::floatCompare(allele_posteriors.at(sample_idx).at(allele_idx), 0)) {

                    assert(allele_filter_flags->at(sample_idx).at(allele_idx) == "P");
                }
            }

            writeAlleleFormatField<string>(output_file, allele_filter_flags->at(sample_idx), false);

        } else {

            *output_file << ":" << empty_variant_sample;        
        }
    } 
}


template <typename ValueType>
void GenotypeWriter::writeAlleleFormatField(ofstream * output_file, const vector<ValueType> & allele_format_field_values, const bool add_missing_allele) {

    assert(!(allele_format_field_values.empty()));

    *output_file << ":" << allele_format_field_values.front();
    
    for (ushort i = 1; i < allele_format_field_values.size(); i++) {

        *output_file << "," << allele_format_field_values.at(i);
    }

    if (add_missing_allele) {

        *output_file << ",-1";
    }
}


void GenotypeWriter::addGenotypes(vector<Genotypes *> * variant_genotypes) {

    genotype_queue->push(variant_genotypes);
}


void GenotypeWriter::completedGenotyping() {

    genotype_queue->pushedLast();

    for (auto &thread: writer_threads) {

        thread.join();
    }

    delete genotype_queue;
}


uint GenotypeWriter::writeGenotypesToVariantFile(const string & vcf_filename, const Regions & chromosome_regions, const vector<Sample> & samples, const string options_header, const string & genome_filename, const uint num_variants) {

    uint num_written_variants = 0;

    if ((vcf_filename.substr(vcf_filename.size()-3,3) == ".gz") or (vcf_filename.substr(vcf_filename.size()-5,5) == ".gzip")) {

        ifstream gz_vcf_file(vcf_filename.c_str(), ifstream::binary);
        assert(gz_vcf_file.is_open());            
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(gz_vcf_file);

        num_written_variants = parseAndWriteGenotypes<boost::iostreams::filtering_istream>(&in, chromosome_regions, samples, options_header, genome_filename, num_variants);

        gz_vcf_file.close();

    } else {

        assert((vcf_filename.substr(vcf_filename.size()-4,4) == ".vcf"));
        ifstream vcf_file(vcf_filename.c_str());
        assert(vcf_file.is_open());

        num_written_variants = parseAndWriteGenotypes<ifstream>(&vcf_file, chromosome_regions, samples, options_header, genome_filename, num_variants);
    
        vcf_file.close();
    }

    return num_written_variants;
}


template <typename FileType>
uint GenotypeWriter::parseAndWriteGenotypes(FileType * vcf_file, const Regions & chromosome_regions, const vector<Sample> & samples, const string & options_header, const string & genome_filename, const uint num_variants) {

    string genotype_file_name = output_prefix;
    genotype_file_name.append("_tmp.txt");

    ifstream genotype_file(genotype_file_name.c_str());
    assert(genotype_file.is_open());

    uint num_written_variants = 0;

    vector<pair<bool, string*> > genotypes = vector<pair<bool, string*> >(num_variants, make_pair(false, nullptr));

    while (genotype_file.good()) {

        string variant_idx = "";
        string has_dependency = "";

        getline(genotype_file, variant_idx, '\t');

        if (variant_idx.size() == 0) {

            genotype_file.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;              
        }

        auto genotype_string = &(genotypes.at(stoi(variant_idx)));

        getline(genotype_file, has_dependency, '\t');

        if (has_dependency == "1") {

            genotype_string->first = true;
        
        } else {

            assert(has_dependency == "0");
            genotype_string->first = false;
        }
            
        genotype_string->second = new string();
        getline(genotype_file, *(genotype_string->second));
    }    

    genotype_file.close();

    string output_file_name = output_prefix;
    output_file_name.append(".vcf");

    ofstream output_file(output_file_name);
    assert(output_file.is_open());

    auto genotypes_it = genotypes.begin();

    string info_header = "";

    bool has_genotypes = true;

    while (vcf_file->good()) {

        if (vcf_file->peek() == '#') {

            string header_line;
            getline(*vcf_file, header_line);

            if ((header_line.substr(0,12) == "##fileformat") or (header_line.substr(0,8) == "##contig")) {

                output_file << header_line << "\n";
            
            } else if (header_line.substr(0,6) == "##INFO") {

                info_header.append(header_line + "\n");
            
            } else if (header_line.substr(0,6) == "#CHROM") {

                vector<string> header_line_split;
                boost::split(header_line_split, header_line, boost::is_any_of("\t"));
        
                assert(header_line_split.size() >= 8);
                output_file << writeHeader(samples, options_header, info_header, genome_filename);

                if (header_line_split.size() == 8) {

                    has_genotypes = false;
                }
            }

            continue;
        }

        string cur_chromosome = "";
        getline(*vcf_file, cur_chromosome, '\t');

        if (cur_chromosome.size() == 0) {

            vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
            continue;              
        }

        string cur_position = "";
        getline(*vcf_file, cur_position, '\t');    

        if (!genotypes_it->second) {

            assert(!genotypes_it->first);

            if (chromosome_regions.isNotIn(cur_chromosome, stoi(cur_position), stoi(cur_position))) {

                vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
                genotypes_it++;

                continue;
            }
        }

        num_written_variants++;

        output_file << cur_chromosome << "\t";
        output_file << cur_position << "\t";    

        string variant_line = "";

        getline(*vcf_file, variant_line, '\t');
        output_file << variant_line << "\t";

        getline(*vcf_file, variant_line, '\t');
        output_file << variant_line << "\t";

        string alternative_alleles = "";
        getline(*vcf_file, alternative_alleles, '\t');
        
        if (alternative_alleles.back() == '*') {

            alternative_alleles.pop_back();
            assert(alternative_alleles.back() == ',');

            alternative_alleles.pop_back();
            assert(alternative_alleles.back() != ',');
        }
        
        assert(alternative_alleles.find("*") == string::npos);

        if (genotypes_it->first) {

            output_file << alternative_alleles << ",*\t";

        } else {

            output_file << alternative_alleles << "\t";
        }

        vcf_file->ignore(numeric_limits<streamsize>::max(), '\t');
        vcf_file->ignore(numeric_limits<streamsize>::max(), '\t');
        
        if (has_genotypes) {

            getline(*vcf_file, variant_line, '\t');
            vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
       
        } else {

            getline(*vcf_file, variant_line);
        }

        if (genotypes_it->second) {

            writeGenotypesVariant(&output_file, genotypes_it->second, variant_line);
            delete genotypes_it->second;

        } else {

            writeUnsupportedVariant(&output_file, variant_line, cur_chromosome, cur_position, count(alternative_alleles.begin(), alternative_alleles.end(), ',') + 1);  
        }

        output_file << "\n";
        genotypes_it++;
    }

    assert(genotypes_it == genotypes.end());

    output_file.close();

    return num_written_variants;
}


string GenotypeWriter::writeHeader(const vector<Sample> & samples, const string & options_header, const string & info_header, const string & genome_filename) {

    stringstream header_stream;

    header_stream << "##reference=file:" << genome_filename << "\n";
    header_stream << options_header;

    header_stream << "##FILTER=<ID=UV,Description=\"Variant filtered due to lack of support for variant type\">" << "\n";
    
    header_stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele counts in called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequencies calculated from called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles in called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=ACP,Number=R,Type=Float,Description=\"Allele call probabilites (used for variant quality calculation)\">" << "\n";

    header_stream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type\">" << "\n";
    header_stream << "##INFO=<ID=VCS,Number=1,Type=Integer,Description=\"Variant cluster size\">" << "\n";
    header_stream << "##INFO=<ID=VCI,Number=1,Type=String,Description=\"Variant cluster id (<chromosome>:<start>-<end>)\">" << "\n";
    header_stream << "##INFO=<ID=VCGS,Number=1,Type=Integer,Description=\"Variant cluster group size (number of variant clusters)\">" << "\n";
    header_stream << "##INFO=<ID=VCGI,Number=1,Type=String,Description=\"Variant cluster group id (<chromosome>:<start>-<end>)\">" << "\n";

    header_stream << "##INFO=<ID=HC,Number=1,Type=Integer,Description=\"Number of haplotype candidates used for inference in variant cluster\">" << "\n";
    header_stream << "##INFO=<ID=HCR,Number=0,Type=Flag,Description=\"Variant cluster has complex region (approximate smallmer counting utilized)\">" << "\n";
    header_stream << "##INFO=<ID=HRS,Number=0,Type=Flag,Description=\"Variant cluster has redundant haplotype sequences\">" << "\n";
    header_stream << "##INFO=<ID=AE,Number=.,Type=String,Description=\"Allele(s) excluded (mutually exclusive with <ANC>). 0: Reference allele\">" << "\n";
    header_stream << "##INFO=<ID=ANC,Number=.,Type=String,Description=\"Allele(s) not covered by a haplotype candidate (mutually exclusive with <AE>). 0: Reference allele\">" << "\n";

    header_stream << info_header;

    header_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n";
    header_stream << "##FORMAT=<ID=GPP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">" << "\n";
    header_stream << "##FORMAT=<ID=MAP,Number=R,Type=Float,Description=\"Marginal allele probabilities\">" << "\n";
    header_stream << "##FORMAT=<ID=NK,Number=R,Type=Float,Description=\"Mean number of informative (unique multiplicity compared to the other alleles in the variant) kmers across gibbs samples. -1: Allele missing \"*\" or not sampled\">" << "\n";
    header_stream << "##FORMAT=<ID=NOK,Number=R,Type=Float,Description=\"Mean number of observed (coverage above zero) informative kmers across gibbs samples. -1: Allele missing \"*\" or not sampled\">" << "\n";
    header_stream << "##FORMAT=<ID=NUK,Number=R,Type=Float,Description=\"Mean number of unique (not intercluster or multicluster) informative kmers across gibbs samples. -1: Allele missing \"*\" or not sampled\">" << "\n";
    header_stream << "##FORMAT=<ID=AFF,Number=R,Type=String,Description=\"Allele filter flag (\"P\": Pass, \"F\": Filtered)\">" << "\n";
        
    header_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for (auto &sit: samples) {

        header_stream << "\t" << sit.name;
    }

    header_stream << "\n";

    return header_stream.str();
}


void GenotypeWriter::writeGenotypesVariant(ofstream * output_file, string * genotypes, const string & info_field_string) {

    auto cur_column = genotypes->find("\t");
    assert(cur_column != string::npos);

    for (ushort i = 0; i < 2; i++) {
        
        cur_column++;
        cur_column = genotypes->find("\t", cur_column);
        assert(cur_column != string::npos);
    }

    stringstream info_and_format;

    assert(!(info_field_string.empty()));

    if (info_field_string != ".") {

        info_and_format << ";" << info_field_string;    
    }

    info_and_format << format_column;

    genotypes->insert(cur_column, info_and_format.str());
    *output_file << *genotypes;
}


void GenotypeWriter::writeUnsupportedVariant(ofstream * output_file, const string & info_field_string, const string & cur_chromosome, const string & cur_position, const uint num_alternative_alleles) {

    *output_file << "0\tUV\t";

    assert(!(cur_chromosome.empty()));
    assert(!(cur_position.empty()));

    stringstream num_allele_zero_stream;
    num_allele_zero_stream << "0";

    for (ushort i = 1; i < num_alternative_alleles; i++) {

        num_allele_zero_stream << ",0";
    }

    *output_file << "AC=" << num_allele_zero_stream.str() << ";AF=" << num_allele_zero_stream.str() << ";AN=0;";
    *output_file << "ACP=0," << num_allele_zero_stream.str() << ";";
    *output_file << "VT=Unsupported;VCS=1;VCI=" << cur_chromosome << ":" << cur_position << "-" << cur_position << ";VCGS=1;VCGI=" << cur_chromosome << ":" << cur_position << "-" << cur_position;
    *output_file << ";HC=0;AE=0";

    for (ushort i = 1; i < (num_alternative_alleles + 1); i++) {

        *output_file << "," << i;
    }

    assert(!(info_field_string.empty()));

    if (info_field_string != ".") {

        *output_file << ";" << info_field_string; 
    }

    *output_file << format_column;

    auto ploidy = chromosome_ploidy.getPloidy(VariantFileParser::classifyGenomeChromosome(VariantFileParser::simplifyChromosomeId(cur_chromosome)));
    assert(ploidy.size() == num_samples);

    for (auto &ploid: ploidy) {

        if (ploid == Utils::Ploidy::Null) {

            *output_file << "\t:";
        
        } else if (ploid == Utils::Ploidy::Haploid) {

            *output_file << "\t.:";
        
        } else {

            *output_file << "\t./.:";            
        }

        *output_file << empty_variant_sample;
    }
}


