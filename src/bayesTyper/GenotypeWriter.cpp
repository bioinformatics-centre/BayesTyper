
/*
GenotypeWriter.cpp - This file is part of BayesTyper (v1.1)


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
#include "KmerStats.hpp"
#include "VariantInfo.hpp"
#include "VariantFileParser.hpp"
#include "OptionsContainer.hpp"
#include "ChromosomePloidy.hpp"

using namespace std;

static const string format_column = "GT:GPP:APP:NAK:FAK:MAC";
static const string empty_variant_sample = ".:.:.:.:.";

static const pair<uint, uint> sample_fak_estimates_scaling_factors(1000 + 1, 1000);
static const pair<uint, uint> sample_mac_estimates_scaling_factors(Utils::uchar_overflow * 10 + 1, 10);

GenotypeWriter::GenotypeWriter(const vector<Sample> & samples, const string & output_prefix_in, const ushort num_threads) : samples(samples), output_prefix(output_prefix_in) {

    genotypes_queue = new ProducerConsumerQueue<vector<Genotypes *> *>(Utils::queue_size_thread_scaling * num_threads);
    writer_threads.push_back(thread(&GenotypeWriter::writeTemporaryFile, this));

    sample_fak_estimates = vector<vector<ulong> >(samples.size(), vector<ulong>(sample_fak_estimates_scaling_factors.first, 0));
    sample_mac_estimates = vector<vector<ulong> >(samples.size(), vector<ulong>(sample_mac_estimates_scaling_factors.first, 0));
}

void GenotypeWriter::writeTemporaryFile() {

    ofstream output_file(output_prefix + "_genotypes_tmp.txt");
    assert(output_file.is_open());

    vector<Genotypes *> * variant_genotypes = nullptr;

    while (genotypes_queue->pop(&variant_genotypes)) {

        for (auto & genotypes: *variant_genotypes) {

            assert(genotypes);
            assert(genotypes->variant_info.num_alleles > 1);
            assert(genotypes->variant_info.num_alleles < Utils::ushort_overflow);

            output_file << genotypes->variant_info.variant_id << "\t" << genotypes->variant_info.has_dependency;

            writeQualityAndFilter(&output_file, genotypes->variant_stats, genotypes->variant_info);
            writeVariantStats(&output_file, genotypes->variant_stats, genotypes->variant_info);

            output_file << ";VCS=" << genotypes->variant_cluster_size << ";VCI=" << genotypes->variant_cluster_id << ";VCGS=" << genotypes->variant_cluster_group_size << ";VCGI=" << genotypes->variant_cluster_group_id << ";HC=" << genotypes->num_candidates;

            if (genotypes->has_redundant_sequence) {

                output_file << ";HRS";
            }
            
            writeAlleleCover(&output_file, &(genotypes->non_covered_alleles), &(genotypes->variant_info));

            output_file << "\t" << format_column;

            writeSamples(&output_file, genotypes->sample_stats, genotypes->variant_info);
            addAlleleKmerEstimates(genotypes->sample_stats);

            output_file << "\n";

            delete genotypes;
        }

        delete variant_genotypes;
    }

    output_file.close();
}

template <typename ValueType>
void GenotypeWriter::writeAlleleField(ofstream * output_file, const vector<ValueType> & allele_field_values) {

    auto allele_field_values_it = allele_field_values.begin();
    assert(allele_field_values_it != allele_field_values.end());
        
    *output_file << *allele_field_values_it;
    allele_field_values_it++;

    while (allele_field_values_it != allele_field_values.end()) {

        *output_file << "," << *allele_field_values_it;
        allele_field_values_it++; 
    }
}

void GenotypeWriter::writeQualityAndFilter(ofstream * output_file, const Genotypes::VariantStats & variant_stats, const VariantInfo & variant_info) {

    assert(variant_stats.allele_call_probabilities.size() == variant_info.num_alleles);

    float max_probability = 0;

    for (ushort i = 1; i < (variant_info.num_alleles - variant_info.has_dependency); i++) {

        max_probability = max(max_probability, variant_stats.allele_call_probabilities.at(i));  
    }

    if (Utils::floatCompare(max_probability, 1)) {

        *output_file << "\t99";
    
    } else if (Utils::floatCompare(max_probability, 0)) {

        *output_file << "\t0";

    } else {

        assert(max_probability > 0);
        assert(max_probability < 1);

        *output_file << "\t" << -10 * log10(1 - max_probability);
    }

    *output_file << "\t.";   
}

void GenotypeWriter::writeVariantStats(ofstream * output_file, const Genotypes::VariantStats & variant_stats, const VariantInfo & variant_info) {

    assert((variant_stats.alt_allele_counts.size() + 1) == variant_info.num_alleles);
    assert((variant_stats.alt_allele_frequency.size() + 1) == variant_info.num_alleles);
    assert(variant_stats.allele_call_probabilities.size() == variant_info.num_alleles);
        
    *output_file << "\tAC=";
    writeAlleleField<uint>(output_file, variant_stats.alt_allele_counts);

    *output_file <<";AF=";
    writeAlleleField<float>(output_file, variant_stats.alt_allele_frequency);

    *output_file << ";AN=" << variant_stats.total_count;

    *output_file << ";ACP=";
    writeAlleleField<float>(output_file, variant_stats.allele_call_probabilities);
}

void GenotypeWriter::writeAlleleCover(ofstream * output_file, vector<ushort> * non_covered_alleles, VariantInfo * variant_info) {

    assert(variant_info->num_alleles > (variant_info->excluded_alt_alleles.size() + variant_info->has_dependency));

    if (!(variant_info->excluded_alt_alleles.empty())) {
        
        sort(variant_info->excluded_alt_alleles.begin(), variant_info->excluded_alt_alleles.end());
        
        *output_file << ";AE=";
        writeAlleleField<ushort>(output_file, variant_info->excluded_alt_alleles);
    }

    assert(non_covered_alleles->size() <= variant_info->num_alleles);

    if (!(non_covered_alleles->empty())) {

        sort(non_covered_alleles->begin(), non_covered_alleles->end());        
        *output_file << ";ANC=";
        writeAlleleField<ushort>(output_file, *non_covered_alleles);
    }
}

void GenotypeWriter::writeSamples(ofstream * output_file, const vector<Genotypes::SampleStats> & sample_stats, const VariantInfo & variant_info) {

    assert(samples.size() == sample_stats.size());

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        *output_file << "\t";

        if (!(sample_stats.at(sample_idx).genotype_estimate.empty())) {      
            
            assert(sample_stats.at(sample_idx).genotype_estimate.size() <= 2);

            for (ushort i = 0; i < sample_stats.at(sample_idx).genotype_estimate.size(); i++) {

                if (i > 0) {

                    *output_file << "/";                    
                }

                if (sample_stats.at(sample_idx).genotype_estimate.at(i) != Utils::ushort_overflow) {

                    assert(sample_stats.at(sample_idx).genotype_estimate.at(i) < variant_info.num_alleles);
                    *output_file << sample_stats.at(sample_idx).genotype_estimate.at(i);

                } else {

                    *output_file << ".";
                }
            }

            const uint num_genotypes = (variant_info.num_alleles * (variant_info.num_alleles - 1)) / 2 + variant_info.num_alleles;

            assert((sample_stats.at(sample_idx).genotype_posteriors.size() == num_genotypes) or (sample_stats.at(sample_idx).genotype_posteriors.size() == variant_info.num_alleles));

            *output_file << ":";
            writeAlleleField<float>(output_file, sample_stats.at(sample_idx).genotype_posteriors);

            assert(sample_stats.at(sample_idx).allele_posteriors.size() == variant_info.num_alleles);

            *output_file << ":";
            writeAlleleField<float>(output_file, sample_stats.at(sample_idx).allele_posteriors);

            assert(sample_stats.at(sample_idx).allele_kmer_stats.count_stats.size() == variant_info.num_alleles);
            
            *output_file << ":";
            writeAlleleKmerStats(output_file, sample_stats.at(sample_idx).allele_kmer_stats);

        } else {

            assert(sample_stats.at(sample_idx).genotype_posteriors.empty());
            assert(sample_stats.at(sample_idx).allele_posteriors.empty());

            *output_file << ":" << empty_variant_sample;        
        }
    } 
}

void GenotypeWriter::writeAlleleKmerStats(ofstream * output_file, const AlleleKmerStats & allele_kmer_stats) {

    assert(allele_kmer_stats.fraction_stats.size() == allele_kmer_stats.count_stats.size());
    assert(allele_kmer_stats.mean_stats.size() == allele_kmer_stats.count_stats.size());

    stringstream allele_kmer_counts;
    stringstream allele_kmer_fractions;
    stringstream allele_kmer_means;

    allele_kmer_counts << allele_kmer_stats.count_stats.front().getMean().first;
    allele_kmer_fractions << allele_kmer_stats.fraction_stats.front().getMean().first;
    allele_kmer_means << allele_kmer_stats.mean_stats.front().getMean().first;

    for (ushort allele_idx = 1; allele_idx < allele_kmer_stats.count_stats.size(); allele_idx++) {

        allele_kmer_counts << "," << allele_kmer_stats.count_stats.at(allele_idx).getMean().first;
        allele_kmer_fractions << "," << allele_kmer_stats.fraction_stats.at(allele_idx).getMean().first;
        allele_kmer_means << "," << allele_kmer_stats.mean_stats.at(allele_idx).getMean().first;
    }

    *output_file << allele_kmer_counts.str() << ":" << allele_kmer_fractions.str() << ":" << allele_kmer_means.str();
}

void GenotypeWriter::addAlleleKmerEstimates(const vector<Genotypes::SampleStats> & sample_stats) {

    assert(samples.size() == sample_stats.size());

    assert(samples.size() == sample_fak_estimates.size());
    assert(samples.size() == sample_mac_estimates.size());

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        if (sample_stats.at(sample_idx).is_allele_kmer_estimate_variant) {

            for (ushort i = 0; i < sample_stats.at(sample_idx).genotype_estimate.size(); i++) {

                auto allele_idx = sample_stats.at(sample_idx).genotype_estimate.at(i);
                assert(allele_idx != Utils::ushort_overflow);

                assert(sample_stats.at(sample_idx).allele_posteriors.at(allele_idx) >= Utils::min_posterior_allele_kmer_estimate);
                assert(sample_stats.at(sample_idx).allele_kmer_stats.count_stats.at(allele_idx).getMean().first > 0);

                sample_fak_estimates.at(sample_idx).at(round(sample_stats.at(sample_idx).allele_kmer_stats.fraction_stats.at(allele_idx).getMean().first * sample_fak_estimates_scaling_factors.second))++;
                sample_mac_estimates.at(sample_idx).at(round(sample_stats.at(sample_idx).allele_kmer_stats.mean_stats.at(allele_idx).getMean().first * sample_mac_estimates_scaling_factors.second))++;
            }
        }
    } 
}

void GenotypeWriter::addGenotypes(vector<Genotypes *> * variant_genotypes) {

    genotypes_queue->push(variant_genotypes);
}

void GenotypeWriter::completedGenotyping() {

    genotypes_queue->pushedLast();

    for (auto &thread: writer_threads) {

        thread.join();
    }

    delete genotypes_queue;
}

void GenotypeWriter::writeSampleAlleleKmerFractionCumDistFunc() {

    writeSampleAttributeCumDistFunc("FAK", sample_fak_estimates, sample_fak_estimates_scaling_factors);
}

void GenotypeWriter::writeSampleAlleleKmerCoverageCumDistFunc() {

    writeSampleAttributeCumDistFunc("MAC", sample_mac_estimates, sample_mac_estimates_scaling_factors);
}

void GenotypeWriter::writeSampleAttributeCumDistFunc(const string & attribute, const vector<vector<ulong> > & allele_kmer_estimates, const pair<uint, uint> & allele_kmer_estimates_scaling_factors) {

    ofstream cdf_file(output_prefix + "_sample_attribute_" + attribute + "_cdf.txt");
    assert(cdf_file.is_open());

    cdf_file << attribute;

    assert(samples.size() == allele_kmer_estimates.size());
    
    vector<double> allele_kmer_estimates_total_sum;
    allele_kmer_estimates_total_sum.reserve(samples.size());

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        assert(allele_kmer_estimates.at(sample_idx).size() == allele_kmer_estimates_scaling_factors.first);
        allele_kmer_estimates_total_sum.emplace_back(accumulate(allele_kmer_estimates.at(sample_idx).begin(), allele_kmer_estimates.at(sample_idx).end(), 0));

        cdf_file << "\t" << samples.at(sample_idx).name;
    } 

    cdf_file << endl;

    assert(allele_kmer_estimates_total_sum.size() == allele_kmer_estimates.size());
    vector<ulong> allele_kmer_estimates_cumulative_sum(allele_kmer_estimates.size(), 0); 

    for (uint value_idx = 0; value_idx < allele_kmer_estimates_scaling_factors.first; value_idx++) {

        bool write_cdf_line = false;

        for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

            allele_kmer_estimates_cumulative_sum.at(sample_idx) += allele_kmer_estimates.at(sample_idx).at(value_idx);

            if (allele_kmer_estimates.at(sample_idx).at(value_idx) > 0) {

                write_cdf_line = true;
            }
        }

        if (write_cdf_line) {

            cdf_file << static_cast<float>(value_idx) / allele_kmer_estimates_scaling_factors.second;

            for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

                cdf_file << "\t" << allele_kmer_estimates_cumulative_sum.at(sample_idx) / allele_kmer_estimates_total_sum.at(sample_idx);
            }

            cdf_file << endl;
        }
    }

    cdf_file.close();
}   

uint GenotypeWriter::writeGenotypesToVariantCallFormat(const string & vcf_filename, const Regions & chromosome_regions, const string & options_header, const string & genome_filename, const uint num_variants) {

    cout << "\n[" << Utils::getLocalTime() << "] Writing genotypes to " << output_prefix << ".vcf ..." << endl;

    uint num_written_variants = 0;

    if ((vcf_filename.substr(vcf_filename.size()-3,3) == ".gz") or (vcf_filename.substr(vcf_filename.size()-5,5) == ".gzip")) {

        ifstream gz_vcf_file(vcf_filename.c_str(), ifstream::binary);
        assert(gz_vcf_file.is_open());            
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(gz_vcf_file);

        num_written_variants = parseAndWriteGenotypes<boost::iostreams::filtering_istream>(&in, chromosome_regions, options_header, genome_filename, num_variants);

        gz_vcf_file.close();

    } else {

        assert((vcf_filename.substr(vcf_filename.size()-4,4) == ".vcf"));
        ifstream vcf_file(vcf_filename.c_str());
        assert(vcf_file.is_open());

        num_written_variants = parseAndWriteGenotypes<ifstream>(&vcf_file, chromosome_regions, options_header, genome_filename, num_variants);
    
        vcf_file.close();
    }

    return num_written_variants;
}

template <typename FileType>
uint GenotypeWriter::parseAndWriteGenotypes(FileType * vcf_file, const Regions & chromosome_regions, const string & options_header, const string & genome_filename, const uint num_variants) {

    string genotypes_tmp_filename = output_prefix + "_genotypes_tmp.txt";

    ifstream genotypes_file(genotypes_tmp_filename);
    assert(genotypes_file.is_open());

    uint num_written_variants = 0;

    vector<pair<bool, string*> > genotypes = vector<pair<bool, string*> >(num_variants, make_pair(false, nullptr));

    while (genotypes_file.good()) {

        string variant_idx = "";
        string has_dependency = "";

        getline(genotypes_file, variant_idx, '\t');

        if (variant_idx.size() == 0) {

            genotypes_file.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;              
        }

        auto genotype_string = &(genotypes.at(stoi(variant_idx)));

        getline(genotypes_file, has_dependency, '\t');

        if (has_dependency == "1") {

            genotype_string->first = true;
        
        } else {

            assert(has_dependency == "0");
            genotype_string->first = false;
        }
            
        genotype_string->second = new string();
        getline(genotypes_file, *(genotype_string->second));
    }    

    genotypes_file.close();
    remove(genotypes_tmp_filename.c_str());

    ChromosomePloidy chromosome_ploidy = ChromosomePloidy(samples);

    ofstream output_file(output_prefix + ".vcf");
    assert(output_file.is_open());

    auto genotypes_it = genotypes.begin();

    while (vcf_file->good()) {

        if (vcf_file->peek() == '#') {

            string header_line;
            getline(*vcf_file, header_line);

            if ((header_line.substr(0,12) == "##fileformat") or (header_line.substr(0,8) == "##contig")) {

                output_file << header_line << "\n";
            
            } else if (header_line.substr(0,6) == "#CHROM") {

                vector<string> header_line_split;
                boost::split(header_line_split, header_line, boost::is_any_of("\t"));
        
                assert(header_line_split.size() >= 8);
                output_file << writeHeader(samples, options_header, genome_filename);
            }

            continue;
        }

        string chromosome = "";
        getline(*vcf_file, chromosome, '\t');

        if (chromosome.size() == 0) {

            vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
            continue;              
        }

        string position = "";
        getline(*vcf_file, position, '\t');    

        if (!genotypes_it->second) {

            assert(!genotypes_it->first);

            if (!(chromosome_regions.overlaps(chromosome, stoi(position), stoi(position)))) {

                vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
                genotypes_it++;

                continue;
            }
        }

        num_written_variants++;

        string id = "";
        getline(*vcf_file, id, '\t');

        string ref_allele = "";
        getline(*vcf_file, ref_allele, '\t');

        string alt_allele = "";
        getline(*vcf_file, alt_allele, '\t');

        vcf_file->ignore(numeric_limits<streamsize>::max(), '\n');
        
        if (alt_allele.back() == '*') {

            alt_allele.pop_back();
            assert(alt_allele.back() == ',');

            alt_allele.pop_back();
            assert(alt_allele.back() != ',');
        }
        
        assert(alt_allele.find("*") == string::npos);

        if (genotypes_it->first) {

            alt_allele.append(",*");
        } 

        output_file << chromosome << "\t" << position << "\t" << id << "\t" << ref_allele << "\t" << alt_allele;

        if (genotypes_it->second) {

            output_file << "\t" << *genotypes_it->second;
            delete genotypes_it->second;

        } else {

            writeUnsupportedVariant(&output_file, chromosome, position, count(alt_allele.begin(), alt_allele.end(), ',') + 2, chromosome_ploidy);  
        }

        output_file << "\n";
        genotypes_it++;
    }

    assert(genotypes_it == genotypes.end());

    output_file.close();

    return num_written_variants;
}

string GenotypeWriter::writeHeader(const vector<Sample> & samples, const string & options_header, const string & genome_filename) {

    stringstream header_stream;

    header_stream << "##reference=file:" << genome_filename << "\n";
    header_stream << options_header;

    header_stream << "##FILTER=<ID=UV,Description=\"Unsupported variant type\">" << "\n";
    
    header_stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele counts in called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequencies calculated from called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles in called genotypes\">" << "\n";
    header_stream << "##INFO=<ID=ACP,Number=R,Type=Float,Description=\"Allele call probabilites (used for variant quality calculation)\">" << "\n";

    header_stream << "##INFO=<ID=VCS,Number=1,Type=Integer,Description=\"Variant cluster size\">" << "\n";
    header_stream << "##INFO=<ID=VCI,Number=1,Type=String,Description=\"Variant cluster id (<chromosome>:<start>-<end>)\">" << "\n";
    header_stream << "##INFO=<ID=VCGS,Number=1,Type=Integer,Description=\"Variant cluster group size (number of variant clusters)\">" << "\n";
    header_stream << "##INFO=<ID=VCGI,Number=1,Type=String,Description=\"Variant cluster group id (<chromosome>:<start>-<end>)\">" << "\n";

    header_stream << "##INFO=<ID=HC,Number=1,Type=Integer,Description=\"Number of haplotype candidates used for inference in variant cluster\">" << "\n";
    header_stream << "##INFO=<ID=HRS,Number=0,Type=Flag,Description=\"Variant cluster has redundant haplotype sequences\">" << "\n";
    header_stream << "##INFO=<ID=AE,Number=.,Type=String,Description=\"Allele(s) excluded (mutually exclusive with <ANC>). 0: Reference allele\">" << "\n";
    header_stream << "##INFO=<ID=ANC,Number=.,Type=String,Description=\"Allele(s) not covered by a haplotype candidate (mutually exclusive with <AE>). 0: Reference allele\">" << "\n";

    header_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n";
    header_stream << "##FORMAT=<ID=GPP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">" << "\n";
    header_stream << "##FORMAT=<ID=APP,Number=R,Type=Float,Description=\"Allele posterior probabilities\">" << "\n";
    header_stream << "##FORMAT=<ID=NAK,Number=R,Type=Float,Description=\"Mean number of allele kmers across gibbs samples. -1: Not sampled\">" << "\n";
    header_stream << "##FORMAT=<ID=FAK,Number=R,Type=Float,Description=\"Mean fraction of observed allele kmers across gibbs samples. -1: Not sampled\">" << "\n";
    header_stream << "##FORMAT=<ID=MAC,Number=R,Type=Float,Description=\"Mean mean allele kmer coverage across gibbs samples. -1: Not sampled\">" << "\n";
        
    header_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for (auto &sit: samples) {

        header_stream << "\t" << sit.name;
    }

    header_stream << "\n";

    return header_stream.str();
}

void GenotypeWriter::writeUnsupportedVariant(ofstream * output_file, const string & chromosome, const string & position, const uint num_alleles, const ChromosomePloidy & chromosome_ploidy) {

    assert(num_alleles >= 2);

    stringstream alt_allele_zero_ss;
    alt_allele_zero_ss << "0";

    for (ushort i = 0; i < (num_alleles - 2); i++) {

        alt_allele_zero_ss << ",0";
    }

    *output_file << "\t0\tUV\tAC=" << alt_allele_zero_ss.str() << ";AF=" << alt_allele_zero_ss.str() << ";AN=0;" << "ACP=0," << alt_allele_zero_ss.str() << ";VCS=1;VCI=" << chromosome << ":" << position << "-" << position << ";VCGS=1;VCGI=" << chromosome << ":" << position << "-" << position << ";HC=0;AE=0";

    for (ushort allele_idx = 1; allele_idx < num_alleles; allele_idx++) {

        *output_file << "," << allele_idx;
    }

    *output_file << "\t" << format_column;

    auto ploidy = chromosome_ploidy.getPloidy(VariantFileParser::classifyGenomeChromosome(VariantFileParser::simplifyChromosomeId(chromosome)));
    assert(ploidy.size() == samples.size());

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


