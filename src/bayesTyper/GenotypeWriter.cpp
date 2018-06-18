
/*
GenotypeWriter.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <stdio.h>

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


static const string format_column = "GT:GPP:APP:NAK:FAK:MAC:SAF";
static const string empty_variant_sample = ".:.:.:.:.:.";

GenotypeWriter::GenotypeWriter(const string & output_prefix, const ushort num_threads, const vector<Sample> & samples_in, const Chromosomes & chromosomes, const Filters & filters) : samples(samples_in), tmp_filename(output_prefix + "_tmp.txt.gz") {

    assert(!(tmp_outfile.is_open()));
    assert(tmp_outfile_fstream.empty());

    assert(tmp_filename.substr(tmp_filename.size() - 3, 3) == ".gz");

    tmp_outfile.open(tmp_filename, ios::binary);
    assert(tmp_outfile.is_open());

    tmp_outfile_fstream.push(boost::iostreams::gzip_compressor());
    tmp_outfile_fstream.push(boost::ref(tmp_outfile));
    
    assert(tmp_outfile_fstream.is_complete());   

    genotypes_queue = new ProducerConsumerQueue<vector<Genotypes *> *>(Utils::queue_size_thread_scaling * num_threads);
    writing_threads.push_back(thread(&GenotypeWriter::writeGenotypes, this, ref(chromosomes), ref(filters)));
}

void GenotypeWriter::writeGenotypes(const Chromosomes & chromosomes, const Filters & filters) {

    vector<Genotypes *> * variant_genotypes = nullptr;

    while (genotypes_queue->pop(&variant_genotypes)) {

        for (auto & genotypes: *variant_genotypes) {

            auto tmp_outfile_chrom_stats_it = tmp_outfile_chrom_stats.emplace(genotypes->chrom_name, 0);
            tmp_outfile_chrom_stats_it.first->second++;

            assert(genotypes);

            const ushort num_alleles = genotypes->variant_info.numberOfAlleles();

            assert(genotypes->variant_info.numberOfAlleles() > 1);
            assert(genotypes->variant_info.numberOfAlleles() < Utils::ushort_overflow);

            tmp_outfile_fstream << genotypes->chrom_name << "\t" << genotypes->variant_info.position << "\t" << genotypes->variant_info.id;

            auto chromosomes_it = chromosomes.find(genotypes->chrom_name);
            assert(chromosomes_it != chromosomes.cend());

            writeAlleleSequences(genotypes->variant_info, chromosomes_it->second);
            writeQualityAndFilter(genotypes->variant_stats, genotypes->num_homozygote_genotypes, filters);
            writeVariantStats(genotypes->variant_stats, num_alleles);

            tmp_outfile_fstream << ";VCS=" << genotypes->variant_cluster_size << ";VCR=" << genotypes->variant_cluster_region << ";VCGS=" << genotypes->variant_cluster_group_size << ";VCGR=" << genotypes->variant_cluster_group_region << ";HC=" << genotypes->num_candidates;
            
            writeAlleleCover(&(genotypes->non_covered_alleles), num_alleles);
            writeAlleleOrigin(genotypes->variant_info);

            tmp_outfile_fstream << "\t" << format_column;

            writeSamples(genotypes->sample_stats, num_alleles);

            tmp_outfile_fstream << "\n";

            delete genotypes;
        }

        delete variant_genotypes;
    }
}

template <typename ValueType>
void GenotypeWriter::writeAlleleField(const vector<ValueType> & allele_field_values) {

    auto allele_field_values_it = allele_field_values.begin();
    assert(allele_field_values_it != allele_field_values.end());
        
    tmp_outfile_fstream << *allele_field_values_it;
    allele_field_values_it++;

    while (allele_field_values_it != allele_field_values.end()) {

        tmp_outfile_fstream << "," << *allele_field_values_it;
        allele_field_values_it++; 
    }
}

void GenotypeWriter::writeAlleleSequences(const VariantInfo & variant_info, const string & chrom_sequence) {

    const uint max_ref_length = variant_info.maxReferenceLength();
    assert(max_ref_length > 0);

    tmp_outfile_fstream << "\t" << max_ref_length;

    assert(!(variant_info.alt_alleles.empty()));

    tmp_outfile_fstream << "\t";

    for (ushort alt_allele_idx = 0; alt_allele_idx < variant_info.alt_alleles.size(); alt_allele_idx++) {

        if (alt_allele_idx > 0) {

            tmp_outfile_fstream << ",";
        }

        assert(variant_info.alt_alleles.at(alt_allele_idx).ref_length <= max_ref_length);

        tmp_outfile_fstream << variant_info.alt_alleles.at(alt_allele_idx).sequence << chrom_sequence.substr(variant_info.position + variant_info.alt_alleles.at(alt_allele_idx).ref_length - 1, max_ref_length - variant_info.alt_alleles.at(alt_allele_idx).ref_length);
    }

    if (variant_info.has_dependency) {

        tmp_outfile_fstream << ",*";
    }  
}

void GenotypeWriter::writeQualityAndFilter(const Genotypes::VariantStats & variant_stats, const ushort num_homozygote_genotypes, const Filters & filters) {

    if (Utils::floatCompare(variant_stats.max_alt_allele_call_probability, 1)) {

        tmp_outfile_fstream << "\t99";
    
    } else if (Utils::floatCompare(variant_stats.max_alt_allele_call_probability, 0)) {

        tmp_outfile_fstream << "\t0";

    } else {

        assert(variant_stats.max_alt_allele_call_probability > 0);
        assert(variant_stats.max_alt_allele_call_probability < 1);

        tmp_outfile_fstream << "\t" << -10 * log10(1 - variant_stats.max_alt_allele_call_probability);
    }

    if (variant_stats.total_count == 0) {

        tmp_outfile_fstream << "\tAN0";

        if ((samples.size() >= filters.minFilterSamples()) and (num_homozygote_genotypes < filters.minHomozygoteGenotypes())) {

            tmp_outfile_fstream << ";HOM";
        } 

    } else if ((samples.size() >= filters.minFilterSamples()) and (num_homozygote_genotypes < filters.minHomozygoteGenotypes())) {

        tmp_outfile_fstream << "\tHOM";

    } else {

        tmp_outfile_fstream << "\tPASS";
    }
}

void GenotypeWriter::writeVariantStats(const Genotypes::VariantStats & variant_stats, const ushort num_alleles) {

    assert((variant_stats.alt_allele_counts.size() + 1) == num_alleles);
    assert((variant_stats.alt_allele_frequency.size() + 1) == num_alleles);
    assert(variant_stats.allele_call_probabilities.size() == num_alleles);
        
    tmp_outfile_fstream << "\tAC=";
    writeAlleleField<uint>(variant_stats.alt_allele_counts);

    tmp_outfile_fstream <<";AF=";
    writeAlleleField<float>(variant_stats.alt_allele_frequency);

    tmp_outfile_fstream << ";AN=" << variant_stats.total_count;

    tmp_outfile_fstream << ";ACP=";
    writeAlleleField<float>(variant_stats.allele_call_probabilities);
}

void GenotypeWriter::writeAlleleCover(vector<ushort> * non_covered_alleles, const ushort num_alleles) {

    assert(non_covered_alleles->size() <= num_alleles);

    if (!(non_covered_alleles->empty())) {

        sort(non_covered_alleles->begin(), non_covered_alleles->end());        
        tmp_outfile_fstream << ";ANC=";
        writeAlleleField<ushort>(*non_covered_alleles);
    }
}

void GenotypeWriter::writeAlleleOrigin(const VariantInfo & variant_info) {
    
    assert(!(variant_info.alt_alleles.empty()));

    tmp_outfile_fstream << ";ACO=";

    for (ushort alt_allele_idx = 0; alt_allele_idx < variant_info.alt_alleles.size(); alt_allele_idx++) {

        if (alt_allele_idx > 0) {

            tmp_outfile_fstream << ",";
        }

        if (variant_info.alt_alleles.at(alt_allele_idx).aco_att.empty()) {

            tmp_outfile_fstream << ".";            
        
        } else {

            tmp_outfile_fstream << variant_info.alt_alleles.at(alt_allele_idx).aco_att;                        
        }
    }

    if (variant_info.has_dependency) {

        tmp_outfile_fstream << ",.";
    }  
}

void GenotypeWriter::writeSamples(const vector<Genotypes::SampleStats> & sample_stats, const ushort num_alleles) {

    assert(samples.size() == sample_stats.size());

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        tmp_outfile_fstream << "\t";

        if (!(sample_stats.at(sample_idx).genotype_estimate.empty())) {      
            
            assert(sample_stats.at(sample_idx).genotype_estimate.size() <= 2);

            for (ushort i = 0; i < sample_stats.at(sample_idx).genotype_estimate.size(); i++) {

                if (i > 0) {

                    tmp_outfile_fstream << "/";                    
                }

                if (sample_stats.at(sample_idx).genotype_estimate.at(i) != Utils::ushort_overflow) {

                    assert(sample_stats.at(sample_idx).genotype_estimate.at(i) < num_alleles);
                    tmp_outfile_fstream << sample_stats.at(sample_idx).genotype_estimate.at(i);

                } else {

                    tmp_outfile_fstream << ".";
                }
            }

            const uint num_genotypes = (num_alleles * (num_alleles - 1)) / 2 + num_alleles;

            assert((sample_stats.at(sample_idx).genotype_posteriors.size() == num_genotypes) or (sample_stats.at(sample_idx).genotype_posteriors.size() == num_alleles));

            tmp_outfile_fstream << ":";
            writeAlleleField<float>(sample_stats.at(sample_idx).genotype_posteriors);

            assert(sample_stats.at(sample_idx).allele_posteriors.size() == num_alleles);

            tmp_outfile_fstream << ":";
            writeAlleleField<float>(sample_stats.at(sample_idx).allele_posteriors);

            assert(sample_stats.at(sample_idx).allele_kmer_stats.count_stats.size() == num_alleles);
            
            tmp_outfile_fstream << ":";
            writeAlleleKmerStats(sample_stats.at(sample_idx).allele_kmer_stats);

            assert(sample_stats.at(sample_idx).allele_filters.size() == num_alleles);

            tmp_outfile_fstream << ":";
            writeAlleleField<ushort>(sample_stats.at(sample_idx).allele_filters);

        } else {

            assert(sample_stats.at(sample_idx).genotype_posteriors.empty());
            assert(sample_stats.at(sample_idx).allele_posteriors.empty());

            tmp_outfile_fstream << ":" << empty_variant_sample;        
        }
    } 
}

void GenotypeWriter::writeAlleleKmerStats(const AlleleKmerStats & allele_kmer_stats) {

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

    tmp_outfile_fstream << allele_kmer_counts.str() << ":" << allele_kmer_fractions.str() << ":" << allele_kmer_means.str();
}

void GenotypeWriter::addGenotypes(vector<Genotypes *> * variant_genotypes) {

    genotypes_queue->push(variant_genotypes);
}

void GenotypeWriter::finalise(const string & output_prefix, const Chromosomes & chromosomes, const string & graph_options_header, const OptionsContainer & options_container, const Filters & filters) {

    genotypes_queue->pushedLast();

    for (auto & writing_thread: writing_threads) {

        writing_thread.join();
    }

    delete genotypes_queue;

    tmp_outfile_fstream.reset();
    assert(!(tmp_outfile.is_open()));


    cout << "[" << Utils::getLocalTime() << "] Sorting genotyped variants ..." << endl;

    assert(tmp_filename.substr(tmp_filename.size() - 3, 3) == ".gz");

    ifstream tmp_infile(tmp_filename, ios::binary);
    assert(tmp_infile.is_open());

    boost::iostreams::filtering_istream tmp_infile_fstream;
    
    tmp_infile_fstream.push(boost::iostreams::gzip_decompressor());
    tmp_infile_fstream.push(boost::ref(tmp_infile));

    assert(tmp_infile_fstream.is_complete());    

    unordered_map<string, vector<GenotypedVariant> > genotyped_variants;

    for (auto & chrom_stat: tmp_outfile_chrom_stats) {

        auto genotyped_variants_it = genotyped_variants.emplace(chrom_stat.first, vector<GenotypedVariant>());
        assert(genotyped_variants_it.second);

        genotyped_variants_it.first->second.reserve(chrom_stat.second);
    }

    string chrom_name = "";

    while (getline(tmp_infile_fstream, chrom_name, '\t')) {

        auto genotyped_variants_it = genotyped_variants.find(chrom_name);
        assert(genotyped_variants_it != genotyped_variants.end());

        genotyped_variants_it->second.emplace_back(GenotypedVariant());

        string position = "";

        assert(getline(tmp_infile_fstream, position, '\t'));
        genotyped_variants_it->second.back().position = stoi(position);

        assert(getline(tmp_infile_fstream, genotyped_variants_it->second.back().variant_id, '\t'));

        string max_ref_length = "";

        assert(getline(tmp_infile_fstream, max_ref_length, '\t'));
        genotyped_variants_it->second.back().max_ref_length = stoi(max_ref_length);

        assert(getline(tmp_infile_fstream, genotyped_variants_it->second.back().genotypes, '\n'));
    }

    tmp_infile_fstream.reset();
    assert(!(tmp_infile.is_open()));

    assert(remove(tmp_filename.c_str()) == 0);


    ofstream variants_outfile;
    boost::iostreams::filtering_ostream variants_outfile_fstream;

    assert(!(variants_outfile.is_open()));
    assert(variants_outfile_fstream.empty());

    string genotype_filename = "";
    
    if (options_container.getValue<bool>("gzip-output")) {

        genotype_filename = output_prefix + ".vcf.gz";

        variants_outfile_fstream.push(boost::iostreams::gzip_compressor());        
        variants_outfile.open(genotype_filename, ios::binary);
    
    } else {

        genotype_filename = output_prefix + ".vcf";

        variants_outfile.open(genotype_filename);
    }

    assert(variants_outfile.is_open());

    variants_outfile_fstream.push(boost::ref(variants_outfile));
    assert(variants_outfile_fstream.is_complete());    

    variants_outfile_fstream << generateHeader(options_container.getValue<string>("genome-file"), chromosomes, graph_options_header, options_container.getHeader(), filters);

    ulong num_genotyped_variants = 0;

    auto chromosomes_it = chromosomes.cbegin();

    while (chromosomes_it != chromosomes.cend()) {

        auto genotyped_variants_it = genotyped_variants.find(chromosomes_it->first);
        
        if (genotyped_variants_it != genotyped_variants.end()) {

            sort(genotyped_variants_it->second.begin(), genotyped_variants_it->second.end(), GenotypeWriter::genotypedVariantCompare);
            assert(genotyped_variants_it->second.size() == tmp_outfile_chrom_stats.at(genotyped_variants_it->first));

            for (auto & genotyped_variant: genotyped_variants_it->second) {

                assert(genotyped_variant.position > 0);
                assert(genotyped_variant.max_ref_length > 0);

                assert(!(genotyped_variant.variant_id.empty()));
                assert(!(genotyped_variant.genotypes.empty()));

                variants_outfile_fstream << genotyped_variants_it->first << "\t" << genotyped_variant.position << "\t" << genotyped_variant.variant_id << "\t" << chromosomes_it->second.substr(genotyped_variant.position - 1, genotyped_variant.max_ref_length) << "\t" << genotyped_variant.genotypes << endl;
            }

            num_genotyped_variants += genotyped_variants_it->second.size();
        }        

        chromosomes_it++;
    }

    variants_outfile_fstream.reset();
    assert(!(variants_outfile.is_open()));

    cout << "[" << Utils::getLocalTime() << "] Wrote " << num_genotyped_variants << " genotyped variants to " << genotype_filename << endl;       
}

string GenotypeWriter::generateHeader(const string & genome_filename, const Chromosomes & chromosomes, const string & graph_options_header, const string & genotype_options_header, const Filters & filters) {

    stringstream header_ss;

    header_ss << "##fileformat=VCFv4.2\n";
    header_ss << "##reference=file:" << genome_filename << "\n";

    auto chromosomes_it = chromosomes.cbegin();

    while (chromosomes_it != chromosomes.cend()) {

        if (!(chromosomes.isDecoy(chromosomes_it->first))) {

            header_ss << "##contig=<ID=" << chromosomes_it->first << ",length=" << chromosomes_it->second.size() << ">\n";
        }
        
        chromosomes_it++;                     
    }

    header_ss << graph_options_header;
    header_ss << genotype_options_header;

    if (samples.size() >= filters.minFilterSamples()) {

        header_ss << "##FILTER=<ID=HOM,Description=\"Less than " + to_string(filters.minHomozygoteGenotypes()) + " homozygote genotypes (calculated before other filters)\">\n";
    }

    header_ss << "##FILTER=<ID=AN0,Description=\"No called genotypes (AN = 0)\">\n";
    
    header_ss << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternative allele counts in called genotypes\">\n";
    header_ss << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternative allele frequencies in called genotypes\">\n";
    header_ss << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
    header_ss << "##INFO=<ID=ACP,Number=R,Type=Float,Description=\"Allele call probabilites (maximum APP across samples)\">\n";

    header_ss << "##INFO=<ID=VCS,Number=1,Type=Integer,Description=\"Variant cluster size\">\n";
    header_ss << "##INFO=<ID=VCR,Number=1,Type=String,Description=\"Variant cluster region (<chromosome>:<start>-<end>)\">\n";
    header_ss << "##INFO=<ID=VCGS,Number=1,Type=Integer,Description=\"Variant cluster group size (number of variant clusters)\">\n";
    header_ss << "##INFO=<ID=VCGR,Number=1,Type=String,Description=\"Variant cluster group region (<chromosome>:<start>-<end>)\">\n";

    header_ss << "##INFO=<ID=HC,Number=1,Type=Integer,Description=\"Number of haplotype candidates used for inference in variant cluster\">\n";
    header_ss << "##INFO=<ID=ANC,Number=.,Type=String,Description=\"Allele(s) not covered by a haplotype candidate ('0': Reference allele)\">\n";
    header_ss << "##INFO=<ID=ACO,Number=A,Type=String,Description=\"Alternative allele call-set origin(s) (<call-set>:...)\">\n";

    header_ss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    header_ss << "##FORMAT=<ID=GPP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">\n";
    header_ss << "##FORMAT=<ID=APP,Number=R,Type=Float,Description=\"Allele posterior probabilities\">\n";
    header_ss << "##FORMAT=<ID=NAK,Number=R,Type=Float,Description=\"Mean number of allele kmers across gibbs samples ('-1': Not sampled)\">\n";
    header_ss << "##FORMAT=<ID=FAK,Number=R,Type=Float,Description=\"Mean fraction of observed allele kmers across gibbs samples ('-1': Not sampled or NAK = 0)\">\n";
    header_ss << "##FORMAT=<ID=MAC,Number=R,Type=Float,Description=\"Mean allele kmer coverage (mean value) across gibbs samples ('-1': Not sampled or NAK = 0)\">\n";
    header_ss << "##FORMAT=<ID=SAF,Number=R,Type=Integer,Description=\"Sample specific allele filter ('0': PASS, '1': NAK, '2': FAK, '3': NAK and FAK)\">\n";
        
    header_ss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for (auto & sample: samples) {

        header_ss << "\t" << sample.name;
    }

    header_ss << "\n";

    return header_ss.str();
}

bool GenotypeWriter::genotypedVariantCompare(const GenotypedVariant & first_genotyped_variant, const GenotypedVariant & second_genotyped_variant) { 

    return (first_genotyped_variant.position < second_genotyped_variant.position);
}
