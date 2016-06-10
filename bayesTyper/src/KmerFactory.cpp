
/*
KmerFactory.cpp - This file is part of BayesTyper (v0.9)


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


#include <map>
#include <vector>
#include <string>
#include <unordered_map>

#include "KmerFactory.hpp"
#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "VariantFileParser.hpp"
#include "VariantClusterGroup.hpp"
#include "KmerCounter.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "NegativeBinomialDistribution.hpp"
#include "PartitionGraph.hpp"


using namespace std;

static const uint smallmer_set_buckets_per_thread = 100000;

KmerFactory::KmerFactory(const OptionsContainer & options_container) : prng_seed(options_container.getValue<uint>("random-seed")), vcf_file(options_container.getValue<string>("vcf-file")), genome_file(options_container.getValue<string>("genome-file")), decoy_file(options_container.getValue<string>("decoy-file")), num_threads(options_container.getValue<ushort>("threads")), max_allele_length(options_container.getValue<uint>("maximum-allele-length")), min_allele_kmers(options_container.getValue<ushort>("minimum-number-of-kmers-per-allele")), num_haplotype_candidates_per_sample(options_container.getValue<ushort>("number-of-haplotype-candidates-per-sample")) {

	number_of_variants = 0;
	max_alternative_alleles = 0;
}


template<uchar kmer_size>
KmerHash * KmerFactory::buildKmerHash(vector<VariantClusterGroup*> * variant_cluster_groups, const vector<Sample> & samples) {

    assert(variant_cluster_groups->empty());

	VariantFileParser vcf_parser = VariantFileParser(genome_file, max_allele_length, min_allele_kmers, num_threads, prng_seed);
    vcf_parser.addDecoys<kmer_size>(decoy_file);

    vcf_parser.readVariantFile<kmer_size>(vcf_file, variant_cluster_groups);
    vcf_parser.writeSizeStatisticsFiles();

    assert(!(variant_cluster_groups->empty()));
    
    number_of_variants = vcf_parser.getNumberOfVariants();
    max_alternative_alleles = vcf_parser.getMaxAlternativeAlleles();

    cout << "[" << Utils::getLocalTime() << "] Sorting variant cluster groups by decreasing complexity (number of variants) ..." << endl;

    sort(variant_cluster_groups->begin(), variant_cluster_groups->end(), VariantClusterGroupCompare);

    cout << "[" << Utils::getLocalTime() << "] Finished sorting\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Counting smallmers ..." << endl;

    KmerHash * kmer_hash = new KmerHashHybrid<kmer_size>(num_threads, samples.size());
    Utils::SmallmerSet * smallmer_set = new Utils::SmallmerSet(smallmer_set_buckets_per_thread * num_threads);
    
    KmerCounter<kmer_size> kmer_counter(kmer_hash, smallmer_set, samples.size(), num_threads, prng_seed);

    ulong num_unique_small_mers = kmer_counter.countVariantClusterSmallmers(variant_cluster_groups);

    cout << "[" << Utils::getLocalTime() << "] Counted " << num_unique_small_mers << " unique smallmers (" << to_string(Utils::small_kmer_size) << " nt)\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    for (ushort i = 0; i < static_cast<ushort>(Utils::ChromosomeClass::CHROMOSOME_CLASS_SIZE); i++) {

        cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in genomic regions between variant clusters on the " << static_cast<Utils::ChromosomeClass>(i) << " chromosome(s) ..." << endl;

        ulong num_intercluster_kmers = kmer_counter.countInterclusterKmers(vcf_parser.getInterclusterIntervals(static_cast<Utils::ChromosomeClass>(i)), static_cast<Utils::ChromosomeClass>(i));

        cout << "[" << Utils::getLocalTime() << "] Counted " << num_intercluster_kmers << " kmers passing the smallmer filter" << endl;
    }

    cout << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Parsing and filtering KMC2 k-mer tables from " << samples.size() << " sample(s) ...\n" << endl;

    ulong num_sample_kmers = kmer_counter.countSampleKmers(samples); 
  
    cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_sample_kmers << " kmers passing the smallmer filter\n" << endl;
 
    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in variant clusters ..." << endl;

    ulong num_unique_kmers = kmer_counter.countVariantClusterKmers(variant_cluster_groups, num_haplotype_candidates_per_sample);

    cout << "[" << Utils::getLocalTime() << "] Counted " << num_unique_kmers << " unique sample kmers\n" << endl; 
     
    delete smallmer_set;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    genomic_count_distributions = kmer_hash->estimateGenomicCountDistributions(samples);
    assert(genomic_count_distributions.size() == samples.size());

    kmer_hash->markKmers();

    return kmer_hash; 
}


ulong KmerFactory::numberOfVariants() {

	return number_of_variants;
}


ushort KmerFactory::maxAlternativeAlleles() {

	return max_alternative_alleles;
}


const vector<NegativeBinomialDistribution> & KmerFactory::genomicCountDistributions() {

    return genomic_count_distributions;
}


template KmerHash * KmerFactory::buildKmerHash<31>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::buildKmerHash<39>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::buildKmerHash<47>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::buildKmerHash<55>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::buildKmerHash<63>(vector<VariantClusterGroup*> *, const vector<Sample> &);




