
/*
KmerFactory.cpp - This file is part of BayesTyper (v1.1)


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
#include <numeric>

#include "KmerFactory.hpp"
#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "VariantFileParser.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantClusterGroup.hpp"
#include "KmerCounter.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "NegativeBinomialDistribution.hpp"


using namespace std;


KmerFactory::KmerFactory(const OptionsContainer & options_container) : prng_seed(options_container.getValue<uint>("random-seed")), vcf_file(options_container.getValue<string>("vcf-file")), genome_file(options_container.getValue<string>("genome-file")), output_prefix(options_container.getValue<string>("output-prefix")), decoy_file(options_container.getValue<string>("decoy-file")), num_threads(options_container.getValue<ushort>("threads")), max_allele_length(options_container.getValue<uint>("max-allele-length")), copy_number_variant_threshold(options_container.getValue<double>("copy-number-variant-threshold")), max_sample_haplotype_candidates(options_container.getValue<ushort>("max-number-of-sample-haplotype-candidates")), num_genomic_rate_gc_bias_bins(options_container.getValue<uchar>("number-of-genomic-rate-gc-bias-bins")) {

	number_of_variants = 0;
    min_sample_kmer_count = 0;
	max_alternative_alleles = 0;
}


template <uchar kmer_size>
KmerHash * KmerFactory::initKmerHash(vector<VariantClusterGroup*> * variant_cluster_groups, const vector<Sample> & samples) {

    assert(variant_cluster_groups->empty());
    assert(!(samples.empty()));

    vector<VariantClusterGraph*> variant_cluster_graphs;

	VariantFileParser vcf_parser = VariantFileParser(genome_file, max_allele_length, copy_number_variant_threshold, num_threads, prng_seed);
    
    vcf_parser.addDecoys<kmer_size>(decoy_file);
    vcf_parser.readVariantFile<kmer_size>(vcf_file, &variant_cluster_graphs, variant_cluster_groups);

    assert(!(variant_cluster_graphs.empty()));
    assert(!(variant_cluster_groups->empty()));
    
    number_of_variants = vcf_parser.getNumberOfVariants();
    max_alternative_alleles = vcf_parser.getMaxAlternativeAlleles();

    cout << "[" << Utils::getLocalTime() << "] Sorting variant clusters by decreasing complexity (number of variants) ..." << endl;

    sort(variant_cluster_graphs.begin(), variant_cluster_graphs.end(), VariantClusterGraphCompare);

    cout << "[" << Utils::getLocalTime() << "] Finished sorting\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    KmerHash * kmer_hash;

    if (samples.size() == 1) {

        kmer_hash = new HybridKmerHash<kmer_size, 1>(samples.size(), num_threads); 

    } else if ((samples.size() > 1) and (samples.size() < 4)) {

        kmer_hash = new HybridKmerHash<kmer_size, 3>(samples.size(), num_threads); 

    } else if ((samples.size() > 3) and (samples.size() < 11)) {

        kmer_hash = new HybridKmerHash<kmer_size, 10>(samples.size(), num_threads); 

    } else if ((samples.size() > 10) and (samples.size() < 21)) {

        kmer_hash = new HybridKmerHash<kmer_size, 20>(samples.size(), num_threads); 

    } else {

        assert(samples.size() <= 30);
        kmer_hash = new HybridKmerHash<kmer_size, 30>(samples.size(), num_threads);
    } 

    cout << "\n[" << Utils::getLocalTime() << "] Counting smallmers ..." << endl;

    Utils::SmallmerSet smallmer_set = Utils::SmallmerSet(Utils::thread_buckets);

    KmerCounter<kmer_size> kmer_counter(kmer_hash, &smallmer_set, num_threads);
    ulong num_unique_small_mers = kmer_counter.countVariantClusterSmallmers(&variant_cluster_graphs);

    cout << "[" << Utils::getLocalTime() << "] Counted " << num_unique_small_mers << " unique smallmers (" << to_string(Utils::small_kmer_size) << " nt)\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in genomic regions between variant clusters including the decoy sequence(s) ..." << endl;

    ulong num_intercluster_kmers = kmer_counter.countInterclusterKmers(vcf_parser.getInterclusterRegions());

    cout << "[" << Utils::getLocalTime() << "] Counted " << num_intercluster_kmers << " unique kmers passing the smallmer filter\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Parsing and filtering KMC k-mer tables from " << samples.size() << " sample(s) ...\n" << endl;

    ulong num_sample_unique_kmers = kmer_counter.countSampleKmers(samples);
    min_sample_kmer_count = kmer_counter.minSampleKmerCount(); 
  
    cout << "\n[" << Utils::getLocalTime() << "] Counted " << num_sample_unique_kmers << " unique kmers passing smallmer filter\n" << endl;
 
    cout << "\n[" << Utils::getLocalTime() << "] " << Utils::getMaxMemoryUsage() << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Sorting variant cluster groups by decreasing complexity (number of variants) ..." << endl;

    sort(variant_cluster_groups->begin(), variant_cluster_groups->end(), VariantClusterGroupCompare);

    cout << "[" << Utils::getLocalTime() << "] Finished sorting" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in variant clusters ..." << endl;

    kmer_counter.countVariantClusterKmers(variant_cluster_groups, prng_seed, samples.size(), max_sample_haplotype_candidates);

    auto diploid_kmer_counts = kmer_hash->calculateKmerStats(num_genomic_rate_gc_bias_bins);

    cout << "\n[" << Utils::getLocalTime() << "] Estimating haploid kmer count distribtion(s) ...\n" << endl;

    estimateGenomicCountDistributions(diploid_kmer_counts, samples);
    
    return kmer_hash; 
}


void KmerFactory::estimateGenomicCountDistributions(const vector<vector<vector<ulong> > > & diploid_kmer_counts, const vector<Sample> & samples) {

    assert(genomic_count_distributions.empty());
    genomic_count_distributions.reserve(diploid_kmer_counts.size());

    ofstream kmer_coverage_estimates_file(output_prefix + "_kmer_coverage_estimates.txt");
    assert(kmer_coverage_estimates_file.is_open());

    kmer_coverage_estimates_file << "Sample\tMean\tVariance" << endl;

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        genomic_count_distributions.emplace_back(vector<NegativeBinomialDistribution>());
        genomic_count_distributions.back().reserve(num_genomic_rate_gc_bias_bins);

        vector<ulong> sample_diploid_kmer_counts(Utils::uchar_overflow + 1, 0);

        for (ushort bias_idx = 0; bias_idx < num_genomic_rate_gc_bias_bins; bias_idx++) {

            ulong num_diploid_kmers = 0;

            assert(diploid_kmer_counts.at(sample_idx).at(bias_idx).size() == (Utils::uchar_overflow + 1));

            for (ushort count_idx = 0; count_idx < sample_diploid_kmer_counts.size(); count_idx++) {

                num_diploid_kmers += diploid_kmer_counts.at(sample_idx).at(bias_idx).at(count_idx);
                sample_diploid_kmer_counts.at(count_idx) += diploid_kmer_counts.at(sample_idx).at(bias_idx).at(count_idx);
            }

            genomic_count_distributions.back().emplace_back(NegativeBinomialDistribution::methodOfMomentsEst(diploid_kmer_counts.at(sample_idx).at(bias_idx)));
            genomic_count_distributions.back().back().size(genomic_count_distributions.back().back().size()/2);

            cout << "[" << Utils::getLocalTime() << "] Estimated fixed negative binomial distribution for sample " << samples.at(sample_idx).name << " and GC bias bin " << bias_idx + 1 << " with mean " << genomic_count_distributions.back().back().mean() << " and variance " << genomic_count_distributions.back().back().var() << " using " << num_diploid_kmers << " kmers" << endl;
        }

        cout << endl;

        NegativeBinomialDistribution sample_genomic_count_distributions = NegativeBinomialDistribution(NegativeBinomialDistribution::methodOfMomentsEst(sample_diploid_kmer_counts));
        sample_genomic_count_distributions.size(sample_genomic_count_distributions.size()/2);

        kmer_coverage_estimates_file << samples.at(sample_idx).name << "\t" << sample_genomic_count_distributions.mean() << "\t" << sample_genomic_count_distributions.var() << endl;
    }

    kmer_coverage_estimates_file.close();
}


ulong KmerFactory::numberOfVariants() {

	return number_of_variants;
}


ushort KmerFactory::maxAlternativeAlleles() {

	return max_alternative_alleles;
}


uchar KmerFactory::minSampleKmerCount() {

    return min_sample_kmer_count;
}


const vector<vector<NegativeBinomialDistribution> > & KmerFactory::genomicCountDistributions() {

    return genomic_count_distributions;
}


template KmerHash * KmerFactory::initKmerHash<31>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::initKmerHash<39>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::initKmerHash<47>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::initKmerHash<55>(vector<VariantClusterGroup*> *, const vector<Sample> &);
template KmerHash * KmerFactory::initKmerHash<63>(vector<VariantClusterGroup*> *, const vector<Sample> &);




