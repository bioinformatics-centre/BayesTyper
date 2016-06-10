
/*
VariantKmerStats.cpp - This file is part of BayesTyper (v0.9)


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


#include "VariantKmerStats.hpp"
#include "KmerCounts.hpp"
#include "KmerCoverage.hpp"

using namespace std;

VariantKmerStats::VariantKmerStats() {}


VariantKmerStats::VariantKmerStats(const ushort num_alleles, const ushort num_samples) {

    mean_num_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, -1));
	mean_num_observed_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, -1));
    mean_num_unique_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, -1));

    sum_weights_num_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, 0));
    sum_weights_num_observed_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, 0));
    sum_weights_num_unique_kmers = vector<vector<float> >(num_samples, vector<float>(num_alleles, 0));
}


void VariantKmerStats::addKmerCoverage(const ushort allele_idx, const KmerCoverage & kmer_coverage, const vector<double> & sample_weights) {

    assert(mean_num_kmers.size() == sample_weights.size());
    assert(mean_num_kmers.size() == kmer_coverage.getNumberOfSamples());

    for (ushort sample_idx = 0; sample_idx < kmer_coverage.getNumberOfSamples(); sample_idx++) {

        if (sample_weights.at(sample_idx) > 0) {

            auto cur_num_kmers = kmer_coverage.numberOfKmers(sample_idx);

            mean_num_kmers.at(sample_idx).at(allele_idx) = updateMean(mean_num_kmers.at(sample_idx).at(allele_idx), sum_weights_num_kmers.at(sample_idx).at(allele_idx), cur_num_kmers, sample_weights.at(sample_idx));
            sum_weights_num_kmers.at(sample_idx).at(allele_idx) += sample_weights.at(sample_idx);           


            auto cur_num_observed_kmers = kmer_coverage.numberOfObservedKmers(sample_idx);

            mean_num_observed_kmers.at(sample_idx).at(allele_idx) = updateMean(mean_num_observed_kmers.at(sample_idx).at(allele_idx), sum_weights_num_observed_kmers.at(sample_idx).at(allele_idx), cur_num_observed_kmers, sample_weights.at(sample_idx));
            sum_weights_num_observed_kmers.at(sample_idx).at(allele_idx) += sample_weights.at(sample_idx);   


            auto cur_mean_num_unique_kmers = kmer_coverage.numberOfUniqueKmers(sample_idx);

            mean_num_unique_kmers.at(sample_idx).at(allele_idx) = updateMean(mean_num_unique_kmers.at(sample_idx).at(allele_idx), sum_weights_num_unique_kmers.at(sample_idx).at(allele_idx), cur_mean_num_unique_kmers, sample_weights.at(sample_idx));
            sum_weights_num_unique_kmers.at(sample_idx).at(allele_idx) += sample_weights.at(sample_idx);                 
        }
    }
}


float VariantKmerStats::updateMean(const float cur_mean, const float cur_sum_weights, const float new_value, const float new_value_weight) {

    return ((cur_mean * cur_sum_weights) + (new_value * new_value_weight))/(cur_sum_weights + new_value_weight);
}

