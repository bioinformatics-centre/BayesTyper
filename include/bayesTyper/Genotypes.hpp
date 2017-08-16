
/*
Genotypes.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef __bayesTyper__Genotypes_hpp
#define __bayesTyper__Genotypes_hpp

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include "boost/functional/hash.hpp"

#include "Utils.hpp"
#include "KmerStats.hpp"
#include "VariantInfo.hpp"

using namespace std;


class Genotypes {

	public:

		Genotypes() {}

		string variant_cluster_id;
		string variant_cluster_group_id;

		ushort num_candidates;
		ushort variant_cluster_size;
		uint variant_cluster_group_size;

		bool has_redundant_sequence;

		VariantInfo variant_info;
		vector<ushort> non_covered_alleles;

		struct VariantStats {

			uint total_count;

			vector<uint> alt_allele_counts;
			vector<float> alt_allele_frequency;
			vector<float> allele_call_probabilities;

			VariantStats() {};
			VariantStats(const ushort num_alleles) : total_count(0), alt_allele_counts(num_alleles - 1, 0), alt_allele_frequency(num_alleles - 1, 0), allele_call_probabilities(num_alleles, 0) {}
		};

	    VariantStats variant_stats;

		struct SampleStats {

		    vector<ushort> genotype_estimate;
		    vector<float> genotype_posteriors;

		    vector<float> allele_posteriors;
			AlleleKmerStats allele_kmer_stats;

			bool is_allele_kmer_estimate_variant;

			SampleStats() {

				is_allele_kmer_estimate_variant = false;
			};

			SampleStats(const ushort num_alleles, const uint num_genotypes) : genotype_posteriors(num_genotypes, 0), allele_posteriors(num_alleles, 0) {

				genotype_estimate.reserve(2);
				is_allele_kmer_estimate_variant = false;
			}
		};

		vector<SampleStats> sample_stats;
};

#endif