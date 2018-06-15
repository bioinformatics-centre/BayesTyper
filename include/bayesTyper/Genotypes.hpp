
/*
Genotypes.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

		string chrom_name;
		VariantInfo variant_info;

		ushort variant_cluster_size;
		string variant_cluster_region;
		
		uint variant_cluster_group_size;
		string variant_cluster_group_region;

		ushort num_candidates;
		vector<ushort> non_covered_alleles;

		ushort num_homozygote_genotypes;

		struct VariantStats {

			uint total_count;

			vector<uint> alt_allele_counts;
			vector<float> alt_allele_frequency;

			float max_alt_allele_call_probability;
			vector<float> allele_call_probabilities;

			VariantStats() {};
			VariantStats(const ushort num_alleles) : total_count(0), alt_allele_counts(num_alleles - 1, 0), alt_allele_frequency(num_alleles - 1, 0), max_alt_allele_call_probability(0), allele_call_probabilities(num_alleles, 0) {}
		};

	    VariantStats variant_stats;

		struct SampleStats {

		    vector<ushort> genotype_estimate;
		    bool is_homozygote;

		    vector<float> genotype_posteriors;

		    vector<float> allele_posteriors;
			AlleleKmerStats allele_kmer_stats;

			vector<string> allele_filters;

			SampleStats(const ushort num_alleles, const uint num_genotypes) : is_homozygote(false), genotype_posteriors(num_genotypes, 0), allele_posteriors(num_alleles, 0), allele_filters(num_alleles, "P") {

				genotype_estimate.reserve(2);
			}
		};

		vector<SampleStats> sample_stats;
};

#endif