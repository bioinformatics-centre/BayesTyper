
/*
VariantClusterHaplotypes.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__VariantClusterHaplotypes_hpp
#define __bayesTyper__VariantClusterHaplotypes_hpp

#include <vector>
#include <unordered_map>
#include <random>

#include "boost/functional/hash.hpp"

#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "KmerStats.hpp"
#include "Sample.hpp"


class VariantClusterHaplotypes {

	public:

		VariantClusterHaplotypes() {};

		struct HaplotypeInfo {

			vector<ushort> variant_allele_indices;
			vector<uint> nested_variant_cluster_indices;

			HaplotypeInfo(const ushort num_variants) : variant_allele_indices(num_variants, Utils::ushort_overflow) {}
		};

		vector<HaplotypeInfo> haplotypes;

		Eigen::MatrixXuchar haplotype_unique_kmer_multiplicities;		
		Eigen::MatrixXuchar haplotype_multicluster_kmer_multiplicities;		

		struct KmerInfo {

			KmerCounts * const counts;
			const uchar bias_idx;
			const vector<pair<ushort, vector<bool> > > variant_haplotype_indices;

			KmerInfo(KmerCounts * const counts_in, const uchar bias_idx_in, const vector<pair<ushort, vector<bool> > > & variant_haplotype_indices_in) : counts(counts_in), bias_idx(bias_idx_in), variant_haplotype_indices(variant_haplotype_indices_in) {}
		};

		vector<KmerInfo> unique_kmers;
		vector<KmerInfo> multicluster_kmers;

		vector<uint> unique_kmer_subset_indices;
		vector<uint> multicluster_kmer_subset_indices;

		Eigen::MatrixXuchar sample_multicluster_kmer_multiplicities;

		struct VariantKmerStatsCache {

			bool update;
			
			vector<KmerStats> allele_kmer_stats_1;
			vector<KmerStats> allele_kmer_stats_2;

			VariantKmerStatsCache(const ushort num_variants) {

				update = true;

				allele_kmer_stats_1 = vector<KmerStats>(num_variants);
				allele_kmer_stats_2 = vector<KmerStats>(num_variants);
			}
		};
		
		vector<VariantKmerStatsCache> variant_kmer_stats_caches;

		bool empty();

		uchar getUniqueKmerMultiplicity(const uint, const pair<ushort,ushort> &, const Utils::Gender);
		uchar getMulticlusterKmerMultiplicity(const uint, const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort, const Utils::Gender);
		uchar getPreviousMulticlusterKmerMultiplicity(const uint, const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort, const Utils::Gender);

		bool sampleKmerSubset(mt19937 *, const uint max_haplotype_variant_kmers);

		bool isMulticlusterKmerUpdated(const uint, const ushort);
		void updateMulticlusterKmerMultiplicities(const pair<ushort,ushort> &, const pair<ushort,ushort> &, const ushort);
		void updateVariantKmerStats(vector<VariantKmerStats> *, const pair<ushort,ushort> &, const ushort, const Utils::Gender);	

	private:

		uchar getDiplotypeUniqueKmerMultiplicity(const uint, const pair<ushort,ushort> &);
		uchar getDiplotypeMulticlusterKmerMultiplicity(const uint, const pair<ushort,ushort> &);

		bool isKmerInSubset(vector<vector<uint> > *, const uint, const vector<pair<ushort, vector<bool> > > &);

		void updateVariantKmerStatsCaches(const pair<ushort,ushort> &, const ushort, const Utils::Gender);
		void updateVariantKmerStatsCache(KmerInfo &, const pair<ushort,ushort> &, const ushort, const float);
};

#endif