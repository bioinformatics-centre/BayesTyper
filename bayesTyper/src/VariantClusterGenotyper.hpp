
/*
VariantClusterGenotyper.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__VariantClusterGenotyper_hpp
#define __bayesTyper__VariantClusterGenotyper_hpp

#include <string>
#include <vector>
#include <random> 
#include <unordered_map>
#include <limits>
#include <bitset>

#include "../Eigen/Dense"

#include "boost/functional/hash.hpp"

#include "Utils.hpp"
#include "VariantClusterGraph.hpp"
#include "KmerHash.hpp"
#include "KmerCounts.hpp"
#include "CountDistribution.hpp"
#include "HaplotypeFrequencyDistribution.hpp"
#include "DiscreteSampler.hpp"
#include "Genotypes.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "LinearMap.hpp"
#include "VariantKmerStats.hpp"


using namespace std;

class CountDistribution;
class MixedFrequencyDistribution;

class VariantClusterGenotyper {

	public: 

		VariantClusterGenotyper(const vector<Sample> &, const uint);
		~VariantClusterGenotyper();

		void initialise(const uint, const uint, KmerHash *, VariantClusterGraph *, const uchar, const ushort, const uint);

		void reset(const uint);

		void sampleGenotypes(const CountDistribution &, const vector<Utils::Ploidy> &, const bool);
		void sampleCountAllocations(const CountDistribution &);
		void sampleHaplotypeFrequencies();

		const vector<vector<unordered_map<uchar,ulong> > > & getNoiseCounts();
		vector<Utils::Ploidy> getVariantClusterPloidy(const uint);
		vector<Genotypes*> getGenotypes(const vector<Utils::Ploidy> &);

		bool hasInterclusterKmer();
		bool hasMulticlusterKmer();
		bool hasExcludedKmer();
		
	private: 

		void reduceSamplePloidy(Utils::Ploidy *);
		double calcDiplotypeLogProb(const CountDistribution &, const ushort, const pair<ushort, ushort> &, const pair<ushort, ushort> &);
		void sampleSingleGenotype(const CountDistribution &, const ushort, const Utils::Ploidy, const bool, const pair<ushort, ushort> &);
		void allocateKmerCounts(const uchar, const uchar, const ushort, const CountDistribution &);
		pair<ushort, ushort> haplotypeToAlleleId(const ushort, const ushort);

		const vector<Sample> & samples;
		const uint variant_cluster_index;

		mt19937 prng;

		uchar has_complex_region : 1, has_intercluster_kmer : 1, has_multicluster_kmer : 1, has_excluded_kmer : 1, has_redundant_sequence : 1;

		VariantClusterHaplotypes variant_cluster_haplotypes;
		vector<VariantInfo> variant_cluster_info;

		vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > > unique_diplotype_log_probabilities;

		bool use_multicluster_kmers;

		struct MulticlusterUpdateInfo {

			uchar update_multicluster_cache : 1, update_multicluster_multiplicity_indices : 1;

			MulticlusterUpdateInfo () {

				update_multicluster_cache = true;
				update_multicluster_multiplicity_indices = true;
			}
		};

		vector<MulticlusterUpdateInfo> multicluster_update_info;
		vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > > multicluster_diplotype_log_probabilities;

		unordered_map<pair<ushort,ushort>, vector<double>, boost::hash<pair<ushort, ushort> > > diplotype_sampling_frequencies;

		vector<pair<ushort,ushort> > current_diplotypes;
		vector<vector<unordered_map<uchar,ulong> > > noise_counts;

		HaplotypeFrequencyDistribution * haplotype_frequency_distribution;
};

#endif