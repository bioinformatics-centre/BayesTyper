
/*
VariantClusterGenotyper.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "Eigen/Dense"

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
#include "KmerStats.hpp"
#include "Filters.hpp"
#include "SparsityEstimator.hpp"


using namespace std;

class CountDistribution;
class MixedFrequencyDistribution;

class VariantClusterGenotyper {

	public: 

		VariantClusterGenotyper(VariantClusterGraph *, KmerCountsHash *, const vector<Sample> &, const uint, const uchar);
		~VariantClusterGenotyper();

		void reset(const float, const uint);
		void clearCache();

		void sampleDiplotypes(const CountDistribution &, const vector<VariantClusterHaplotypes::NestedVariantClusterInfo> &, const bool);
		void getNoiseCounts(CountAllocation *, const CountDistribution &);
		void sampleHaplotypeFrequencies();

		void updateNestedVariantClusterInfo(vector<VariantClusterHaplotypes::NestedVariantClusterInfo> *, const uint);
		vector<Genotypes*> getGenotypes(const string &, const vector<Utils::Ploidy> &, const Filters & filters);
		
	private: 

		void updateNestedPloidy(Utils::Ploidy *);
		void addNestedKmerStats(vector<KmerStats> *, const uint, const ushort, const vector<KmerStats> &);

		void updateMulticlusterDiplotypeLogProb(const CountDistribution &, const ushort);
		double calcDiplotypeLogProb(const CountDistribution &, const ushort, const pair<ushort, ushort> &);

		void sampleDiplotype(const vector<ushort> &, const CountDistribution &, const ushort, const Utils::Ploidy);

		ushort haplotypeToAlleleIndex(const ushort, const ushort);
		vector<ushort> getNonCoveredAlleles(const ushort);
		vector<Genotypes::SampleStats> getGenotypeSampleStats(const uint, const vector<Utils::Ploidy> &, const Filters & filters);
		Genotypes::VariantStats getGenotypeVariantStats(const uint, const vector<Genotypes::SampleStats> &);

		const vector<Sample> & samples;

		mt19937 prng;

		bool use_multicluster_kmers; 

		VariantClusterHaplotypes variant_cluster_haplotypes;

		vector<VariantInfo> variant_cluster_info;
		vector<vector<AlleleKmerStats> > allele_kmer_stats;

		vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > > unique_diplotype_log_probabilities;
		vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > > multicluster_diplotype_log_probabilities;

		vector<pair<ushort,ushort> > diplotypes;
		unordered_map<pair<ushort,ushort>, vector<uint>, boost::hash<pair<ushort, ushort> > > diplotype_sampling_frequencies;

		HaplotypeFrequencyDistribution * haplotype_frequency_distribution;
};

#endif