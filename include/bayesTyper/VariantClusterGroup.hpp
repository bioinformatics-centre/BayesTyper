
/*
VariantClusterGroup.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__VariantClusterGroup_hpp
#define __bayesTyper__VariantClusterGroup_hpp

#include <unordered_map>
#include <unordered_set>
#include <random>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include "Utils.hpp"
#include "VariantClusterGenotyper.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "CountAllocation.hpp"
#include "CountDistribution.hpp"
#include "KmerHash.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "VariantInfo.hpp"
#include "Regions.hpp"
#include "VariantCluster.hpp"
#include "Filters.hpp"


class VariantClusterGroup {

	private:

		struct VariantClusterVertex {

			uint variant_cluster_idx;

			VariantClusterGraph * graph;
			VariantClusterGenotyper * genotyper;

			VariantClusterVertex() : genotyper(nullptr) {}
			VariantClusterVertex(const uint variant_cluster_idx_in) : variant_cluster_idx(variant_cluster_idx_in), graph(nullptr), genotyper(nullptr) {}

			friend class boost::serialization::access;
			template<class Archive>
	    	void serialize(Archive & ar, const unsigned int version) {

	    		ar & variant_cluster_idx;
	    		ar & graph;
	    	}
		};

		string chrom_name;

		uint start_position;
		uint end_position;

		uint num_variants;
		
		vector<VariantClusterVertex> vertices;
		vector<vector<uint> > out_edges; 

		vector<uint> source_vertices;

		friend class boost::serialization::access;
		template<class Archive>
    	void serialize(Archive & ar, const unsigned int version) {

    		ar & chrom_name;
    		ar & start_position;
    		ar & end_position;
    		ar & num_variants;
    		ar & vertices;
    		ar & out_edges;
    		ar & source_vertices;
    	}

		void runGibbsSample(const uint, const CountDistribution &, const vector<VariantClusterHaplotypes::NestedVariantClusterInfo> &, const bool);

	public:

		VariantClusterGroup() {}
		VariantClusterGroup(const vector<VariantCluster *> &, const vector<VariantClusterGraph *> &, const unordered_map<uint, uint> &);
		~VariantClusterGroup();

		uint numberOfVariants() const;
		uint numberOfVariantClusters() const;
		uint numberOfVariantClusterGroupTrees() const;

		string region() const;

		void findSamplePaths(KmerBloom<Utils::kmer_size> *, const uint, const ushort);
		
		void countPathKmers(unordered_set<bitset<Utils::kmer_size * 2> > *);
		void classifyPathKmers(KmerCountsHash *, KmerBloom<Utils::kmer_size> *);

		bool hasExcludedKmers() const;

		void initGenotyper(KmerCountsHash *, const vector<Sample> &, const uint, const uchar, const float, const uint);
		void clearGenotyperCache();
		void resetGroup();

		void shuffleBranchOrdering(const uint);

		void estimateGenotypes(const CountDistribution &, const ChromosomePloidy &, const bool);
		void getNoiseCounts(CountAllocation *, const CountDistribution &);

		void collectGenotypes(vector<Genotypes*> *, const ChromosomePloidy &, const Filters & filters);
};

bool VariantClusterGroupCompare(VariantClusterGroup *, VariantClusterGroup *);

#endif