
/*
VariantClusterGroup.hpp - This file is part of BayesTyper (v0.9)


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
#include <random>

#include "Utils.hpp"
#include "VariantClusterGenotyper.hpp"
#include "VariantClusterGraph.hpp"
#include "CountAllocation.hpp"
#include "CountDistribution.hpp"
#include "KmerHash.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "LinearMap.hpp"
#include "VariantInfo.hpp"
#include "Regions.hpp"
#include "VariantCluster.hpp"
#include "PartitionGraph.hpp"


class VariantClusterGroup {

	private:

		struct VariantClusterVertex {

			const uint variant_cluster_idx;
			const uint start_position;
			const uint end_position;

			VariantClusterGraph * graph;
			VariantClusterGenotyper * genotyper;

			VariantClusterVertex(const uint variant_cluster_idx_in, const uint start_position_in, const uint end_position_in) : variant_cluster_idx(variant_cluster_idx_in), start_position(start_position_in), end_position(end_position_in), graph(nullptr), genotyper(nullptr) {}
		};
		
		const Utils::ChromosomeClass chromosome_class;
		const string chromosome_id;

		uint start_position;
		uint end_position;

		double complexity; 
		uint number_of_variants;

		bool is_single_nucleotide_polymorphism;

		vector<VariantClusterVertex> vertices;
		vector<uint> source_vertices;
		unordered_map<uint, vector<uint> > out_edges;

		void runGibbsSample(const uint, const CountDistribution &, const vector<Utils::Ploidy> &, const bool);

	public:

		VariantClusterGroup(const vector<VariantCluster *> &, const vector<VariantClusterGraph *> &, const unordered_map<uint, uint> & , const Utils::ChromosomeClass, const string &);
		~VariantClusterGroup();

		double getComplexity() const;

		uint numberOfVariants() const;
		uint numberOfVariantClusters() const;
		uint numberOfVariantClusterGroupTrees() const;

		bool isAutosomal() const;
		bool isSingleNucleotidePolymorphism() const;
		bool hasAmbiguousNucleotide() const;
		bool hasComplexRegion() const;
		bool hasRedundantSequence() const;
		bool hasNonUniqueKmer() const;		
		bool hasInterclusterKmer() const;
		bool hasMulticlusterKmer() const;
		bool hasExcludedKmer() const;
		
		ulong countSmallmers(Utils::SmallmerSet *);
		ulong countKmers(KmerHash *, const uint, const ushort, const uint);

		bool isInChromosomeRegions(const Regions &);

		void initialise(KmerHash *, const vector<Sample> &, const uint, const uchar, const ushort, const uint);
		void shuffleBranches(mt19937 *);

		void estimateGenotypes(const CountDistribution &, const ChromosomePloidy &, const bool);
		void sampleCountAllocations(const CountDistribution &, CountAllocation *);

		void collectGenotypes(vector<Genotypes*> *, const ChromosomePloidy &);
};

bool VariantClusterGroupCompare(VariantClusterGroup *, VariantClusterGroup *);

#endif