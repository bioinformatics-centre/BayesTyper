
/*
VariantClusterGraph.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__VariantClusterGraph_hpp
#define __bayesTyper__VariantClusterGraph_hpp

#include <vector>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <list>
#include <atomic>
#include <random>
#include <deque>

#include "boost/graph/adjacency_list.hpp"
#include "../Eigen/Dense"

#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantKmerStats.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "KmerInfo.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "VariantClusterGraphKmerPath.hpp"
#include "PartitionGraph.hpp"


using namespace std;

class VariantClusterGraph {

	protected:

		typedef pair<string::const_iterator, string::const_iterator> StringItPair;

		typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS, VariantClusterGraphVertex *> Graph;
		
		typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
		typedef boost::graph_traits<Graph>::in_edge_iterator in_edge_it;

		Graph graph;

		vector<VariantInfo> variant_cluster_info;

		uchar has_ambiguous_nucleotide : 1, has_complex_region : 1, has_redundant_sequence : 1, has_non_unique_kmer : 1, has_excluded_kmer : 1;

		void addVertex(vertex_t, vector<StringItPair>, vector<uint>, pair<ushort, ushort>, unordered_set<ushort> &, const uchar);
		void addSequence(vector<vector<bool> > *, string::const_iterator, string::const_iterator);

		ulong calculateVertexSmallmerComplexity(const vertex_t, uint, uchar, const ushort);
		ulong countVertexSmallmers(Utils::SmallmerSet *, KmerPair<Utils::small_kmer_size>, uint, const vertex_t, set<vertex_t> *, ushort, const ushort, ushort);

	public:

		VariantClusterGraph(VariantCluster *, string &, const uchar);
		virtual ~VariantClusterGraph();

		ulong countSmallmers(Utils::SmallmerSet *);
		const vector<VariantInfo> & getInfo();

		bool hasAmbiguousNucleotide();
		bool hasComplexRegion();
		bool hasRedundantSequence();
		bool hasNonUniqueKmer();
		bool hasExcludedKmer();

		virtual ulong countKmers(KmerHash *, const uint, const ushort, const uint, double *) = 0;
		virtual void getBestHaplotypeCandidates(KmerHash *, VariantClusterHaplotypes *, const uint, const ushort) = 0;

		// DEBUG
		// set<uint> variant_ids;
		// void printInfoAndGraph();
};


template <uchar kmer_size>
class VariantClusterGraphLockedKmerSize : public VariantClusterGraph {

	private:

		typedef vector<pair<bitset<kmer_size * 2>, KmerInfo> > KmerMultiplicities;
		typedef unordered_map<bitset<kmer_size * 2>, PathsKmerInfo<kmer_size> > KmerMultiplicitiesIndex;

		bool updateKmerCounts(SampleKmerCounts *, const ushort, const bool);

		template <typename PathType>
		vector<vector<PathType *> > * findBestPathsPerSample(KmerHashHybrid<kmer_size> *, const ushort, vector<KmerMultiplicities *> *, mt19937 *, const ushort);

		template <typename PathType>
		void mergePaths(const in_edge_it, const vertex_t, vector<PathType *> *, vector<PathType *> *);

		bool hasEqualSourceSequence(in_edge_it, const in_edge_it);

		template <typename PathType>
		bool swapRedundantPath(PathType &, PathType &);

		template <typename PathType>
		void addKmer(PathType *, KmerHashHybrid<kmer_size> *, const ushort, const bool);

		template <typename PathType>
		void filterBestPaths(vector<PathType *> *, unordered_set<const VariantClusterGraphVertex *> *, const uint, const bool);

		template <typename PathType>
		void deleteVisitedVertexPaths(vertex_t, unordered_map<vertex_t, vector<vector<PathType *> > * > *, unordered_map<vertex_t, vector<vertex_t> > *);

		template <typename PathType>
		vector<PathType *> collapsePaths(vector<vector<PathType *> > *, const uint);

		template <typename PathType>
		KmerMultiplicitiesIndex createKmerMultiplicitiesIndex(vector<PathType *> *, unordered_set<KmerMultiplicities *> *, const bool);

	public:

		VariantClusterGraphLockedKmerSize(VariantCluster *, string &);

		ulong countKmers(KmerHash *, const uint, const ushort, const uint, double *);
		void getBestHaplotypeCandidates(KmerHash *, VariantClusterHaplotypes *, const uint, const ushort);
};


#endif