
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
#include "boost/graph/compressed_sparse_row_graph.hpp"

#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "VariantClusterGraphPath.hpp"


using namespace std;

class VariantClusterGraph {

	public:

		VariantClusterGraph(VariantCluster *, const string &, const uchar);
		virtual ~VariantClusterGraph() {};

		ulong countSmallmers(Utils::SmallmerSet *);
		const vector<VariantInfo> & getInfo();

		bool hasAmbiguousNucleotide();
		bool hasComplexRegion();
		bool hasRedundantSequence();
		bool hasInterclusterKmer();
		bool hasExcludedKmer();

		virtual double countKmers(KmerHash *, const uint, const uint, const ushort, const ushort) = 0;
		virtual void getBestHaplotypeCandidates(KmerHash *, VariantClusterHaplotypes *, const uint, const ushort, const ushort, const uchar) = 0;

	protected:

		typedef pair<string::const_iterator, string::const_iterator> StringItPair;

		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VariantClusterGraphVertex> Graph;
		
		typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
		typedef boost::graph_traits<Graph>::in_edge_iterator in_edge_it;

		Graph graph;

		vector<VariantInfo> variant_cluster_info;

		uchar has_ambiguous_nucleotide : 1, has_complex_region : 1, has_redundant_sequence : 1, has_intercluster_kmer : 1, has_excluded_kmer : 1;

		void addVertices(vertex_t *, const vector<StringItPair> &, const pair<ushort, ushort> &, const unordered_set<ushort> &, const vector<uint> &);
		void initVertex(vertex_t *, StringItPair, const pair<ushort, ushort> &, const vector<ushort> &, const uint);

		ulong calculateVertexSmallmerComplexity(const vertex_t, uint, uchar, const ushort);
		ulong countVertexSmallmers(Utils::SmallmerSet *, KmerPair<Utils::small_kmer_size>, uint, const vertex_t, uchar);
};


template <uchar kmer_size>
class VariantClusterKmerGraph : public VariantClusterGraph {

	public:

		VariantClusterKmerGraph(VariantCluster *, const string &);

		double countKmers(KmerHash *, const uint, const uint, const ushort, const ushort);
		void getBestHaplotypeCandidates(KmerHash *, VariantClusterHaplotypes *, const uint, const ushort, const ushort, const uchar);

	private:

		struct KmerPathInfo {

			KmerCounts * kmer_counts;
			vector<uchar> multiplicities;

			KmerPathInfo(KmerCounts * kmer_counts_in, const ushort num_paths) : kmer_counts(kmer_counts_in), multiplicities(num_paths, 0) {}
		};

		struct VariantKmerPathInfo : public KmerPathInfo {

			vector<pair<ushort, vector<bool> > > variant_path_indices;

			VariantKmerPathInfo(KmerCounts * kmer_counts, const ushort num_paths) : KmerPathInfo(kmer_counts, num_paths) {}
		};

		struct KmerMultiplicitiesIndex {

        	HybridHash<KmerPathInfo, 8, kmer_size * 2 - 8> index;
			
			uint num_kmers;
		
			KmerMultiplicitiesIndex() {

				num_kmers = 0;
			}
		};

		struct VariantKmerMultiplicitiesIndex {

        	HybridHash<VariantKmerPathInfo, 8, kmer_size * 2 - 8> index;
			
			uint num_unique_kmers;
			uint num_multicluster_kmers;
		
			VariantKmerMultiplicitiesIndex() {

				num_unique_kmers = 0;
				num_multicluster_kmers = 0;
			}
		};

		vector<VariantClusterGraphPath<kmer_size> *> findBestPaths(KmerHash *, mt19937 *, const ushort, const ushort);
		
		void mergePaths(vector<VariantClusterGraphPath<kmer_size> *> *, vector<VariantClusterGraphPath<kmer_size> *> *, const bool);
		bool swapRedundantPath(VariantClusterGraphPath<kmer_size> &, VariantClusterGraphPath<kmer_size> &);
		void filterBestPaths(vector<VariantClusterGraphPath<kmer_size> *> *, const uint);
		
		vector<VariantClusterGraphPath<kmer_size> *> collapsePaths(vector<vector<VariantClusterGraphPath<kmer_size> *> > *);
		
		KmerMultiplicitiesIndex indexKmerMultiplicities(KmerHash *, const vector<VariantClusterGraphPath<kmer_size> *> &);
		VariantKmerMultiplicitiesIndex indexVariantKmerMultiplicities(KmerHash *, const vector<VariantClusterGraphPath<kmer_size> *> &);
};


#endif