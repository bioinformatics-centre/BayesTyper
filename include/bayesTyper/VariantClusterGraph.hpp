
/*
VariantClusterGraph.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include "boost/graph/adj_list_serialize.hpp"

#include "boost/functional/hash.hpp"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include "KmerBloom.hpp"
#include "Utils.hpp"
#include "VariantCluster.hpp"
#include "KmerHash.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "VariantInfo.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "VariantClusterGraphPath.hpp"


using namespace std;

class VariantClusterGraph {

	public:

		VariantClusterGraph() {};
		VariantClusterGraph(VariantCluster *, const string &);
		virtual ~VariantClusterGraph() {};

		const vector<VariantInfo> & getInfo() const;
		bool isSimpleCluster();

		void findSamplePaths(KmerBloom<Utils::kmer_size> *, const uint, const ushort);
		
		void countPathKmers(unordered_set<bitset<Utils::kmer_size * 2> > *);
		void classifyPathKmers(KmerCountsHash *, KmerBloom<Utils::kmer_size> *);
		
		VariantClusterHaplotypes getHaplotypeCandidates(KmerCountsHash *, const uchar);

	private:

		typedef pair<string::const_iterator, string::const_iterator> StringItPair;

		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VariantClusterGraphVertex> Graph;
		typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;

		Graph graph;

		ulong num_path_kmers;
		bool is_simple_cluster;

		vector<VariantInfo> variant_cluster_info;
		vector<bool> best_paths_indices;

		void addVertices(vertex_t *, const vector<StringItPair> &, const pair<ushort, ushort> &, const unordered_set<ushort> &, const vector<uint> &, const bool);
		void initVertex(vertex_t *, StringItPair, const pair<ushort, ushort> &, const vector<ushort> &, const uint, const bool);

		void mergePaths(vector<VariantClusterGraphPath *> *, vector<VariantClusterGraphPath *> *, const bool);
		bool isPathsRedundant(const vector<typename VariantClusterGraphPath::PathVertexInfo> &, const vector<typename VariantClusterGraphPath::PathVertexInfo> &) const;
		void filterPaths(vector<VariantClusterGraphPath *> *, const uint, const bool);		
		void addPathIndices(vector<VariantClusterGraphPath *> *);

		void updateVariantPathIndices(vector<pair<ushort, vector<bool> > > *, unordered_map<pair<ushort, ushort>, pair<uint, uint>, boost::hash<pair<ushort, ushort> > > *, const uint, const ushort, const ushort);

		friend class boost::serialization::access;

		template<class Archive>
    	void serialize(Archive & ar, const unsigned int version) {

    		ar & graph;
    		ar & num_path_kmers;
    		ar & is_simple_cluster;
    		ar & variant_cluster_info;
    		ar & best_paths_indices;
    	}

};

bool VariantClusterGraphCompare(VariantClusterGraph *, VariantClusterGraph *);


#endif