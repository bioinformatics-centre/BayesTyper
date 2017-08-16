
/*
VariantClusterGraphPath.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef __bayesTyper__VariantClusterGraphPath_hpp
#define __bayesTyper__VariantClusterGraphPath_hpp

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <algorithm>

#include "boost/functional/hash.hpp"
#include "boost/graph/adjacency_list.hpp"

#include "Utils.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterGraphVertex.hpp"

using namespace std;


template <uchar kmer_size>
class VariantClusterGraphPath {

	private:

		vector<pair<uint, const VariantClusterGraphVertex *> > path;
		pair<uint, uint> path_score;

	public:

		VariantClusterGraphPath(const ushort);

		void addVertex(const uint, const VariantClusterGraphVertex &);
		void updateScore(KmerCounts *, const ushort);

		const vector<pair<uint, const VariantClusterGraphVertex *> > & getPath() const;

		double getPathKmerScore() const;
		uint getPathVertexScore(const unordered_set<uint> &) const;

		void updateCoveredVertices(unordered_set<uint> *) const;

		bool operator == (const VariantClusterGraphPath<kmer_size> &) const;

		KmerPair<kmer_size> kmer_pair;		
};


#endif