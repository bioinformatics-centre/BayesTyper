
/*
VariantClusterGraphPath.hpp - This file is part of BayesTyper (v0.9)


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

#include "Utils.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "VariantClusterGraphVertex.hpp"

using namespace std;


template <uchar kmer_size>
class VariantClusterGraphPath {
	
	public:

		VariantClusterGraphPath(const ushort);

		bool operator == (const VariantClusterGraphPath<kmer_size> & rhs) const;

		void addVertex(const VariantClusterGraphVertex &);
		void updateScore(KmerCounts *, const ushort);

		const vector<const VariantClusterGraphVertex *> & getPath() const;

		double getPathKmerScore() const;
		uint getPathVertexScore(const unordered_set<const VariantClusterGraphVertex *> &) const;

		void updateCoveredVertices(unordered_set<const VariantClusterGraphVertex *> *) const;

		KmerPair<kmer_size> kmer_pair;

	private:

		vector<const VariantClusterGraphVertex *> path;
		pair<uint, uint> path_score;
};


#endif