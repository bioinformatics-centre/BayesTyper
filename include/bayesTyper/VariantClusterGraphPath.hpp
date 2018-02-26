
/*
VariantClusterGraphPath.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include <unordered_set>
#include <vector>

#include "Utils.hpp"
#include "Kmer.hpp"
#include "VariantClusterGraphVertex.hpp"
#include "KmerBloom.hpp"

using namespace std;


template <uchar kmer_size>
class VariantClusterGraphPath {

	public:

		struct PathVertexInfo {

			const uint index;
			const VariantClusterGraphVertex * const vertex;
			uchar num_observed_kmers;

			PathVertexInfo(const uint index_in, const VariantClusterGraphVertex * const vertex_in) : index(index_in), vertex(vertex_in) {

				num_observed_kmers = 0;
			}
		};

		VariantClusterGraphPath(const ushort);

		void addVertex(const uint, const VariantClusterGraphVertex &, KmerBloom *);

		const vector<PathVertexInfo> & getPath() const;

		double getKmerScore() const;
		uint getVertexScore(const vector<bool> &, const bool) const;

		void updateObservedCoveredVertices(vector<bool> *, const bool) const;

	private:

		KmerPair<kmer_size> kmer_pair;

		vector<PathVertexInfo> path;
		pair<uint, uint> kmer_score;

		void updateScore(const bool, uint);
};


#endif