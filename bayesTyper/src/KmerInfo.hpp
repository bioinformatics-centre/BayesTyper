
/*
KmerInfo.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__KmerInfo_hpp
#define __bayesTyper__KmerInfo_hpp

#include <vector>
#include <unordered_map>

#include "Utils.hpp"
#include "KmerCounts.hpp"

using namespace std;

struct PathVariants {

	const ushort path_idx;
	vector<ushort> variant_ids;

	PathVariants(const ushort & path_idx_in) : path_idx(path_idx_in) {}
};


class KmerInfo {

	public:

		KmerCounts * kmer_counts;		
		vector<ushort> variant_ids;
		
		KmerInfo(KmerCounts * kmer_counts_in) : kmer_counts(kmer_counts_in) {}		
};


template <uchar kmer_size>
class PathsKmerInfo {

	private:

		typedef vector<pair<bitset<kmer_size * 2>, KmerInfo> > KmerMultiplicities;

	public:

		KmerCounts * kmer_counts;		
		vector<uchar> multiplicities;
		vector<PathVariants> path_variants;

		bool has_max_multiplcity;

		PathsKmerInfo(KmerCounts * kmer_counts_in, const ushort num_elements) : kmer_counts(kmer_counts_in), multiplicities(num_elements, 0), has_max_multiplcity(false) {}
};

#endif