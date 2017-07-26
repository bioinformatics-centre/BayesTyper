
/*
KmerHash.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__KmerHash_hpp
#define __bayesTyper__KmerHash_hpp

#include <vector>
#include <unordered_map>
#include <string>
#include <mutex>

#include "Utils.hpp"
#include "HybridHash.hpp"
#include "KmerCounts.hpp"
#include "Sample.hpp"
#include "NegativeBinomialDistribution.hpp"

using namespace std;


class KmerHash {

	public:

    	KmerHash() {};
    	virtual ~KmerHash() {};

        virtual vector<vector<vector<ulong> > > calculateKmerStats(const uchar) = 0;
};


template<uchar kmer_size>
class BasicKmerHash : public KmerHash {

    public:

        BasicKmerHash() {};
        ~BasicKmerHash() {};

        virtual unique_lock<mutex> getKmerLock(const bitset<kmer_size * 2> &) = 0;
        virtual pair<KmerCounts *, bool> addKmer(const bitset<kmer_size * 2> &, const bool) = 0;
        virtual KmerCounts * findKmer(const bitset<kmer_size * 2> &) = 0;
        virtual void sortKmers() = 0;
};


template<uchar kmer_size, uchar sample_bin>
class HybridKmerHash : public BasicKmerHash<kmer_size> {

	private:

		static const uint hash_root_size = 24;
        const ushort num_samples;

        ThreadedHybridHash<ObservedKmerCounts<sample_bin>, hash_root_size, kmer_size * 2 - hash_root_size> * _hash;

    public:

    	HybridKmerHash(const ushort, const ushort);
    	~HybridKmerHash();

        unique_lock<mutex> getKmerLock(const bitset<kmer_size * 2> &);
        pair<KmerCounts *, bool> addKmer(const bitset<kmer_size * 2> &, const bool);
        KmerCounts * findKmer(const bitset<kmer_size * 2> &);
        void sortKmers();

        vector<vector<vector<ulong> > > calculateKmerStats(const uchar);
};


#endif
