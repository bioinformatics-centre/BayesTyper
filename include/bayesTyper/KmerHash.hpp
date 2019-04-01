
/*
KmerCountsHash.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__KmerCountsHash_hpp
#define __bayesTyper__KmerCountsHash_hpp

#include <vector>
#include <unordered_map>
#include <string>
#include <mutex>

#include "KmerBloom.hpp"

#include "Utils.hpp"
#include "HybridHash.hpp"
#include "KmerCounts.hpp"
#include "KmerStats.hpp"
#include "Sample.hpp"


using namespace std;

template<class T>
class KmerHash {

    private:

        ThreadedHybridHash<T, Utils::kmer_size * 2> * _hash;

    public:

        KmerHash(const ulong, const ushort);
        ~KmerHash();

        unique_lock<mutex> getKmerLock(const bitset<Utils::kmer_size * 2> &);
        pair<T *, bool> addKmer(const bitset<Utils::kmer_size * 2> &);
        T * findKmer(const bitset<Utils::kmer_size * 2> &);

        void shuffle(const uint);
        ulong size();

        void writeRootSizeDistribution(const string &);
        ulong writeKmersToFasta(const string &, bool (*)(T), const uint);
        ulong addKmersToBloomFilter(KmerBloom<Utils::kmer_size> *, bool (*)(T));
};

class KmerCountsHash {

	public:

    	KmerCountsHash() {};
    	virtual ~KmerCountsHash() {};

        virtual unique_lock<mutex> getKmerLock(const bitset<Utils::kmer_size * 2> &) = 0;
        virtual pair<KmerCounts *, bool> addKmer(const bitset<Utils::kmer_size * 2> &, const bool) = 0;
        virtual KmerCounts * findKmer(const bitset<Utils::kmer_size * 2> &) = 0;
        virtual void sortKmers() = 0;

        virtual void writeRootSizeDistribution(const string &) = 0;
        virtual vector<vector<vector<KmerStats> > > calculateKmerStats(const vector<Sample> &) = 0;
};

template<uchar sample_bin>
class ObservedKmerCountsHash : public KmerCountsHash {

    private:

        ThreadedHybridHash<ObservedKmerCounts<sample_bin>, Utils::kmer_size * 2> * _hash;

    public:

        ObservedKmerCountsHash(const ulong, const ushort);
        ~ObservedKmerCountsHash();

        unique_lock<mutex> getKmerLock(const bitset<Utils::kmer_size * 2> &);
        pair<KmerCounts *, bool> addKmer(const bitset<Utils::kmer_size * 2> &, const bool);
        KmerCounts * findKmer(const bitset<Utils::kmer_size * 2> &);
        void sortKmers();

        void writeRootSizeDistribution(const string &);
        vector<vector<vector<KmerStats> > > calculateKmerStats(const vector<Sample> &);
};

#endif
