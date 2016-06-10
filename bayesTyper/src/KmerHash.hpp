
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

        virtual vector<NegativeBinomialDistribution> estimateGenomicCountDistributions(const vector<Sample> &) = 0;

        virtual void setMinKmerCount(const uchar) = 0;
        virtual ushort getNumberOfSamples() = 0;

        virtual void markKmers() = 0;


};


template<uchar kmer_size>
class KmerHashHybrid : public KmerHash {

	private:

		static const uint hash_root_size = 28;

        const ushort num_samples;
        uchar min_kmer_count;

    public:

        typedef typename ThreadedHybridHash<KmerCounts *, hash_root_size, kmer_size * 2 - hash_root_size>::iterator HashIterator;

    	KmerHashHybrid(const ushort, const ushort);
    	~KmerHashHybrid();

    	ThreadedHybridHash<KmerCounts *, hash_root_size, kmer_size * 2 - hash_root_size> * _hash;

        virtual vector<NegativeBinomialDistribution> estimateGenomicCountDistributions(const vector<Sample> &);

        void setMinKmerCount(const uchar);
        ushort getNumberOfSamples();

        void markKmers();

};


#endif
