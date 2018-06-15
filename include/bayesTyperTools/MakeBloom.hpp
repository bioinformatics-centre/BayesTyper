
/*
MakeBloom.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyperTools__MakeBloom_hpp
#define __bayesTyperTools__MakeBloom_hpp

#include <iostream>
#include <string>
#include <assert.h>
#include <bitset>
#include <random>
#include <unordered_set>

#include "kmc_api/kmc_file.h"

#include "KmerBloom.hpp"
#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"

static const uint kmer_size = BT_KMER_SIZE;
typedef KmerBloom<kmer_size> kmer_bloom_t;
typedef vector<CKmerAPI*> kmer_batch_t;

class MakeBloom {

	public:

		void kmc2bloomFile(const string &, const float, const bool, const uint);

	private:

		kmer_bloom_t * kmc2bloom(const string &, const float);
		kmer_bloom_t * kmc2bloomThreaded(const string &, const float, const uint);
		void bloomInsertCallBack(kmer_bloom_t *, ProducerConsumerQueue<kmer_batch_t*> *, ProducerConsumerQueue<kmer_batch_t*> *);

		std::bitset<kmer_size*2> kmc2bits(CKmerAPI &);
		std::bitset<kmer_size*2> random_kmer();
		void testbloom(const kmer_bloom_t &, const unordered_set<bitset<kmer_size * 2> > &);
};

#endif
