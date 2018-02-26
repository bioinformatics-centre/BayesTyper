
/*
bloomMaker.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef BLOOM_MAKER
#define BLOOM_MAKER

#include <iostream>
#include <string>
#include <assert.h>
#include <bitset>
#include <random>
#include <unordered_set>

#include "kmc_api/kmc_file.h"

#include "KmerBloom.hpp"

#include "Utils.hpp"

template <uchar kmer_size>
class BloomMaker {

public:

	static void kmc2bloomFile(const string &, float, bool test = false);
	static BasicKmerBloom<kmer_size> kmc2bloom(const string &, float);
	static std::bitset<kmer_size*2> kmc2bits(CKmerAPI &);
	static std::bitset<kmer_size*2> random_kmer();



private:

	static void testbloom(BasicKmerBloom<kmer_size> &, const unordered_set<bitset<kmer_size * 2> > &);
};

#endif
