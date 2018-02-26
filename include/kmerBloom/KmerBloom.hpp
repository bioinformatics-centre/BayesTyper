
/*
KmerBloom.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__KmerBloom_hpp
#define __bayesTyper__KmerBloom_hpp

#include <string>
#include <iostream>
#include <bitset>
#include <array>
#include <assert.h>
#include <mutex>


#include "libbf/bf/bloom_filter/basic.hpp"

typedef unsigned char uchar;
typedef unsigned short int ushort;
typedef unsigned int uint;

class KmerBloom {

public:

	KmerBloom() {};
	virtual ~KmerBloom() {};
};

template <uchar kmer_size>
class BasicKmerBloom : public KmerBloom {

public:

	BasicKmerBloom(uint64_t, float);
	BasicKmerBloom(const std::string &);
	~BasicKmerBloom();

	void addKmer(const std::bitset<kmer_size*2> &);
	void addKmer(const std::string &);

	bool lookup(const std::bitset<kmer_size*2> &) const;
	bool lookup(const std::string &) const;

	void save(const std::string &) const;

	static uint64_t ntPrehash(const std::bitset<kmer_size*2> &);
	static uint64_t ntPrehash(const std::string &);

private:

	static void writeBloomBitvector(const bf::bitvector &, std::ofstream &);
	static bf::bitvector readBloomBitvector(std::ifstream &, uint64_t);

	static void processReadBuffer(bf::bitvector *, uint64_t, char *, uint, uchar);

	uint64_t num_kmers;
	const uint seed = 0;
	const bool double_hashing = true;
	const bool partition = false;
	
	bf::basic_bloom_filter * bloom;
};

template <uchar kmer_size>
class ThreadedBasicKmerBloom : public KmerBloom {

public:

	ThreadedBasicKmerBloom(const uint, const uint64_t, const float);
	~ThreadedBasicKmerBloom();

	void addKmer(const std::bitset<kmer_size*2> &);
	bool lookup(const std::bitset<kmer_size*2> &) const;

	std::unique_lock<std::mutex> lockBloom(const std::bitset<kmer_size*2> & kmer);

private:

	uint64_t rootIndex(const std::string &) const;

	std::vector<BasicKmerBloom<kmer_size> *> blooms;
	std::vector<std::mutex> mutexes;
};

#endif
