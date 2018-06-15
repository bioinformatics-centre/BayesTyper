
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


#ifndef __KmerBloom__KmerBloom_hpp
#define __KmerBloom__KmerBloom_hpp

#include <string>
#include <iostream>
#include <bitset>
#include <array>
#include <assert.h>
#include <mutex>
#include <vector>

#include "ntHash/BloomFilter.hpp"


typedef unsigned char uchar;
typedef unsigned short int ushort;
typedef unsigned int uint;

template <uchar kmer_size>
class KmerBloom {

	public:

		KmerBloom(const uint64_t, const float);
		KmerBloom(const std::string &);
		~KmerBloom();

		void addKmer(const char *);
		void addKmer(const std::string &);
		void addKmer(const std::bitset<kmer_size*2> &);

		bool lookup(const char *) const;
		bool lookup(const std::bitset<kmer_size*2> &) const;
		bool lookup(const std::string &) const;

		void save(const std::string &) const;

		static std::string bitToNt(const std::bitset<kmer_size * 2> &);
		static uint64_t calcOptNumBloomBits(const float, const uint64_t);
		static uint calcOptNumHashes(const uint64_t, const uint64_t);

	private:

		uint64_t num_kmers;
		uint64_t num_bloom_bits;
		
		BloomFilter * bloom;
};

template <uchar kmer_size>
class ThreadedKmerBloom {

	public:

		ThreadedKmerBloom(const uint64_t, const float);
		// ThreadedKmerBloom(const std::string &);
		~ThreadedKmerBloom();

		void addKmer(const char *);
		void addKmer(const std::string &);
		void addKmer(const std::bitset<kmer_size*2> &);

		bool lookup(const char *) const;
		bool lookup(const std::string &) const;
		bool lookup(const std::bitset<kmer_size*2> &) const;

		// void save(const std::string &) const;

		std::unique_lock<std::mutex> getKmerLock(const string & kmer);
		std::unique_lock<std::mutex> getKmerLock(const std::bitset<kmer_size*2> & kmer);

	private:

		uint rootIndex(const char *) const;
		uint rootIndex(const string &) const;

		std::vector<KmerBloom<kmer_size> *> blooms;
		std::vector<std::mutex> mutexes;
};

#endif
