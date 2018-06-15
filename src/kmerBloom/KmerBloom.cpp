
/*
KmerBloom<kmer_size>.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>
#include <cmath>

#include "KmerBloom.hpp"

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

template <uchar kmer_size>
KmerBloom<kmer_size>::KmerBloom(const uint64_t num_kmers_in, const float fpr) : num_kmers(num_kmers_in) {

	num_bloom_bits = calcOptNumBloomBits(fpr, num_kmers);
	uint opt_num_hashes = calcOptNumHashes(num_bloom_bits, num_kmers);

	bloom = new BloomFilter(num_bloom_bits, opt_num_hashes, kmer_size);
}

template <uchar kmer_size>
KmerBloom<kmer_size>::KmerBloom(const std::string & prefix) {

	std::ifstream meta_infile(prefix + ".bloomMeta");
	assert(meta_infile.is_open());

	std::string bloom_meta_str;
	
	std::getline(meta_infile, bloom_meta_str);
	meta_infile.close();

	auto bloom_meta = split(bloom_meta_str, '\t');
	assert(bloom_meta.size() == 3);
	
	num_kmers = std::stol(bloom_meta.front());
	num_bloom_bits = std::stol(bloom_meta.at(1));
	
	uint opt_num_hashes = calcOptNumHashes(num_bloom_bits, num_kmers);
	assert(std::stoi(bloom_meta.back()) == kmer_size);

	string kmer_bloom_data_fn = prefix + ".bloomData";
	bloom = new BloomFilter(num_bloom_bits, opt_num_hashes, kmer_size, kmer_bloom_data_fn.c_str());
}

template <uchar kmer_size>
KmerBloom<kmer_size>::~KmerBloom() {

	delete bloom;
}

template <uchar kmer_size>
std::string KmerBloom<kmer_size>::bitToNt(const std::bitset<kmer_size * 2> & bit_seq) {

	std::string nt_seq;
	nt_seq.reserve(kmer_size);

	for (ushort i = 0; i < (kmer_size * 2); i += 2) {

		if (bit_seq[i] == bit_seq[i + 1]) {

			if (bit_seq[i] == 0) {

				nt_seq += "A";
			
			} else {

				nt_seq += "T";
			}				
		
		} else {

			if (bit_seq[i] == 0) {

				nt_seq += "G";
			
			} else {

				nt_seq += "C";
			}					
		}
	}

    return nt_seq;
}

// COPIED VERBATIM FROM https://github.com/mavam/libbf/blob/master/src/bloom_filter/basic.cpp
template <uchar kmer_size>
uint64_t KmerBloom<kmer_size>::calcOptNumBloomBits(const float fpr, const uint64_t num_kmers) {
	
  auto ln2 = std::log(2);
  return std::ceil(-(num_kmers * std::log(fpr) / ln2 / ln2));
}

// COPIED VERBATIM FROM https://github.com/mavam/libbf/blob/master/src/bloom_filter/basic.cpp
template <uchar kmer_size>
uint KmerBloom<kmer_size>::calcOptNumHashes(const uint64_t num_bloom_bits, const uint64_t num_kmers) {
	
	auto frac = static_cast<double>(num_bloom_bits) / static_cast<double>(num_kmers);
  	return std::ceil(frac * std::log(2));
}

template <uchar kmer_size>
void KmerBloom<kmer_size>::save(const std::string & prefix) const {

	std::ofstream meta_outfile(prefix + ".bloomMeta");
	assert(meta_outfile.is_open());
	
	meta_outfile << std::to_string(num_kmers) << "\t" << std::to_string(num_bloom_bits) << "\t" << std::to_string(kmer_size) << std::endl;
	meta_outfile.close();

	string kmer_bloom_data_fn = prefix + ".bloomData";
	bloom->storeFilter(kmer_bloom_data_fn.c_str());
}

template <uchar kmer_size>
inline void KmerBloom<kmer_size>::addKmer(const char * kmer) {

	bloom->insertF(kmer);
}

template <uchar kmer_size>
inline void KmerBloom<kmer_size>::addKmer(const std::string & kmer) {

	addKmer(kmer.c_str());
}

template <uchar kmer_size>
inline void KmerBloom<kmer_size>::addKmer(const std::bitset<kmer_size*2> & kmer) {

	addKmer(bitToNt(kmer));
}

template <uchar kmer_size>
inline bool KmerBloom<kmer_size>::lookup(const char * kmer) const {

	return bloom->containsF(kmer);
}

template <uchar kmer_size>
inline bool KmerBloom<kmer_size>::lookup(const std::string & kmer) const {

	return lookup(kmer.c_str());
}

template <uchar kmer_size>
inline bool KmerBloom<kmer_size>::lookup(const std::bitset<kmer_size*2> & kmer) const {

	return lookup(bitToNt(kmer));
}


template <uchar kmer_size>
ThreadedKmerBloom<kmer_size>::ThreadedKmerBloom(const uint64_t num_kmers, const float fpr) {
		
	uint root_size = 65536;
	
	blooms = std::vector<KmerBloom<kmer_size> *>(root_size);
	mutexes = std::vector<std::mutex>(blooms.size());

	for (auto & bloom: blooms) {

		bloom = new KmerBloom<kmer_size>(std::ceil(num_kmers / static_cast<float>(root_size)), fpr);
	}
}

template <uchar kmer_size>
ThreadedKmerBloom<kmer_size>::~ThreadedKmerBloom() {

	for (auto & bloom: blooms) {

		delete bloom;
	}
}

template <uchar kmer_size>
void ThreadedKmerBloom<kmer_size>::addKmer(const char * kmer) {

	blooms.at(rootIndex(kmer))->addKmer(kmer);
}

template <uchar kmer_size>
void ThreadedKmerBloom<kmer_size>::addKmer(const std::string & kmer) {

	addKmer(kmer.c_str());
}

template <uchar kmer_size>
void ThreadedKmerBloom<kmer_size>::addKmer(const std::bitset<kmer_size*2> & kmer) {

	addKmer(KmerBloom<kmer_size>::bitToNt(kmer));
}

template <uchar kmer_size>
bool ThreadedKmerBloom<kmer_size>::lookup(const char * kmer) const {

	return blooms.at(rootIndex(kmer))->lookup(kmer);
}

template <uchar kmer_size>
bool ThreadedKmerBloom<kmer_size>::lookup(const std::string & kmer) const {

	return lookup(kmer.c_str());
}

template <uchar kmer_size>
bool ThreadedKmerBloom<kmer_size>::lookup(const std::bitset<kmer_size*2> & kmer) const {

	return lookup(KmerBloom<kmer_size>::bitToNt(kmer));
}

template <uchar kmer_size>
std::unique_lock<std::mutex> ThreadedKmerBloom<kmer_size>::getKmerLock(const std::string & kmer) {

	std::unique_lock<std::mutex> bloom_lock(mutexes.at(rootIndex(kmer)));
	return std::move(bloom_lock);
}

template <uchar kmer_size>
std::unique_lock<std::mutex> ThreadedKmerBloom<kmer_size>::getKmerLock(const std::bitset<kmer_size*2> & kmer) {

	std::unique_lock<std::mutex> bloom_lock(mutexes.at(rootIndex(KmerBloom<kmer_size>::bitToNt(kmer))));
	return std::move(bloom_lock);
}

template <uchar kmer_size>
uint ThreadedKmerBloom<kmer_size>::rootIndex(const char * kmer) const {

	return NTP64(kmer, kmer_size, 1029283129) % blooms.size(); 
}

template <uchar kmer_size>
uint ThreadedKmerBloom<kmer_size>::rootIndex(const string & kmer) const {

	return rootIndex(kmer.c_str()); 
}


template class KmerBloom<BT_KMER_SIZE>;
template class ThreadedKmerBloom<BT_KMER_SIZE>;