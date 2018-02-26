
/*
KmerBloom.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "nthash.hpp"

#include "KmerBloom.hpp"

// typedef unsigned short int ushort;

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

template<uchar kmer_size>
std::string bitToNt(const std::bitset<kmer_size * 2> & bit_seq) {

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

template <uchar kmer_size>
BasicKmerBloom<kmer_size>::BasicKmerBloom(const uint64_t num_kmers_in, const float fpr) : num_kmers(num_kmers_in) {

	assert(kmer_size <= 64);

	uint64_t opt_num_cells = bf::basic_bloom_filter::m(fpr, num_kmers);
	uint64_t opt_num_hashes = bf::basic_bloom_filter::k(opt_num_cells, num_kmers);
	opt_num_cells += opt_num_hashes - (opt_num_cells % opt_num_hashes); // Fixes known libbf issue: https://github.com/mavam/libbf/issues/18

	assert(opt_num_cells % opt_num_hashes == 0);
	auto opt_hasher = bf::make_hasher(opt_num_hashes, seed, double_hashing);
	bloom = new bf::basic_bloom_filter(opt_hasher, opt_num_cells, partition);
}

template <uchar kmer_size>
BasicKmerBloom<kmer_size>::BasicKmerBloom(const std::string & fn) {

	assert(kmer_size <= 64);

	std::ifstream kmer_bloom_file_in(fn, std::ios::binary);

	std::string bloom_meta_str;
	std::getline(kmer_bloom_file_in, bloom_meta_str);
	auto bloom_meta = split(bloom_meta_str, '\t');
	assert(bloom_meta.size() == 3);
	num_kmers = std::stol(bloom_meta.front());

	uint64_t num_bv_bits = std::stol(bloom_meta.at(1));
	uint num_hash_functions = bf::basic_bloom_filter::k(num_bv_bits, num_kmers);

	assert(std::stoi(bloom_meta.back()) == kmer_size);

	auto bloom_bv = readBloomBitvector(kmer_bloom_file_in, num_bv_bits);

	bloom = new bf::basic_bloom_filter(bf::make_hasher(num_hash_functions, seed, double_hashing), bloom_bv, partition);

	kmer_bloom_file_in.close();
}

template <uchar kmer_size>
BasicKmerBloom<kmer_size>::~BasicKmerBloom() {

	delete bloom;
}

template <uchar kmer_size>
void BasicKmerBloom<kmer_size>::writeBloomBitvector(const bf::bitvector & bits, std::ofstream & file) {

	std::bitset<8> byte;
	
	for (uint64_t i = 0; i < bits.size(); ++i) {

		if (i > 0 and i % 8 == 0) {

			file << static_cast<char>(byte.to_ulong());
			byte.reset();
		}

		byte[i % 8] = bits[i];
	}

	file << static_cast<char>(byte.to_ulong());
}

template <uchar kmer_size>
inline void BasicKmerBloom<kmer_size>::processReadBuffer(bf::bitvector * bv, uint64_t bv_offset, char * buffer, uint buf_size, uchar num_last_byte_bits) {

	uint64_t buf_idx = 0;

	while (buf_idx < (buf_size - 1)) {

		for (uchar bit_idx = 0; bit_idx < 8; ++bit_idx) {

			(*bv)[bv_offset + buf_idx*8 + static_cast<uint64_t>(bit_idx)] = (1 << bit_idx) & buffer[buf_idx];
		}

		++buf_idx;
	}

	// Process last byte separately
	for (uchar bit_idx = 0; bit_idx < num_last_byte_bits; ++bit_idx) {

		(*bv)[bv_offset + buf_idx*8 + static_cast<uint64_t>(bit_idx)] = (1 << bit_idx) & buffer[buf_idx];
	}
}

template <uchar kmer_size>
bf::bitvector BasicKmerBloom<kmer_size>::readBloomBitvector(std::ifstream & file, uint64_t bv_size) {

	bf::bitvector bv(bv_size);
	uint64_t bv_offset = 0;

	// while (file.good() and bv_offset < bv_size) {

	// 	std::bitset<8> byte = file.get();

	// 	for (uint i = 0; i < 8 and (bv_offset + i) < bv_size; ++i) {

	// 		bv[bv_offset + i] = byte[i];
	// 	}

	// 	bv_offset += 8;
	// }

	// NEW SOLUTION START
	uint buf_size = 1024; 
	char * buffer = new char[buf_size];

	while (file.read(buffer, buf_size)) {

		processReadBuffer(&bv, bv_offset, buffer, buf_size, 8);
		bv_offset += buf_size * 8;
	}

	processReadBuffer(&bv, bv_offset, buffer, file.gcount(), (bv.size() % 8 == 0) ? 8 : bv.size() % 8);		
	// NEW SOLUTION END

	delete[] buffer;

	return bv;
}

template <uchar kmer_size>
void BasicKmerBloom<kmer_size>::save(const std::string & fn) const {

	std::ofstream kmer_bloom_file_out(fn, std::ios::binary);

	kmer_bloom_file_out << std::to_string(num_kmers) << "\t" << std::to_string(bloom->storage().size()) << "\t" << std::to_string(kmer_size) << std::endl;
	writeBloomBitvector(bloom->storage(), kmer_bloom_file_out);
	kmer_bloom_file_out.close();
}

template <uchar kmer_size>
inline uint64_t BasicKmerBloom<kmer_size>::ntPrehash(const std::bitset<kmer_size*2> & kmer) {

	auto kmer_str = bitToNt<kmer_size>(kmer);
	return NT64(kmer_str.c_str(), kmer_size);
}

template <uchar kmer_size>
inline uint64_t BasicKmerBloom<kmer_size>::ntPrehash(const std::string & kmer) {

	return NT64(kmer.c_str(), kmer_size);
}

template <uchar kmer_size>
void BasicKmerBloom<kmer_size>::addKmer(const std::bitset<kmer_size*2> & kmer) {

	bloom->add(ntPrehash(kmer));
}

template <uchar kmer_size>
void BasicKmerBloom<kmer_size>::addKmer(const std::string & kmer) {

	bloom->add(ntPrehash(kmer));
}

template <uchar kmer_size>
bool BasicKmerBloom<kmer_size>::lookup(const std::bitset<kmer_size*2> & kmer) const {

	return bloom->lookup(ntPrehash(kmer));
}

template <uchar kmer_size>
bool BasicKmerBloom<kmer_size>::lookup(const std::string & kmer) const {

	return bloom->lookup(ntPrehash(kmer));
}


template <uchar kmer_size>
ThreadedBasicKmerBloom<kmer_size>::ThreadedBasicKmerBloom(const uint root_size, const uint64_t num_kmers, const float fpr) {

	assert(root_size > 0);

	blooms = std::vector<BasicKmerBloom<kmer_size> *>(root_size);
	mutexes = std::vector<std::mutex>(blooms.size());

	for (auto & bloom: blooms) {

		bloom = new BasicKmerBloom<kmer_size>(std::ceil(num_kmers / static_cast<float>(root_size)), fpr);
	}
}

template <uchar kmer_size>
ThreadedBasicKmerBloom<kmer_size>::~ThreadedBasicKmerBloom() {

	for (auto & bloom: blooms) {

		delete bloom;
	}
}

template <uchar kmer_size>
void ThreadedBasicKmerBloom<kmer_size>::addKmer(const std::bitset<kmer_size*2> & kmer) {

	auto kmer_str = bitToNt<kmer_size>(kmer);
	auto root_idx = rootIndex(kmer_str);

	blooms.at(root_idx)->addKmer(kmer_str);
}

template <uchar kmer_size>
bool ThreadedBasicKmerBloom<kmer_size>::lookup(const std::bitset<kmer_size*2> & kmer) const {

	auto kmer_str = bitToNt<kmer_size>(kmer);

	return blooms.at(rootIndex(kmer_str))->lookup(kmer_str);
}

template <uchar kmer_size>
std::unique_lock<std::mutex> ThreadedBasicKmerBloom<kmer_size>::lockBloom(const std::bitset<kmer_size*2> & kmer) {

	std::unique_lock<std::mutex> bloom_lock(mutexes.at(rootIndex(bitToNt<kmer_size>(kmer))));
	return std::move(bloom_lock);
}

template <uchar kmer_size>
uint64_t ThreadedBasicKmerBloom<kmer_size>::rootIndex(const std::string & kmer) const {

	return (BasicKmerBloom<kmer_size>::ntPrehash(kmer) % blooms.size());
}


template class BasicKmerBloom<31>;
template class BasicKmerBloom<39>;
template class BasicKmerBloom<47>;
template class BasicKmerBloom<55>;
template class BasicKmerBloom<63>;


template class ThreadedBasicKmerBloom<31>;
template class ThreadedBasicKmerBloom<39>;
template class ThreadedBasicKmerBloom<47>;
template class ThreadedBasicKmerBloom<55>;
template class ThreadedBasicKmerBloom<63>;
