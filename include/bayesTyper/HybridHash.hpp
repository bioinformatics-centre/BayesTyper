
/*
HybridHash.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__HybridHash_hpp
#define __bayesTyper__HybridHash_hpp

#include <vector>
#include <bitset>
#include <mutex>
#include <memory>
#include <map>

#include "LinearMap.hpp"
#include "BitsetCompare.hpp"


typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

template<typename ValueType, uchar key_size> class HybridHashIterator;

template<typename ValueType, uchar key_size>
class HybridHash {
	
	friend class HybridHashIterator<ValueType, key_size>;
	
	public:

		typedef HybridHashIterator<ValueType, key_size> iterator;

		HybridHash(const uint, const ulong);
		virtual ~HybridHash() {};

		iterator begin();
		iterator end();

		std::pair<iterator, bool> insert(const std::bitset<key_size> &, ValueType, const bool);
		// iterator erase(iterator);

		iterator find(const std::bitset<key_size> &);
		
		void sort();
		void shuffle(const uint);
		ulong size();

	 	map<ulong, ulong> getRootSizeDistribution();

	protected:
		
		typedef PartialSortedLinearMap<std::bitset<key_size>, ValueType, BitsetLess<key_size> > LeafHash;

		std::vector<LeafHash> root_hash;
		std::hash<std::bitset<key_size> > key_hasher;

		ulong rootHashIndex(const std::bitset<key_size> &);
};


template<typename ValueType, uchar key_size>
class ThreadedHybridHash : public HybridHash<ValueType, key_size> {
	
	public:

		ThreadedHybridHash(const uint, const ulong, const ushort);
		std::unique_lock<std::mutex> lockKey(const std::bitset<key_size> &);
		void sort();		

	private:

		void sortCallback(const ushort);

		const ushort num_threads;
		std::vector<std::mutex> mutexes;
};


template<typename ValueType, uchar key_size>
class HybridHashIterator {

	friend class HybridHash<ValueType, key_size>;

	public:
		
		typedef typename HybridHash<ValueType, key_size>::LeafHash LeafHash;
		typedef typename std::vector<LeafHash> RootHash;

		HybridHashIterator(HybridHash<ValueType, key_size> *, typename RootHash::iterator, typename LeafHash::iterator);
		HybridHashIterator(const HybridHashIterator &);
		HybridHashIterator & operator=(HybridHashIterator);
		std::pair<std::bitset<key_size>, ValueType> & operator*();
		HybridHashIterator & operator++();
		HybridHashIterator operator++(int);
		bool operator==(HybridHashIterator<ValueType, key_size> const &);
		bool operator!=(HybridHashIterator<ValueType, key_size> const &);

	private:

		void swap(HybridHashIterator & first, HybridHashIterator & second);

		HybridHash<ValueType, key_size> * hybrid_hash;
		typename RootHash::iterator root_element_iter;
		typename LeafHash::iterator leaf_element_iter;  
};

#include "HybridHash.tpp"

#endif