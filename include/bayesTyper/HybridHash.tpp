
/*
HybridHash.tpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <map>
#include <thread>
#include <limits>
#include <functional>
#include <random>


template<typename ValueType, uchar key_size>
HybridHash<ValueType, key_size>::HybridHash(const uint root_hash_size, const ulong expected_total_size) {

	assert(key_size <= std::numeric_limits<uchar>::max());

	root_hash = std::vector<LeafHash>(root_hash_size);

	const uint expected_size_leaf_hash_size = std::ceil(expected_total_size / static_cast<float>(root_hash_size));

	for (auto & leaf_hash: root_hash) {

		leaf_hash.reserve(expected_size_leaf_hash_size);
	}
}

template<typename ValueType, uchar key_size>
typename HybridHash<ValueType, key_size>::iterator HybridHash<ValueType, key_size>::begin() {

	for (auto root_hash_iter = root_hash.begin(); root_hash_iter != root_hash.end(); root_hash_iter++) {

		if (root_hash_iter->begin() != root_hash_iter->end()) {

			return iterator(this, root_hash_iter, root_hash_iter->begin());
		}					
	}

	return end();
}

template<typename ValueType, uchar key_size>
typename HybridHash<ValueType, key_size>::iterator HybridHash<ValueType, key_size>::end() {

	typename LeafHash::iterator return_dummy_iterator;
	return iterator(this, root_hash.end(), return_dummy_iterator);
}

template<typename ValueType, uchar key_size>
ulong HybridHash<ValueType, key_size>::rootHashIndex(const std::bitset<key_size> & key) {

	return key_hasher(key) % root_hash.size();
}

template<typename ValueType, uchar key_size>
std::pair<typename HybridHash<ValueType, key_size>::iterator, bool> HybridHash<ValueType, key_size>::insert(const std::bitset<key_size> & key, ValueType value, const bool add_sorted) {

	ulong root_hash_idx = rootHashIndex(key);

	auto leaf_hash_insert_result = root_hash.at(root_hash_idx).insert(std::make_pair(key, value), add_sorted);
	auto return_iterator = iterator(this, root_hash.begin() + root_hash_idx, leaf_hash_insert_result.first);

	return std::pair<iterator, bool> (return_iterator, leaf_hash_insert_result.second);
}

template<typename ValueType, uchar key_size>
typename HybridHash<ValueType, key_size>::iterator HybridHash<ValueType, key_size>::find(const std::bitset<key_size> & key) {

	ulong root_hash_idx = rootHashIndex(key);

	auto leaf_hash_find_result = root_hash.at(root_hash_idx).find(key);

	if (leaf_hash_find_result == root_hash.at(root_hash_idx).end()){

		return end();			
	
	} else {

		return iterator(this, root_hash.begin() + root_hash_idx, leaf_hash_find_result);
	}
}	

template<typename ValueType, uchar key_size>
void HybridHash<ValueType, key_size>::sort() {

	for (auto & leaf_hash: root_hash) {

		leaf_hash.sort();
	}
}

template<typename ValueType, uchar key_size>
void HybridHash<ValueType, key_size>::shuffle(const uint prng_seed) {

    mt19937 prng = mt19937(prng_seed);

	for (auto & leaf_hash: root_hash) {

		leaf_hash.shuffle(&prng);
	}    
}

template<typename ValueType, uchar key_size>
ulong HybridHash<ValueType, key_size>::size() {

	ulong _size = 0;

	for (auto & leaf_hash: root_hash) {

		_size += leaf_hash.size();
	}

	return _size;
}

template<typename ValueType, uchar key_size>
map<ulong, ulong> HybridHash<ValueType, key_size>::getRootSizeDistribution() {

	map<ulong, ulong> root_size_distribution;

	for (auto & leaf_hash: root_hash) {

		auto root_size_distribution_it = root_size_distribution.emplace(leaf_hash.size(), 0);
		root_size_distribution_it.first->second++;
	}

	return root_size_distribution;
}


template<typename ValueType, uchar key_size>
ThreadedHybridHash<ValueType, key_size>::ThreadedHybridHash(const uint root_hash_size, const ulong expected_total_size, const ushort num_threads_in) : HybridHash<ValueType, key_size>(root_hash_size, expected_total_size), num_threads(num_threads_in) {

	mutexes = std::vector<std::mutex>(HybridHash<ValueType, key_size>::root_hash.size());
}


template<typename ValueType, uchar key_size>
std::unique_lock<std::mutex> ThreadedHybridHash<ValueType, key_size>::lockKey(const std::bitset<key_size> & key) {

	std::unique_lock<std::mutex> current_lock(mutexes.at(HybridHash<ValueType, key_size>::rootHashIndex(key)));

	return std::move(current_lock);
}


template<typename ValueType, uchar key_size>
void ThreadedHybridHash<ValueType, key_size>::sortCallback(const ushort thread_idx) {

    ulong root_hash_index = thread_idx;

	while (root_hash_index < HybridHash<ValueType, key_size>::root_hash.size()) {

		HybridHash<ValueType, key_size>::root_hash.at(root_hash_index).sort();
        root_hash_index += num_threads;
	}
}

template<typename ValueType, uchar key_size>
void ThreadedHybridHash<ValueType, key_size>::sort() {

	std::vector<std::thread> sorting_threads;
	sorting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    sorting_threads.push_back(std::thread(&ThreadedHybridHash<ValueType, key_size>::sortCallback, this, thread_idx));
    }

    for (auto & sorting_thread: sorting_threads) {
        	
       	sorting_thread.join();
	}
}


/*
	Iterator for HybridHash
*/

template<typename ValueType, uchar key_size>
HybridHashIterator<ValueType, key_size>::HybridHashIterator(HybridHash<ValueType, key_size> * hybrid_hash_in, typename RootHash::iterator root_element_iter_in, typename LeafHash::iterator leaf_element_iter_in) {

	hybrid_hash = hybrid_hash_in;
	root_element_iter = root_element_iter_in;
	leaf_element_iter = leaf_element_iter_in;
}

template<typename ValueType, uchar key_size>
HybridHashIterator<ValueType, key_size>::HybridHashIterator(const HybridHashIterator<ValueType, key_size> & other) {

	hybrid_hash = other.hybrid_hash;
	root_element_iter = other.root_element_iter;
	leaf_element_iter = other.leaf_element_iter;
}

template<typename ValueType, uchar key_size>
HybridHashIterator<ValueType, key_size> & HybridHashIterator<ValueType, key_size>::operator=(HybridHashIterator<ValueType, key_size> other) {

    swap(*this, other);

    return *this;
}

template<typename ValueType, uchar key_size>
void HybridHashIterator<ValueType, key_size>::swap(HybridHashIterator<ValueType, key_size> & first, HybridHashIterator<ValueType, key_size> & second) {

    std::swap(first.hybrid_hash, second.hybrid_hash);
    std::swap(first.root_element_iter, second.root_element_iter);
    std::swap(first.leaf_element_iter, second.leaf_element_iter);
}

template<typename ValueType, uchar key_size>
std::pair<std::bitset<key_size>, ValueType> & HybridHashIterator<ValueType, key_size>::operator*() {

	assert(root_element_iter != hybrid_hash->root_hash.end());
	return *leaf_element_iter; 
}

template<typename ValueType, uchar key_size>
HybridHashIterator<ValueType, key_size> & HybridHashIterator<ValueType, key_size>::operator++() {

	leaf_element_iter++;

	if (leaf_element_iter == root_element_iter->end()) {

		root_element_iter++;

		while (root_element_iter != hybrid_hash->root_hash.end()) {

			leaf_element_iter = root_element_iter->begin();

			if (leaf_element_iter != root_element_iter->end()) {

				break;
			}			
	
			root_element_iter++;
		}
	}	

	return *this;
}

template<typename ValueType, uchar key_size>
HybridHashIterator<ValueType, key_size> HybridHashIterator<ValueType, key_size>::operator++(int) {

	HybridHashIterator<ValueType, key_size> iterator_copy(*this);
	++(*this);

	return iterator_copy;
}

template<typename ValueType, uchar key_size>
bool HybridHashIterator<ValueType, key_size>::operator==(HybridHashIterator<ValueType, key_size> const & other) {

	if (root_element_iter == other.root_element_iter) {

		if (root_element_iter == hybrid_hash->root_hash.end()) {

			return true;
		
		} else {

			return (leaf_element_iter == other.leaf_element_iter);
		}
	}

	return false;
} 

template<typename ValueType, uchar key_size>
bool HybridHashIterator<ValueType, key_size>::operator!=(HybridHashIterator<ValueType, key_size> const & other) {

	return !(*this == other);
}
