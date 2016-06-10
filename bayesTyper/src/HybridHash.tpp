#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <map>

static const uchar root_split_size = 14;

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHash<ValueType, root_hash_size, leaf_hash_size>::HybridHash() {

	assert(root_split_size <= root_hash_size);
	assert((root_hash_size + leaf_hash_size) <= 128);
	assert(root_hash_size <= 32);

	root_hash_bit_values.reserve(root_hash_size);

	for (auto i = 0; i < root_hash_size; i++) {

		root_hash_bit_values.push_back(pow(2, i));
	}

	root_hash = std::vector<LeafHash*>((ulong)pow(2, root_hash_size));

	for (auto & hash : root_hash) {

		hash = nullptr;
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHash<ValueType, root_hash_size, leaf_hash_size>::~HybridHash() {

	for (auto & hash : root_hash) {

		if (hash) {

			delete hash;
		}
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator HybridHash<ValueType, root_hash_size, leaf_hash_size>::begin() {

	for (auto root_hash_iter = root_hash.begin(); root_hash_iter != root_hash.end(); root_hash_iter++) {

		if (*root_hash_iter) {

			return iterator(this, root_hash_iter, (*root_hash_iter)->begin());
		}					
	}

	return end();
 }


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator HybridHash<ValueType, root_hash_size, leaf_hash_size>::end() {

	typename LeafHash::iterator return_dummy_iterator;
	return iterator(this, root_hash.end(), return_dummy_iterator);
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ulong HybridHash<ValueType, root_hash_size, leaf_hash_size>::bitsubsetToUlong(const std::bitset<root_hash_size + leaf_hash_size> & key, uchar start_idx, uchar end_idx) {

	ulong hash_value = 0;
	uchar bit = 0;

	while (start_idx <= end_idx) {

		if (key[start_idx]) {

			hash_value += root_hash_bit_values.at(bit);
		}

		bit++;
		start_idx++;
	}

	return hash_value;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ulong HybridHash<ValueType, root_hash_size, leaf_hash_size>::bitsubsetToUlongReverse(const std::bitset<root_hash_size + leaf_hash_size> & key, uchar end_idx, uchar start_idx) {

	ulong hash_value = 0;
	uchar bit = 0;

	while (start_idx >= end_idx) {

		if (key[start_idx]) {

			hash_value += root_hash_bit_values.at(bit);
		}

		bit++;
		start_idx--;
	}

	return hash_value;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ulong HybridHash<ValueType, root_hash_size, leaf_hash_size>::bitsubsetToUlongSplit(const std::bitset<root_hash_size + leaf_hash_size> & key, uchar start_first_idx, uchar end_first_idx, uchar start_second_idx, uchar end_second_idx) {

	ulong hash_value = 0;
	uchar bit = 0;

	while (start_first_idx <= end_first_idx) {

		if (key[start_first_idx]) {

			hash_value += root_hash_bit_values.at(bit);
		}

		bit++;
		start_first_idx++;
	}

	while (start_second_idx <= end_second_idx) {

		if (key[start_second_idx]) {

			hash_value += root_hash_bit_values.at(bit);
		}

		bit++;
		start_second_idx++;
	}

	return hash_value;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::bitset<leaf_hash_size> HybridHash<ValueType, root_hash_size, leaf_hash_size>::bitsetToLeafBitsubset(const std::bitset<root_hash_size + leaf_hash_size> & key, uchar start_idx, uchar end_idx) {

	std::bitset<leaf_hash_size> leaf_bitset;
	uchar bit = 0;

	while (start_idx <= end_idx) {

		leaf_bitset[bit] = key[start_idx];

		bit++;
		start_idx++;
	}

	return leaf_bitset;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::pair<typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator, bool> HybridHash<ValueType, root_hash_size, leaf_hash_size>::insert(const std::bitset<root_hash_size + leaf_hash_size> & key, ValueType value) {

	assert(key.size() == (root_hash_size + leaf_hash_size));

	ulong root_hash_value = bitsubsetToUlongSplit(key, 0, root_split_size - 1, leaf_hash_size + root_split_size, root_hash_size + leaf_hash_size - 1); 
	std::bitset<leaf_hash_size> leaf_bitsubset = bitsetToLeafBitsubset(key, root_split_size, leaf_hash_size + root_split_size - 1);

	if (root_hash.at(root_hash_value)) {

		auto leaf_hash_insert_result = root_hash.at(root_hash_value)->insert(std::make_pair(leaf_bitsubset, value));
		auto return_iterator = iterator(this, root_hash.begin() + root_hash_value, leaf_hash_insert_result.first);

		if (leaf_hash_insert_result.second) {

			return std::pair<iterator,bool> (return_iterator, true);
		
		} else {

			return std::pair<iterator,bool> (return_iterator, false);
		}

	} else {

		LeafHash * new_leaf_hash = new LeafHash;
		root_hash.at(root_hash_value) = new_leaf_hash;
		auto leaf_hash_insert_result = new_leaf_hash->insert(std::pair<std::bitset<leaf_hash_size>, ValueType>(leaf_bitsubset, value));
		assert(leaf_hash_insert_result.second);

		return std::pair<iterator,bool> (iterator(this, root_hash.begin() + root_hash_value, leaf_hash_insert_result.first), true);
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator HybridHash<ValueType, root_hash_size, leaf_hash_size>::find(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	assert(key.size() == (root_hash_size + leaf_hash_size));

	ulong root_hash_value = bitsubsetToUlongSplit(key, 0, root_split_size - 1, leaf_hash_size + root_split_size, root_hash_size + leaf_hash_size - 1);

	if (root_hash.at(root_hash_value)) {

		std::bitset<leaf_hash_size> leaf_bitsubset = bitsetToLeafBitsubset(key, root_split_size, leaf_hash_size + root_split_size - 1);

		LeafHash * current_leaf_hash = root_hash.at(root_hash_value);

		auto leaf_hash_find_result = current_leaf_hash->find(leaf_bitsubset);

		if (leaf_hash_find_result == current_leaf_hash->end()){

			return end();			
		
		} else {

			return iterator(this, root_hash.begin() + root_hash_value, leaf_hash_find_result);
		}

	} else {

		return end();
	}
}	


// template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
// typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator HybridHash<ValueType, root_hash_size, leaf_hash_size>::erase(HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator erase_iterator) {

// 	assert(erase_iterator.hybrid_hash == this);

// 	iterator erase_iterator_copy(erase_iterator);
	
// 	erase_iterator++;

// 	(*erase_iterator_copy.root_element_iter)->erase(erase_iterator_copy.leaf_element_iter);

// 	if ((*erase_iterator_copy.root_element_iter)->empty()) {

// 		delete *erase_iterator_copy.root_element_iter;
// 		*erase_iterator_copy.root_element_iter = nullptr;
// 	}

// 	return erase_iterator;
// }


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
void HybridHash<ValueType, root_hash_size, leaf_hash_size>::writeBucketSizeDistribution(std::string filename) {

	std::map<uint, uint> size_distribution;

	for (auto &hit: root_hash) {

		if (hit) {

			assert(!hit->empty());
			auto size_distribution_insert = size_distribution.insert({hit->size(), 0});
			size_distribution_insert.first->second++;
		
		} else {

			auto size_distribution_insert = size_distribution.insert({0, 0});
			size_distribution_insert.first->second++;
		}
	}

	std::ofstream size_distribution_file(filename);
    assert(size_distribution_file.is_open());

    size_distribution_file << "Size\tCount" << std::endl;
    
    for (auto &cit: size_distribution) {

        size_distribution_file << cit.first << "\t" << cit.second  << std::endl;
    } 

    size_distribution_file.close();
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::ThreadedHybridHash(uint num_lock_bins_in) : HybridHash<ValueType,root_hash_size,leaf_hash_size>() {

	ulong root_bin_size = pow(2, root_hash_size);

	if (root_bin_size > num_lock_bins_in) {

		num_lock_bins = num_lock_bins_in;
		hash_scaling_factor = (double) num_lock_bins/root_bin_size;

	} else {

		num_lock_bins = root_bin_size;
		hash_scaling_factor = 1;
	}

	mutexes = std::vector<std::mutex> (num_lock_bins);
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::mutex * ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::getKeyMutex(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	ulong root_hash_value = HybridHash<ValueType,root_hash_size,leaf_hash_size>::bitsubsetToUlongSplit(key, 0, root_split_size - 1, leaf_hash_size + root_split_size, root_hash_size + leaf_hash_size - 1);
	uint mutex_bin_idx = floor(root_hash_value * hash_scaling_factor);
		
	return &mutexes.at(mutex_bin_idx);
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::pair<std::unique_lock<std::mutex>, std::unique_lock<std::mutex> > ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::lockKeyPair(const std::bitset<root_hash_size + leaf_hash_size> & key_1, const std::bitset<root_hash_size + leaf_hash_size> & key_2) {

	std::mutex * current_mutex_1 = getKeyMutex(key_1);
	std::mutex * current_mutex_2 = getKeyMutex(key_2);

	if (current_mutex_1 == current_mutex_2) {
		
		std::unique_lock<std::mutex> current_lock_1(*current_mutex_1);
		std::unique_lock<std::mutex> current_lock_2(*current_mutex_1, std::defer_lock);
		
		return std::pair<std::unique_lock<std::mutex>, std::unique_lock<std::mutex> >(std::move(current_lock_1), std::move(current_lock_2));

	} else {

		std::unique_lock<std::mutex> current_lock_1(*current_mutex_1, std::defer_lock);
		std::unique_lock<std::mutex> current_lock_2(*current_mutex_2, std::defer_lock);
		std::lock(current_lock_1, current_lock_2);

		return std::pair<std::unique_lock<std::mutex>, std::unique_lock<std::mutex> >(std::move(current_lock_1), std::move(current_lock_2));
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::unique_lock<std::mutex> ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::lockKey(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	std::mutex * current_mutex = getKeyMutex(key);
	std::unique_lock<std::mutex> current_lock(*current_mutex);

	return std::move(current_lock);
}

/*
	Iterator for HybridHash
*/

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::HybridHashIterator(HybridHash<ValueType,root_hash_size,leaf_hash_size> * hybrid_hash_in, typename RootHash::iterator root_element_iter_in, typename LeafHash::iterator leaf_element_iter_in) {

	hybrid_hash = hybrid_hash_in;
	root_element_iter = root_element_iter_in;
	leaf_element_iter = leaf_element_iter_in;
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::HybridHashIterator(const HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> & other) {

	hybrid_hash = other.hybrid_hash;
	root_element_iter = other.root_element_iter;
	leaf_element_iter = other.leaf_element_iter;
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> & HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator=(HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> other) {

    swap(*this, other);

    return *this;
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
void HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::swap(HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> & first, HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> & second) {

    std::swap(first.hybrid_hash, second.hybrid_hash);
    std::swap(first.root_element_iter, second.root_element_iter);
    std::swap(first.leaf_element_iter, second.leaf_element_iter);
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ValueType & HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator*() {

	// Sanity check - to be improve
	assert(root_element_iter != hybrid_hash->root_hash.end());

	return leaf_element_iter->second; 
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> & HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator++() {

	leaf_element_iter++;

	if (leaf_element_iter == (*root_element_iter)->end()) {

		root_element_iter++;

		while (root_element_iter != hybrid_hash->root_hash.end()) {

			if (*root_element_iter) {

				leaf_element_iter = (*root_element_iter)->begin();
				break;
			}			
	
			root_element_iter++;
		}
	}	

	return *this;
}

template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator++(int) {

	HybridHashIterator<ValueType,root_hash_size,leaf_hash_size> iterator_copy(*this);
	++(*this);

	return iterator_copy;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
bool HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator==(HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> const & other) {

	if (root_element_iter == other.root_element_iter) {

		if (root_element_iter == hybrid_hash->root_hash.end()) {

			return true;
		
		} else {

			return (leaf_element_iter == other.leaf_element_iter);
		}
	
	}

	return false;
} 


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
bool HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator!=(HybridHashIterator<ValueType, root_hash_size, leaf_hash_size> const & other) {

	return !(*this == other);
}
