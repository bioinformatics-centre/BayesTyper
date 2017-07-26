#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <map>
#include <thread>


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
HybridHash<ValueType, root_hash_size, leaf_hash_size>::HybridHash() {

	assert(root_hash_size > 0);
	assert(leaf_hash_size > 0);

	assert(root_hash_size <= 32);
	assert((root_hash_size % 2) == 0);

	root_hash = std::vector<LeafHash*>(static_cast<ulong>(pow(2, root_hash_size)));

	for (auto & hash : root_hash) {

		hash = nullptr;
	}

	root_hash_bit_values.reserve(root_hash_size);
	
	ushort root_hash_key_index_shift = std::floor(static_cast<float>(root_hash_size + leaf_hash_size)/root_hash_size);
	assert(root_hash_key_index_shift > 0);

	for (ushort i = 0; i < (root_hash_size / 2); i++) {

		root_hash_bit_values.emplace_back(i * root_hash_key_index_shift, pow(2, i * 2));
		root_hash_bit_values.emplace_back(root_hash_size + leaf_hash_size - i * root_hash_key_index_shift - 1, pow(2, i * 2 + 1));
	}

	assert(root_hash_bit_values.size() == root_hash_size);
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
ulong HybridHash<ValueType, root_hash_size, leaf_hash_size>::rootHashIndex(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	ulong root_hash_idx = 0;

	for (auto & root_hash_bit_value: root_hash_bit_values) {

		if (key[root_hash_bit_value.first]) {

			root_hash_idx += root_hash_bit_value.second;
		}
	}

	return root_hash_idx;
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::pair<typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator, bool> HybridHash<ValueType, root_hash_size, leaf_hash_size>::insert(const std::bitset<root_hash_size + leaf_hash_size> & key, ValueType value, const bool add_sorted) {

	ulong root_hash_idx = rootHashIndex(key);

	if (root_hash.at(root_hash_idx)) {

		auto leaf_hash_insert_result = root_hash.at(root_hash_idx)->insert(std::make_pair(key, value), add_sorted);

		auto return_iterator = iterator(this, root_hash.begin() + root_hash_idx, leaf_hash_insert_result.first);

		if (leaf_hash_insert_result.second) {

			return std::pair<iterator,bool> (return_iterator, true);
		
		} else {

			return std::pair<iterator,bool> (return_iterator, false);
		}

	} else {

		LeafHash * new_leaf_hash = new LeafHash();
		root_hash.at(root_hash_idx) = new_leaf_hash;

		auto leaf_hash_insert_result = new_leaf_hash->insert(std::make_pair(key, value), add_sorted);

		assert(leaf_hash_insert_result.second);

		return std::pair<iterator,bool> (iterator(this, root_hash.begin() + root_hash_idx, leaf_hash_insert_result.first), true);
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
typename HybridHash<ValueType, root_hash_size, leaf_hash_size>::iterator HybridHash<ValueType, root_hash_size, leaf_hash_size>::find(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	ulong root_hash_idx = rootHashIndex(key);

	if (root_hash.at(root_hash_idx)) {

		auto leaf_hash_find_result = root_hash.at(root_hash_idx)->find(key);

		if (leaf_hash_find_result == root_hash.at(root_hash_idx)->end()){

			return end();			
		
		} else {

			return iterator(this, root_hash.begin() + root_hash_idx, leaf_hash_find_result);
		}

	} else {

		return end();
	}
}	


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
void HybridHash<ValueType, root_hash_size, leaf_hash_size>::sort() {

	for (auto & leaf_hash: root_hash) {

		if (leaf_hash) {

			assert(!(leaf_hash->empty()));
			leaf_hash->sort();
		} 
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::ThreadedHybridHash(const ushort num_threads_in) : HybridHash<ValueType, root_hash_size, leaf_hash_size>(), num_threads(num_threads_in) {

	mutexes = std::vector<std::mutex>(HybridHash<ValueType, root_hash_size, leaf_hash_size>::root_hash.size());
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
std::unique_lock<std::mutex> ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::lockKey(const std::bitset<root_hash_size + leaf_hash_size> & key) {

	std::unique_lock<std::mutex> current_lock(mutexes.at(HybridHash<ValueType, root_hash_size, leaf_hash_size>::rootHashIndex(key)));

	return std::move(current_lock);
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
void ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::sortCallback(const ushort thread_index) {

    uint root_hash_index = thread_index;

	while (root_hash_index < HybridHash<ValueType, root_hash_size, leaf_hash_size>::root_hash.size()) {

		if (HybridHash<ValueType, root_hash_size, leaf_hash_size>::root_hash.at(root_hash_index)) {

			assert(!(HybridHash<ValueType, root_hash_size, leaf_hash_size>::root_hash.at(root_hash_index)->empty()));
			HybridHash<ValueType, root_hash_size, leaf_hash_size>::root_hash.at(root_hash_index)->sort();
		} 

        root_hash_index += num_threads;
	}
}


template<typename ValueType, uchar root_hash_size, uchar leaf_hash_size>
void ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::sort() {

	std::vector<std::thread> sorting_threads;
	sorting_threads.reserve(num_threads);

	for (uint i = 0; i < num_threads; i++) {

   	    sorting_threads.push_back(std::thread(&ThreadedHybridHash<ValueType, root_hash_size, leaf_hash_size>::sortCallback, this, i));
    }

    for (auto & sorting_thread: sorting_threads) {
        	
       	sorting_thread.join();
	}
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
std::pair<std::bitset<root_hash_size + leaf_hash_size>, ValueType> & HybridHashIterator<ValueType, root_hash_size, leaf_hash_size>::operator*() {

	assert(root_element_iter != hybrid_hash->root_hash.end());
	return *leaf_element_iter; 
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
