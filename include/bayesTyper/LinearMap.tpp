
/*
LinearMap.tpp - This file is part of BayesTyper (v1.1)


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


template<typename KeyType, typename ValueType, typename KeyComparator>
LinearMap<KeyType, ValueType, KeyComparator>::LinearMap() {}

template<typename KeyType, typename ValueType, typename KeyComparator>
void LinearMap<KeyType, ValueType, KeyComparator>::reserve(const uint size) {

	map.reserve(size);
}

template<typename KeyType, typename ValueType, typename KeyComparator>
bool LinearMap<KeyType, ValueType, KeyComparator>::empty() {

	return map.empty();
}

template<typename KeyType, typename ValueType, typename KeyComparator>
uint LinearMap<KeyType, ValueType, KeyComparator>::size() {

	return map.size();
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator LinearMap<KeyType, ValueType, KeyComparator>::begin() {

	return map.begin();
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator LinearMap<KeyType, ValueType, KeyComparator>::end() {

	return map.end();
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator LinearMap<KeyType, ValueType, KeyComparator>::erase(LinearMap<KeyType, ValueType, KeyComparator>::iterator erase_iter) {
	
	return map.erase(erase_iter);
}


template<typename KeyType, typename ValueType, typename KeyComparator>
SortedLinearMap<KeyType, ValueType, KeyComparator>::SortedLinearMap() : LinearMap<KeyType, ValueType, KeyComparator>() {}

template<typename KeyType, typename ValueType, typename KeyComparator>
std::pair<typename LinearMap<KeyType, ValueType, KeyComparator>::iterator, bool> SortedLinearMap<KeyType, ValueType, KeyComparator>::insert(std::pair<KeyType, ValueType> key_value_pair) {

	auto lower_bound_iter = std::lower_bound(LinearMap<KeyType, ValueType, KeyComparator>::map.begin(), LinearMap<KeyType, ValueType, KeyComparator>::map.end(), key_value_pair.first, typename LinearMap<KeyType, ValueType, KeyComparator>::PairComparator());

	if (lower_bound_iter == LinearMap<KeyType, ValueType, KeyComparator>::map.end()) {

		LinearMap<KeyType, ValueType, KeyComparator>::map.emplace_back(key_value_pair);
		return std::make_pair(--LinearMap<KeyType, ValueType, KeyComparator>::map.end(), true);
	
	} else {

		if (lower_bound_iter->first == key_value_pair.first) {

			return std::make_pair(lower_bound_iter, false);

		} else {

			return std::make_pair(LinearMap<KeyType, ValueType, KeyComparator>::map.insert(lower_bound_iter, key_value_pair), true);
		}
	}
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator SortedLinearMap<KeyType, ValueType, KeyComparator>::find(KeyType key) {

	auto lower_bound_iter = std::lower_bound(LinearMap<KeyType, ValueType, KeyComparator>::map.begin(), LinearMap<KeyType, ValueType, KeyComparator>::map.end(), key, typename LinearMap<KeyType, ValueType, KeyComparator>::PairComparator());

	if (lower_bound_iter == LinearMap<KeyType, ValueType, KeyComparator>::map.end()) {

		return LinearMap<KeyType, ValueType, KeyComparator>::map.end();

	} else if (lower_bound_iter->first == key) {

		return lower_bound_iter;

	} else {

		return LinearMap<KeyType, ValueType, KeyComparator>::map.end();
	}
}


template<typename KeyType, typename ValueType, typename KeyComparator>
PartialSortedLinearMap<KeyType, ValueType, KeyComparator>::PartialSortedLinearMap() : LinearMap<KeyType, ValueType, KeyComparator>() {

	sorted_segment_length = 0;
}

template<typename KeyType, typename ValueType, typename KeyComparator>
std::pair<typename LinearMap<KeyType, ValueType, KeyComparator>::iterator, bool> PartialSortedLinearMap<KeyType, ValueType, KeyComparator>::insert(std::pair<KeyType, ValueType> key_value_pair, const bool add_sorted) {

	assert(sorted_segment_length <= (LinearMap<KeyType, ValueType, KeyComparator>::map.size()));
	
	auto sorted_segment_end_it = LinearMap<KeyType, ValueType, KeyComparator>::map.begin();
	std::advance(sorted_segment_end_it, sorted_segment_length);

	auto lower_bound_iter = std::lower_bound(LinearMap<KeyType, ValueType, KeyComparator>::map.begin(), sorted_segment_end_it, key_value_pair.first, typename LinearMap<KeyType, ValueType, KeyComparator>::PairComparator());	
	
	if (lower_bound_iter == LinearMap<KeyType, ValueType, KeyComparator>::map.end()) {

		assert(sorted_segment_length == (LinearMap<KeyType, ValueType, KeyComparator>::map.size()));
		sorted_segment_length += 1;

		return std::make_pair(LinearMap<KeyType, ValueType, KeyComparator>::map.insert(lower_bound_iter, key_value_pair), true);
	
	} else {

		if (lower_bound_iter->first == key_value_pair.first) {

			return std::make_pair(lower_bound_iter, false);

		} else if (add_sorted) {
			
			sorted_segment_length += 1;

			return std::make_pair(LinearMap<KeyType, ValueType, KeyComparator>::map.insert(lower_bound_iter, key_value_pair), true);

		} else {

			return std::make_pair(LinearMap<KeyType, ValueType, KeyComparator>::map.insert(LinearMap<KeyType, ValueType, KeyComparator>::map.end(), key_value_pair), true);
		}
	}
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator PartialSortedLinearMap<KeyType, ValueType, KeyComparator>::find(KeyType key) {

	assert(sorted_segment_length <= (LinearMap<KeyType, ValueType, KeyComparator>::map.size()));

	auto sorted_segment_end_it = LinearMap<KeyType, ValueType, KeyComparator>::map.begin();
	std::advance(sorted_segment_end_it, sorted_segment_length);

	auto lower_bound_iter = std::lower_bound(LinearMap<KeyType, ValueType, KeyComparator>::map.begin(), sorted_segment_end_it, key, typename LinearMap<KeyType, ValueType, KeyComparator>::PairComparator());

	if (lower_bound_iter == LinearMap<KeyType, ValueType, KeyComparator>::map.end()) {

		assert(sorted_segment_length == (LinearMap<KeyType, ValueType, KeyComparator>::map.size()));

		return lower_bound_iter;

	} else if (lower_bound_iter->first == key) {

		return lower_bound_iter;

	} else {

		auto map_unsorted_it = sorted_segment_end_it;

		while (map_unsorted_it != LinearMap<KeyType, ValueType, KeyComparator>::map.end()) {

			if (map_unsorted_it->first == key) {

				return map_unsorted_it;
			}

			map_unsorted_it++;
		}

		return LinearMap<KeyType, ValueType, KeyComparator>::map.end();
	}
}

template<typename KeyType, typename ValueType, typename KeyComparator>
void PartialSortedLinearMap<KeyType, ValueType, KeyComparator>::sort() {

	std::sort(LinearMap<KeyType, ValueType, KeyComparator>::map.begin(), LinearMap<KeyType, ValueType, KeyComparator>::map.end(), typename LinearMap<KeyType, ValueType, KeyComparator>::PairComparator());
	sorted_segment_length = LinearMap<KeyType, ValueType, KeyComparator>::map.size();
}


