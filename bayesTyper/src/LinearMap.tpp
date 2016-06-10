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
std::pair<typename LinearMap<KeyType, ValueType, KeyComparator>::iterator, bool> LinearMap<KeyType, ValueType, KeyComparator>::insert(std::pair<KeyType, ValueType> key_value_pair) {

	auto lower_bound_iter = std::lower_bound(map.begin(), map.end(), key_value_pair.first, PairComparator());

	// Insert at end
	if (lower_bound_iter == map.end()) {

		map.emplace_back(key_value_pair);
		return std::make_pair(--map.end(), true);
	
	} else {

		if (lower_bound_iter->first == key_value_pair.first) {

			return std::make_pair(lower_bound_iter, false);

		} else {

			return std::make_pair(map.insert(lower_bound_iter, key_value_pair), true);
		}
	}
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator LinearMap<KeyType, ValueType, KeyComparator>::find(KeyType key) {

	auto equal_range_iters = std::equal_range(map.begin(), map.end(), key, PairComparator());

	if (equal_range_iters.first == equal_range_iters.second) {

		return map.end();

	} else {

		return equal_range_iters.first;
	}
}

template<typename KeyType, typename ValueType, typename KeyComparator>
typename LinearMap<KeyType, ValueType, KeyComparator>::iterator LinearMap<KeyType, ValueType, KeyComparator>::erase(LinearMap<KeyType, ValueType, KeyComparator>::iterator erase_iter) {
	
	return map.erase(erase_iter);
}