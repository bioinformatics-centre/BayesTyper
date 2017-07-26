#include <iostream>
#include <math.h>

template<uchar bitset_size>
PerfectSet<bitset_size>::PerfectSet() {

	_set = std::vector<bool>(static_cast<ulong>(pow(2, bitset_size)), false);
}


template<uchar bitset_size>
PerfectSet<bitset_size>::PerfectSet(const ulong & set_size) {
	
	_set = std::vector<bool>(set_size, false);
}


template<uchar bitset_size>
bool PerfectSet<bitset_size>::insert(const std::bitset<bitset_size> & key) {

	return insert(key.to_ullong());
}


template<uchar bitset_size>
bool PerfectSet<bitset_size>::insert(const ulong & idx) {

	if (_set.at(idx)) {

		return false;
	
	} else {

		_set.at(idx) = true;
		return true;
	}
}


template<uchar bitset_size>
bool PerfectSet<bitset_size>::count(const std::bitset<bitset_size> & key) {

	return count(key.to_ullong());
}


template<uchar bitset_size>
bool PerfectSet<bitset_size>::count(const ulong & idx) {

	return _set.at(idx);
}


template<uchar bitset_size>
ThreadedPerfectSet<bitset_size>::ThreadedPerfectSet(const ulong & num_lock_bins) {

	ulong set_size = static_cast<ulong>(pow(2, bitset_size));
	_set_bins.reserve(num_lock_bins);

	max_bin_size = set_size / num_lock_bins;

	if ((set_size % num_lock_bins) > 0) {

		max_bin_size++;
	}

	assert(max_bin_size > 0);

	ulong cum_num_values = 0;

	while (cum_num_values + max_bin_size < set_size) {

		_set_bins.emplace_back(max_bin_size);
		cum_num_values += max_bin_size;
	}

	assert(set_size > cum_num_values);
	assert((set_size - cum_num_values) <= max_bin_size);
	
	_set_bins.emplace_back(set_size - cum_num_values);
	_set_bins.shrink_to_fit();

	_mutex_bins = std::vector<std::mutex>(_set_bins.size());
}


template<uchar bitset_size>
bool ThreadedPerfectSet<bitset_size>::insert(const std::bitset<bitset_size> & key) {

	ulong idx = key.to_ullong();
	std::unique_lock<std::mutex> current_lock(_mutex_bins.at(idx / max_bin_size));

	return _set_bins.at(idx / max_bin_size).insert(idx % max_bin_size);
}


template<uchar bitset_size>
bool ThreadedPerfectSet<bitset_size>::count(const std::bitset<bitset_size> & key) {

	ulong idx = key.to_ullong();

	return _set_bins.at(idx / max_bin_size).count(idx % max_bin_size);
}

