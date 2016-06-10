
/*
CountAllocation.cpp - This file is part of BayesTyper (v0.9)


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


#include "CountAllocation.hpp"

#include "Utils.hpp"

CountAllocation::CountAllocation(const uchar num_noise_sources_in, const ushort num_samples_in) : num_noise_sources(num_noise_sources_in), num_samples(num_samples_in) {

	noise_counts = vector<vector<unordered_map<uchar,ulong> > > (num_samples, vector<unordered_map<uchar,ulong> >(num_noise_sources));
}

const vector<vector<unordered_map<uchar,ulong> > > & CountAllocation::noiseCounts() const {

	return noise_counts;
}

void CountAllocation::mergeCountContainers(unordered_map<uchar,ulong> * merge_recipent_counts, const unordered_map<uchar,ulong> & merge_donor_counts) {

	for (auto merge_donor_count : merge_donor_counts) {
	
		auto merge_recipent_counts_insert_result = merge_recipent_counts->insert(merge_donor_count);

		if (!merge_recipent_counts_insert_result.second) {

			assert(merge_recipent_counts_insert_result.first->second > 0);
			merge_recipent_counts_insert_result.first->second += merge_donor_count.second;
		}	
	}
}

void CountAllocation::addNoiseCounts(const vector<vector<unordered_map<uchar,ulong> > > & noise_counts_in) {

	assert(noise_counts_in.size() == num_samples);

	for (uint i = 0; i < num_samples; i++) {

		assert(noise_counts_in.at(i).size() == num_noise_sources);

		for (uint j = 0; j < num_noise_sources; j++) {

			mergeCountContainers(&noise_counts.at(i).at(j), noise_counts_in.at(i).at(j));
		}
	}
}

void CountAllocation::mergeInCountAllocations(const CountAllocation & count_allocation_in) {

	assert(count_allocation_in.noise_counts.size() == num_samples);
	
	for (uint i = 0; i < count_allocation_in.noiseCounts().size(); i++) {

		assert(count_allocation_in.noiseCounts().at(i).size() == num_noise_sources);

		for (uint j = 0; j < count_allocation_in.noiseCounts().at(i).size(); j++) {

			for (auto noise_counts_local_iter = count_allocation_in.noiseCounts().at(i).at(j).begin(); noise_counts_local_iter != count_allocation_in.noiseCounts().at(i).at(j).end(); noise_counts_local_iter++) {

				auto count_emplace_result = noise_counts.at(i).at(j).emplace(noise_counts_local_iter->first, 0);
				count_emplace_result.first->second += noise_counts_local_iter->second;
			}
		}
	}
} 

