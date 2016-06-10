
/*
KmerCounter.cpp - This file is part of BayesTyper (v0.9)


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
#include <fstream>
#include <thread>
#include <mutex>
#include <limits>
#include <assert.h>
#include <algorithm>

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "kmc_file.h"

#include "KmerCounter.hpp"
#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Kmer.hpp"
#include "KmerCounts.hpp"
#include "Sample.hpp"

static const uint precounted_kmer_batch_size = 1000000;

enum class FastqLineType : uchar {header = 0, sequence, seperator, quality};

template <uchar kmer_size>
KmerCounter<kmer_size>::KmerCounter(KmerHash * const kmer_hash_in, Utils::SmallmerSet * const smallmer_set_in, const ushort num_samples_in, const ushort num_threads_in, const uint prng_seed_in) : smallmer_set(smallmer_set_in), num_samples(num_samples_in), num_threads(num_threads_in), prng_seed(prng_seed_in) {

	kmer_hash_hybrid = static_cast<KmerHashHybrid<kmer_size> * >(kmer_hash_in);
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterSmallmersCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const ushort thread_index, mutex * counting_mutex, ulong * num_unique_smallmers) {

    ulong local_num_unique_smallmers = 0;
    uint variant_cluster_group_index = thread_index;

	while (variant_cluster_group_index < variant_cluster_groups->size()) {

		local_num_unique_smallmers += variant_cluster_groups->at(variant_cluster_group_index)->countSmallmers(smallmer_set);
        variant_cluster_group_index += num_threads;
	} 

    unique_lock<mutex> counting_lock(*counting_mutex);
    *num_unique_smallmers += local_num_unique_smallmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countVariantClusterSmallmers(vector<VariantClusterGroup *> * variant_cluster_groups) {

    ulong num_unique_smallmers = 0;

	vector<thread> counting_threads;
    mutex counting_mutex;

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countVariantClusterSmallmersCallback, this, variant_cluster_groups, i, &counting_mutex, &num_unique_smallmers));
    }

    for (auto & counting_thread: counting_threads) {
        	
       	counting_thread.join();
	}

    return num_unique_smallmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countInterclusterKmersCallback(vector<pair<string::const_iterator, string::const_iterator> >::iterator current_interval, vector<pair<string::const_iterator, string::const_iterator> >::iterator end_interval, const Utils::ChromosomeClass chromosome_class, mutex * counting_mutex, ulong * num_intercluster_kmers) {

    ulong local_num_intercluster_kmers = 0;

    KmerPair<kmer_size> kmer_pair = KmerPair<kmer_size>();
    KmerForward<Utils::small_kmer_size> small_mer = KmerForward<Utils::small_kmer_size>();

    while (current_interval != end_interval) {

    	auto sit = current_interval->first;
    	auto eit = current_interval->second;

    	kmer_pair.reset();
    	small_mer.reset();

		while (sit != eit) {

			auto kmer_pair_move = kmer_pair.move(*sit);
			auto small_mer_move = small_mer.move(*sit);

			if (small_mer_move.second) {

				if (smallmer_set->count(small_mer.getKmer())) {

					if (kmer_pair_move.second) {

						auto hash_lock = kmer_hash_hybrid->_hash->lockKey(kmer_pair.getLexicographicalLowestKmer());
						auto intercluster_kmer_insert = kmer_hash_hybrid->_hash->insert(kmer_pair.getLexicographicalLowestKmer(), nullptr);
					
						if (intercluster_kmer_insert.second) {

							*intercluster_kmer_insert.first = new SampleKmerCounts(num_samples);				
						} 

						local_num_intercluster_kmers++;

						static_cast<SampleKmerCounts*>(*intercluster_kmer_insert.first)->addInterclusterMultiplicity(chromosome_class);
					}

				} else {
					
					kmer_pair.shrink(Utils::small_kmer_size);
				}
			}

			sit++;
		} 

		current_interval++;
	} 

	lock_guard<mutex> counting_lock(*counting_mutex);
    *num_intercluster_kmers += local_num_intercluster_kmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countInterclusterKmers(vector<pair<string::const_iterator, string::const_iterator> > & intercluster_intervals, const Utils::ChromosomeClass chromosome_class) {

    ulong num_intercluster_kmers = 0;

	auto thread_intercluster_allocations = Utils::allocateToThreads(intercluster_intervals, num_threads);

	vector<thread> counting_threads;
	mutex counting_mutex;

	for (uint i = 0; i < thread_intercluster_allocations.size(); i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countInterclusterKmersCallback, this, thread_intercluster_allocations.at(i).first, thread_intercluster_allocations.at(i).second, chromosome_class, &counting_mutex, &num_intercluster_kmers));
    }

    for (auto & counting_thread : counting_threads) {
        	
       	counting_thread.join();
	}

	return num_intercluster_kmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countSampleKmersCallBack(ProducerConsumerQueue<pair<uint, PrecountedKmerBatch*> > * precounted_kmer_queue, mutex * counting_mutex, uint64 * total_processed_kmers, uint64 * total_included_kmers) {

	uint64 thread_processed_kmers = 0;
	uint64 thread_included_kmers = 0;

	KmerForward<kmer_size> kmer;
    KmerForward<Utils::small_kmer_size> small_mer;

	pair<uint, PrecountedKmerBatch*> precounted_kmer_batch;
	while (precounted_kmer_queue->pop(&precounted_kmer_batch)) {

		thread_processed_kmers += precounted_kmer_batch.second->size();

		for (auto precounted_kmer_iter = precounted_kmer_batch.second->begin(); precounted_kmer_iter != precounted_kmer_batch.second->end(); precounted_kmer_iter++) {

	    	kmer.reset();
			small_mer.reset();

			bool found_all_small_mers = true;

			for (uchar kmer_idx = 0; kmer_idx < kmer_size; kmer_idx++) {

				char current_kmer_char = precounted_kmer_iter->first.get_asci_symbol(kmer_idx);
	    	    auto kmer_move = kmer.move(current_kmer_char);
				auto small_mer_move = small_mer.move(current_kmer_char);
   				
   				assert(kmer_move.first);
   				assert(small_mer_move.first);

   				if (small_mer_move.second) {

					if (!smallmer_set->count(small_mer.getKmer())) {

						found_all_small_mers = false;
						break;
					}
				}
	    	}

	    	if (found_all_small_mers) {

				auto hash_lock = kmer_hash_hybrid->_hash->lockKey(kmer.getKmer());
				auto sample_kmer_insert = kmer_hash_hybrid->_hash->insert(kmer.getKmer(), nullptr);

				if (sample_kmer_insert.second) {

					*sample_kmer_insert.first = new SampleKmerCounts(num_samples);
				}

				static_cast<SampleKmerCounts*>(*sample_kmer_insert.first)->addCount(precounted_kmer_batch.first, precounted_kmer_iter->second);

				thread_included_kmers++;
	    	}
		}

		delete precounted_kmer_batch.second;
	}

	lock_guard<mutex> counting_lock(*counting_mutex);
	(*total_processed_kmers) += thread_processed_kmers;
	(*total_included_kmers) += thread_included_kmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countSampleKmers(const vector<Sample> & samples) {

	mutex counting_mutex;
    uint64 expected_processed_kmers = 0;
	uint64 total_processed_kmers = 0;
	uint64 total_included_kmers = 0;
    uint32 min_sample_count_global = 1;
    uint32 max_sample_count_global = Utils::uchar_overflow;

	ProducerConsumerQueue<pair<uint, PrecountedKmerBatch*> > precounted_kmer_queue(Utils::queue_size_scaling * num_threads);

	// Spawn worker threads 
	vector<thread> filter_threads(num_threads);

	for (int i=0; i < num_threads; i++) {

        filter_threads.at(i) = thread(&KmerCounter::countSampleKmersCallBack, this, &precounted_kmer_queue, &counting_mutex, &total_processed_kmers, &total_included_kmers);
    }  

	for (uint sample_idx = 0; sample_idx < num_samples; sample_idx++) {

	    // Read KMC2 database
		CKMCFile kmer_database;

		if (!kmer_database.OpenForListing(samples.at(sample_idx).file)) {

			cout << "ERROR: Could not open KMC database " << samples.at(sample_idx).file<< endl;
			exit(1);
		}

		uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_count, max_count;
		uint64 total_kmers;
		assert(kmer_database.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_count, max_count, total_kmers));

		assert(db_kmer_size == kmer_size);
		expected_processed_kmers += total_kmers;

		if (sample_idx == 0) {

			min_sample_count_global = min_count;
			max_sample_count_global = max_count;
		
		} else {

			assert(min_sample_count_global == min_count);
			assert(max_sample_count_global == max_count);
		}

		cout << "[" << Utils::getLocalTime() << "] Reading kmer table with " << total_kmers << " kmers for sample " << samples.at(sample_idx).name << ". Only kmers with " << min_count << " <= multiplicity <= " << max_count << " were included by KMC2 ..." << endl;

		PrecountedKmerBatch * current_kmer_batch = new PrecountedKmerBatch(precounted_kmer_batch_size, make_pair(CKmerAPI(kmer_size),0));

		uint current_kmer_batch_idx = 0;
		uint num_kmers_from_file = 0;

		while (kmer_database.ReadNextKmer(current_kmer_batch->at(current_kmer_batch_idx).first, current_kmer_batch->at(current_kmer_batch_idx).second)) {
			
			current_kmer_batch_idx++;
			num_kmers_from_file++;

			if (current_kmer_batch_idx == precounted_kmer_batch_size) {

				precounted_kmer_queue.push(std::make_pair(sample_idx, current_kmer_batch));
				current_kmer_batch = new PrecountedKmerBatch(precounted_kmer_batch_size, make_pair(CKmerAPI(kmer_size),0));
				current_kmer_batch_idx = 0;
			}
		}

		// Shrink and push last vector
		if (current_kmer_batch_idx > 0) {

			current_kmer_batch->resize(current_kmer_batch_idx);
			precounted_kmer_queue.push(std::make_pair(sample_idx, current_kmer_batch));

		} else {

			delete current_kmer_batch;
		}
	}

	precounted_kmer_queue.pushedLast();

	for (auto & thread : filter_threads) {

    	thread.join();
    } 

    assert(expected_processed_kmers == total_processed_kmers);
    kmer_hash_hybrid->setMinKmerCount(min_sample_count_global);

	return total_included_kmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterKmersCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const ushort num_haplotype_candidates_per_sample, const ushort thread_index, mutex * thread_mutex, ulong * num_unique_kmers) {

    ulong local_num_unique_kmers = 0;
    uint variant_cluster_group_index = thread_index;

	while (variant_cluster_group_index < variant_cluster_groups->size()) {

		local_num_unique_kmers += variant_cluster_groups->at(variant_cluster_group_index)->countKmers(kmer_hash_hybrid, prng_seed, num_haplotype_candidates_per_sample, variant_cluster_group_index);
        variant_cluster_group_index += num_threads;
	} 

    unique_lock<mutex> counting_lock(*thread_mutex);
    *num_unique_kmers += local_num_unique_kmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countVariantClusterKmers(vector<VariantClusterGroup *> * variant_cluster_groups, const ushort num_haplotype_candidates_per_sample) {

    ulong num_unique_kmers = 0;

	vector<thread> counting_threads;
    mutex thread_mutex;

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countVariantClusterKmersCallback, this, variant_cluster_groups, num_haplotype_candidates_per_sample, i, &thread_mutex, &num_unique_kmers));
    }

    for (auto & counting_thread: counting_threads) {
        	
       	counting_thread.join();
	}

    return num_unique_kmers;
}


template class KmerCounter<31>;
template class KmerCounter<39>;
template class KmerCounter<47>;
template class KmerCounter<55>;
template class KmerCounter<63>;


