
/*
KmerCounter.cpp - This file is part of BayesTyper (v1.1)


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

#include "kmc_api/kmc_file.h"

#include "KmerCounter.hpp"
#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Kmer.hpp"
#include "Sequence.hpp"
#include "KmerCounts.hpp"
#include "Sample.hpp"
#include "VariantFileParser.hpp"


static const uint kmer_batch_size = 500000;

template <uchar kmer_size>
KmerCounter<kmer_size>::KmerCounter(KmerHash * const kmer_hash_in, Utils::SmallmerSet * const smallmer_set_in, const ushort num_threads_in) : kmer_hash(kmer_hash_in), smallmer_set(smallmer_set_in), num_threads(num_threads_in) {

	min_sample_kmer_count = 0;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterSmallmersCallback(vector<VariantClusterGraph *> * variant_cluster_graphs, const ushort thread_index, mutex * counting_mutex, ulong * num_unique_smallmers) {

    ulong local_num_unique_smallmers = 0;
    uint variant_cluster_graph_index = thread_index;

	while (variant_cluster_graph_index < variant_cluster_graphs->size()) {

		local_num_unique_smallmers += variant_cluster_graphs->at(variant_cluster_graph_index)->countSmallmers(smallmer_set);
        variant_cluster_graph_index += num_threads;
	} 

    unique_lock<mutex> counting_lock(*counting_mutex);
    *num_unique_smallmers += local_num_unique_smallmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countVariantClusterSmallmers(vector<VariantClusterGraph *> * variant_cluster_graphs) {

    mutex counting_mutex;
    ulong num_unique_smallmers = 0;

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countVariantClusterSmallmersCallback, this, variant_cluster_graphs, i, &counting_mutex, &num_unique_smallmers));
    }

    for (auto & counting_thread: counting_threads) {
        	
       	counting_thread.join();
	}

    return num_unique_smallmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countInterclusterKmersCallback(const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const ushort thread_index, mutex * counting_mutex, ulong * num_intercluster_kmers) {

    ulong local_num_intercluster_kmers = 0;
    uint intercluster_regions_index = thread_index;

    KmerPair<kmer_size> kmer_pair = KmerPair<kmer_size>();
    KmerForward<Utils::small_kmer_size> small_mer = KmerForward<Utils::small_kmer_size>();

	while (intercluster_regions_index < intercluster_regions.size()) {

    	auto sit = intercluster_regions.at(intercluster_regions_index).start;

    	kmer_pair.reset();
    	small_mer.reset();

		while (sit != intercluster_regions.at(intercluster_regions_index).end) {

			auto nt_bits = Sequence::ntToBit<1>(*sit);

			if (nt_bits.second) {

				auto kmer_pair_move = kmer_pair.move(nt_bits.first);

				if (small_mer.move(nt_bits.first)) {

					if (smallmer_set->count(small_mer.getKmer())) {

						if (kmer_pair_move) {

							auto kmer_lock = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->getKmerLock(kmer_pair.getLexicographicalLowestKmer());
							auto added_kmer = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->addKmer(kmer_pair.getLexicographicalLowestKmer(), true);

							if (added_kmer.second) {

								local_num_intercluster_kmers++;
							}

							assert(added_kmer.first);
							added_kmer.first->addInterclusterMultiplicity(intercluster_regions.at(intercluster_regions_index).chromosome_class);
						}

					} else {
						
						kmer_pair.shrink(Utils::small_kmer_size);
					}
				}

			} else {

				kmer_pair.reset();
    			small_mer.reset();
			}

			sit++;
		} 

        intercluster_regions_index += num_threads;
	} 

	lock_guard<mutex> counting_lock(*counting_mutex);
    *num_intercluster_kmers += local_num_intercluster_kmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countInterclusterKmers(const vector<VariantFileParser::InterClusterRegion> & intercluster_regions) {

	mutex counting_mutex;
    ulong num_intercluster_kmers = 0;

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countInterclusterKmersCallback, this, intercluster_regions, i, &counting_mutex, &num_intercluster_kmers));
    }

    for (auto & counting_thread: counting_threads) {
        	
       	counting_thread.join();
	}

    static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->sortKmers();

	return num_intercluster_kmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countSampleKmersCallBack(ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_forward, ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_backward, mutex * counting_mutex, ulong * total_included_kmers) {

	ulong thread_included_kmers = 0;

	KmerForward<kmer_size> kmer;
    KmerForward<Utils::small_kmer_size> small_mer;

	KmerBatchInfo * kmer_batch_info;

	while (kmer_batch_queue_forward->pop(&kmer_batch_info)) {

		auto kmer_it = kmer_batch_info->begin_it;

		while (kmer_it != kmer_batch_info->end_it) {

			assert(kmer_it->second <= Utils::uchar_overflow);

	    	kmer.reset();
			small_mer.reset();

			bool found_all_small_mers = true;

			for (uchar kmer_idx = 0; kmer_idx < kmer_size; kmer_idx++) {

				auto nt_bits = Sequence::ntToBit<1>(kmer_it->first.get_asci_symbol(kmer_idx));
				assert(nt_bits.second);

   				if (small_mer.move(nt_bits.first)) {

					if (!smallmer_set->count(small_mer.getKmer())) {

						found_all_small_mers = false;
						break;
					}
				}

	    	    kmer.move(nt_bits.first);
	    	}

	    	if (found_all_small_mers) {

				auto kmer_lock = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->getKmerLock(kmer.getKmer());
				auto added_kmer = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->addKmer(kmer.getKmer(), false);

				if (added_kmer.second) {

					thread_included_kmers++;
				}

				assert(added_kmer.first);
				added_kmer.first->addSampleCount(kmer_batch_info->sample_idx, kmer_it->second);
	    	}

	    	kmer_it++;
		}

		kmer_batch_queue_backward->push(kmer_batch_info);
	}

	lock_guard<mutex> counting_lock(*counting_mutex);
	(*total_included_kmers) += thread_included_kmers;
}

template <uchar kmer_size>
ulong KmerCounter<kmer_size>::countSampleKmers(const vector<Sample> & samples) {
	
    uint min_kmer_count_all = 0;
    ulong max_kmer_count_all = 0;

	mutex counting_mutex;
	ulong total_included_kmers = 0;
    
	vector<KmerBatch> kmer_batches(Utils::queue_size_thread_scaling * num_threads, KmerBatch(kmer_batch_size, make_pair(CKmerAPI(kmer_size), 0)));

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		CKMCFile kmer_database;

		if (!kmer_database.OpenForListing(samples.at(sample_idx).file)) {

			cout << "\nERROR: Could not open KMC table " << samples.at(sample_idx).file << "\n" << endl;
			exit(1);
		}

		uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count;
		uint64 max_kmer_count, total_kmers;

		assert(kmer_database.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count, max_kmer_count, total_kmers));

		cout << "[" << Utils::getLocalTime() << "] Parsing kmer table with " << total_kmers << " kmers for sample " << samples.at(sample_idx).name << " (only kmers observed at least " << min_kmer_count << " time(s) were included by KMC) ..." << endl;

		if (mode) {

			cout << "\nERROR: Does not support Quake's compatible counting mode (-q) in KMC\n" << endl;
			exit(1);
		}

		assert(db_kmer_size == kmer_size);

		if (sample_idx == 0) {

			min_kmer_count_all = min_kmer_count;
			max_kmer_count_all = max_kmer_count;
		
		} else {

			assert(min_kmer_count_all == min_kmer_count);
			assert(max_kmer_count_all == max_kmer_count);
		}

		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_forward(Utils::queue_size_thread_scaling * num_threads);
		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_backward(Utils::queue_size_thread_scaling * num_threads);

		for (uint kmer_batch_idx = 0; kmer_batch_idx < (Utils::queue_size_thread_scaling * num_threads); kmer_batch_idx++) {

			kmer_batch_queue_backward.push(new KmerBatchInfo(sample_idx, kmer_batches.at(kmer_batch_idx).begin(), kmer_batches.at(kmer_batch_idx).end()));
		}

		vector<thread> counting_threads;
		counting_threads.reserve(num_threads);

		for (uint i = 0; i < num_threads; i++) {

	        counting_threads.push_back(thread(&KmerCounter::countSampleKmersCallBack, this, &kmer_batch_queue_forward, &kmer_batch_queue_backward, &counting_mutex, &total_included_kmers));
	    }  

	    KmerBatchInfo * kmer_batch_info;

	    assert(kmer_batch_queue_backward.pop(&kmer_batch_info));
		assert(kmer_batch_info->begin_it != kmer_batch_info->end_it);

		auto kmer_batch_it = kmer_batch_info->begin_it;

		bool is_kmer = kmer_database.ReadNextKmer(kmer_batch_it->first, kmer_batch_it->second);
		assert(is_kmer);

		while (true) {

			while (is_kmer) {

				kmer_batch_it++;

				if (kmer_batch_it == kmer_batch_info->end_it) {

					kmer_batch_queue_forward.push(kmer_batch_info);
					break;
				}

				is_kmer = kmer_database.ReadNextKmer(kmer_batch_it->first, kmer_batch_it->second);
			}

			if (!is_kmer) {

				kmer_batch_info->end_it = kmer_batch_it;
				kmer_batch_queue_forward.push(kmer_batch_info);
				break;
			}		

		    assert(kmer_batch_queue_backward.pop(&kmer_batch_info));
			assert(kmer_batch_info->begin_it != kmer_batch_info->end_it);

			kmer_batch_it = kmer_batch_info->begin_it;
		}

		kmer_batch_queue_forward.pushedLast();

		for (auto & counting_thread: counting_threads) {

    		counting_thread.join();
    	} 

		kmer_batch_queue_backward.pushedLast();

		while (kmer_batch_queue_backward.pop(&kmer_batch_info)) {

			delete kmer_batch_info;
		}

    	static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->sortKmers();
	}

	assert(min_kmer_count_all <= Utils::uchar_overflow);
	min_sample_kmer_count = min_kmer_count_all;

	return total_included_kmers;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterKmersCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const uint prng_seed, const ushort num_samples, const ushort max_sample_haplotype_candidates, const ushort thread_index) {

    uint variant_cluster_group_index = thread_index;

	while (variant_cluster_group_index < variant_cluster_groups->size()) {

		variant_cluster_groups->at(variant_cluster_group_index)->countKmers(kmer_hash, variant_cluster_group_index, prng_seed, num_samples, max_sample_haplotype_candidates);
        variant_cluster_group_index += num_threads;
	} 
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterKmers(vector<VariantClusterGroup *> * variant_cluster_groups, const uint prng_seed, const ushort num_samples, const ushort max_sample_haplotype_candidates) {

	vector<thread> counting_threads;

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countVariantClusterKmersCallback, this, variant_cluster_groups, prng_seed, num_samples, max_sample_haplotype_candidates, i));
    }

    for (auto & counting_thread: counting_threads) {
        	
       	counting_thread.join();
	}
}

template <uchar kmer_size>
uchar KmerCounter<kmer_size>::minSampleKmerCount() {

	return min_sample_kmer_count;
}


template class KmerCounter<31>;
template class KmerCounter<39>;
template class KmerCounter<47>;
template class KmerCounter<55>;
template class KmerCounter<63>;



