
/*
KmerCounter.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

static const uint max_intercluster_kmers_per_bias_bin = 10000000;
static const uint kmer_batch_size = 500000;

template <uchar kmer_size>
KmerCounter<kmer_size>::KmerCounter(const ushort num_threads_in, const ulong expected_num_path_kmers, const uchar num_genomic_rate_gc_bias_bins) : num_threads(num_threads_in), max_intercluster_kmers(max_intercluster_kmers_per_bias_bin * num_genomic_rate_gc_bias_bins) {

	path_bloom = new ThreadedBasicKmerBloom<kmer_size>(1024, expected_num_path_kmers + max_intercluster_kmers, 0.001);
}

template <uchar kmer_size>
KmerCounter<kmer_size>::~KmerCounter() {

	delete path_bloom;
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::findVariantClusterPathsCallback(vector<VariantClusterGraph *> * variant_cluster_graphs, KmerBloom * kmer_bloom, const ushort sample_idx, const uint prng_seed, const ushort max_sample_haplotype_candidates, const ushort thread_idx) {

    uint variant_cluster_graph_idx = thread_idx;

	while (variant_cluster_graph_idx < variant_cluster_graphs->size()) {

		variant_cluster_graphs->at(variant_cluster_graph_idx)->findSamplePaths(path_bloom, kmer_bloom, prng_seed + variant_cluster_graph_idx + sample_idx, max_sample_haplotype_candidates);
        variant_cluster_graph_idx += num_threads;
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::findVariantClusterPaths(vector<VariantClusterGraph *> * variant_cluster_graphs, const vector<Sample> & samples, const uint prng_seed, const ushort max_sample_haplotype_candidates) {

	KmerBloom * cur_kmer_bloom = nullptr;
	KmerBloom * next_kmer_bloom = nullptr;

	cout << endl;

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		cout << "[" << Utils::getLocalTime() << "] Finding variant cluster paths for sample " << samples.at(sample_idx).name << " ..." << endl;

		if (sample_idx == 0) {

			cur_kmer_bloom = new BasicKmerBloom<kmer_size>(samples.at(sample_idx).file + ".bloom");

		} else {

			cur_kmer_bloom = next_kmer_bloom;
		}

		vector<thread> finding_threads;

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	   	    finding_threads.push_back(thread(&KmerCounter<kmer_size>::findVariantClusterPathsCallback, this, variant_cluster_graphs, cur_kmer_bloom, sample_idx, prng_seed, max_sample_haplotype_candidates, thread_idx));
	    }

	    if (static_cast<uint>(sample_idx + 1) < samples.size()) {

	    	next_kmer_bloom = new BasicKmerBloom<kmer_size>(samples.at(sample_idx + 1).file + ".bloom");
	    }

	    for (auto & finding_thread: finding_threads) {

	       	finding_thread.join();
		}

		delete cur_kmer_bloom;
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::selectParameterInterclusterKmersCallback(const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const ushort thread_idx, const uint max_thread_intercluster_kmers) {

    uint intercluster_regions_idx = thread_idx;
    uint num_thread_intercluster_kmers = 0;

    KmerPair<kmer_size> kmer_pair = KmerPair<kmer_size>();

	while (intercluster_regions_idx < intercluster_regions.size()) {

		assert((intercluster_regions.at(intercluster_regions_idx).end - intercluster_regions.at(intercluster_regions_idx).start) >= kmer_size);
    	auto sit = intercluster_regions.at(intercluster_regions_idx).start;

    	kmer_pair.reset();

		while (sit != intercluster_regions.at(intercluster_regions_idx).end) {

			auto nt_bits = Sequence::ntToBit<1>(*sit);

			if (nt_bits.second) {

				if (kmer_pair.move(nt_bits.first)) {

					auto lowest_kmer = kmer_pair.getLexicographicalLowestKmer();
                    auto bloom_lock = path_bloom->lockBloom(lowest_kmer);

					if (!(path_bloom->lookup(lowest_kmer))) {

						num_thread_intercluster_kmers++;
						path_bloom->addKmer(lowest_kmer);
					}

					if (num_thread_intercluster_kmers == max_thread_intercluster_kmers) {

						break;
					}
				}

			} else {

				kmer_pair.reset();
			}

			sit++;
		}

		if (num_thread_intercluster_kmers == max_thread_intercluster_kmers) {

			break;
		}

        intercluster_regions_idx += num_threads;
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countInterclusterKmersCallback(KmerHash * kmer_hash, const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const ushort thread_idx) {

    uint intercluster_regions_idx = thread_idx;

    KmerPair<kmer_size> kmer_pair = KmerPair<kmer_size>();

	while (intercluster_regions_idx < intercluster_regions.size()) {

		assert((intercluster_regions.at(intercluster_regions_idx).end - intercluster_regions.at(intercluster_regions_idx).start) >= kmer_size);
    	auto sit = intercluster_regions.at(intercluster_regions_idx).start;

    	kmer_pair.reset();

		while (sit != intercluster_regions.at(intercluster_regions_idx).end) {

			auto nt_bits = Sequence::ntToBit<1>(*sit);

			if (nt_bits.second) {

				if (kmer_pair.move(nt_bits.first)) {

					auto lowest_kmer = kmer_pair.getLexicographicalLowestKmer();

					if (path_bloom->lookup(lowest_kmer)) {

						auto kmer_lock = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->getKmerLock(lowest_kmer);

						auto kmer_counts = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->addKmer(lowest_kmer, true);
						assert(kmer_counts.first);

						kmer_counts.first->addInterclusterMultiplicity(intercluster_regions.at(intercluster_regions_idx).chromosome_class);
					} 
				}

			} else {

				kmer_pair.reset();
			}

			sit++;
		}

        intercluster_regions_idx += num_threads;
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countInterclusterKmers(KmerHash * kmer_hash, const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const uchar num_genomic_rate_gc_bias_bins) {

    cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in genomic regions between variant clusters and decoy sequence(s) ..." << endl;

	vector<thread> selection_threads;
	selection_threads.reserve(num_threads);

	for (uint i = 0; i < num_threads; i++) {

   	    selection_threads.push_back(thread(&KmerCounter<kmer_size>::selectParameterInterclusterKmersCallback, this, intercluster_regions, i, ceil((max_intercluster_kmers)/static_cast<float>(num_threads))));
    }

    for (auto & selection_thread: selection_threads) {

       	selection_thread.join();
	}

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (uint i = 0; i < num_threads; i++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countInterclusterKmersCallback, this, kmer_hash, intercluster_regions, i));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::parseSampleKmersCallBack(KmerHash * kmer_hash, ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_forward, ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_backward) {

	KmerForward<kmer_size> kmer;
	KmerBatchInfo * kmer_batch_info;

	while (kmer_batch_queue_forward->pop(&kmer_batch_info)) {

		auto kmer_it = kmer_batch_info->begin_it;

		while (kmer_it != kmer_batch_info->end_it) {

			assert(kmer_it->second <= Utils::uchar_overflow);

	    	kmer.reset();

			for (uchar kmer_idx = 0; kmer_idx < kmer_size; kmer_idx++) {

				bitset<2> nt_bits(kmer_it->first.get_num_symbol(kmer_idx));
	    	    kmer.move(nt_bits);
	    	}

			if (path_bloom->lookup(kmer.getKmer())) {

				auto kmer_lock = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->getKmerLock(kmer.getKmer());

				auto kmer_added = static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->addKmer(kmer.getKmer(), false);
				assert(kmer_added.first);
					
				kmer_added.first->addSampleCount(kmer_batch_info->sample_idx, kmer_it->second);
			}

	    	kmer_it++;
		}

		kmer_batch_queue_backward->push(kmer_batch_info);
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::parseSampleKmers(KmerHash * kmer_hash, const vector<Sample> & samples) {

	vector<KmerBatch> kmer_batches(Utils::queue_size_thread_scaling * num_threads, KmerBatch(kmer_batch_size, make_pair(CKmerAPI(kmer_size), 0)));

	cout << endl;

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		CKMCFile kmer_database;

		if (!kmer_database.OpenForListing(samples.at(sample_idx).file)) {

			cout << "\nERROR: Could not open KMC table " << samples.at(sample_idx).file << "\n" << endl;
			exit(1);
		}

		uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count;
		uint64 max_kmer_count, total_kmers;

		assert(kmer_database.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count, max_kmer_count, total_kmers));

		cout << "[" << Utils::getLocalTime() << "] Parsing KMC table containing " << total_kmers << " kmers for sample " << samples.at(sample_idx).name << " (kmers observed less than " << min_kmer_count << " time(s) were excluded by KMC) ..." << endl;

		assert(db_kmer_size == kmer_size);
		assert(!mode);

		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_forward(Utils::queue_size_thread_scaling * num_threads);
		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_backward(Utils::queue_size_thread_scaling * num_threads);

		for (uint kmer_batch_idx = 0; kmer_batch_idx < (Utils::queue_size_thread_scaling * num_threads); kmer_batch_idx++) {

			kmer_batch_queue_backward.push(new KmerBatchInfo(sample_idx, kmer_batches.at(kmer_batch_idx).begin(), kmer_batches.at(kmer_batch_idx).end()));
		}

		vector<thread> parsing_threads;
		parsing_threads.reserve(num_threads);

		for (uint i = 0; i < num_threads; i++) {

	        parsing_threads.push_back(thread(&KmerCounter::parseSampleKmersCallBack, this, kmer_hash, &kmer_batch_queue_forward, &kmer_batch_queue_backward));
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

		for (auto & parsing_thread: parsing_threads) {

    		parsing_thread.join();
    	}

		kmer_batch_queue_backward.pushedLast();

		while (kmer_batch_queue_backward.pop(&kmer_batch_info)) {

			delete kmer_batch_info;
		}

		static_cast<BasicKmerHash<kmer_size> *>(kmer_hash)->sortKmers();
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterKmersCallback(KmerHash * kmer_hash, vector<VariantClusterGroup *> * variant_cluster_groups, const ushort thread_idx) {

    uint variant_cluster_group_idx = thread_idx;

	while (variant_cluster_group_idx < variant_cluster_groups->size()) {

		variant_cluster_groups->at(variant_cluster_group_idx)->countKmers(kmer_hash, variant_cluster_group_idx);
        variant_cluster_group_idx += num_threads;
	}
}

template <uchar kmer_size>
void KmerCounter<kmer_size>::countVariantClusterKmers(KmerHash * kmer_hash, vector<VariantClusterGroup *> * variant_cluster_groups) {

    cout << "\n[" << Utils::getLocalTime() << "] Counting kmers in variant cluster paths ..." << endl;

	vector<thread> counting_threads;

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter<kmer_size>::countVariantClusterKmersCallback, this, kmer_hash, variant_cluster_groups, thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}

template class KmerCounter<31>;
template class KmerCounter<39>;
template class KmerCounter<47>;
template class KmerCounter<55>;
template class KmerCounter<63>;
