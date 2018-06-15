
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
#include <unordered_set>

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "kmc_api/kmc_file.h"

#include "KmerCounter.hpp"
#include "Utils.hpp"
#include "KmerHash.hpp"
#include "Kmer.hpp"
#include "Nucleotide.hpp"
#include "KmerCounts.hpp"
#include "Sample.hpp"
#include "VariantFileParser.hpp"


static const uint kmer_batch_size = 1000000;


KmerCounter::KmerCounter(const vector<Sample> & samples_in, const OptionsContainer & options_container) : samples(samples_in), num_threads(options_container.getValue<ushort>("threads")), prng_seed(options_container.getValue<uint>("random-seed")) {}

void KmerCounter::findVariantClusterPathsCallback(vector<VariantClusterGroup *> * variant_cluster_groups, KmerBloom<Utils::kmer_size> * sample_kmer_bloom, const ushort max_sample_haplotype_candidates, const ushort sample_idx, const ushort thread_idx) {

    uint variant_cluster_group_idx = thread_idx;

	while (variant_cluster_group_idx < variant_cluster_groups->size()) {

		variant_cluster_groups->at(variant_cluster_group_idx)->findSamplePaths(sample_kmer_bloom, prng_seed * (sample_idx + 1) + variant_cluster_group_idx, max_sample_haplotype_candidates);
        variant_cluster_group_idx += num_threads;
	}
}

void KmerCounter::findVariantClusterPaths(InferenceUnit * inference_unit, const ushort max_sample_haplotype_candidates) {
 
 	cout << "[" << Utils::getLocalTime() << "] Finding variant cluster paths for " << samples.size() << " sample(s) ..." << endl;

	assert(!(samples.empty()));

	KmerBloom<Utils::kmer_size> * cur_sample_kmer_bloom = nullptr;
	KmerBloom<Utils::kmer_size> * next_sample_kmer_bloom = new KmerBloom<Utils::kmer_size>(samples.front().file);

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		cur_sample_kmer_bloom = next_sample_kmer_bloom;

		vector<thread> finding_threads;
		finding_threads.reserve(num_threads);

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	   	    finding_threads.push_back(thread(&KmerCounter::findVariantClusterPathsCallback, this, &(inference_unit->variant_cluster_groups), cur_sample_kmer_bloom, max_sample_haplotype_candidates, sample_idx, thread_idx));
	    }

	    if (static_cast<uint>(sample_idx + 1) < samples.size()) {

	    	next_sample_kmer_bloom = new KmerBloom<Utils::kmer_size>(samples.at(sample_idx + 1).file);
	    }
	    
	    for (auto & finding_thread: finding_threads) {

	       	finding_thread.join();
		}

		delete cur_sample_kmer_bloom;
	}
}

void KmerCounter::countPathMultigroupKmersCallback(BooleanKmerHash * multigroup_kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, vector<VariantClusterGroup *> * variant_cluster_groups, mutex * counting_mutex, ulong * num_kmers_global, const ushort thread_idx) {

	ulong num_kmers_local = 0;

    uint variant_cluster_group_idx = thread_idx;

    unordered_set<bitset<Utils::kmer_size * 2> > variant_cluster_group_kmers;

	while (variant_cluster_group_idx < variant_cluster_groups->size()) {

		variant_cluster_group_kmers.clear();
		variant_cluster_groups->at(variant_cluster_group_idx)->countPathKmers(&variant_cluster_group_kmers);

        num_kmers_local += variant_cluster_group_kmers.size();
	
		for (auto & kmer: variant_cluster_group_kmers) {

			const string kmer_str = Nucleotide::bitToNt<Utils::kmer_size>(kmer);

			auto bloom_lock = path_kmer_bloom->getKmerLock(kmer_str);

			if (path_kmer_bloom->lookup(kmer_str)) {

				auto hash_lock = multigroup_kmer_hash->getKmerLock(kmer);
	
				auto multigroup_kmer = multigroup_kmer_hash->addKmer(kmer);
				assert(multigroup_kmer.first);
			
			} else {

				path_kmer_bloom->addKmer(kmer_str);
			}
		}

		variant_cluster_group_idx += num_threads;
	}

	lock_guard<mutex> counting_lock(*counting_mutex);
	*num_kmers_global += num_kmers_local;
}

void KmerCounter::countPathMultigroupKmers(BooleanKmerHash * multigroup_kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, InferenceUnit * inference_unit) {

    cout << "[" << Utils::getLocalTime() << "] Counting multigroup kmers in variant cluster paths ..." << endl;

    assert(inference_unit->num_path_kmers == 0);

    mutex counting_mutex;
    ulong num_kmers = 0;

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter::countPathMultigroupKmersCallback, this, multigroup_kmer_hash, path_kmer_bloom, &(inference_unit->variant_cluster_groups), &counting_mutex, &num_kmers, thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}

	inference_unit->num_path_kmers = num_kmers;
}

void KmerCounter::countInterclusterParameterKmersCallback(BooleanKmerHash * parameter_kmer_hash, const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const float parameter_kmer_fraction, const Chromosomes & chromosomes, const ThreadedKmerBloom<Utils::kmer_size> & path_kmer_bloom, const ushort thread_idx) {

    uint intercluster_regions_idx = thread_idx;

    KmerPair<Utils::kmer_size> kmer_pair;

    mt19937 prng = mt19937(prng_seed * (thread_idx + 1));
    bernoulli_distribution bernoulli_dist(parameter_kmer_fraction);

	while (intercluster_regions_idx < intercluster_regions.size()) {

    	kmer_pair.reset();

		auto chromosomes_it = chromosomes.find(intercluster_regions.at(intercluster_regions_idx).chrom_name);
		assert(chromosomes_it != chromosomes.cend());

    	auto sequence_cit = chromosomes_it->second.cbegin() + intercluster_regions.at(intercluster_regions_idx).start;
    	auto sequence_ceit = chromosomes_it->second.cbegin() + intercluster_regions.at(intercluster_regions_idx).end + 1;

		while (sequence_cit != sequence_ceit) {

			if (kmer_pair.move(Nucleotide::ntToBit<1>(*sequence_cit))) {

				auto lowest_kmer = kmer_pair.getLexicographicalLowestKmer();
				
				if (!(path_kmer_bloom.lookup(lowest_kmer))) {

					if (intercluster_regions.at(intercluster_regions_idx).chrom_class == Utils::ChromClass::Decoy) {

						auto hash_lock = parameter_kmer_hash->getKmerLock(lowest_kmer);

						auto new_parameter_kmer = parameter_kmer_hash->addKmer(lowest_kmer);
						assert(new_parameter_kmer.first);

						*(new_parameter_kmer.first) = false;

					} else if (bernoulli_dist(prng)) {

						auto hash_lock = parameter_kmer_hash->getKmerLock(lowest_kmer);

						auto new_parameter_kmer = parameter_kmer_hash->addKmer(lowest_kmer);
						assert(new_parameter_kmer.first);
						
						if (new_parameter_kmer.second) {

							*(new_parameter_kmer.first) = true;
						} 
					}
				}
			}

			sequence_cit++;
		}

        intercluster_regions_idx += num_threads;
	}
}

void KmerCounter::countInterclusterParameterKmers(BooleanKmerHash * parameter_kmer_hash, const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const float parameter_kmer_fraction, const Chromosomes & chromosomes, const ThreadedKmerBloom<Utils::kmer_size> & path_kmer_bloom) {

    cout << "[" << Utils::getLocalTime() << "] Counting parameter kmers in inter-cluster regions ..." << endl;

	assert(parameter_kmer_fraction > 0);
	assert(parameter_kmer_fraction <= 1);

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter::countInterclusterParameterKmersCallback, this, parameter_kmer_hash, ref(intercluster_regions), parameter_kmer_fraction, ref(chromosomes), ref(path_kmer_bloom), thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}

void KmerCounter::countPathKmersCallback(ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, vector<VariantClusterGroup *> * variant_cluster_groups, const ushort thread_idx) {

    uint variant_cluster_group_idx = thread_idx;

    unordered_set<bitset<Utils::kmer_size * 2> > variant_cluster_group_kmers;

	while (variant_cluster_group_idx < variant_cluster_groups->size()) {

		variant_cluster_group_kmers.clear();

		variant_cluster_groups->at(variant_cluster_group_idx)->countPathKmers(&variant_cluster_group_kmers);
        variant_cluster_group_idx += num_threads;
	
		for (auto & kmer: variant_cluster_group_kmers) {

			auto bloom_lock = path_kmer_bloom->getKmerLock(kmer);
			path_kmer_bloom->addKmer(kmer);
		}
	}
}

void KmerCounter::countPathKmers(ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, InferenceUnit * inference_unit) {

	cout << "[" << Utils::getLocalTime() << "] Counting kmers in variant cluster paths ..." << endl;

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter::countPathKmersCallback, this, path_kmer_bloom, &(inference_unit->variant_cluster_groups), thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}

void KmerCounter::countInterclusterKmersCallback(KmerCountsHash * kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, const vector<VariantFileParser::InterClusterRegion> & intercluster_regions, const Chromosomes & chromosomes, const ushort thread_idx) {

    uint intercluster_regions_idx = thread_idx;

    KmerPair<Utils::kmer_size> kmer_pair;

	while (intercluster_regions_idx < intercluster_regions.size()) {

    	kmer_pair.reset();

		auto chromosomes_it = chromosomes.find(intercluster_regions.at(intercluster_regions_idx).chrom_name);
		assert(chromosomes_it != chromosomes.cend());

    	auto sequence_cit = chromosomes_it->second.cbegin() + intercluster_regions.at(intercluster_regions_idx).start;
    	auto sequence_ceit = chromosomes_it->second.cbegin() + intercluster_regions.at(intercluster_regions_idx).end + 1;

		while (sequence_cit != sequence_ceit) {

			if (kmer_pair.move(Nucleotide::ntToBit<1>(*sequence_cit))) {

				auto lowest_kmer = kmer_pair.getLexicographicalLowestKmer();

				if (path_kmer_bloom->lookup(lowest_kmer)) {

					auto hash_lock = kmer_hash->getKmerLock(lowest_kmer);

					auto kmer_counts = kmer_hash->addKmer(lowest_kmer, true);
					assert(kmer_counts.first);

					kmer_counts.first->addInterclusterMultiplicity(intercluster_regions.at(intercluster_regions_idx).chrom_class);
				}
			}

			sequence_cit++;
		}

        intercluster_regions_idx += num_threads;
	}
}

void KmerCounter::countInterclusterKmers(KmerCountsHash * kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, const string & intercluster_regions_prefix, const Chromosomes & chromosomes) {

	cout << "[" << Utils::getLocalTime() << "] Counting kmers in inter-cluster regions and decoy sequence(s) ..." << endl;

    vector<VariantFileParser::InterClusterRegion> intercluster_regions;

    {
	    ifstream regions_infile(intercluster_regions_prefix + ".txt.gz", std::ios::binary);
	    assert(regions_infile.is_open());

		boost::iostreams::filtering_istream regions_infile_fstream;
		
		regions_infile_fstream.push(boost::iostreams::gzip_decompressor());
    	regions_infile_fstream.push(boost::ref(regions_infile));

		assert(regions_infile_fstream.is_complete());    

	    vector<string> line(4, "");

	    while (regions_infile_fstream.good()) {

	    	for (ushort i = 0; i < 3; i++) {

		        getline(regions_infile_fstream, line.at(i), '\t');
		    }

		    getline(regions_infile_fstream, line.at(3), '\n');

		    if (line.at(0).empty() or (line.at(0) == "Chrom")) {

		    	continue;
		    }
	        
	        intercluster_regions.emplace_back(line.at(0), static_cast<Utils::ChromClass>(stoi(line.at(1))), stoi(line.at(2)), stoi(line.at(3)));
		}
	}

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter::countInterclusterKmersCallback, this, kmer_hash, path_kmer_bloom, ref(intercluster_regions), ref(chromosomes), thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}

void KmerCounter::parseSampleKmersCallBack(KmerCountsHash * kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom, ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_forward, ProducerConsumerQueue<KmerBatchInfo *> * kmer_batch_queue_backward) {

	KmerBatchInfo * kmer_batch_info;

	char * kmer_seq = new char [Utils::kmer_size + 1];
	bitset<Utils::kmer_size * 2> kmer_bitset;

	while (kmer_batch_queue_forward->pop(&kmer_batch_info)) {

		auto kmer_it = kmer_batch_info->begin_it;

		while (kmer_it != kmer_batch_info->end_it) {

			assert(kmer_it->second <= Utils::uchar_overflow);

			for (uchar kmer_idx = 0; kmer_idx < Utils::kmer_size; kmer_idx++) {

				kmer_seq[kmer_idx] = kmer_it->first.get_asci_symbol(kmer_idx);
				
				bitset<2> nt_bits = kmer_it->first.get_num_symbol(kmer_idx);
	    	    kmer_bitset[kmer_idx * 2] = nt_bits[0];
	    	    kmer_bitset[kmer_idx * 2 + 1] = nt_bits[1];
	    	}

			if (path_kmer_bloom->lookup(kmer_seq)) {

				auto hash_lock = kmer_hash->getKmerLock(kmer_bitset);

				auto kmer_counts = kmer_hash->addKmer(kmer_bitset, false);
				assert(kmer_counts.first);
					
				kmer_counts.first->addSampleCount(kmer_batch_info->sample_idx, kmer_it->second);
			}

	    	kmer_it++;
		}

		kmer_batch_queue_backward->push(kmer_batch_info);
	}

	delete[] kmer_seq;
}

void KmerCounter::parseSampleKmers(KmerCountsHash * kmer_hash, ThreadedKmerBloom<Utils::kmer_size> * path_kmer_bloom) {

	vector<KmerBatch> kmer_batches(Utils::queue_size_thread_scaling * num_threads, KmerBatch(kmer_batch_size, make_pair(CKmerAPI(Utils::kmer_size), 0)));

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		CKMCFile kmc_table;

		if (!kmc_table.OpenForListing(samples.at(sample_idx).file)) {

			cout << "\nERROR: Could not open KMC table " << samples.at(sample_idx).file << "\n" << endl;
			exit(1);
		}

		CKMCFileInfo kmc_table_info;
		assert(kmc_table.Info(kmc_table_info));

		assert(kmc_table_info.kmer_length == Utils::kmer_size);
		assert(!(kmc_table_info.mode));

		cout << "[" << Utils::getLocalTime() << "] Parsing KMC table containing " << kmc_table_info.total_kmers << " kmers for sample " << samples.at(sample_idx).name << " ..." << endl;

		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_forward(Utils::queue_size_thread_scaling * num_threads);
		ProducerConsumerQueue<KmerBatchInfo *> kmer_batch_queue_backward(Utils::queue_size_thread_scaling * num_threads);

		for (uint kmer_batch_idx = 0; kmer_batch_idx < (Utils::queue_size_thread_scaling * num_threads); kmer_batch_idx++) {

			kmer_batch_queue_backward.push(new KmerBatchInfo(sample_idx, kmer_batches.at(kmer_batch_idx).begin(), kmer_batches.at(kmer_batch_idx).end()));
		}

		vector<thread> parsing_threads;
		parsing_threads.reserve(num_threads);

		for (uint i = 0; i < num_threads; i++) {

	        parsing_threads.push_back(thread(&KmerCounter::parseSampleKmersCallBack, this, kmer_hash, path_kmer_bloom, &kmer_batch_queue_forward, &kmer_batch_queue_backward));
	    }

	    KmerBatchInfo * kmer_batch_info;

	    assert(kmer_batch_queue_backward.pop(&kmer_batch_info));
		assert(kmer_batch_info->begin_it != kmer_batch_info->end_it);

		auto kmer_batch_it = kmer_batch_info->begin_it;

		bool is_kmer = kmc_table.ReadNextKmer(kmer_batch_it->first, kmer_batch_it->second);
		assert(is_kmer);

		while (true) {

			while (is_kmer) {

				kmer_batch_it++;

				if (kmer_batch_it == kmer_batch_info->end_it) {

					kmer_batch_queue_forward.push(kmer_batch_info);
					break;
				}

				is_kmer = kmc_table.ReadNextKmer(kmer_batch_it->first, kmer_batch_it->second);
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

		kmer_hash->sortKmers();
	}
}

void KmerCounter::classifyPathKmersCallback(KmerCountsHash * kmer_hash, KmerBloom<Utils::kmer_size> * multigroup_kmers_bloom, vector<VariantClusterGroup *> * variant_cluster_groups, const ushort thread_idx) {

    uint variant_cluster_group_idx = thread_idx;

	while (variant_cluster_group_idx < variant_cluster_groups->size()) {

		variant_cluster_groups->at(variant_cluster_group_idx)->classifyPathKmers(kmer_hash, multigroup_kmers_bloom);
        variant_cluster_group_idx += num_threads;
	}
}

void KmerCounter::classifyPathKmers(KmerCountsHash * kmer_hash, InferenceUnit * inference_unit, const string & multigroup_kmers_bloom_prefix) {

	cout << "[" << Utils::getLocalTime() << "] Classifying kmers in variant cluster paths ..." << endl;

    KmerBloom<Utils::kmer_size> multigroup_kmers_bloom(multigroup_kmers_bloom_prefix);

	vector<thread> counting_threads;
	counting_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

   	    counting_threads.push_back(thread(&KmerCounter::classifyPathKmersCallback, this, kmer_hash, &multigroup_kmers_bloom, &(inference_unit->variant_cluster_groups), thread_idx));
    }

    for (auto & counting_thread: counting_threads) {

       	counting_thread.join();
	}
}
