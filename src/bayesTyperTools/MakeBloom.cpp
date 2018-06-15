
/*
MakeBloom.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <unordered_set>
#include <thread>

#include "MakeBloom.hpp"

using namespace std;

static const ulong max_tests = 10000000;

void MakeBloom::kmc2bloomFile(const string & kmc_table_prefix, const float false_positive_rate, const bool test, const uint num_threads) {

    cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") makeBloom ...\n" << endl;

	assert(kmer_size <= static_cast<uint>(Utils::uchar_overflow));

	auto kmc_bloom = kmc2bloomThreaded(kmc_table_prefix, false_positive_rate, num_threads);
	std::cout << "\n[" << Utils::getLocalTime() << "] Saving bloom filter to " << kmc_table_prefix << " ..." << std::endl;
	kmc_bloom->save(kmc_table_prefix);
	delete kmc_bloom;

	std::cout << "[" << Utils::getLocalTime() << "] Completed saving bloom filter\n" << std::endl;

	if (test) {

		std::cout << "\n[" << Utils::getLocalTime() << "] Constructing test set based on the KMC kmers ...\n" << std::endl;

		CKMCFile kmc_table;

		if (!kmc_table.OpenForListing(kmc_table_prefix)) {

			cout << "ERROR: Could not open KMC table " << kmc_table_prefix << endl;
			exit(1);
		}

		CKMCFileInfo kmc_table_info;
		assert(kmc_table.Info(kmc_table_info));

		assert(kmc_table_info.kmer_length == kmer_size);
		assert(!(kmc_table_info.mode));

		CKmerAPI kmer(kmer_size);
		uint32 count;

		unordered_set<bitset<kmer_size * 2> > kmer_test_set;
		uint64_t kmers_parsed = 0;

		while (kmc_table.ReadNextKmer(kmer, count)) {

			assert(kmer_test_set.insert(kmc2bits(kmer)).second);
			++kmers_parsed;

			if (kmers_parsed % 100000000 == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << kmers_parsed << " kmers" << std::endl;
			}
		}

		std::cout << "\n[" << Utils::getLocalTime() << "] Reading bloom filter from file..." << std::endl;
		kmer_bloom_t kmc_bloom_serialised(kmc_table_prefix);
		std::cout << "[" << Utils::getLocalTime() << "] Done reading bloom filter from file!" << std::endl;

		testbloom(kmc_bloom_serialised, kmer_test_set);

		kmer_test_set.clear();
		const uint num_random_tests = 10000000;

		std::cout << "\n[" << Utils::getLocalTime() << "] Constructing test set based on " << num_random_tests << " random kmers ..." << std::endl;

		std::bitset<kmer_size*2> cur_kmer;

	    std::random_device rd;
	    std::mt19937 gen(rd());
	    std::bernoulli_distribution d(0.5);

		kmer_bloom_t random_kmer_bloom(num_random_tests, false_positive_rate);

		while (kmer_test_set.size() < num_random_tests) {

			cur_kmer<<=2;
			cur_kmer[0] = d(gen);
			cur_kmer[1] = d(gen);

			if (kmer_test_set.insert(cur_kmer).second) {

				random_kmer_bloom.addKmer(cur_kmer);
			}
		}

		testbloom(random_kmer_bloom, kmer_test_set);
	}
}

std::bitset<kmer_size*2> MakeBloom::kmc2bits(CKmerAPI & kmer) {

	std::bitset<kmer_size * 2> bits;

	for (unsigned int i = 0; i < kmer_size; ++i) {

		std::bitset<2> nt(kmer.get_num_symbol(i));
		bits[i*2] = nt[0];
		bits[i*2 + 1] = nt[1];
	}

	return bits;
}

kmer_bloom_t * MakeBloom::kmc2bloom(const string & kmc_table_fn, const float false_positive_rate) {

	CKMCFile kmc_table;

	if (!kmc_table.OpenForListing(kmc_table_fn)) {

		cout << "ERROR: Could not open KMC table " << kmc_table_fn << endl;
		exit(1);
	}

	CKMCFileInfo kmc_table_info;
	assert(kmc_table.Info(kmc_table_info));

	assert(kmc_table_info.kmer_length == kmer_size);
	assert(!(kmc_table_info.mode));

	std::cout << "[" << Utils::getLocalTime() << "] Making bloom filter of " << kmc_table_info.total_kmers << " kmers with a false positive rate of " << false_positive_rate << " ...\n" << std::endl;

	kmer_bloom_t * kmer_bloom = new kmer_bloom_t(kmc_table_info.total_kmers, false_positive_rate);

	ulong kmers_parsed = 0;

	CKmerAPI kmer(kmer_size);
	uint32 count;

	char * kmer_seq = new char [kmer_size + 1];

	while (kmc_table.ReadNextKmer(kmer, count)) {

		kmer.to_string(kmer_seq);
		kmer_bloom->addKmer(kmer_seq);

		kmers_parsed++;

		if (kmers_parsed % 100000000 == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << kmers_parsed << " kmers" << std::endl;
		}
	}

	delete[] kmer_seq;

	return kmer_bloom;
}

void MakeBloom::bloomInsertCallBack(kmer_bloom_t * kmer_bloom, ProducerConsumerQueue<kmer_batch_t*> * kmer_queue, ProducerConsumerQueue<kmer_batch_t*> * used_kmer_queue) {

  char * kmer_seq = new char [kmer_size + 1];
  kmer_batch_t * cur_read_batch;

	while (kmer_queue->pop(&cur_read_batch)) {

    for (uint kmer_idx = 0; kmer_idx < cur_read_batch->size(); ++kmer_idx) {

      cur_read_batch->at(kmer_idx)->to_string(kmer_seq);
      kmer_bloom->addKmer(kmer_seq);
    }

    used_kmer_queue->push(cur_read_batch);
	}

  delete[] kmer_seq;
}

kmer_bloom_t * MakeBloom::kmc2bloomThreaded(const string & kmc_table_fn, const float false_positive_rate, const uint num_threads) {

	CKMCFile kmc_table;

	if (!kmc_table.OpenForListing(kmc_table_fn)) {

		cout << "ERROR: Could not open KMC table " << kmc_table_fn << endl;
		exit(1);
	}

	CKMCFileInfo kmc_table_info;
	assert(kmc_table.Info(kmc_table_info));

	assert(kmc_table_info.kmer_length == kmer_size);
	assert(!(kmc_table_info.mode));

	std::cout << "[" << Utils::getLocalTime() << "] Making bloom filter of " << kmc_table_info.total_kmers << " kmers with a false positive rate of " << false_positive_rate << " ...\n" << std::endl;

	kmer_bloom_t * kmer_bloom = new kmer_bloom_t(kmc_table_info.total_kmers, false_positive_rate);

  uint num_batches = 2 * num_threads;
  uint batch_size = 1000000;

  ProducerConsumerQueue<kmer_batch_t*> kmer_queue(15000);
  ProducerConsumerQueue<kmer_batch_t*> used_kmer_queue(15000);

  for (uint batch_idx = 0; batch_idx < num_batches; ++batch_idx) {

    kmer_batch_t * kmer_batch = new kmer_batch_t();
    for (uint batch_kmer_idx = 0; batch_kmer_idx < batch_size; ++batch_kmer_idx) {

      kmer_batch->push_back(new CKmerAPI(kmer_size));
    }

    used_kmer_queue.push(kmer_batch);
  }

  vector<thread> bloom_insert_threads;
  for (uint thread_idx = 0; thread_idx < num_threads; ++thread_idx) {

    bloom_insert_threads.push_back(thread(&MakeBloom::bloomInsertCallBack, this, kmer_bloom, &kmer_queue, &used_kmer_queue));
  }

  kmer_batch_t * cur_write_batch;
  assert(used_kmer_queue.pop(&cur_write_batch));

  uint cur_write_batch_kmer_idx = 0;

  uint32 count;
	ulong kmers_parsed = 0;

  while (kmc_table.ReadNextKmer(*cur_write_batch->at(cur_write_batch_kmer_idx), count)) {

    ++cur_write_batch_kmer_idx;

    if (cur_write_batch_kmer_idx >= batch_size) {

      kmer_queue.push(cur_write_batch);
      assert(used_kmer_queue.pop(&cur_write_batch));
      cur_write_batch_kmer_idx = 0;
    }

    kmers_parsed++;
    if (kmers_parsed % 100000000 == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << kmers_parsed << " kmers" << std::endl;
		}
	}

  for (uint erase_idx = batch_size - 1; erase_idx >= cur_write_batch_kmer_idx; --erase_idx) {

    delete cur_write_batch->at(erase_idx);
    cur_write_batch->erase(cur_write_batch->begin() + erase_idx);
  }

  kmer_queue.push(cur_write_batch);
  kmer_queue.pushedLast();

  for (auto & thread : bloom_insert_threads) {

    thread.join();
  }

  used_kmer_queue.pushedLast();
  while (used_kmer_queue.pop(&cur_write_batch)) {

    for (uint kmer_idx = 0; kmer_idx < cur_write_batch->size(); ++kmer_idx) {

      delete cur_write_batch->at(kmer_idx);
    }

    delete cur_write_batch;
  }

	return kmer_bloom;
}

std::bitset<kmer_size*2> MakeBloom::random_kmer() {

	std::bitset<kmer_size*2> kmer;
    std::random_device rd;
    std::mt19937 gen( rd());
    std::bernoulli_distribution d(0.5);
    for (uint i = 0; i < kmer_size*2; ++i) {

		kmer[i] = d(gen);
    }

    return kmer;
}

void MakeBloom::testbloom(const kmer_bloom_t & kmer_bloom, const unordered_set<bitset<kmer_size * 2> > & kmer_test_set) {

	ulong num_tests = min(static_cast<ulong>(kmer_test_set.size()), max_tests);

	std::cout << "\n[" << Utils::getLocalTime() << "] Testing bloom filter using " << num_tests << " kmers ...\n" << std::endl;

    std::random_device rd;
    std::uniform_int_distribution<int> int_dist(0,kmer_size * 2 - 1);

	uint kmers_tested = 0;
	uint kmers_found_original = 0;
	uint kmers_found_perturbed_tp = 0;
	uint kmers_found_perturbed_fp = 0;
	uint kmers_found_random_tp = 0;
	uint kmers_found_random_fp = 0;

	for (auto & kmer: kmer_test_set) {

		kmers_tested++;

		auto kmer_bits = kmer;
		kmers_found_original += static_cast<unsigned int>(kmer_bloom.lookup(kmer_bits));

		auto bit_idx = int_dist(rd);
		kmer_bits[bit_idx] = ~kmer_bits[bit_idx];

		if (kmer_test_set.count(kmer_bits) > 0) {

			kmers_found_perturbed_tp++;

		} else {

			kmers_found_perturbed_fp += static_cast<unsigned int>(kmer_bloom.lookup(kmer_bits));
		}

		auto kmer_random = random_kmer();

		if (kmer_test_set.count(kmer_random) > 0) {

			kmers_found_random_tp++;

		} else {

			kmers_found_random_fp += static_cast<unsigned int>(kmer_bloom.lookup(kmer_random));
		}

		if (kmers_tested % 1000000 == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Tested " << kmers_tested << " kmers" << std::endl;
		}

		if (kmers_tested == num_tests) {

			break;
		}
	}

	std::cout << "\n[" << Utils::getLocalTime() << "] Out of " << num_tests << " kmers tested: " << endl;
	cout << "\n\t- Number of original kmers found: " << kmers_found_original << " (" << kmers_found_original/static_cast<float>(num_tests) << ")" << endl;
	cout << "\t- Number of perturbed kmers found (TP): " << kmers_found_perturbed_tp << " (" << kmers_found_perturbed_tp/static_cast<float>(num_tests) << ")" << endl;
	cout << "\t- Number of perturbed kmers found (FP): " << kmers_found_perturbed_fp << " (" << kmers_found_perturbed_fp/static_cast<float>(num_tests) << ")" << endl;
	cout << "\t- Number of random kmers found (TP): " << kmers_found_random_tp << " (" << kmers_found_random_tp/static_cast<float>(num_tests) << ")" << endl;
	cout << "\t- Number of random kmers found (FP): " << kmers_found_random_fp << " (" << kmers_found_random_fp/static_cast<float>(num_tests) << ")" << endl;
	cout << endl;
}
