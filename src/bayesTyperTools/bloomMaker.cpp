
/*
BloomMaker.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "bloomMaker.hpp"

using namespace std;

static const ulong max_tests = 10000000;

template <uchar kmer_size>
std::bitset<kmer_size*2> BloomMaker<kmer_size>::kmc2bits(CKmerAPI & kmer) {

	std::bitset<kmer_size * 2> bits;

	for (unsigned int i = 0; i < kmer_size; ++i) {

		std::bitset<2> nt(kmer.get_num_symbol(i));
		bits[i*2] = nt[0];
		bits[i*2 + 1] = nt[1];
	}

	return bits;
}

template <uchar kmer_size>
BasicKmerBloom<kmer_size> BloomMaker<kmer_size>::kmc2bloom(const string & kmc_table_fn, float fpr) {

	CKMCFile kmc_table;

	if (!kmc_table.OpenForListing(kmc_table_fn)) {

		cout << "ERROR: Could not open KMC table " << kmc_table_fn << endl;
		exit(1);
	}

	uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count;
	uint64 max_kmer_count, total_kmers;

	assert(kmc_table.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count, max_kmer_count, total_kmers));
	assert(db_kmer_size == kmer_size);

	std::cout << "[" << Utils::getLocalTime() << "] Making bloom filter with " << total_kmers << " kmers with fpr " << fpr << ".\n" << std::endl;

	BasicKmerBloom<kmer_size> kmer_bloom(total_kmers, fpr);

	CKmerAPI kmer(db_kmer_size);
	uint32 count;
	ulong kmers_parsed = 0;

	while (kmc_table.ReadNextKmer(kmer, count)) {

		kmer_bloom.addKmer(kmer.to_string());
		kmers_parsed++;

		if (kmers_parsed % 10000000 == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << kmers_parsed << " kmers" << std::endl;
		}
	}

	return kmer_bloom;
}

template<uchar kmer_size>
std::bitset<kmer_size*2> BloomMaker<kmer_size>::random_kmer() {

	std::bitset<kmer_size*2> kmer;
    std::random_device rd;
    std::mt19937 gen( rd());
    std::bernoulli_distribution d(0.5);
    for (uint i = 0; i < kmer_size*2; ++i) {

		kmer[i] = d(gen);
    }

    return kmer;
}

template <uchar kmer_size>
void BloomMaker<kmer_size>::testbloom(BasicKmerBloom<kmer_size> & kmer_bloom, const unordered_set<bitset<kmer_size * 2> > & kmer_test_set) {

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

template <uchar kmer_size>
void BloomMaker<kmer_size>::kmc2bloomFile(const string & kmc_table_prefix, float fpr, bool test) {

    cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") makeBloom ...\n" << endl;

	auto kmc_bloom = kmc2bloom(kmc_table_prefix, fpr);
	string bloom_fn = kmc_table_prefix + ".bloom";
	std::cout << "\n[" << Utils::getLocalTime() << "] Saving bloom filter to " << bloom_fn << " ..." << std::endl;
	kmc_bloom.save(bloom_fn);

	std::cout << "[" << Utils::getLocalTime() << "] Completed saving bloom filter\n" << std::endl;
	
	if (test) {

		std::cout << "\n[" << Utils::getLocalTime() << "] Constructing test set based on the KMC kmers ..." << std::endl;
	
		CKMCFile kmc_table;

		if (!kmc_table.OpenForListing(kmc_table_prefix)) {

			cout << "ERROR: Could not open KMC table " << kmc_table_prefix << endl;
			exit(1);
		}

		uint32 db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count;
		uint64 max_kmer_count, total_kmers;

		assert(kmc_table.Info(db_kmer_size, mode, counter_size, lut_prefix_length, signature_len, min_kmer_count, max_kmer_count, total_kmers));
		assert(db_kmer_size == kmer_size);

		CKmerAPI kmer(db_kmer_size);
		uint32 count;

		unordered_set<bitset<kmer_size * 2> > kmer_test_set;

		while (kmc_table.ReadNextKmer(kmer, count)) {

			assert(kmer_test_set.insert(kmc2bits(kmer)).second);
		}

		std::cout << "[" << Utils::getLocalTime() << "] Reading bloom filter from file...\n" << std::endl;
		BasicKmerBloom<kmer_size> kmc_bloom_serialised(kmc_table_prefix + ".bloom");
		std::cout << "[" << Utils::getLocalTime() << "] Done reading bloom filter from file!\n" << std::endl;

		testbloom(kmc_bloom_serialised, kmer_test_set);


		// std::cout << "\n[" << Utils::getLocalTime() << "] Constructing test set based on random kmers ..." << std::endl;

		// kmer_test_set.clear();
		// uint num_random_tests = 10000000;

		// std::bitset<kmer_size*2> cur_kmer;

	 //    std::random_device rd;
	 //    std::mt19937 gen( rd());
	 //    std::bernoulli_distribution d(0.5);

		// BasicKmerBloom<kmer_size> random_kmer_bloom(num_random_tests, fpr);

		// while (kmer_test_set.size() < num_random_tests) {

		// 	cur_kmer<<=2;
		// 	cur_kmer[0] = d(gen);
		// 	cur_kmer[1] = d(gen);		

		// 	if (kmer_test_set.insert(cur_kmer).second) {

		// 		random_kmer_bloom.addKmer(cur_kmer);
		// 	}
		// }

		// testbloom(random_kmer_bloom, kmer_test_set);
	}
}




template class BloomMaker<31>;
template class BloomMaker<39>;
template class BloomMaker<47>;
template class BloomMaker<55>;
template class BloomMaker<63>;