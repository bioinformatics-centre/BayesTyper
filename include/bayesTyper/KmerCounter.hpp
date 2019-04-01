
/*
KmerCounter.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__KmerCounter_hpp
#define __bayesTyper__KmerCounter_hpp

#include <map>
#include <mutex>
#include <unordered_set>

#include "kmc_api/kmer_api.h"
#include "boost/graph/adjacency_list.hpp"

#include "KmerBloom.hpp"

#include "Utils.hpp"
#include "KmerHash.hpp"
#include "Sample.hpp"
#include "VariantClusterGraph.hpp"
#include "VariantFileParser.hpp"
#include "InferenceUnit.hpp"
#include "OptionsContainer.hpp"
#include "Chromosomes.hpp"
#include "ChromosomePloidy.hpp"


class KmerCounter {

	public:

		KmerCounter(const vector<Sample> &, const OptionsContainer &);

		void findVariantClusterPaths(InferenceUnit *, const ushort);	
		void countPathMultigroupKmers(KmerHash<bool> *, ThreadedKmerBloom<Utils::kmer_size> *, InferenceUnit *);
		void countInterclusterParameterKmers(KmerHash<bool> *, const vector<VariantFileParser::InterClusterRegion> &, const Chromosomes &, const ThreadedKmerBloom<Utils::kmer_size> &, const float);

		void countPathKmers(ThreadedKmerBloom<Utils::kmer_size> *, InferenceUnit *);
		void countInterclusterKmers(KmerCountsHash *, ThreadedKmerBloom<Utils::kmer_size> *, const string &, const Chromosomes &, const ChromosomePloidy &);
		
		void parseSampleKmers(KmerCountsHash *, ThreadedKmerBloom<Utils::kmer_size> *);
		void classifyPathKmers(KmerCountsHash *, InferenceUnit *, const string & );

	private:

		const vector<Sample> samples;

		const ushort num_threads;
		const uint prng_seed;

		typedef vector<pair<CKmerAPI, uint32> > KmerBatch;

		struct KmerBatchInfo {

			const ushort sample_idx;
			const KmerBatch::iterator begin_it;
			KmerBatch::iterator end_it;

			KmerBatchInfo(const ushort sample_idx_in, const KmerBatch::iterator begin_it_in, const KmerBatch::iterator end_it_in) : sample_idx(sample_idx_in), begin_it(begin_it_in), end_it(end_it_in) {}
		};

		void findVariantClusterPathsCallback(vector<VariantClusterGroup *> *, KmerBloom<Utils::kmer_size> *, const ushort, const ushort, const ushort);
		void countPathMultigroupKmersCallback(KmerHash<bool> *, ThreadedKmerBloom<Utils::kmer_size> *, vector<VariantClusterGroup *> *, mutex *, ulong *, const ushort);
		void countInterclusterParameterKmersCallback(KmerHash<bool> *, const vector<VariantFileParser::InterClusterRegion> &, const Chromosomes &, const ThreadedKmerBloom<Utils::kmer_size> &, const float, const ushort);

		void countPathKmersCallback(ThreadedKmerBloom<Utils::kmer_size> *, vector<VariantClusterGroup *> *, const ushort);
		void countInterclusterKmersCallback(KmerCountsHash *, ThreadedKmerBloom<Utils::kmer_size> *, const vector<VariantFileParser::InterClusterRegion> &, const Chromosomes &, const ChromosomePloidy &, const ushort);

		void parseSampleKmersCallBack(KmerCountsHash *, ThreadedKmerBloom<Utils::kmer_size> *, ProducerConsumerQueue<KmerBatchInfo *> *, ProducerConsumerQueue<KmerBatchInfo *> *);
		void classifyPathKmersCallback(KmerCountsHash *, KmerBloom<Utils::kmer_size> *, vector<VariantClusterGroup *> *, const ushort);
};


#endif
