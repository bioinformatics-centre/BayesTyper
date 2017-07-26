
/*
KmerCounter.hpp - This file is part of BayesTyper (v0.9)


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

#include "kmc_api/kmer_api.h"
#include "boost/graph/adjacency_list.hpp"

#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Sample.hpp"
#include "VariantClusterGroup.hpp"
#include "VariantFileParser.hpp"


template <uchar kmer_size>
class KmerCounter {

	public:

		KmerCounter(KmerHash * const, Utils::SmallmerSet * const, const ushort);

		ulong countVariantClusterSmallmers(vector<VariantClusterGroup *> *);
		ulong countInterclusterKmers(const vector<VariantFileParser::InterClusterRegion> &);
		ulong countSampleKmers(const vector<Sample> &);
		void countVariantClusterKmers(vector<VariantClusterGroup *> *, const uint, const ushort, const ushort);

		uchar minSampleKmerCount();
		
	private:

		typedef vector<pair<CKmerAPI, uint32> > KmerBatch;

		struct KmerBatchInfo {

			const ushort sample_idx;
			const KmerBatch::iterator begin_it;
			KmerBatch::iterator end_it;

			KmerBatchInfo(const ushort sample_idx_in, const KmerBatch::iterator begin_it_in, const KmerBatch::iterator end_it_in) : sample_idx(sample_idx_in), begin_it(begin_it_in), end_it(end_it_in) {}
		};

		KmerHash * const kmer_hash;
		Utils::SmallmerSet * const smallmer_set;
		const ushort num_threads;

		uchar min_sample_kmer_count;

		void countVariantClusterSmallmersCallback(vector<VariantClusterGroup *> *, const ushort, mutex *, ulong *);
		void countInterclusterKmersCallback(const vector<VariantFileParser::InterClusterRegion> &, const ushort, mutex *, ulong *);
		void countSampleKmersCallBack(ProducerConsumerQueue<KmerBatchInfo *> *, ProducerConsumerQueue<KmerBatchInfo *> *, mutex *, ulong *);
		void countVariantClusterKmersCallback(vector<VariantClusterGroup *> *, const uint, const ushort, const ushort, const ushort);
};


#endif
