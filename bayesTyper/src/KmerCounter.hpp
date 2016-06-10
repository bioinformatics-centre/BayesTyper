
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

#include "kmer_api.h"
#include "boost/graph/adjacency_list.hpp"

#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"
#include "KmerHash.hpp"
#include "PerfectSet.hpp"
#include "Sample.hpp"
#include "VariantClusterGroup.hpp"
#include "PartitionGraph.hpp"


template <uchar kmer_size>
class KmerCounter {

	public:

		KmerCounter(KmerHash * const, Utils::SmallmerSet * const, const ushort, const ushort, const uint);

		ulong countVariantClusterSmallmers(vector<VariantClusterGroup *> *);
		ulong countInterclusterKmers(vector<pair<string::const_iterator, string::const_iterator> > &, const Utils::ChromosomeClass);
		ulong countSampleKmers(const vector<Sample> &);
		ulong countVariantClusterKmers(vector<VariantClusterGroup *> *, const ushort);
		

	private:

		typedef pair<string *, string *> ReadSequences;
		typedef std::vector<std::pair<CKmerAPI,uint32> > PrecountedKmerBatch;

		KmerHashHybrid<kmer_size> * kmer_hash_hybrid;
		Utils::SmallmerSet * const smallmer_set;

		const ushort num_samples;
		const ushort num_threads;
		const uint prng_seed;

		void countVariantClusterSmallmersCallback(vector<VariantClusterGroup *> *, const ushort, mutex *, ulong *);
		void countInterclusterKmersCallback(vector<pair<string::const_iterator, string::const_iterator> >::iterator, vector<pair<string::const_iterator, string::const_iterator> >::iterator, const Utils::ChromosomeClass, mutex *, ulong *);
		void countSampleKmersCallBack(ProducerConsumerQueue<pair<uint, PrecountedKmerBatch*> > *, mutex *, uint64 *, uint64 *);
		void countVariantClusterKmersCallback(vector<VariantClusterGroup *> *, const ushort, const ushort, mutex *, ulong *);

};


#endif
