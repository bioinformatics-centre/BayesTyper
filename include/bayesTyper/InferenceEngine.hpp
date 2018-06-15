
/*
InferenceEngine.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__InferenceEngine_hpp
#define __bayesTyper__InferenceEngine_hpp

#include <vector>
#include <random> 
#include <mutex>

#include "ProducerConsumerQueue.hpp"

#include "Utils.hpp"
#include "VariantClusterGroup.hpp" 
#include "CountDistribution.hpp"
#include "CountAllocation.hpp"
#include "KmerHash.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "GenotypeWriter.hpp"
#include "OptionsContainer.hpp"
#include "Regions.hpp"
#include "InferenceUnit.hpp"
#include "Filters.hpp"


using namespace std;

class InferenceEngine {

	public:

		InferenceEngine(const vector<Sample> &, const OptionsContainer &);		

		void estimateNoiseParameters(CountDistribution *, InferenceUnit *, KmerCountsHash *, const string &);
		void genotypeVariantClusterGroups(InferenceUnit *, KmerCountsHash *, const CountDistribution &, const Filters & filters, GenotypeWriter *);	

	private: 

		const vector<Sample> samples;

		const ChromosomePloidy chrom_ploidy;

		const ushort num_threads;
		const uint prng_seed;
		const ushort num_gibbs_burn;
		const ushort num_gibbs_samples;
		const ushort num_gibbs_chains;
		const float kmer_subsampling_rate;
		const uint max_haplotype_variant_kmers;

		struct VariantClusterGroupBatch {

			uint first_variant_cluster_group_idx;
			uint num_variants;

			vector<VariantClusterGroup*>::iterator start_it; 
			vector<VariantClusterGroup*>::iterator end_it;

			VariantClusterGroupBatch() {}
			VariantClusterGroupBatch(const uint first_variant_cluster_group_idx_in, const uint num_variants_in, const vector<VariantClusterGroup*>::iterator start_it_in, const vector<VariantClusterGroup*>::iterator end_it_in) : first_variant_cluster_group_idx(first_variant_cluster_group_idx_in), num_variants(num_variants_in), start_it(start_it_in), end_it(end_it_in) {}
		};

		void initNoiseEstimationGroupsCallback(vector<VariantClusterGroup *> *, const vector<uint> &, KmerCountsHash *, const ushort, const ushort);
		void sampleNoiseCountsCallback(vector<VariantClusterGroup *> *, const vector<uint> &, CountAllocation *, mutex *, const CountDistribution &, const ushort);
		void resetNoiseEstimationGroups(vector<VariantClusterGroup *> *, const vector<uint> &);

		void genotypeVariantClusterGroupsCallback(ProducerConsumerQueue<VariantClusterGroupBatch> *, KmerCountsHash *, const CountDistribution &, const Filters & filters, GenotypeWriter *, uint *, mutex *);
};

#endif