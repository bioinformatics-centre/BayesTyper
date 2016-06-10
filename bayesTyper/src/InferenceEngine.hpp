
/*
InferenceEngine.hpp - This file is part of BayesTyper (v0.9)


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


using namespace std;

class InferenceEngine {

	public:

		InferenceEngine(const vector<Sample> &, const uchar, const OptionsContainer &);		

		void estimateNoiseParameters(CountDistribution *, vector<VariantClusterGroup*> *, KmerHash *, const vector<Sample> &);
		void genotypeVariantClusterGroups(vector<VariantClusterGroup*> *, KmerHash *, const CountDistribution &, const vector<Sample> &, GenotypeWriter *);	

	private: 

		const ushort num_samples;
		const uchar num_noise_sources;

		const ChromosomePloidy chromosome_ploidy;
		const Regions chromosome_regions;

		const uint prng_seed;
		const ushort num_threads;
		const ushort num_haplotype_candidates_per_sample;
		const ushort gibbs_burn;
		const ushort gibbs_samples;
		const ushort num_gibbs_sampling_chains;
		const ushort num_parameter_estimation_samples;
		const uint num_parameter_estimation_variants;
		const uint max_multicluster_kmers;

		struct VariantClusterGroupCounts {

			uint num_unique;	
			uint num_intercluster;
			uint num_multicluster;
			uint num_skipped;

			VariantClusterGroupCounts() : num_unique(0), num_intercluster(0), num_multicluster(0), num_skipped(0) {}

			VariantClusterGroupCounts& operator+=(const VariantClusterGroupCounts& rhs) {
			    
				this->num_unique += rhs.num_unique;
				this->num_intercluster += rhs.num_intercluster;
				this->num_multicluster += rhs.num_multicluster;
				this->num_skipped += rhs.num_skipped;

			    return *this;
			}
		};

		struct VariantClusterGroupBatch {

			uint number_of_variants;
			vector<VariantClusterGroup*>::iterator start_it; 
			vector<VariantClusterGroup*>::iterator end_it;

			VariantClusterGroupBatch() {}
			VariantClusterGroupBatch(const uint number_of_variants_in, const vector<VariantClusterGroup*>::iterator start_it_in, const vector<VariantClusterGroup*>::iterator end_it_in) : number_of_variants(number_of_variants_in), start_it(start_it_in), end_it(end_it_in) {}
		};

		void allocateShuffledIndicesToThreads(vector<vector<uint> > *, const uint);
		void selectVariantsClusterGroupsForParameterEstimationCallback(vector<VariantClusterGroup*> *, const vector<uint> &, vector<VariantClusterGroup*> *, KmerHash *, const vector<Sample> &, const uint);
		void allocateCountsForParameterEstimationCallback(vector<VariantClusterGroup*> *, const CountDistribution &, CountAllocation *, mutex *);
		vector<vector<double> > meanParameterEstimationSamples(const vector<vector<vector<double> > > &);
		
		void genotypeVariantClusterGroupsCallback(ProducerConsumerQueue<VariantClusterGroupBatch> *, KmerHash *, const CountDistribution &, const vector<Sample> &, GenotypeWriter *, VariantClusterGroupCounts *, mutex *);
};

#endif