
/*
InferenceEngine.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <vector>
#include <thread>
#include <mutex>
#include <random>

#include "ProducerConsumerQueue.hpp"

#include "InferenceEngine.hpp"
#include "Utils.hpp"
#include "VariantClusterGroup.hpp" 
#include "CountDistribution.hpp"
#include "Genotypes.hpp"
#include "GenotypeWriter.hpp"
#include "KmerHash.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "OptionsContainer.hpp"
#include "Regions.hpp"


static const uint noise_estimation_batch_size = 100000;

static const uint variant_cluster_groups_batch_size = 1000;
static const double genotyping_stdout_frequency = 100000;

static const uchar num_genomic_rate_gc_bias_bins = 1;


InferenceEngine::InferenceEngine(const vector<Sample> & samples_in, const OptionsContainer & options_container) : samples(samples_in), chrom_ploidy(samples), num_threads(options_container.getValue<ushort>("threads")), prng_seed(options_container.getValue<uint>("random-seed")), num_gibbs_burn(options_container.getValue<ushort>("gibbs-burn-in")), num_gibbs_samples(options_container.getValue<ushort>("gibbs-samples")), num_gibbs_chains(options_container.getValue<ushort>("number-of-gibbs-chains")), kmer_subsampling_rate(options_container.getValue<float>("kmer-subsampling-rate")), max_haplotype_variant_kmers(options_container.getValue<uint>("max-haplotype-variant-kmers")) {}


void InferenceEngine::initNoiseEstimationGroupsCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & noise_estimation_group_indices, KmerCountsHash * kmer_hash, const ushort gibbs_chain_idx, const ushort thread_idx) {

    const uint num_noise_estimation_groups = min(noise_estimation_batch_size, static_cast<uint>(noise_estimation_group_indices.size()));

    uint noise_estimation_group_indices_idx = thread_idx;

	while (noise_estimation_group_indices_idx < num_noise_estimation_groups) {

		const uint variant_cluster_group_idx = noise_estimation_group_indices.at(noise_estimation_group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));

		assert(variant_cluster_groups->at(variant_cluster_group_idx)->numberOfVariants() == 1);
		assert(variant_cluster_groups->at(variant_cluster_group_idx)->numberOfVariantClusters() == 1);
		assert(variant_cluster_groups->at(variant_cluster_group_idx)->numberOfVariantClusterGroupTrees() == 1);

		variant_cluster_groups->at(variant_cluster_group_idx)->initGenotyper(kmer_hash, samples, prng_seed * (gibbs_chain_idx + 1) + variant_cluster_group_idx, num_genomic_rate_gc_bias_bins, kmer_subsampling_rate, max_haplotype_variant_kmers);
        
        noise_estimation_group_indices_idx += num_threads;
	}
}

void InferenceEngine::sampleNoiseCountsCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & noise_estimation_group_indices, CountAllocation * noise_counts_global, mutex * noise_counts_lock, const CountDistribution & count_distribution, const ushort thread_idx) {

	CountAllocation noise_counts_local(samples.size());

    const uint num_noise_estimation_groups = min(noise_estimation_batch_size, static_cast<uint>(noise_estimation_group_indices.size()));

    uint noise_estimation_group_indices_idx = thread_idx;

	while (noise_estimation_group_indices_idx < num_noise_estimation_groups) {

		const uint variant_cluster_group_idx = noise_estimation_group_indices.at(noise_estimation_group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));

		variant_cluster_groups->at(variant_cluster_group_idx)->estimateGenotypes(count_distribution, chrom_ploidy, false);
		variant_cluster_groups->at(variant_cluster_group_idx)->getNoiseCounts(&noise_counts_local, count_distribution);	
        
        noise_estimation_group_indices_idx += num_threads;
	}
	
	lock_guard<mutex> counts_locker(*noise_counts_lock);
	noise_counts_global->mergeInCountAllocations(noise_counts_local);	
}

void InferenceEngine::resetNoiseEstimationGroups(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & noise_estimation_group_indices) {

    const uint num_noise_estimation_groups = min(noise_estimation_batch_size, static_cast<uint>(noise_estimation_group_indices.size()));

    uint noise_estimation_group_indices_idx = 0;

	while (noise_estimation_group_indices_idx < num_noise_estimation_groups) {

		const uint variant_cluster_group_idx = noise_estimation_group_indices.at(noise_estimation_group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));
		variant_cluster_groups->at(variant_cluster_group_idx)->resetGroup();
        
        noise_estimation_group_indices_idx++;;
	}
}

void InferenceEngine::estimateNoiseParameters(CountDistribution * count_distribution, InferenceUnit * inference_unit, KmerCountsHash * kmer_hash, const string & output_prefix) {

	cout << "[" << Utils::getLocalTime() << "] Estimating noise model parameters using " << num_gibbs_chains << " parallel gibbs sampling chains each with " << num_gibbs_burn + num_gibbs_samples << " iterations (" << num_gibbs_burn << " burn-in) ..." << endl;

	vector<uint> noise_estimation_group_indices;
	noise_estimation_group_indices.reserve(inference_unit->variant_cluster_groups.size());

	for (uint variant_cluster_group_idx = 0; variant_cluster_group_idx < inference_unit->variant_cluster_groups.size(); variant_cluster_group_idx++) {

		if (inference_unit->variant_cluster_groups.at(variant_cluster_group_idx)->isAutosomalSimpleSNV()) {
		
			noise_estimation_group_indices.push_back(variant_cluster_group_idx);
		}		
	}

	noise_estimation_group_indices.shrink_to_fit();

	vector<double> mean_noise_rates(samples.size(), 0);

    ofstream noise_outfile(output_prefix + ".txt");
    assert(noise_outfile.is_open());

    noise_outfile << "Chain\tIteration";

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

	  	noise_outfile << "\t" << samples.at(sample_idx).name;
    }

    noise_outfile << endl;

    mt19937 prng = mt19937(prng_seed);

	for (ushort gibbs_chain_idx = 0; gibbs_chain_idx < num_gibbs_chains; gibbs_chain_idx++) {

	    shuffle(noise_estimation_group_indices.begin(), noise_estimation_group_indices.end(), prng);

	    vector<thread> initing_threads;
	    initing_threads.reserve(num_threads);

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

		    initing_threads.push_back(thread(&InferenceEngine::initNoiseEstimationGroupsCallback, this, &(inference_unit->variant_cluster_groups), ref(noise_estimation_group_indices), kmer_hash, gibbs_chain_idx, thread_idx));
	    }

		for (auto & initing_thread: initing_threads) {
	    	
	    	initing_thread.join();
		}

	    assert(count_distribution->getNoiseRates().size() == samples.size());
	    noise_outfile << gibbs_chain_idx + 1 << "\t0\t" << count_distribution->getNoiseRates() << endl;

		mutex noise_counts_lock;

		for (uint iteration = 1; iteration <= (num_gibbs_burn + num_gibbs_samples); iteration++) {
		
			CountAllocation noise_counts(samples.size());

			vector<thread> estimation_threads;
			estimation_threads.reserve(num_threads);

			for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	    	    estimation_threads.push_back(thread(&InferenceEngine::sampleNoiseCountsCallback, this, &(inference_unit->variant_cluster_groups), ref(noise_estimation_group_indices), &noise_counts, &noise_counts_lock, ref(*count_distribution), thread_idx));
		    }

	    	for (auto & estimation_thread: estimation_threads) {
	        	
	        	estimation_thread.join();
			}

			count_distribution->sampleNoiseParameters(noise_counts);

		    assert(count_distribution->getNoiseRates().size() == samples.size());
	    	noise_outfile << gibbs_chain_idx + 1 << "\t" << iteration << "\t" << count_distribution->getNoiseRates() << endl;

	    	if (num_gibbs_burn < iteration) {

			    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

				  	mean_noise_rates.at(sample_idx) += count_distribution->getNoiseRates().at(sample_idx);
			    }
			}
		}

		resetNoiseEstimationGroups(&(inference_unit->variant_cluster_groups), noise_estimation_group_indices);
		count_distribution->resetNoiseRates();
	}

	cout << "[" << Utils::getLocalTime() << "] Calculated final noise model parameters by averaging " << num_gibbs_samples * num_gibbs_chains << " parameter estimates (" << num_gibbs_samples << " per gibbs sampling chain)" << endl;

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

	  	mean_noise_rates.at(sample_idx) /= (num_gibbs_samples * num_gibbs_chains);
    }

	count_distribution->setNoiseRates(mean_noise_rates);

	noise_outfile << "0\t0\t" << count_distribution->getNoiseRates() << endl;
	noise_outfile.close();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote parameters to " << output_prefix << ".txt" << endl;
}

void InferenceEngine::genotypeVariantClusterGroupsCallback(ProducerConsumerQueue<VariantClusterGroupBatch> * variant_cluster_group_batch_queue, KmerCountsHash * kmer_hash, const CountDistribution & count_distribution, const Filters & filters, GenotypeWriter * genotype_writer, uint * num_genotyped_variants, mutex * thread_lock) {

	assert(samples.size() == samples.size());

	VariantClusterGroupBatch variant_cluster_group_batch;

	while (variant_cluster_group_batch_queue->pop(&variant_cluster_group_batch)) {

		vector<Genotypes*> * variant_genotypes = new vector<Genotypes*>();
		variant_genotypes->reserve(variant_cluster_group_batch.num_variants);

		uint variant_cluster_group_idx = variant_cluster_group_batch.first_variant_cluster_group_idx;

		auto variant_cluster_group_it = variant_cluster_group_batch.start_it;

		while (variant_cluster_group_it != variant_cluster_group_batch.end_it) {

			(*variant_cluster_group_it)->initGenotyper(kmer_hash, samples, prng_seed + variant_cluster_group_idx, num_genomic_rate_gc_bias_bins, kmer_subsampling_rate, max_haplotype_variant_kmers);
			
			for (ushort gibbs_chain_idx = 0; gibbs_chain_idx < num_gibbs_chains; gibbs_chain_idx++) {

				if (gibbs_chain_idx > 0) {

					(*variant_cluster_group_it)->resetGenotyper(kmer_subsampling_rate, max_haplotype_variant_kmers);
				}

				(*variant_cluster_group_it)->shuffleBranchOrdering(prng_seed * (gibbs_chain_idx + 1) + variant_cluster_group_idx);

				for (ushort i = 0; i < num_gibbs_burn; i++) {

					(*variant_cluster_group_it)->estimateGenotypes(count_distribution, chrom_ploidy, false);
				}

				for (ushort i = 0; i < num_gibbs_samples; i++) {

					(*variant_cluster_group_it)->estimateGenotypes(count_distribution, chrom_ploidy, true);
				}					
			}

			(*variant_cluster_group_it)->collectGenotypes(variant_genotypes, chrom_ploidy, filters);

			variant_cluster_group_idx++;
	
			delete *variant_cluster_group_it;			
			variant_cluster_group_it++;
		}

		thread_lock->lock();

		uint genotyping_stdout_idx_1 = floor(*num_genotyped_variants / genotyping_stdout_frequency);
		*num_genotyped_variants += variant_genotypes->size();

		uint genotyping_stdout_idx_2 = floor(*num_genotyped_variants / genotyping_stdout_frequency);
		assert(genotyping_stdout_idx_1 <= genotyping_stdout_idx_2);

		while (genotyping_stdout_idx_1 < genotyping_stdout_idx_2) {

			cout << "[" << Utils::getLocalTime() << "] Genotyped " << static_cast<uint>(genotyping_stdout_frequency * (genotyping_stdout_idx_1 + 1)) << " variants" << endl;	
			genotyping_stdout_idx_1++;
		}

		thread_lock->unlock();

		genotype_writer->addGenotypes(variant_genotypes);
	}
}

void InferenceEngine::genotypeVariantClusterGroups(InferenceUnit * inference_unit, KmerCountsHash * kmer_hash, const CountDistribution & count_distribution, const Filters & filters, GenotypeWriter * genotype_writer) {

	cout << "[" << Utils::getLocalTime() << "] Running " << num_gibbs_chains << " parallel gibbs sampling chains each with " << num_gibbs_burn + num_gibbs_samples << " iterations (" << num_gibbs_burn << " burn-in) on " << inference_unit->num_variants << " variants ...\n" << endl;

    ProducerConsumerQueue<VariantClusterGroupBatch> variant_cluster_group_batch_queue(Utils::queue_size_thread_scaling * num_threads);

    uint num_genotyped_variants = 0;
    mutex thread_lock;

	vector<thread> genotyping_threads;
	genotyping_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        genotyping_threads.push_back(thread(&InferenceEngine::genotypeVariantClusterGroupsCallback, this, &variant_cluster_group_batch_queue, kmer_hash, ref(count_distribution), ref(filters), genotype_writer, &num_genotyped_variants, &thread_lock));
    }  

    auto variant_cluster_group_it = inference_unit->variant_cluster_groups.begin();
    auto first_variant_cluster_group_it = variant_cluster_group_it;

	uint num_batch_variants = 0;
	
	while (variant_cluster_group_it != inference_unit->variant_cluster_groups.end()) {

		assert(*variant_cluster_group_it);
		num_batch_variants += (*variant_cluster_group_it)->numberOfVariants();

		variant_cluster_group_it++;

		if (num_batch_variants >= variant_cluster_groups_batch_size) {

			variant_cluster_group_batch_queue.push(VariantClusterGroupBatch(first_variant_cluster_group_it - inference_unit->variant_cluster_groups.begin(), num_batch_variants, first_variant_cluster_group_it, variant_cluster_group_it));	
			
			first_variant_cluster_group_it = variant_cluster_group_it;
			num_batch_variants = 0;
		}
 	}

 	variant_cluster_group_batch_queue.push(VariantClusterGroupBatch(first_variant_cluster_group_it - inference_unit->variant_cluster_groups.begin(), num_batch_variants, first_variant_cluster_group_it, variant_cluster_group_it));
	variant_cluster_group_batch_queue.pushedLast();

	for(auto & genotyping_thread: genotyping_threads) {
    	
    	genotyping_thread.join();
	}

	assert(num_genotyped_variants <= inference_unit->num_variants);

	cout << "\n[" << Utils::getLocalTime() << "] Out of " << inference_unit->num_variants << " variants:\n" << endl; 
    cout << "\t- " << num_genotyped_variants << " were genotyped" << endl;
    cout << "\t- " << inference_unit->num_variants - num_genotyped_variants << " were skipped (unsupported)" << endl;
}


