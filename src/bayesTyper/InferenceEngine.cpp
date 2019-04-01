
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


static const uint noise_variants_batch_size = 100000;
static const uint variant_cluster_groups_batch_size = 1000;

static const double genotyping_stdout_frequency = 100000;

static const uchar num_genomic_rate_gc_bias_bins = 1;


InferenceEngine::InferenceEngine(const vector<Sample> & samples_in, const ChromosomePloidy & chrom_ploidy_in, const OptionsContainer & options_container) : samples(samples_in), chrom_ploidy(chrom_ploidy_in), num_threads(options_container.getValue<ushort>("threads")), prng_seed(options_container.getValue<uint>("random-seed")), num_gibbs_burn(options_container.getValue<ushort>("gibbs-burn-in")), num_gibbs_samples(options_container.getValue<ushort>("gibbs-samples")), num_gibbs_chains(options_container.getValue<ushort>("number-of-gibbs-chains")), kmer_subsampling_rate(options_container.getValue<float>("kmer-subsampling-rate")), max_haplotype_variant_kmers(options_container.getValue<uint>("max-haplotype-variant-kmers")) {}

void InferenceEngine::initGenotypersCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & group_indices, const uint group_indices_idx_end, KmerCountsHash * kmer_hash, const ushort gibbs_chain_idx, const ushort thread_idx) {

    uint group_indices_idx = thread_idx;

	while (group_indices_idx < group_indices_idx_end) {

		const uint variant_cluster_group_idx = group_indices.at(group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));

		variant_cluster_groups->at(variant_cluster_group_idx)->initGenotyper(kmer_hash, samples, prng_seed + (variant_cluster_group_idx + 1) * (gibbs_chain_idx + 1), num_genomic_rate_gc_bias_bins, kmer_subsampling_rate, max_haplotype_variant_kmers);
 		variant_cluster_groups->at(variant_cluster_group_idx)->shuffleBranchOrdering(prng_seed + (variant_cluster_group_idx + 1) * (gibbs_chain_idx + 1));

        group_indices_idx += num_threads;
	}
}

void InferenceEngine::sampleGenotypesCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & group_indices, const uint group_indices_idx_end, const CountDistribution & count_distribution, const bool collect_samples, CountAllocation * noise_counts_global, mutex * counts_mutex, const ushort thread_idx) {

	CountAllocation noise_counts_local(samples.size());

    uint group_indices_idx = thread_idx;

	while (group_indices_idx < group_indices_idx_end) {

		const uint variant_cluster_group_idx = group_indices.at(group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));

		variant_cluster_groups->at(variant_cluster_group_idx)->estimateGenotypes(count_distribution, chrom_ploidy, collect_samples);
		variant_cluster_groups->at(variant_cluster_group_idx)->getNoiseCounts(&noise_counts_local, count_distribution);	

		variant_cluster_groups->at(variant_cluster_group_idx)->clearGenotyperCache();
        group_indices_idx += num_threads;
	}	

	lock_guard<mutex> counts_locker(*counts_mutex);
	noise_counts_global->mergeInCountAllocations(noise_counts_local);	
}

void InferenceEngine::resetGroupsCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & group_indices, const uint group_indices_idx_end, const ushort thread_idx) {

    uint group_indices_idx = thread_idx;

	while (group_indices_idx < group_indices_idx_end) {

		const uint variant_cluster_group_idx = group_indices.at(group_indices_idx);

		assert(variant_cluster_groups->at(variant_cluster_group_idx));
		variant_cluster_groups->at(variant_cluster_group_idx)->resetGroup();
        
        group_indices_idx += num_threads;
	}
}

void InferenceEngine::collectGenotypesCallback(vector<VariantClusterGroup *> * variant_cluster_groups, const vector<uint> & group_indices, const uint group_indices_idx_end, const Filters & filters, GenotypeWriter * genotype_writer, const ushort thread_idx) {

	vector<Genotypes*> * variant_genotypes = new vector<Genotypes*>();

    uint group_indices_idx = thread_idx;

	while (group_indices_idx < group_indices_idx_end) {

		const uint variant_cluster_group_idx = group_indices.at(group_indices_idx);
		assert(variant_cluster_groups->at(variant_cluster_group_idx));

		variant_cluster_groups->at(variant_cluster_group_idx)->collectGenotypes(variant_genotypes, chrom_ploidy, filters);   
        delete variant_cluster_groups->at(variant_cluster_group_idx);	

        group_indices_idx += num_threads;
	}

	genotype_writer->addGenotypes(variant_genotypes);
}

void InferenceEngine::estimateNoise(CountDistribution * count_distribution, InferenceUnit * inference_unit, KmerCountsHash * kmer_hash, const string & output_prefix) {

	cout << "[" << Utils::getLocalTime() << "] Estimating noise model parameters using " << num_gibbs_chains << " independent gibbs sampling chains each with " << num_gibbs_burn + num_gibbs_samples << " iterations (" << num_gibbs_burn << " burn-in) ..." << endl;

	uint num_noise_variants = 0;

	vector<uint> noise_group_indices;
	noise_group_indices.reserve(inference_unit->variant_cluster_groups.size());

	for (uint variant_cluster_group_idx = 0; variant_cluster_group_idx < inference_unit->variant_cluster_groups.size(); variant_cluster_group_idx++) {

		if (inference_unit->variant_cluster_groups.at(variant_cluster_group_idx)->numberOfVariantClusters() == 1) {

			num_noise_variants += inference_unit->variant_cluster_groups.at(variant_cluster_group_idx)->numberOfVariants();
			noise_group_indices.push_back(variant_cluster_group_idx);
		}		
	}

	noise_group_indices.shrink_to_fit();

	vector<double> mean_noise_rates(samples.size(), 0);

    ofstream noise_outfile(output_prefix + ".txt");

    if (!noise_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << output_prefix + ".txt" << "\n" << endl;
        exit(1);
    }

    noise_outfile << "Chain\tIteration";

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

	  	noise_outfile << "\t" << samples.at(sample_idx).name;
    }

    noise_outfile << endl;

    mt19937 prng = mt19937(prng_seed);

	for (ushort gibbs_chain_idx = 0; gibbs_chain_idx < num_gibbs_chains; gibbs_chain_idx++) {

	    shuffle(noise_group_indices.begin(), noise_group_indices.end(), prng);
	    
	    uint noise_group_indices_idx_end = 0;
	    num_noise_variants = 0;

		while ((num_noise_variants < noise_variants_batch_size) and (noise_group_indices_idx_end < noise_group_indices.size())) {

			num_noise_variants += inference_unit->variant_cluster_groups.at(noise_group_indices.at(noise_group_indices_idx_end))->numberOfVariants();
			noise_group_indices_idx_end++;	
		}

		sort(noise_group_indices.begin(), noise_group_indices.begin() += noise_group_indices_idx_end);

	    vector<thread> initing_threads;
	    initing_threads.reserve(num_threads);

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

		    initing_threads.push_back(thread(&InferenceEngine::initGenotypersCallback, this, &(inference_unit->variant_cluster_groups), ref(noise_group_indices), noise_group_indices_idx_end, kmer_hash, gibbs_chain_idx, thread_idx));
	    }

		for (auto & initing_thread: initing_threads) {
	    	
	    	initing_thread.join();
		}

	    assert(count_distribution->getNoiseRates().size() == samples.size());
	    noise_outfile << gibbs_chain_idx + 1 << "\t0\t" << count_distribution->getNoiseRates() << endl;

		mutex counts_mutex;

		for (uint iteration = 1; iteration <= (num_gibbs_burn + num_gibbs_samples); iteration++) {
		
			CountAllocation noise_counts(samples.size());

			vector<thread> sampling_threads;
			sampling_threads.reserve(num_threads);

			for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	    	    sampling_threads.push_back(thread(&InferenceEngine::sampleGenotypesCallback, this, &(inference_unit->variant_cluster_groups), ref(noise_group_indices), noise_group_indices_idx_end, ref(*count_distribution), false, &noise_counts, &counts_mutex, thread_idx));
		    }

	    	for (auto & sampling_thread: sampling_threads) {
	        	
	        	sampling_thread.join();
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

	    vector<thread> reset_threads;
	    reset_threads.reserve(num_threads);

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

		    reset_threads.push_back(thread(&InferenceEngine::resetGroupsCallback, this, &(inference_unit->variant_cluster_groups), ref(noise_group_indices), noise_group_indices_idx_end, thread_idx));
	    }

		for (auto & reset_thread: reset_threads) {
	    	
	    	reset_thread.join();
		}

		count_distribution->resetNoiseRates();
	}

	cout << "[" << Utils::getLocalTime() << "] Calculated final noise model parameters by averaging " << num_gibbs_samples * num_gibbs_chains << " parameter estimates (" << num_gibbs_samples << " per gibbs sampling chain)" << endl;


    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

	  	mean_noise_rates.at(sample_idx) /= num_gibbs_samples * num_gibbs_chains;
    }

	count_distribution->setNoiseRates(mean_noise_rates);

	noise_outfile << "0\t0\t" << count_distribution->getNoiseRates() << endl;
	noise_outfile.close();

	if (num_noise_variants < noise_variants_batch_size) {

		cout << "\nWARNING: Low number of variants used for Poisson parameter estimation (" << num_noise_variants << " < " << noise_variants_batch_size << ")" << endl;
		cout << "WARNING: Try using the noise genotyping mode (--noise-genotyping) if noise rate estimates are not stable between chains" << endl;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Wrote noise parameters to " << output_prefix << ".txt" << endl;
}

void InferenceEngine::estimateGenotypesCallback(ProducerConsumerQueue<VariantClusterGroupBatch> * variant_cluster_group_batch_queue, KmerCountsHash * kmer_hash, const CountDistribution & count_distribution, const Filters & filters, GenotypeWriter * genotype_writer, uint * num_genotyped_variants, mutex * stdout_mutex) {

	VariantClusterGroupBatch variant_cluster_group_batch;

	while (variant_cluster_group_batch_queue->pop(&variant_cluster_group_batch)) {

		vector<Genotypes*> * variant_genotypes = new vector<Genotypes*>();
		variant_genotypes->reserve(variant_cluster_group_batch.num_variants);

		uint variant_cluster_group_idx = variant_cluster_group_batch.first_variant_cluster_group_idx;
		auto variant_cluster_group_it = variant_cluster_group_batch.start_it;

		while (variant_cluster_group_it != variant_cluster_group_batch.end_it) {
			
			for (ushort gibbs_chain_idx = 0; gibbs_chain_idx < num_gibbs_chains; gibbs_chain_idx++) {

				(*variant_cluster_group_it)->initGenotyper(kmer_hash, samples, prng_seed + (variant_cluster_group_idx + 1), num_genomic_rate_gc_bias_bins, kmer_subsampling_rate, max_haplotype_variant_kmers);
				(*variant_cluster_group_it)->shuffleBranchOrdering(prng_seed + (variant_cluster_group_idx + 1) * (gibbs_chain_idx + 1));

				for (ushort i = 0; i < num_gibbs_burn; i++) {

					(*variant_cluster_group_it)->estimateGenotypes(count_distribution, chrom_ploidy, false);
				}

				for (ushort i = 0; i < num_gibbs_samples; i++) {

					(*variant_cluster_group_it)->estimateGenotypes(count_distribution, chrom_ploidy, true);
				}					
			}

			(*variant_cluster_group_it)->collectGenotypes(variant_genotypes, chrom_ploidy, filters);	
			delete *variant_cluster_group_it;			

			variant_cluster_group_idx++;
			variant_cluster_group_it++;
		}

		stdout_mutex->lock();

		uint genotyping_stdout_idx_1 = floor(*num_genotyped_variants / genotyping_stdout_frequency);
		*num_genotyped_variants += variant_genotypes->size();

		uint genotyping_stdout_idx_2 = floor(*num_genotyped_variants / genotyping_stdout_frequency);
		assert(genotyping_stdout_idx_1 <= genotyping_stdout_idx_2);

		while (genotyping_stdout_idx_1 < genotyping_stdout_idx_2) {

			cout << "[" << Utils::getLocalTime() << "] Genotyped " << static_cast<uint>(genotyping_stdout_frequency * (genotyping_stdout_idx_1 + 1)) << " variants" << endl;	
			genotyping_stdout_idx_1++;
		}

		stdout_mutex->unlock();

		genotype_writer->addGenotypes(variant_genotypes);
	}
}

void InferenceEngine::estimateGenotypes(InferenceUnit * inference_unit, KmerCountsHash * kmer_hash, const CountDistribution & count_distribution, const Filters & filters, GenotypeWriter * genotype_writer) {

	cout << "[" << Utils::getLocalTime() << "] Estimating genotypes using " << num_gibbs_chains << " independent gibbs sampling chains each with " << num_gibbs_burn + num_gibbs_samples << " iterations (" << num_gibbs_burn << " burn-in) ...\n" << endl;

    ProducerConsumerQueue<VariantClusterGroupBatch> variant_cluster_group_batch_queue(Utils::queue_size_thread_scaling * num_threads);

    uint num_genotyped_variants = 0;
    mutex stdout_mutex;

	vector<thread> genotyping_threads;
	genotyping_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        genotyping_threads.push_back(thread(&InferenceEngine::estimateGenotypesCallback, this, &variant_cluster_group_batch_queue, kmer_hash, ref(count_distribution), ref(filters), genotype_writer, &num_genotyped_variants, &stdout_mutex));
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
}

void InferenceEngine::estimateNoiseAndGenotypes(InferenceUnit * inference_unit, CountDistribution * count_distribution, KmerCountsHash * kmer_hash, const Filters & filters, GenotypeWriter * genotype_writer, const string & output_prefix) {

	cout << "[" << Utils::getLocalTime() << "] Estimating noise model parameters and genotypes using " << num_gibbs_chains << " independent gibbs sampling chains each with " << num_gibbs_burn + num_gibbs_samples << " iterations (" << num_gibbs_burn << " burn-in) ...\n" << endl;

    ofstream noise_outfile(output_prefix + ".txt");

    if (!noise_outfile.is_open()) {

        cerr << "\nERROR: Unable to write file " << output_prefix + ".txt" << "\n" << endl;
        exit(1);
    }

    noise_outfile << "Chain\tIteration";

    for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

	  	noise_outfile << "\t" << samples.at(sample_idx).name;
    }

    noise_outfile << endl;

    vector<uint> group_indices(inference_unit->variant_cluster_groups.size());
	iota(group_indices.begin(), group_indices.end(), 0);

	for (ushort gibbs_chain_idx = 0; gibbs_chain_idx < num_gibbs_chains; gibbs_chain_idx++) {

	    vector<thread> initing_threads;
	    initing_threads.reserve(num_threads);

		for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

		    initing_threads.push_back(thread(&InferenceEngine::initGenotypersCallback, this, &(inference_unit->variant_cluster_groups), ref(group_indices), group_indices.size(), kmer_hash, gibbs_chain_idx, thread_idx));
	    }

		for (auto & initing_thread: initing_threads) {
	    	
	    	initing_thread.join();
		}

	    assert(count_distribution->getNoiseRates().size() == samples.size());
	    noise_outfile << gibbs_chain_idx + 1 << "\t0\t" << count_distribution->getNoiseRates() << endl;

		mutex counts_mutex;

		for (uint iteration = 1; iteration <= (num_gibbs_burn + num_gibbs_samples); iteration++) {
		
			CountAllocation noise_counts(samples.size());

			vector<thread> sampling_threads;
			sampling_threads.reserve(num_threads);

			for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	    	    sampling_threads.push_back(thread(&InferenceEngine::sampleGenotypesCallback, this, &(inference_unit->variant_cluster_groups), ref(group_indices), group_indices.size(), ref(*count_distribution), iteration > num_gibbs_burn, &noise_counts, &counts_mutex, thread_idx));
		    }

	    	for (auto & sampling_thread: sampling_threads) {
	        	
	        	sampling_thread.join();
			}

			count_distribution->sampleNoiseParameters(noise_counts);

		    assert(count_distribution->getNoiseRates().size() == samples.size());
	    	noise_outfile << gibbs_chain_idx + 1 << "\t" << iteration << "\t" << count_distribution->getNoiseRates() << endl;
		}

		count_distribution->resetNoiseRates();
		cout << "[" << Utils::getLocalTime() << "] Finished " << gibbs_chain_idx + 1 << " gibbs sampling chain" << (gibbs_chain_idx > 0 ? "s" : "") << endl;
	}

	noise_outfile.close();

	cout << "\n[" << Utils::getLocalTime() << "] Wrote noise parameters to " << output_prefix << ".txt" << endl;
	cout << "\n[" << Utils::getLocalTime() << "] Collecting genotypes ..." << endl;

	vector<thread> collection_threads;
	collection_threads.reserve(num_threads);

	for (ushort thread_idx = 0; thread_idx < num_threads; thread_idx++) {

	    collection_threads.push_back(thread(&InferenceEngine::collectGenotypesCallback, this, &(inference_unit->variant_cluster_groups), ref(group_indices), group_indices.size(), ref(filters), genotype_writer, thread_idx));
    }

	for (auto & collection_thread: collection_threads) {
    	
    	collection_thread.join();
	}
}




