
/*
VariantClusterGenotyper.cpp - This file is part of BayesTyper (v1.1)


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


#include <string>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <assert.h>
#include <math.h>
#include <algorithm>

#include "boost/functional/hash.hpp"

#include "VariantClusterGenotyper.hpp"
#include "Utils.hpp"
#include "VariantClusterGraph.hpp"
#include "KmerHash.hpp"
#include "KmerCounts.hpp"
#include "CountDistribution.hpp"
#include "SparsityEstimator.hpp"
#include "Genotypes.hpp"
#include "FrequencyDistribution.hpp"
#include "HaplotypeFrequencyDistribution.hpp"
#include "VariantClusterHaplotypes.hpp"
#include "ChromosomePloidy.hpp"
#include "Sample.hpp"
#include "KmerStats.hpp"
#include "VariantInfo.hpp"

using namespace std;

VariantClusterGenotyper::VariantClusterGenotyper(const vector<Sample> & samples_in, const float kmer_subsampling_rate) : samples(samples_in), variant_cluster_haplotypes(kmer_subsampling_rate) {

	is_parameter_estimation_cluster = false;
	use_multicluster_kmers = false;

	has_redundant_sequence = false;	
	has_excluded_kmer = false;
}


VariantClusterGenotyper::~VariantClusterGenotyper() {

	delete haplotype_frequency_distribution;
}


void VariantClusterGenotyper::initialise(KmerHash * kmer_hash, VariantClusterGraph * variant_cluster_graph, const uint prng_seed, const ushort max_sample_haplotype_candidates, const uchar num_genomic_rate_gc_bias_bins) {

	prng = mt19937(prng_seed);

	unique_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());
	multicluster_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());

	diplotype_sample = vector<pair<ushort,ushort> >(samples.size(), make_pair(Utils::ushort_overflow, Utils::ushort_overflow));

	variant_cluster_graph->getBestHaplotypeCandidates(kmer_hash, &variant_cluster_haplotypes, samples.size(), max_sample_haplotype_candidates, num_genomic_rate_gc_bias_bins);
	
	variant_cluster_info = variant_cluster_graph->getInfo();
	allele_kmer_stats.reserve(variant_cluster_info.size());

	for (auto & variant_info: variant_cluster_info) {

		allele_kmer_stats.emplace_back(vector<AlleleKmerStats>(samples.size(), AlleleKmerStats(variant_info.num_alleles)));
	}

	has_redundant_sequence = variant_cluster_graph->hasRedundantSequence();
	has_excluded_kmer = variant_cluster_graph->hasExcludedKmer();

	assert(variant_cluster_haplotypes.unique_kmers.size() == uint(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.rows()));
    assert(variant_cluster_haplotypes.multicluster_kmers.size() == uint(variant_cluster_haplotypes.haplotype_multicluster_kmer_multiplicities.rows()));

	assert(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols() > 0);
	assert(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols() < Utils::ushort_overflow);

	assert(variant_cluster_haplotypes.haplotype_multicluster_kmer_multiplicities.cols() > 0);
	assert(variant_cluster_haplotypes.haplotype_multicluster_kmer_multiplicities.cols() < Utils::ushort_overflow);

	assert(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols() == variant_cluster_haplotypes.haplotype_multicluster_kmer_multiplicities.cols());
	assert(variant_cluster_haplotypes.haplotypes.size() == uint(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols()));

	assert(variant_cluster_haplotypes.unique_kmer_subset_indices.empty());
	assert(variant_cluster_haplotypes.multicluster_kmer_subset_indices.empty());

    variant_cluster_haplotypes.kmer_stats_cache = vector<VariantClusterHaplotypes::KmerStatsCache>(samples.size(), variant_cluster_info.size());

	Eigen::RowVectorXbool non_zero_unqiue_kmer_counts;
	non_zero_unqiue_kmer_counts = Eigen::RowVectorXbool(variant_cluster_haplotypes.unique_kmers.size());

	for (uint kmer_idx = 0; kmer_idx < variant_cluster_haplotypes.unique_kmers.size(); kmer_idx++) {

		if (variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts) {

			assert(!(variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts->hasMulticlusterOccurrence()));

			if (variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts->getInterclusterMultiplicity(Utils::Gender::Male) > 0) {

				non_zero_unqiue_kmer_counts(kmer_idx) = false;
		
			} else {

				non_zero_unqiue_kmer_counts(kmer_idx) = true;			
			} 
		
		} else {

			non_zero_unqiue_kmer_counts(kmer_idx) = false;
		}
	}
	
	for (uint kmer_idx = 0; kmer_idx < variant_cluster_haplotypes.multicluster_kmers.size(); kmer_idx++) {

		assert(variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).counts);
		assert(variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).counts->hasMulticlusterOccurrence());
	}
		
	SparsityEstimator sparsity_estimator(prng_seed);
	uint minimum_set_cover = 0;

	if (non_zero_unqiue_kmer_counts.sum() > 0) {

		minimum_set_cover = sparsity_estimator.estimateMinimumSetCover(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities, &non_zero_unqiue_kmer_counts);

	} else {

		minimum_set_cover = variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols();
	}

	assert(minimum_set_cover > 0);
	assert(minimum_set_cover <= variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols());
	assert(non_zero_unqiue_kmer_counts.sum() == 0);

	if (minimum_set_cover == variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols()) {

		haplotype_frequency_distribution = new HaplotypeFrequencyDistribution(new FrequencyDistribution(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols(), prng_seed));

	} else {

		double sparsity = minimum_set_cover / static_cast<double>(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols());
		haplotype_frequency_distribution = new HaplotypeFrequencyDistribution(new SparseFrequencyDistribution(variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols(), prng_seed, sparsity));
	}
}


void VariantClusterGenotyper::restart(const uint max_haplotype_variant_kmers) {

	use_multicluster_kmers = false;

	assert(unique_diplotype_log_probabilities.size() == samples.size());
	assert(multicluster_diplotype_log_probabilities.size() == samples.size());

	variant_cluster_haplotypes.sampleKmerSubset(&prng, max_haplotype_variant_kmers, samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		unique_diplotype_log_probabilities.at(sample_idx).clear();
		multicluster_diplotype_log_probabilities.at(sample_idx).clear();
	}

	haplotype_frequency_distribution->reset();
}


bool VariantClusterGenotyper::hasMulticlusterKmer() {

	return (variant_cluster_haplotypes.multicluster_kmers.size() > 0);
}


bool VariantClusterGenotyper::hasExcludedKmer() {

	return has_excluded_kmer;
}


void VariantClusterGenotyper::reduceSamplePloidy(Utils::Ploidy * sample_ploidy) {

	assert(*sample_ploidy != Utils::Ploidy::Null);

	if (*sample_ploidy == Utils::Ploidy::Diploid) {

		*sample_ploidy = Utils::Ploidy::Haploid;
	
	} else {

		*sample_ploidy = Utils::Ploidy::Null;
	}
}


vector<Utils::Ploidy> VariantClusterGenotyper::getVariantClusterPloidy(const uint variant_cluster_idx) {

	vector<Utils::Ploidy> variant_cluster_ploidy;
	variant_cluster_ploidy.reserve(samples.size());

	assert(diplotype_sample.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		Utils::Ploidy sample_ploidy = Utils::Ploidy::Diploid;

		if (diplotype_sample.at(sample_idx).first != Utils::ushort_overflow) {

			if (!(binary_search(variant_cluster_haplotypes.haplotypes.at(diplotype_sample.at(sample_idx).first).nested_variant_cluster_indices.begin(), variant_cluster_haplotypes.haplotypes.at(diplotype_sample.at(sample_idx).first).nested_variant_cluster_indices.end(), variant_cluster_idx))) {

				reduceSamplePloidy(&sample_ploidy);				
			}
		
		} else {

			reduceSamplePloidy(&sample_ploidy);			
		} 

		if (diplotype_sample.at(sample_idx).second != Utils::ushort_overflow) {

			if (!(binary_search(variant_cluster_haplotypes.haplotypes.at(diplotype_sample.at(sample_idx).second).nested_variant_cluster_indices.begin(), variant_cluster_haplotypes.haplotypes.at(diplotype_sample.at(sample_idx).second).nested_variant_cluster_indices.end(), variant_cluster_idx))) {

				reduceSamplePloidy(&sample_ploidy);			
			}
		
		} else {

			reduceSamplePloidy(&sample_ploidy);
		} 

		variant_cluster_ploidy.push_back(sample_ploidy);
	}

	return variant_cluster_ploidy;
}


ushort VariantClusterGenotyper::haplotypeToAlleleIndex(const ushort haplotype_idx, const ushort variant_idx) {

	if (haplotype_idx < Utils::ushort_overflow) {

		assert(variant_cluster_haplotypes.haplotypes.at(haplotype_idx).variant_allele_indices.at(variant_idx) < variant_cluster_info.at(variant_idx).num_alleles);
		return variant_cluster_haplotypes.haplotypes.at(haplotype_idx).variant_allele_indices.at(variant_idx);

	} else {

		return variant_cluster_info.at(variant_idx).num_alleles - 1;
	}
}


vector<ushort> VariantClusterGenotyper::getNonCoveredAlleles(const ushort variant_idx) {

	vector<bool> is_allele_covered(variant_cluster_info.at(variant_idx).num_alleles, false);

	for (auto & haplotype_info: variant_cluster_haplotypes.haplotypes) {

		assert(haplotype_info.variant_allele_indices.at(variant_idx) < variant_cluster_info.at(variant_idx).num_alleles);
		is_allele_covered.at(haplotype_info.variant_allele_indices.at(variant_idx)) = true;
	}

	if (variant_cluster_info.at(variant_idx).has_dependency) {

		is_allele_covered.back() = true;
	}

	vector<ushort> non_covered_alleles;

	for (ushort allele_idx = 0; allele_idx < is_allele_covered.size(); allele_idx++) {

		if (!(is_allele_covered.at(allele_idx))) {

			non_covered_alleles.push_back(allele_idx);
		}
	}

	return non_covered_alleles;
}


vector<Genotypes::SampleStats> VariantClusterGenotyper::getGenotypeSampleStats(const uint variant_idx, const vector<Utils::Ploidy> & chromosome_ploidy) {

	vector<Genotypes::SampleStats> sample_stats;
	sample_stats.reserve(samples.size());

	assert(allele_kmer_stats.at(variant_idx).size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

            sample_stats.emplace_back(variant_cluster_info.at(variant_idx).num_alleles, (variant_cluster_info.at(variant_idx).num_alleles * (variant_cluster_info.at(variant_idx).num_alleles - 1)) / 2 + variant_cluster_info.at(variant_idx).num_alleles);

        } else if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

            sample_stats.emplace_back(variant_cluster_info.at(variant_idx).num_alleles, variant_cluster_info.at(variant_idx).num_alleles);
        
        } else {

        	assert(chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Null);
            sample_stats.emplace_back(0, 0);
        }

		uint num_iterations = 0;

	    pair<vector<pair<ushort, ushort> >, float> max_posterior_genotypes;
	    max_posterior_genotypes.second = 0;

		for (auto & diplotype_sampling_frequency: diplotype_sampling_frequencies) {

			assert(diplotype_sampling_frequency.first.first <= diplotype_sampling_frequency.first.second);
			assert(diplotype_sampling_frequency.second.size() == samples.size());

			if (diplotype_sampling_frequency.second.at(sample_idx) > 0) {

				auto allele_idx_1 = haplotypeToAlleleIndex(diplotype_sampling_frequency.first.first, variant_idx); 
				auto allele_idx_2 = haplotypeToAlleleIndex(diplotype_sampling_frequency.first.second, variant_idx);

				pair<ushort, ushort> genotype_estimate(allele_idx_1, allele_idx_2);

				if (allele_idx_1 > allele_idx_2) {

					genotype_estimate.first = allele_idx_2;
					genotype_estimate.second = allele_idx_1;
				}

				auto genotype_idx = genotype_estimate.first;

	            if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

	            	assert((diplotype_sampling_frequency.first.first != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);
	            	assert((diplotype_sampling_frequency.first.second != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);

	                genotype_idx = (genotype_estimate.second * (genotype_estimate.second + 1)) / 2 + genotype_estimate.first;
	                
	                sample_stats.back().genotype_posteriors.at(genotype_idx) += diplotype_sampling_frequency.second.at(sample_idx);
	                sample_stats.back().allele_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);

	                if (genotype_estimate.first != genotype_estimate.second) {
	                    
	                    sample_stats.back().allele_posteriors.at(genotype_estimate.second) += diplotype_sampling_frequency.second.at(sample_idx);
	                }

	            } else if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

	            	assert((diplotype_sampling_frequency.first.first != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);
	            	assert(diplotype_sampling_frequency.first.second == Utils::ushort_overflow);

	                sample_stats.back().genotype_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);
	                sample_stats.back().allele_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);
	            
	            } else {

	            	assert(chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Null);
	           		assert(diplotype_sampling_frequency.first.first == Utils::ushort_overflow);
	           		assert(diplotype_sampling_frequency.first.second == Utils::ushort_overflow);
	            }

	            num_iterations += diplotype_sampling_frequency.second.at(sample_idx);

	            if (chromosome_ploidy.at(sample_idx) != Utils::Ploidy::Null) {

					if (Utils::floatCompare(max_posterior_genotypes.second, sample_stats.back().genotype_posteriors.at(genotype_idx))) {
			 		
			            max_posterior_genotypes.first.push_back(genotype_estimate);

			 		} else if (max_posterior_genotypes.second < sample_stats.back().genotype_posteriors.at(genotype_idx)) {

			            max_posterior_genotypes.first.clear();
			            max_posterior_genotypes.first.push_back(genotype_estimate);

			            max_posterior_genotypes.second = sample_stats.back().genotype_posteriors.at(genotype_idx);
			        } 
			    }
			}
        }

        max_posterior_genotypes.second /= num_iterations;

       	for (auto &posterior: sample_stats.back().genotype_posteriors) {

            posterior /= num_iterations;
        }

        for (auto &posterior: sample_stats.back().allele_posteriors) {

            posterior /= num_iterations;
        }

        if (is_parameter_estimation_cluster and (max_posterior_genotypes.second >= Utils::min_posterior_allele_kmer_estimate)) {

        	assert(chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Diploid);
        	assert(max_posterior_genotypes.first.size() == 1);

	        sample_stats.back().is_allele_kmer_estimate_variant = true;
	    }

        assert(sample_stats.back().genotype_estimate.empty());

        if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

	        assert(!(max_posterior_genotypes.first.empty()));

            if (max_posterior_genotypes.first.size() == 1) {

			    assert(max_posterior_genotypes.first.front().first <= max_posterior_genotypes.first.front().second);
			    assert(max_posterior_genotypes.first.front().second < variant_cluster_info.at(variant_idx).num_alleles);

                sample_stats.back().genotype_estimate.push_back(max_posterior_genotypes.first.front().first);
                sample_stats.back().genotype_estimate.push_back(max_posterior_genotypes.first.front().second);
 
            } else {

                sample_stats.back().genotype_estimate.push_back(Utils::ushort_overflow);
                sample_stats.back().genotype_estimate.push_back(Utils::ushort_overflow);            
            }
        
        } else if (chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

	        assert(!(max_posterior_genotypes.first.empty()));

            if (max_posterior_genotypes.first.size() == 1) {

			    assert(max_posterior_genotypes.first.front().first < variant_cluster_info.at(variant_idx).num_alleles);
                sample_stats.back().genotype_estimate.push_back(max_posterior_genotypes.first.front().first);
                
            } else {

                sample_stats.back().genotype_estimate.push_back(Utils::ushort_overflow);
            }    
        
        } else {
			
			assert(chromosome_ploidy.at(sample_idx) == Utils::Ploidy::Null);
	        assert(max_posterior_genotypes.first.empty());
        } 
		
		sample_stats.back().allele_kmer_stats = allele_kmer_stats.at(variant_idx).at(sample_idx);
	}

    return sample_stats;
}


Genotypes::VariantStats VariantClusterGenotyper::getGenotypeVariantStats(const uint variant_idx, const vector<Genotypes::SampleStats> & sample_stats) {

	assert(sample_stats.size() == samples.size());

	Genotypes::VariantStats variant_stat(variant_cluster_info.at(variant_idx).num_alleles);

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		for (auto & allele_idx: sample_stats.at(sample_idx).genotype_estimate) {

			if (allele_idx != Utils::ushort_overflow) {

				assert(allele_idx < variant_cluster_info.at(variant_idx).num_alleles);

				variant_stat.total_count++;

				if (allele_idx > 0) {

					variant_stat.alt_allele_counts.at(allele_idx - 1)++;
				}
			}
		}

		assert((variant_stat.allele_call_probabilities.size() == sample_stats.at(sample_idx).allele_posteriors.size()) or sample_stats.at(sample_idx).allele_posteriors.empty());

		for (ushort allele_idx = 0; allele_idx < sample_stats.at(sample_idx).allele_posteriors.size(); allele_idx++) {

			variant_stat.allele_call_probabilities.at(allele_idx) = max(variant_stat.allele_call_probabilities.at(allele_idx), sample_stats.at(sample_idx).allele_posteriors.at(allele_idx));
		}
	}

	assert(variant_stat.alt_allele_counts.size() == variant_stat.alt_allele_frequency.size());

	if (variant_stat.total_count > 0) {

		for (ushort alt_allele_idx = 0; alt_allele_idx < variant_stat.alt_allele_counts.size(); alt_allele_idx++) {

			variant_stat.alt_allele_frequency.at(alt_allele_idx) = variant_stat.alt_allele_counts.at(alt_allele_idx)/static_cast<float>(variant_stat.total_count);
		}
	}

    return variant_stat;
}


vector<Genotypes*> VariantClusterGenotyper::getGenotypes(const vector<Utils::Ploidy> & chromosome_ploidy) {

	assert(samples.size() == chromosome_ploidy.size());
	assert(variant_cluster_info.size() == allele_kmer_stats.size());

	vector<Genotypes*> variant_genotypes;
	variant_genotypes.reserve(variant_cluster_info.size());

	for (ushort variant_idx = 0; variant_idx < variant_cluster_info.size(); variant_idx++) {

		auto genotypes = new Genotypes();

		genotypes->num_candidates = variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols();
		genotypes->variant_cluster_size = variant_cluster_info.size();

		genotypes->variant_info = variant_cluster_info.at(variant_idx);

		genotypes->has_redundant_sequence = is_parameter_estimation_cluster;
		genotypes->has_redundant_sequence = has_redundant_sequence;

	    genotypes->non_covered_alleles = getNonCoveredAlleles(variant_idx);
	    genotypes->sample_stats = getGenotypeSampleStats(variant_idx, chromosome_ploidy);
	    genotypes->variant_stats = getGenotypeVariantStats(variant_idx, genotypes->sample_stats);

		variant_genotypes.push_back(genotypes);
	}

	return variant_genotypes;
}


void VariantClusterGenotyper::updateMulticlusterDiplotypeLogProb(const CountDistribution & count_distribution, const ushort sample_idx) {

	for (uint kmer_subset_idx = 0; kmer_subset_idx < variant_cluster_haplotypes.multicluster_kmer_subset_indices.size(); kmer_subset_idx++) {

		if (variant_cluster_haplotypes.isMulticlusterKmerUpdated(kmer_subset_idx, sample_idx)) {

			for (auto & multicluster_diplotype_log_probability: multicluster_diplotype_log_probabilities.at(sample_idx)) {

				const uint kmer_idx = variant_cluster_haplotypes.multicluster_kmer_subset_indices.at(kmer_subset_idx);

				auto prev_kmer_multiplicity = variant_cluster_haplotypes.getPreviousMulticlusterKmerMultiplicity(kmer_subset_idx, multicluster_diplotype_log_probability.first, diplotype_sample.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

				multicluster_diplotype_log_probability.second -= count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).bias_idx, prev_kmer_multiplicity, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx));

				auto kmer_multiplicity = variant_cluster_haplotypes.getMulticlusterKmerMultiplicity(kmer_idx, multicluster_diplotype_log_probability.first, diplotype_sample.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

				multicluster_diplotype_log_probability.second += count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).bias_idx, kmer_multiplicity, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx));
			}
		}
	}
}


double VariantClusterGenotyper::calcDiplotypeLogProb(const CountDistribution & count_distribution, const ushort sample_idx, const pair<ushort, ushort> & diplotype) {

	double diplotype_log_prob = 0;

	assert(diplotype.first != Utils::ushort_overflow);

	if (diplotype.second == Utils::ushort_overflow) {

		diplotype_log_prob += log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second);

	} else {

		if (diplotype.first == diplotype.second) {

			diplotype_log_prob += 2 * log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second);

		} else {

			diplotype_log_prob += log(2) + log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second) + log(haplotype_frequency_distribution->getElementFrequency(diplotype.second).second);
		}
	}

	auto unique_diplotype_log_probabilities_emplace = unique_diplotype_log_probabilities.at(sample_idx).emplace(diplotype, 0);

	if (unique_diplotype_log_probabilities_emplace.second) {

		for (auto & kmer_idx: variant_cluster_haplotypes.unique_kmer_subset_indices) {

			auto bias_idx = variant_cluster_haplotypes.unique_kmers.at(kmer_idx).bias_idx;
			auto kmer_multiplicity = variant_cluster_haplotypes.getUniqueKmerMultiplicity(kmer_idx, diplotype, samples.at(sample_idx).gender);

			if (variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts) {

				unique_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCountLogProb(sample_idx, bias_idx, kmer_multiplicity, variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts->getSampleCount(sample_idx));
			
			} else {

				unique_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCountLogProb(sample_idx, bias_idx, kmer_multiplicity, 0);				
			}
		}
	}

	diplotype_log_prob += unique_diplotype_log_probabilities_emplace.first->second;

	if (use_multicluster_kmers) {

		auto multicluster_diplotype_log_probabilities_emplace = multicluster_diplotype_log_probabilities.at(sample_idx).emplace(diplotype, 0);

		if (multicluster_diplotype_log_probabilities_emplace.second) {

			for (auto & kmer_idx: variant_cluster_haplotypes.multicluster_kmer_subset_indices) {

				auto kmer_multiplicity = variant_cluster_haplotypes.getMulticlusterKmerMultiplicity(kmer_idx, diplotype, diplotype_sample.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

				multicluster_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).bias_idx, kmer_multiplicity, variant_cluster_haplotypes.multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx));
			}
		}

		diplotype_log_prob += multicluster_diplotype_log_probabilities_emplace.first->second;
	}

	assert(isfinite(diplotype_log_prob));		
	
	if (!((diplotype_log_prob < 0) or Utils::doubleCompare(diplotype_log_prob, 0))) {

		cout << "\n### DEBUG: " << diplotype_log_prob << " diplotype log probability (" << diplotype.first << ", " << diplotype.second << ")\n" << endl;
	}

	return diplotype_log_prob;
}


void VariantClusterGenotyper::sampleGenotypes(const CountDistribution & count_distribution, const vector<Utils::Ploidy> & variant_cluster_ploidy, const bool collect_samples) {

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		auto prev_diplotype = diplotype_sample.at(sample_idx);

		updateMulticlusterDiplotypeLogProb(count_distribution, sample_idx);		
		sampleGenotype(count_distribution, sample_idx, variant_cluster_ploidy.at(sample_idx));

		variant_cluster_haplotypes.updateMulticlusterKmerMultiplicities(diplotype_sample.at(sample_idx), prev_diplotype, sample_idx);

		if (collect_samples) {

			auto diplotype_sampling_frequencies_emplace = diplotype_sampling_frequencies.emplace(diplotype_sample.at(sample_idx), vector<uint>(samples.size(), 0));
			diplotype_sampling_frequencies_emplace.first->second.at(sample_idx)++;

			variant_cluster_haplotypes.updateAlleleKmerStats(&allele_kmer_stats, diplotype_sample.at(sample_idx), sample_idx, samples.at(sample_idx).gender);
		}
	}

	use_multicluster_kmers = !(variant_cluster_haplotypes.multicluster_kmer_subset_indices.empty());
}


void VariantClusterGenotyper::sampleGenotype(const CountDistribution & count_distribution, const ushort sample_idx, const Utils::Ploidy sample_variant_cluster_ploidy) {

	const ushort num_haplotypes = variant_cluster_haplotypes.haplotype_unique_kmer_multiplicities.cols();

	uint num_expected_diplotypes = (num_haplotypes * num_haplotypes)/2;
	LogDiscreteSampler diplotype_sampler(num_expected_diplotypes);
	
	vector<pair<ushort, ushort> > diplotypes;
	diplotypes.reserve(num_expected_diplotypes);

	if (sample_variant_cluster_ploidy == Utils::Ploidy::Diploid) {

		for (ushort i = 0; i < num_haplotypes; i++) {

			for (ushort j = i; j < num_haplotypes; j++) {

				if (haplotype_frequency_distribution->getElementFrequency(i).first) {

					if (haplotype_frequency_distribution->getElementFrequency(j).first) {

						diplotype_sampler.addOutcome(calcDiplotypeLogProb(count_distribution, sample_idx, make_pair(i, j)));
						diplotypes.emplace_back(i, j);
					}

				} else {

					break;
				}
			}
		}

	} else if (sample_variant_cluster_ploidy == Utils::Ploidy::Haploid) {

		for (ushort i = 0; i < num_haplotypes; i++) {

			if (haplotype_frequency_distribution->getElementFrequency(i).first) {

				diplotype_sampler.addOutcome(calcDiplotypeLogProb(count_distribution, sample_idx, make_pair(i, Utils::ushort_overflow)));
				diplotypes.emplace_back(i, Utils::ushort_overflow);
			}
		}

	} else {

		assert(sample_variant_cluster_ploidy == Utils::Ploidy::Null);

		diplotype_sampler.addOutcome(0);
		diplotypes.emplace_back(Utils::ushort_overflow, Utils::ushort_overflow);
	}

	assert(diplotypes.size() > 0);
	diplotype_sample.at(sample_idx) = diplotypes.at(diplotype_sampler.sample(&prng));

	haplotype_frequency_distribution->incrementObservationCount(diplotype_sample.at(sample_idx).first);
	haplotype_frequency_distribution->incrementObservationCount(diplotype_sample.at(sample_idx).second);
}


void VariantClusterGenotyper::getCountAllocations(CountAllocation * sample_noise_counts, const CountDistribution & count_distribution) {

	is_parameter_estimation_cluster = true;

	assert(!has_excluded_kmer);
	assert(!use_multicluster_kmers);

	assert(variant_cluster_info.size() == 1);

	assert(variant_cluster_haplotypes.multicluster_kmers.empty());
	assert(variant_cluster_haplotypes.multicluster_kmer_subset_indices.empty());
	assert(variant_cluster_haplotypes.haplotype_multicluster_kmer_multiplicities.rows() == 0);

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		assert(diplotype_sample.at(sample_idx).first != Utils::ushort_overflow);
		assert(diplotype_sample.at(sample_idx).second != Utils::ushort_overflow);

		unique_diplotype_log_probabilities.at(sample_idx).clear();
		assert(multicluster_diplotype_log_probabilities.at(sample_idx).empty());

		for (auto & kmer_idx: variant_cluster_haplotypes.unique_kmer_subset_indices) {

			auto kmer_multiplicity = variant_cluster_haplotypes.getUniqueKmerMultiplicity(kmer_idx, diplotype_sample.at(sample_idx), samples.at(sample_idx).gender);

			if (kmer_multiplicity == 0) {

				if (variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts) {

					sample_noise_counts->addCount(sample_idx, variant_cluster_haplotypes.unique_kmers.at(kmer_idx).counts->getSampleCount(sample_idx));
				
				} else {

					sample_noise_counts->addCount(sample_idx, 0);
				}
			}
		}
	}
}


void VariantClusterGenotyper::sampleHaplotypeFrequencies() {

  	assert((haplotype_frequency_distribution->sumHaplotypeCount() + haplotype_frequency_distribution->sumMissingCount()) == (2 * samples.size()));
  	haplotype_frequency_distribution->sampleFrequencies();
}
