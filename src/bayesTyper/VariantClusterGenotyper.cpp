
/*
VariantClusterGenotyper.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


VariantClusterGenotyper::VariantClusterGenotyper(VariantClusterGraph * variant_cluster_graph, KmerCountsHash * kmer_hash, const vector<Sample> & samples_in, const uint prng_seed, const uchar num_genomic_rate_gc_bias_bins) : samples(samples_in) {

	prng = mt19937(prng_seed);
	
	use_multicluster_kmers = false;

	unique_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());
	multicluster_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());

	diplotypes = vector<pair<ushort,ushort> >(samples.size(), make_pair(Utils::ushort_overflow, Utils::ushort_overflow));
	
	variant_cluster_info = variant_cluster_graph->getInfo();
	allele_kmer_stats.reserve(variant_cluster_info.size());

	for (auto & variant_info: variant_cluster_info) {

		allele_kmer_stats.emplace_back(vector<AlleleKmerStats>(samples.size(), AlleleKmerStats(variant_info.numberOfAlleles())));
	}

	variant_cluster_haplotypes = variant_cluster_graph->getHaplotypeCandidates(kmer_hash, num_genomic_rate_gc_bias_bins);

	assert(variant_cluster_haplotypes.haplotype_kmer_multiplicities.cols() > 0);
	assert(variant_cluster_haplotypes.haplotype_kmer_multiplicities.cols() < Utils::ushort_overflow);

	assert(variant_cluster_haplotypes.haplotypes.size() == static_cast<uint>(variant_cluster_haplotypes.haplotype_kmer_multiplicities.cols()));

	assert(variant_cluster_haplotypes.unique_kmer_subset_indices.empty());
	assert(variant_cluster_haplotypes.multicluster_kmer_subset_indices.empty());

    variant_cluster_haplotypes.kmer_stats_cache = vector<VariantClusterHaplotypes::KmerStatsCache>(samples.size(), variant_cluster_info.size());

	Utils::RowVectorXbool non_zero_kmer_counts = Utils::RowVectorXbool::Zero(variant_cluster_haplotypes.kmers.size());

	for (uint kmer_idx = 0; kmer_idx < variant_cluster_haplotypes.kmers.size(); kmer_idx++) {

		if (variant_cluster_haplotypes.kmers.at(kmer_idx).counts) {

			non_zero_kmer_counts(kmer_idx) = true;	
		} 
	}

	SparsityEstimator sparsity_estimator = SparsityEstimator(prng_seed);

	auto minimum_haplotype_cover = sparsity_estimator.estimateMinimumColumnCover(variant_cluster_haplotypes.haplotype_kmer_multiplicities, non_zero_kmer_counts, false);
	assert(minimum_haplotype_cover.size() <= variant_cluster_haplotypes.haplotypes.size());

	haplotype_frequency_distribution = new SparseHaplotypeFrequencyDistribution(minimum_haplotype_cover, variant_cluster_haplotypes.haplotypes.size(), prng_seed);
}

VariantClusterGenotyper::~VariantClusterGenotyper() {

	delete haplotype_frequency_distribution;
}

void VariantClusterGenotyper::reset(const float kmer_subsampling_rate, const uint max_haplotype_variant_kmers) {

	use_multicluster_kmers = false;

	assert(unique_diplotype_log_probabilities.size() == samples.size());
	assert(multicluster_diplotype_log_probabilities.size() == samples.size());

	variant_cluster_haplotypes.sampleKmerSubset(&prng, kmer_subsampling_rate, max_haplotype_variant_kmers, samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		unique_diplotype_log_probabilities.at(sample_idx).clear();
		multicluster_diplotype_log_probabilities.at(sample_idx).clear();
	}

	haplotype_frequency_distribution->reset();
}

void VariantClusterGenotyper::clearCache() {

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		unique_diplotype_log_probabilities.at(sample_idx).clear();
		multicluster_diplotype_log_probabilities.at(sample_idx).clear();
	}
}

void VariantClusterGenotyper::updateNestedPloidy(Utils::Ploidy * nested_ploidy) {

	assert(*nested_ploidy != Utils::Ploidy::Null);

	if (*nested_ploidy == Utils::Ploidy::Diploid) {

		*nested_ploidy = Utils::Ploidy::Haploid;
	
	} else {

		*nested_ploidy = Utils::Ploidy::Null;
	}
}

void VariantClusterGenotyper::addNestedKmerStats(vector<KmerStats> * nested_kmer_stats, const uint variant_cluster_idx, const ushort haplotype_idx, const vector<KmerStats> & haplotype_kmer_stats_cache) {

	assert(nested_kmer_stats->size() < 2);
	assert(haplotype_idx != Utils::ushort_overflow);

	auto nested_variant_cluster_dependency_it = variant_cluster_haplotypes.nested_variant_cluster_dependency.find(variant_cluster_idx);
	assert(nested_variant_cluster_dependency_it != variant_cluster_haplotypes.nested_variant_cluster_dependency.end());

	ushort variant_idx = Utils::ushort_overflow; 

	for (auto & nested_variant_idx: nested_variant_cluster_dependency_it->second) {

		auto allele_idx = variant_cluster_haplotypes.haplotypes.at(haplotype_idx).variant_allele_indices.at(nested_variant_idx);
		assert(allele_idx < variant_cluster_info.at(nested_variant_idx).numberOfAlleles());

		if (!variant_cluster_info.at(nested_variant_idx).isMissing(allele_idx)) {

			assert(allele_idx > 0);
			variant_idx = nested_variant_idx;

			break;
		}
	}

	assert(variant_idx != Utils::ushort_overflow);
	nested_kmer_stats->push_back(haplotype_kmer_stats_cache.at(variant_idx));
}

void VariantClusterGenotyper::updateNestedVariantClusterInfo(vector<VariantClusterHaplotypes::NestedVariantClusterInfo> * nested_variant_cluster_info, const uint variant_cluster_idx) {

	assert(nested_variant_cluster_info->size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		if (diplotypes.at(sample_idx).first != Utils::ushort_overflow) {

			if (!binary_search(variant_cluster_haplotypes.haplotypes.at(diplotypes.at(sample_idx).first).nested_variant_cluster_indices.begin(), variant_cluster_haplotypes.haplotypes.at(diplotypes.at(sample_idx).first).nested_variant_cluster_indices.end(), variant_cluster_idx)) {

				updateNestedPloidy(&(nested_variant_cluster_info->at(sample_idx).nested_ploidy));				
				addNestedKmerStats(&(nested_variant_cluster_info->at(sample_idx).nested_kmer_stats), variant_cluster_idx, diplotypes.at(sample_idx).first, variant_cluster_haplotypes.kmer_stats_cache.at(sample_idx).haplotype_1);
			}
		} 

		if (diplotypes.at(sample_idx).second != Utils::ushort_overflow) {

			if (!binary_search(variant_cluster_haplotypes.haplotypes.at(diplotypes.at(sample_idx).second).nested_variant_cluster_indices.begin(), variant_cluster_haplotypes.haplotypes.at(diplotypes.at(sample_idx).second).nested_variant_cluster_indices.end(), variant_cluster_idx)) {

				updateNestedPloidy(&(nested_variant_cluster_info->at(sample_idx).nested_ploidy));			
				addNestedKmerStats(&(nested_variant_cluster_info->at(sample_idx).nested_kmer_stats), variant_cluster_idx, diplotypes.at(sample_idx).second, variant_cluster_haplotypes.kmer_stats_cache.at(sample_idx).haplotype_2);
			}	
		}
	}
}

ushort VariantClusterGenotyper::haplotypeToAlleleIndex(const ushort haplotype_idx, const ushort variant_idx) {

	if (haplotype_idx != Utils::ushort_overflow) {

		assert(variant_cluster_haplotypes.haplotypes.at(haplotype_idx).variant_allele_indices.at(variant_idx) < variant_cluster_info.at(variant_idx).numberOfAlleles());
		return variant_cluster_haplotypes.haplotypes.at(haplotype_idx).variant_allele_indices.at(variant_idx);

	} else {

		return variant_cluster_info.at(variant_idx).numberOfAlleles() - 1;
	}
}

vector<ushort> VariantClusterGenotyper::getNonCoveredAlleles(const ushort variant_idx) {

	vector<bool> is_allele_covered(variant_cluster_info.at(variant_idx).numberOfAlleles(), false);

	for (auto & haplotype_info: variant_cluster_haplotypes.haplotypes) {

		assert(haplotype_info.variant_allele_indices.at(variant_idx) < variant_cluster_info.at(variant_idx).numberOfAlleles());
		is_allele_covered.at(haplotype_info.variant_allele_indices.at(variant_idx)) = true;
	}

	if (variant_cluster_info.at(variant_idx).has_dependency) {

		is_allele_covered.back() = true;
	}

	vector<ushort> non_covered_alleles;

	for (ushort allele_idx = 0; allele_idx < is_allele_covered.size(); allele_idx++) {

		if (!is_allele_covered.at(allele_idx)) {

			non_covered_alleles.push_back(allele_idx);
		}
	}

	return non_covered_alleles;
}

vector<Genotypes::SampleStats> VariantClusterGenotyper::getGenotypeSampleStats(const uint variant_idx, const vector<Utils::Ploidy> & chrom_ploidy, const Filters & filters) {

	vector<Genotypes::SampleStats> sample_stats;
	sample_stats.reserve(samples.size());

	assert(allele_kmer_stats.at(variant_idx).size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

        if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

            sample_stats.emplace_back(variant_cluster_info.at(variant_idx).numberOfAlleles(), (variant_cluster_info.at(variant_idx).numberOfAlleles() * (variant_cluster_info.at(variant_idx).numberOfAlleles() - 1)) / 2 + variant_cluster_info.at(variant_idx).numberOfAlleles());

        } else if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

            sample_stats.emplace_back(variant_cluster_info.at(variant_idx).numberOfAlleles(), variant_cluster_info.at(variant_idx).numberOfAlleles());
        
        } else {

        	assert(chrom_ploidy.at(sample_idx) == Utils::Ploidy::Null);
            sample_stats.emplace_back(0, 0);
        }

		uint num_iterations = 0;

	    pair<vector<pair<ushort, ushort> >, float> max_posterior_genotypes;
	    max_posterior_genotypes.second = 0;

		for (auto & diplotype_sampling_frequency: diplotype_sampling_frequencies) {

			assert(diplotype_sampling_frequency.first.first <= diplotype_sampling_frequency.first.second);
			assert(diplotype_sampling_frequency.second.size() == samples.size());

			if (diplotype_sampling_frequency.second.at(sample_idx) > 0) {

				pair<ushort, ushort> genotype_estimate(Utils::ushort_overflow, Utils::ushort_overflow);
				auto genotype_idx = genotype_estimate.first;

	            if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

	            	assert((diplotype_sampling_frequency.first.first != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);
	            	assert((diplotype_sampling_frequency.first.second != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);

					genotype_estimate.first = haplotypeToAlleleIndex(diplotype_sampling_frequency.first.first, variant_idx); 
					genotype_estimate.second = haplotypeToAlleleIndex(diplotype_sampling_frequency.first.second, variant_idx);

					if (genotype_estimate.first > genotype_estimate.second) {

						auto allele_idx_tmp = genotype_estimate.first;

						genotype_estimate.first = genotype_estimate.second;
						genotype_estimate.second = allele_idx_tmp;
					}

	                genotype_idx = (genotype_estimate.second * (genotype_estimate.second + 1)) / 2 + genotype_estimate.first;
	                
	                sample_stats.back().genotype_posteriors.at(genotype_idx) += diplotype_sampling_frequency.second.at(sample_idx);
	                sample_stats.back().allele_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);

	                if (genotype_estimate.first != genotype_estimate.second) {
	                    
	                    sample_stats.back().allele_posteriors.at(genotype_estimate.second) += diplotype_sampling_frequency.second.at(sample_idx);
	                }

	            } else if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

	            	assert((diplotype_sampling_frequency.first.first != Utils::ushort_overflow) or variant_cluster_info.at(variant_idx).has_dependency);
	            	assert(diplotype_sampling_frequency.first.second == Utils::ushort_overflow);

					genotype_estimate.first = haplotypeToAlleleIndex(diplotype_sampling_frequency.first.first, variant_idx); 
				 	genotype_idx = genotype_estimate.first;

	                sample_stats.back().genotype_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);
	                sample_stats.back().allele_posteriors.at(genotype_estimate.first) += diplotype_sampling_frequency.second.at(sample_idx);
	            
	            } else {

	            	assert(chrom_ploidy.at(sample_idx) == Utils::Ploidy::Null);
	           		assert(diplotype_sampling_frequency.first.first == Utils::ushort_overflow);
	           		assert(diplotype_sampling_frequency.first.second == Utils::ushort_overflow);
	            }

	            num_iterations += diplotype_sampling_frequency.second.at(sample_idx);

	            if (chrom_ploidy.at(sample_idx) != Utils::Ploidy::Null) {

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

       	for (auto & posterior: sample_stats.back().genotype_posteriors) {

            posterior /= num_iterations;
        }

        for (auto & posterior: sample_stats.back().allele_posteriors) {

            posterior /= num_iterations;
        }

		sample_stats.back().allele_kmer_stats = allele_kmer_stats.at(variant_idx).at(sample_idx);

    	assert(sample_stats.back().allele_kmer_stats.count_stats.size() == variant_cluster_info.at(variant_idx).numberOfAlleles());
    	assert(sample_stats.back().allele_kmer_stats.fraction_stats.size() == variant_cluster_info.at(variant_idx).numberOfAlleles());
		assert(sample_stats.back().allele_kmer_stats.mean_stats.size() == variant_cluster_info.at(variant_idx).numberOfAlleles());

    	assert((sample_stats.back().allele_posteriors.size() == sample_stats.back().allele_kmer_stats.count_stats.size()) or sample_stats.back().allele_posteriors.empty());
    	assert(sample_stats.back().allele_posteriors.size() == sample_stats.back().allele_filters.size());

	    for (ushort allele_idx = 0; allele_idx < sample_stats.back().allele_posteriors.size(); allele_idx++) {

	        if (!Utils::floatCompare(sample_stats.back().allele_posteriors.at(allele_idx), 0)) {

            	assert(sample_stats.back().allele_posteriors.at(allele_idx) > 0);

            	auto allele_count_stat = sample_stats.back().allele_kmer_stats.count_stats.at(allele_idx).getMean();

            	assert(allele_count_stat.first >= 0);
            	assert(allele_count_stat.second);

		        if (Utils::floatLess(allele_count_stat.first, filters.minNumberOfKmers())) {

		        	sample_stats.back().allele_filters.at(allele_idx) += 1;      
		        } 

            	auto allele_fraction_stat = sample_stats.back().allele_kmer_stats.fraction_stats.at(allele_idx).getMean();

                assert(Utils::floatCompare(allele_count_stat.first, 0) == !allele_fraction_stat.second);

                if (!Utils::floatCompare(allele_count_stat.first, 0)) {

            		assert(allele_fraction_stat.first >= 0);
            		assert(allele_fraction_stat.second);

			        if (Utils::floatLess(allele_fraction_stat.first, filters.minFractionObservedKmers(sample_idx))) {

			        	sample_stats.back().allele_filters.at(allele_idx) += 2;      
			        }
                }

		        assert(sample_stats.back().allele_filters.at(allele_idx) >= 0);
                assert(sample_stats.back().allele_filters.at(allele_idx) <= 3);
		    }
	    }

        assert(sample_stats.back().genotype_estimate.empty());

	    if (Utils::floatCompare(max_posterior_genotypes.second, 1)) {

	        sample_stats.back().genotype_quality = 99;
	    
	   	} else if (Utils::floatCompare(max_posterior_genotypes.second, 0)) {

	        sample_stats.back().genotype_quality = 0;

	    } else {

	        assert(max_posterior_genotypes.second > 0);
	        assert(max_posterior_genotypes.second < 1);

	        sample_stats.back().genotype_quality = -10 * log10(1 - max_posterior_genotypes.second);
	    }

        if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Diploid) {

			sample_stats.back().genotype_estimate = vector<ushort>(2, Utils::ushort_overflow);
	        assert(!max_posterior_genotypes.first.empty());

            if (max_posterior_genotypes.first.size() == 1) {

            	if (!Utils::floatLess(max_posterior_genotypes.second, filters.minGenotypePosterior())) {

				    assert(max_posterior_genotypes.first.front().first <= max_posterior_genotypes.first.front().second);
				    assert(max_posterior_genotypes.first.front().second < variant_cluster_info.at(variant_idx).numberOfAlleles());

				    if ((sample_stats.back().allele_filters.at(max_posterior_genotypes.first.front().first) == 0) and (sample_stats.back().allele_filters.at(max_posterior_genotypes.first.front().second) == 0)) {

				    	sample_stats.back().genotype_estimate.front() = max_posterior_genotypes.first.front().first;
				    	sample_stats.back().genotype_estimate.back() = max_posterior_genotypes.first.front().second;
				    }
				}
            }
        
        } else if (chrom_ploidy.at(sample_idx) == Utils::Ploidy::Haploid) {

			sample_stats.back().genotype_estimate = vector<ushort>(1, Utils::ushort_overflow);
	        assert(!max_posterior_genotypes.first.empty());

            if ((max_posterior_genotypes.first.size() == 1) and !Utils::floatLess(max_posterior_genotypes.second, filters.minGenotypePosterior())) {

			    assert(max_posterior_genotypes.first.front().first < variant_cluster_info.at(variant_idx).numberOfAlleles());

			    if (sample_stats.back().allele_filters.at(max_posterior_genotypes.first.front().first) == 0) {

			    	sample_stats.back().genotype_estimate.front() = max_posterior_genotypes.first.front().first;
                }
            }    
        
        } else {
			
			assert(chrom_ploidy.at(sample_idx) == Utils::Ploidy::Null);
	        assert(max_posterior_genotypes.first.empty());
        } 
	}

    return sample_stats;
}

Genotypes::VariantStats VariantClusterGenotyper::getGenotypeVariantStats(const uint variant_idx, const vector<Genotypes::SampleStats> & sample_stats) {

	assert(sample_stats.size() == samples.size());

	Genotypes::VariantStats variant_stat(variant_cluster_info.at(variant_idx).numberOfAlleles());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		for (auto & allele_idx: sample_stats.at(sample_idx).genotype_estimate) {

			if (allele_idx != Utils::ushort_overflow) {

				assert(allele_idx < variant_cluster_info.at(variant_idx).numberOfAlleles());
				assert(sample_stats.at(sample_idx).allele_filters.at(allele_idx) == 0);

				variant_stat.total_count++;

				if (allele_idx > 0) {

					variant_stat.alt_allele_counts.at(allele_idx - 1)++;
				}
			}
		}

		assert((sample_stats.at(sample_idx).allele_posteriors.size() == variant_stat.allele_call_probabilities.size()) or sample_stats.at(sample_idx).allele_posteriors.empty());

		for (ushort allele_idx = 0; allele_idx < sample_stats.at(sample_idx).allele_posteriors.size(); allele_idx++) {

			assert(sample_stats.at(sample_idx).allele_filters.at(allele_idx) <= 3);

			if (sample_stats.at(sample_idx).allele_filters.at(allele_idx) == 0) {

				variant_stat.allele_call_probabilities.at(allele_idx) = max(variant_stat.allele_call_probabilities.at(allele_idx), sample_stats.at(sample_idx).allele_posteriors.at(allele_idx));
			}
		}
	}

	assert(variant_cluster_info.at(variant_idx).alt_alleles.size() < variant_stat.allele_call_probabilities.size());

    for (ushort alt_allele_idx = 0; alt_allele_idx < variant_cluster_info.at(variant_idx).alt_alleles.size(); alt_allele_idx++) {

    	assert(!variant_cluster_info.at(variant_idx).isMissing(alt_allele_idx + 1));
        variant_stat.max_alt_allele_call_probability = max(variant_stat.max_alt_allele_call_probability, variant_stat.allele_call_probabilities.at(alt_allele_idx + 1));  
    }

	assert(variant_stat.alt_allele_counts.size() == variant_stat.alt_allele_frequency.size());

	if (variant_stat.total_count > 0) {

		for (ushort alt_allele_idx = 0; alt_allele_idx < variant_stat.alt_allele_counts.size(); alt_allele_idx++) {

			variant_stat.alt_allele_frequency.at(alt_allele_idx) = variant_stat.alt_allele_counts.at(alt_allele_idx)/static_cast<float>(variant_stat.total_count);
		}
	}

    return variant_stat;
}

vector<Genotypes*> VariantClusterGenotyper::getGenotypes(const string & chrom_name, const vector<Utils::Ploidy> & chrom_ploidy, const Filters & filters) {

	assert(samples.size() == chrom_ploidy.size());
	assert(variant_cluster_info.size() == allele_kmer_stats.size());

	vector<Genotypes*> variant_genotypes;
	variant_genotypes.reserve(variant_cluster_info.size());

	uint start_position = variant_cluster_info.front().position;
	uint end_position = 0;

	for (auto & variant_info: variant_cluster_info) {

		assert(start_position <= variant_info.position);
        end_position = max(end_position, static_cast<uint>(variant_info.position + variant_info.maxReferenceLength() - 1));
    }

    assert(start_position <= end_position);

	for (ushort variant_idx = 0; variant_idx < variant_cluster_info.size(); variant_idx++) {

		auto genotypes = new Genotypes();

		genotypes->chrom_name = chrom_name;
		genotypes->variant_info = variant_cluster_info.at(variant_idx);

		genotypes->variant_cluster_size = variant_cluster_info.size();
		genotypes->variant_cluster_region = chrom_name + ":" + to_string(start_position) + "-" + to_string(end_position);

		genotypes->num_candidates = variant_cluster_haplotypes.haplotype_kmer_multiplicities.cols();
	    genotypes->non_covered_alleles = getNonCoveredAlleles(variant_idx);

	    genotypes->sample_stats = getGenotypeSampleStats(variant_idx, chrom_ploidy, filters);
	    genotypes->variant_stats = getGenotypeVariantStats(variant_idx, genotypes->sample_stats);

		variant_genotypes.push_back(genotypes);
	}

	return variant_genotypes;
}

void VariantClusterGenotyper::updateMulticlusterDiplotypeLogProb(const CountDistribution & count_distribution, const ushort sample_idx) {

	if (!multicluster_diplotype_log_probabilities.at(sample_idx).empty()) {

		for (uint kmer_subset_idx = 0; kmer_subset_idx < variant_cluster_haplotypes.multicluster_kmer_subset_indices.size(); kmer_subset_idx++) {

			if (variant_cluster_haplotypes.isMulticlusterKmerUpdated(kmer_subset_idx, sample_idx)) {

				for (auto & multicluster_diplotype_log_probability: multicluster_diplotype_log_probabilities.at(sample_idx)) {

					const uint kmer_idx = variant_cluster_haplotypes.multicluster_kmer_subset_indices.at(kmer_subset_idx);

					auto kmer_counts = variant_cluster_haplotypes.kmers.at(kmer_idx).counts;
					assert(kmer_counts);

					auto prev_kmer_multiplicity = variant_cluster_haplotypes.getPreviousMulticlusterKmerMultiplicity(kmer_subset_idx, multicluster_diplotype_log_probability.first, diplotypes.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

					multicluster_diplotype_log_probability.second -= count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.kmers.at(kmer_idx).bias_idx, prev_kmer_multiplicity, kmer_counts->getSampleCount(sample_idx));

					auto kmer_multiplicity = variant_cluster_haplotypes.getMulticlusterKmerMultiplicity(kmer_idx, multicluster_diplotype_log_probability.first, diplotypes.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

					multicluster_diplotype_log_probability.second += count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.kmers.at(kmer_idx).bias_idx, kmer_multiplicity, kmer_counts->getSampleCount(sample_idx));
				}
			}
		}
	}
}

double VariantClusterGenotyper::calcDiplotypeLogProb(const CountDistribution & count_distribution, const ushort sample_idx, const pair<ushort, ushort> & diplotype) {

	double diplotype_log_prob = 0;

	assert(diplotype.first != Utils::ushort_overflow);

	if (diplotype.second == Utils::ushort_overflow) {

		diplotype_log_prob += log(haplotype_frequency_distribution->getFrequency(diplotype.first).second);

	} else {

		if (diplotype.first == diplotype.second) {

			diplotype_log_prob += 2 * log(haplotype_frequency_distribution->getFrequency(diplotype.first).second);

		} else {

			diplotype_log_prob += log(2) + log(haplotype_frequency_distribution->getFrequency(diplotype.first).second) + log(haplotype_frequency_distribution->getFrequency(diplotype.second).second);
		}
	}

	auto unique_diplotype_log_probabilities_emplace = unique_diplotype_log_probabilities.at(sample_idx).emplace(diplotype, 0);

	if (unique_diplotype_log_probabilities_emplace.second) {

		for (auto & kmer_idx: variant_cluster_haplotypes.unique_kmer_subset_indices) {

			auto bias_idx = variant_cluster_haplotypes.kmers.at(kmer_idx).bias_idx;
			auto kmer_multiplicity = variant_cluster_haplotypes.getUniqueKmerMultiplicity(kmer_idx, diplotype, samples.at(sample_idx).gender);

			auto kmer_counts = variant_cluster_haplotypes.kmers.at(kmer_idx).counts;

			if (kmer_counts) {

				unique_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCountLogProb(sample_idx, bias_idx, kmer_multiplicity, kmer_counts->getSampleCount(sample_idx));

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

				auto kmer_counts = variant_cluster_haplotypes.kmers.at(kmer_idx).counts;
				assert(kmer_counts);

				auto kmer_multiplicity = variant_cluster_haplotypes.getMulticlusterKmerMultiplicity(kmer_idx, diplotype, diplotypes.at(sample_idx), sample_idx, samples.at(sample_idx).gender);

				multicluster_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCountLogProb(sample_idx, variant_cluster_haplotypes.kmers.at(kmer_idx).bias_idx, kmer_multiplicity, kmer_counts->getSampleCount(sample_idx));
			}
		}

		diplotype_log_prob += multicluster_diplotype_log_probabilities_emplace.first->second;
	}

	assert(isfinite(diplotype_log_prob));

	return diplotype_log_prob;
}

void VariantClusterGenotyper::sampleDiplotypes(const CountDistribution & count_distribution, const vector<VariantClusterHaplotypes::NestedVariantClusterInfo> & nested_variant_cluster_info, const bool collect_samples) {

	assert(nested_variant_cluster_info.size() == samples.size());

	vector<ushort> non_zero_haplotypes;
	non_zero_haplotypes.reserve(variant_cluster_haplotypes.haplotypes.size());

	for (ushort haplotype_idx = 0; haplotype_idx < variant_cluster_haplotypes.haplotypes.size(); haplotype_idx++) {

		if (haplotype_frequency_distribution->getFrequency(haplotype_idx).first) {

			non_zero_haplotypes.emplace_back(haplotype_idx);
		}
	}

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		auto prev_diplotype = diplotypes.at(sample_idx);
		
		updateMulticlusterDiplotypeLogProb(count_distribution, sample_idx);		
		sampleDiplotype(non_zero_haplotypes, count_distribution, sample_idx, nested_variant_cluster_info.at(sample_idx).nested_ploidy);

		variant_cluster_haplotypes.updateMulticlusterKmerMultiplicities(diplotypes.at(sample_idx), prev_diplotype, sample_idx);

		if (collect_samples) {

			auto diplotype_sampling_frequencies_emplace = diplotype_sampling_frequencies.emplace(diplotypes.at(sample_idx), vector<uint>(samples.size(), 0));
			diplotype_sampling_frequencies_emplace.first->second.at(sample_idx)++;
		}
	}

	if (collect_samples) {

		variant_cluster_haplotypes.updateAlleleKmerStats(&allele_kmer_stats, samples, variant_cluster_info, nested_variant_cluster_info, diplotypes);
	}

	use_multicluster_kmers = !variant_cluster_haplotypes.multicluster_kmer_subset_indices.empty();
}

void VariantClusterGenotyper::sampleDiplotype(const vector<ushort> & non_zero_haplotypes, const CountDistribution & count_distribution, const ushort sample_idx, const Utils::Ploidy ploidy) {

	const uint num_expected_diplotypes = (non_zero_haplotypes.size() * (non_zero_haplotypes.size() - 1)) / 2 + non_zero_haplotypes.size();
	LogDiscreteSampler diplotype_sampler(num_expected_diplotypes);
	
	vector<pair<ushort, ushort> > diplotype_samples;
	diplotype_samples.reserve(num_expected_diplotypes);

	if (ploidy == Utils::Ploidy::Diploid) {

		auto non_zero_haplotypes_it_1 = non_zero_haplotypes.begin();

		while (non_zero_haplotypes_it_1 != non_zero_haplotypes.end()) {

			auto non_zero_haplotypes_it_2 = non_zero_haplotypes_it_1;

			while (non_zero_haplotypes_it_2 != non_zero_haplotypes.end()) {

				diplotype_sampler.addOutcome(calcDiplotypeLogProb(count_distribution, sample_idx, make_pair(*non_zero_haplotypes_it_1, *non_zero_haplotypes_it_2)));
				diplotype_samples.emplace_back(*non_zero_haplotypes_it_1, *non_zero_haplotypes_it_2);
				non_zero_haplotypes_it_2++;
			}

			non_zero_haplotypes_it_1++;
		}

	} else if (ploidy == Utils::Ploidy::Haploid) {

		for (auto & haplotype_idx: non_zero_haplotypes) {

			diplotype_sampler.addOutcome(calcDiplotypeLogProb(count_distribution, sample_idx, make_pair(haplotype_idx, Utils::ushort_overflow)));
			diplotype_samples.emplace_back(haplotype_idx, Utils::ushort_overflow);
		}

	} else {

		assert(ploidy == Utils::Ploidy::Null);

		diplotype_sampler.addOutcome(0);
		diplotype_samples.emplace_back(Utils::ushort_overflow, Utils::ushort_overflow);
	}

	assert(diplotype_samples.size() > 0);

	diplotypes.at(sample_idx) = diplotype_samples.at(diplotype_sampler.sample(&prng));

	haplotype_frequency_distribution->incrementCount(diplotypes.at(sample_idx).first);
	haplotype_frequency_distribution->incrementCount(diplotypes.at(sample_idx).second);
}

void VariantClusterGenotyper::getNoiseCounts(CountAllocation * noise_counts, const CountDistribution & count_distribution) {

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		for (auto & kmer_idx: variant_cluster_haplotypes.unique_kmer_subset_indices) {

			if (variant_cluster_haplotypes.getUniqueKmerMultiplicity(kmer_idx, diplotypes.at(sample_idx), samples.at(sample_idx).gender) == 0) {

				auto kmer_counts = variant_cluster_haplotypes.kmers.at(kmer_idx).counts;

				if (kmer_counts) {

					assert(!kmer_counts->hasMulticlusterOccurrence());
					noise_counts->addCount(sample_idx, kmer_counts->getSampleCount(sample_idx));

				} else {

					noise_counts->addCount(sample_idx, 0);
				}
			}
		}
	}
}

void VariantClusterGenotyper::sampleHaplotypeFrequencies() {

  	assert((haplotype_frequency_distribution->numHaplotypeCount() + haplotype_frequency_distribution->numMissingCount()) == (2 * samples.size()));
  	haplotype_frequency_distribution->sampleFrequencies();
}
