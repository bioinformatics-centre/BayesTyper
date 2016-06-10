
/*
VariantClusterGenotyper.cpp - This file is part of BayesTyper (v0.9)


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
#include "LinearMap.hpp"
#include "VariantKmerStats.hpp"
#include "VariantInfo.hpp"

using namespace std;

VariantClusterGenotyper::VariantClusterGenotyper(const vector<Sample> & samples_in, const uint variant_cluster_index_in) : samples(samples_in), variant_cluster_index(variant_cluster_index_in) {

	has_complex_region = false;
	has_intercluster_kmer = false;
	has_multicluster_kmer = false;
	has_excluded_kmer = false;
	has_redundant_sequence = false;

	use_multicluster_kmers = false;
	multicluster_update_info = vector<MulticlusterUpdateInfo>(samples.size());
}


VariantClusterGenotyper::~VariantClusterGenotyper() {

	delete haplotype_frequency_distribution;
}


void VariantClusterGenotyper::initialise(const uint prng_seed, const uint seed_offset, KmerHash * kmer_hash, VariantClusterGraph * variant_cluster_graph, const uchar num_noise_sources, const ushort num_haplotype_candidates_per_sample, const uint max_multicluster_kmers) {

	prng = mt19937(prng_seed + seed_offset);

	noise_counts = vector<vector<unordered_map<uchar,ulong> > >(samples.size(), vector<unordered_map<uchar,ulong> >(num_noise_sources));

	unique_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());
	multicluster_diplotype_log_probabilities = vector<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > > >(samples.size());

	current_diplotypes = vector<pair<ushort,ushort> >(samples.size(), make_pair(Utils::ushort_overflow, Utils::ushort_overflow));

	variant_cluster_graph->getBestHaplotypeCandidates(kmer_hash, &variant_cluster_haplotypes, prng_seed, num_haplotype_candidates_per_sample);
	
	variant_cluster_info = variant_cluster_graph->getInfo();
	has_complex_region = variant_cluster_graph->hasComplexRegion();
	has_excluded_kmer = variant_cluster_graph->hasExcludedKmer();
	has_redundant_sequence = variant_cluster_graph->hasRedundantSequence();

	assert(variant_cluster_haplotypes.kmers.size() == uint(variant_cluster_haplotypes.kmer_haplotype_multiplicities.rows()));
	assert(variant_cluster_haplotypes.kmers.size() == (variant_cluster_haplotypes.unique_kmer_indices.size() + variant_cluster_haplotypes.multicluster_kmer_indices.size()));

	assert(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols() > 0);
	assert(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols() < Utils::ushort_overflow);

	assert(variant_cluster_haplotypes.variants.size() == uint(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols()));
	assert(variant_cluster_haplotypes.variants.size() == variant_cluster_haplotypes.nested_variant_cluster_indices.size());
	assert(variant_cluster_haplotypes.variants.size() == variant_cluster_haplotypes.haplotype_variant_kmer_indices.size());
	assert(variant_cluster_haplotypes.variants.size() == variant_cluster_haplotypes.haplotype_multicluster_kmer_indices.size());
	assert(variant_cluster_haplotypes.variants.size() == variant_cluster_haplotypes.redundant_multicluster_haplotypes.size());

	Eigen::RowVectorXbool non_zero_counts;
	non_zero_counts = Eigen::RowVectorXbool(variant_cluster_haplotypes.kmers.size());

	for (uint i = 0; i < variant_cluster_haplotypes.kmers.size(); i++) {

		non_zero_counts(i) = true;

		if (variant_cluster_haplotypes.kmers.at(i)->getInterclusterMultiplicity(Utils::Sex::Male) > 0) {

			has_intercluster_kmer = true;
			non_zero_counts(i) = false;
		} 

		if (variant_cluster_haplotypes.kmers.at(i)->isMulti()) {

			has_multicluster_kmer = true;
			non_zero_counts(i) = false;

		} 

		if (variant_cluster_haplotypes.kmers.at(i)->isEmpty()) {

			non_zero_counts(i) = false;
		}
	}

	if (has_intercluster_kmer or has_multicluster_kmer) {

		assert(variant_cluster_graph->hasNonUniqueKmer());
	}

	assert(has_multicluster_kmer == !(variant_cluster_haplotypes.multicluster_kmer_indices.empty()));
	variant_cluster_haplotypes.resetMulticlusterKmers(max_multicluster_kmers, &prng, samples.size());

	SparsityEstimator sparsity_estimator(prng_seed + seed_offset);
	uint minimum_set_cover = 0;

	if (non_zero_counts.sum() > 0) {

		minimum_set_cover = sparsity_estimator.estimateMinimumSetCover(variant_cluster_haplotypes.kmer_haplotype_multiplicities, &non_zero_counts);

	} else {

		minimum_set_cover = variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols();
	}

	assert(minimum_set_cover > 0);
	assert(minimum_set_cover <= variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols());
	assert(non_zero_counts.sum() == 0);

	if (minimum_set_cover == variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols()) {

		haplotype_frequency_distribution = new HaplotypeFrequencyDistribution(new FrequencyDistribution(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols(), prng_seed + seed_offset));

	} else {

		double sparsity = minimum_set_cover/static_cast<double>(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols());
		haplotype_frequency_distribution = new HaplotypeFrequencyDistribution(new SparseFrequencyDistribution(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols(), prng_seed + seed_offset, sparsity));
	}
}


void VariantClusterGenotyper::reset(const uint max_multicluster_kmers) {

	assert(current_diplotypes.size() == samples.size());
	assert(multicluster_update_info.size() == multicluster_update_info.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		variant_cluster_haplotypes.removeDiplotypeMulticlusterMultiplicities(current_diplotypes.at(sample_idx), sample_idx, Utils::uint_overflow);

		current_diplotypes.at(sample_idx).first = Utils::ushort_overflow;
		current_diplotypes.at(sample_idx).second = Utils::ushort_overflow;

		multicluster_update_info.at(sample_idx).update_multicluster_cache = true;
		multicluster_update_info.at(sample_idx).update_multicluster_multiplicity_indices = true;
	}

	use_multicluster_kmers = false;
	variant_cluster_haplotypes.resetMulticlusterKmers(max_multicluster_kmers, &prng, samples.size());

	for (auto & log_probabilites: multicluster_diplotype_log_probabilities) {

		log_probabilites.clear();
	}

	haplotype_frequency_distribution->reset();
}


const vector<vector<unordered_map<uchar,ulong> > > & VariantClusterGenotyper::getNoiseCounts() {

	return noise_counts;
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

	assert(current_diplotypes.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		Utils::Ploidy sample_ploidy = Utils::Ploidy::Diploid;

		if (current_diplotypes.at(sample_idx).first != Utils::ushort_overflow) {

			if (variant_cluster_haplotypes.nested_variant_cluster_indices.at(current_diplotypes.at(sample_idx).first).count(variant_cluster_idx) < 1) {

				reduceSamplePloidy(&sample_ploidy);
			}
		
		} else {

			reduceSamplePloidy(&sample_ploidy);			
		} 

		if (current_diplotypes.at(sample_idx).second != Utils::ushort_overflow) {

			if (variant_cluster_haplotypes.nested_variant_cluster_indices.at(current_diplotypes.at(sample_idx).second).count(variant_cluster_idx) < 1) {

				reduceSamplePloidy(&sample_ploidy);			
			}
		
		} else {

			reduceSamplePloidy(&sample_ploidy);
		} 

		variant_cluster_ploidy.push_back(sample_ploidy);
	}

	return variant_cluster_ploidy;
}


pair<ushort, ushort> VariantClusterGenotyper::haplotypeToAlleleId(const ushort haplotype_idx, const ushort variant_idx) {

	if (haplotype_idx < Utils::ushort_overflow) {

		auto variant_it = variant_cluster_haplotypes.variants.at(haplotype_idx).find(variant_idx);

		if (variant_it != variant_cluster_haplotypes.variants.at(haplotype_idx).end()) {

			assert(variant_it->second < Utils::ushort_overflow);
			return *variant_it;
		
		} else {

			return make_pair(Utils::ushort_overflow, Utils::ushort_overflow);
		}
	
	} else {

		return make_pair(Utils::ushort_overflow, Utils::ushort_overflow);
	}
}


vector<Genotypes*> VariantClusterGenotyper::getGenotypes(const vector<Utils::Ploidy> & sample_chromosome_ploidy) {

	assert(samples.size() == sample_chromosome_ploidy.size());
	assert(variant_cluster_haplotypes.variants.size() == uint(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols()));
	assert(variant_cluster_haplotypes.haplotype_variant_kmer_indices.size() == uint(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols()));

	vector<vector<pair<ushort, ushort> > > variant_haplotype_allele_cover(variant_cluster_info.size());

	for (auto &variant: variant_haplotype_allele_cover) {

		variant.reserve(variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols() / 2);
	}

	for (ushort haplotype_idx = 0; haplotype_idx < variant_cluster_haplotypes.variants.size(); haplotype_idx++) {

		for (auto &variant: variant_cluster_haplotypes.variants.at(haplotype_idx)) {

			variant_haplotype_allele_cover.at(variant.first).emplace_back(haplotype_idx, variant.second);
		}
	}

	vector<Genotypes*> estimated_genotypes(variant_cluster_info.size());

	for (ushort variant_idx = 0; variant_idx < variant_cluster_info.size(); variant_idx++) {

		auto git = new Genotypes(sample_chromosome_ploidy);

		git->filter_status = Utils::FilterStatus::PASS;

		git->num_candidates = variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols();
		git->variant_cluster_size = variant_cluster_info.size();

		git->variant_info = variant_cluster_info.at(variant_idx);
		git->variant_kmer_stats = VariantKmerStats(variant_cluster_info.at(variant_idx).num_alleles, samples.size());
	
		for (auto & haplotype_allele: variant_haplotype_allele_cover.at(variant_idx)) {

			git->covered_alleles.insert(haplotype_allele.second);
		}
	
		git->has_complex_region = has_complex_region;
		git->has_redundant_sequence = has_redundant_sequence;

		for (auto & diplotype_sampling_frequency: diplotype_sampling_frequencies) {

			assert(diplotype_sampling_frequency.first.first <= diplotype_sampling_frequency.first.second);

			auto allele_id_1 = haplotypeToAlleleId(diplotype_sampling_frequency.first.first, variant_idx); 
			auto allele_id_2 = haplotypeToAlleleId(diplotype_sampling_frequency.first.second, variant_idx);

			if (allele_id_1.first != Utils::ushort_overflow) {

				assert(allele_id_1.first == variant_idx);
				assert(allele_id_1.second != Utils::ushort_overflow);

				auto allele_kmer_coverage = variant_cluster_haplotypes.getAlleleKmerCoverage(allele_id_1, diplotype_sampling_frequency.first.first, variant_haplotype_allele_cover.at(variant_idx), samples.size());

				git->variant_kmer_stats.addKmerCoverage(allele_id_1.second, allele_kmer_coverage, diplotype_sampling_frequency.second);

				if (diplotype_sampling_frequency.first.first == diplotype_sampling_frequency.first.second) {

					assert(allele_id_1 == allele_id_2);

					git->variant_kmer_stats.addKmerCoverage(allele_id_2.second, allele_kmer_coverage, diplotype_sampling_frequency.second);
				}
			}

			if ((allele_id_2.first != Utils::ushort_overflow) and (diplotype_sampling_frequency.first.first != diplotype_sampling_frequency.first.second)) {

				assert(allele_id_2.first == variant_idx);
				assert(allele_id_2.second != Utils::ushort_overflow);

				auto allele_kmer_coverage = variant_cluster_haplotypes.getAlleleKmerCoverage(allele_id_2, diplotype_sampling_frequency.first.second, variant_haplotype_allele_cover.at(variant_idx), samples.size());

				git->variant_kmer_stats.addKmerCoverage(allele_id_2.second, allele_kmer_coverage, diplotype_sampling_frequency.second);
			}

			for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

				if (diplotype_sampling_frequency.second.at(sample_idx) > 0) {

					pair<unordered_map<pair<ushort, ushort>, double, boost::hash<pair<ushort, ushort> > >::iterator, bool> estimated_genotypes_insert;

					if (allele_id_1.second <= allele_id_2.second) {

						estimated_genotypes_insert = git->posteriors.at(sample_idx).second.emplace(pair<ushort, ushort>(allele_id_1.second, allele_id_2.second), 0);

					} else {

						estimated_genotypes_insert = git->posteriors.at(sample_idx).second.emplace(pair<ushort, ushort>(allele_id_2.second, allele_id_1.second), 0);
					}

					estimated_genotypes_insert.first->second += diplotype_sampling_frequency.second.at(sample_idx);
				}
			}
		}

		estimated_genotypes.at(variant_idx) = git;
	}

	return estimated_genotypes;
}


bool VariantClusterGenotyper::hasInterclusterKmer() {

	return has_intercluster_kmer;
}


bool VariantClusterGenotyper::hasMulticlusterKmer() {

	return has_multicluster_kmer;
}


bool VariantClusterGenotyper::hasExcludedKmer() {

	return has_excluded_kmer;
}


double VariantClusterGenotyper::calcDiplotypeLogProb(const CountDistribution & count_distribution, const ushort sample_idx, const pair<ushort, ushort> & diplotype, const pair<ushort, ushort> & last_diplotype) {

	double diplotype_log_prob = 0;

	if (diplotype.second == Utils::ushort_overflow) {

		diplotype_log_prob += log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second);

	} else {

		assert(diplotype.first != Utils::ushort_overflow);

		// Homozygote
		if (diplotype.first == diplotype.second) {

			diplotype_log_prob += 2 * log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second);

		// Heterozygote
		} else {

			diplotype_log_prob += log(2) + log(haplotype_frequency_distribution->getElementFrequency(diplotype.first).second) + log(haplotype_frequency_distribution->getElementFrequency(diplotype.second).second);
		}
	}

	auto unique_diplotype_log_probabilities_emplace = unique_diplotype_log_probabilities.at(sample_idx).emplace(diplotype, 0);

	if (unique_diplotype_log_probabilities_emplace.second) {

		assert(unique_diplotype_log_probabilities_emplace.first->second == 0);

		for (uint & kmer_idx: variant_cluster_haplotypes.unique_kmer_indices) {

			unique_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCombinedGenotypeNoiseCountLogProb(sample_idx, variant_cluster_haplotypes.getDiplotypeInterclusterKmerMultiplicity(kmer_idx, diplotype, samples.at(sample_idx).sex), variant_cluster_haplotypes.kmers.at(kmer_idx)->getCount(sample_idx));
		}
	}

	diplotype_log_prob += unique_diplotype_log_probabilities_emplace.first->second;

	if ((!(variant_cluster_haplotypes.multicluster_kmer_indices_subset.empty())) and use_multicluster_kmers) {

		auto multicluster_diplotype_log_probabilities_emplace = multicluster_diplotype_log_probabilities.at(sample_idx).emplace(diplotype, 0);

		if (multicluster_diplotype_log_probabilities_emplace.second or multicluster_update_info.at(sample_idx).update_multicluster_cache) {

			multicluster_diplotype_log_probabilities_emplace.first->second = 0;

			for (uint & kmer_idx: variant_cluster_haplotypes.multicluster_kmer_indices_subset) {

				multicluster_diplotype_log_probabilities_emplace.first->second += count_distribution.calcCombinedGenotypeNoiseCountLogProb(sample_idx, variant_cluster_haplotypes.getDiplotypeMultilusterKmerMultiplicity(kmer_idx, diplotype, last_diplotype, sample_idx, samples.at(sample_idx).sex), variant_cluster_haplotypes.kmers.at(kmer_idx)->getCount(sample_idx));
			}
		}

		diplotype_log_prob += multicluster_diplotype_log_probabilities_emplace.first->second;
	}

	return diplotype_log_prob;
}


void VariantClusterGenotyper::sampleGenotypes(const CountDistribution & count_distribution, const vector<Utils::Ploidy> & variant_cluster_ploidy, const bool collect_samples) {

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		multicluster_update_info.at(sample_idx).update_multicluster_cache = variant_cluster_haplotypes.hasUpdatedMulticlusterMultiplicity(sample_idx);

		if (multicluster_update_info.at(sample_idx).update_multicluster_multiplicity_indices) {

			variant_cluster_haplotypes.updateMulticlusterMultiplicityIndices(variant_cluster_index, sample_idx);
		}

		multicluster_update_info.at(sample_idx).update_multicluster_multiplicity_indices = false;

		auto last_sample_diplotype = current_diplotypes.at(sample_idx);

		sampleSingleGenotype(count_distribution, sample_idx, variant_cluster_ploidy.at(sample_idx), collect_samples, last_sample_diplotype);

		if (!(variant_cluster_haplotypes.multicluster_kmer_indices_subset.empty())) {

			bool is_multicluster_multiplicity_constant = variant_cluster_haplotypes.isMulticlusterMultiplicityConstant(last_sample_diplotype, current_diplotypes.at(sample_idx));

			if (!is_multicluster_multiplicity_constant) {

				multicluster_update_info.at(sample_idx).update_multicluster_multiplicity_indices = true;	

				variant_cluster_haplotypes.removeDiplotypeMulticlusterMultiplicities(last_sample_diplotype, sample_idx, variant_cluster_index);
				variant_cluster_haplotypes.addDiplotypeMulticlusterMultiplicities(current_diplotypes.at(sample_idx), last_sample_diplotype, sample_idx, variant_cluster_index);			
			}
		}
	}

	use_multicluster_kmers = true;
}


void VariantClusterGenotyper::sampleSingleGenotype(const CountDistribution & count_distribution, const ushort sample_idx, const Utils::Ploidy sample_variant_cluster_ploidy, const bool collect_samples, const pair<ushort, ushort> & last_diplotype) {

	const ushort num_alleles = variant_cluster_haplotypes.kmer_haplotype_multiplicities.cols();

	uint num_expected_diplotypes = (num_alleles * num_alleles)/2;
	LogDiscreteSampler diplotype_sampler(num_expected_diplotypes);
	
	vector<pair<ushort,ushort> > diplotypes;
	diplotypes.reserve(num_expected_diplotypes);

	if (sample_variant_cluster_ploidy == Utils::Ploidy::Diploid) {

		for (ushort i = 0; i < num_alleles; i++) {

			for (ushort j = i; j < num_alleles; j++) {

				// If any of the allele frequencies are zero, skip calculation
				if (haplotype_frequency_distribution->getElementFrequency(i).first) {

					if (haplotype_frequency_distribution->getElementFrequency(j).first) {

						double diplotype_log_prob = calcDiplotypeLogProb(count_distribution, sample_idx, pair<ushort, ushort>(i, j), last_diplotype);

						diplotype_sampler.addOutcome(diplotype_log_prob);
						diplotypes.emplace_back(i, j);
					}

				} else {

					break;
				}
			}
		}

	// Only one allelel is exclusive
	} else if (sample_variant_cluster_ploidy == Utils::Ploidy::Haploid) {

		for (ushort i = 0; i < num_alleles; i++) {

			// If any of the allele frequencies are zero, skip calculation
			if (haplotype_frequency_distribution->getElementFrequency(i).first) {

				double diplotype_log_prob = calcDiplotypeLogProb(count_distribution, sample_idx, pair<ushort, ushort>(i, Utils::ushort_overflow), last_diplotype);

				diplotype_sampler.addOutcome(diplotype_log_prob);
				diplotypes.emplace_back(i, Utils::ushort_overflow);
			}
		}

	// Both allelels are exclusive
	} else {

		assert(sample_variant_cluster_ploidy == Utils::Ploidy::Null);

		diplotype_sampler.addOutcome(0);
		diplotypes.emplace_back(Utils::ushort_overflow, Utils::ushort_overflow);
	}

	assert(diplotypes.size() > 0);
	auto sampled_diplotype = diplotypes.at(diplotype_sampler.sample(prng));

	if (collect_samples) {

		auto diplotype_sampling_frequencies_emplace_result = diplotype_sampling_frequencies.insert({sampled_diplotype, vector<double>(samples.size(), 0)});
		diplotype_sampling_frequencies_emplace_result.first->second.at(sample_idx)++;
	}

	current_diplotypes.at(sample_idx) = sampled_diplotype;

	haplotype_frequency_distribution->incrementObservationCount(sampled_diplotype.first);
	haplotype_frequency_distribution->incrementObservationCount(sampled_diplotype.second);
}


void VariantClusterGenotyper::sampleCountAllocations(const CountDistribution & count_distribution) {

	assert(!has_intercluster_kmer);
	assert(!has_multicluster_kmer);
	assert(!has_excluded_kmer);

	assert(variant_cluster_haplotypes.multicluster_kmer_indices.empty());

	assert(unique_diplotype_log_probabilities.size() == samples.size());
	assert(multicluster_diplotype_log_probabilities.size() == samples.size());
	assert(noise_counts.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		unique_diplotype_log_probabilities.at(sample_idx).clear();
		multicluster_diplotype_log_probabilities.at(sample_idx).clear();

		for (auto & sample_noise_counts: noise_counts.at(sample_idx)) {

			sample_noise_counts.clear();
		}

		for (uint kmer_idx = 0; kmer_idx < variant_cluster_haplotypes.kmers.size(); kmer_idx++) {

			allocateKmerCounts(variant_cluster_haplotypes.kmers.at(kmer_idx)->getCount(sample_idx), variant_cluster_haplotypes.getDiplotypeInterclusterKmerMultiplicity(kmer_idx, current_diplotypes.at(sample_idx), samples.at(sample_idx).sex), sample_idx, count_distribution);
		}
	}
}


void VariantClusterGenotyper::allocateKmerCounts(const uchar kmer_count, const uchar kmer_multiplicity, const ushort sample_idx, const CountDistribution & count_distribution) {

	assert(kmer_multiplicity <= 2);

	uchar remaining_total_count = kmer_count;

	if (kmer_multiplicity > 0) {

		unique_ptr<LogDiscreteSampler> genotype_count_log_posterior = count_distribution.calcGenotypeCountLogPosterior(sample_idx, kmer_multiplicity, remaining_total_count);

		remaining_total_count -= genotype_count_log_posterior->sample(prng);
	} 

	vector<uchar> noise_indices(noise_counts.at(sample_idx).size()); 
	iota(noise_indices.begin(), noise_indices.end(), 0);

	auto noise_indices_it = noise_indices.begin();

	uchar remaining_noise_count = remaining_total_count;

	for (ushort i = 0; i < noise_indices.size() - 1; i++) {

		unique_ptr<LogDiscreteSampler> noise_split_count_log_posterior = count_distribution.calcNoiseSplitCountLogPosterior(sample_idx, *noise_indices_it, remaining_noise_count);
		uchar noise_split_count_sample = noise_split_count_log_posterior->sample(prng);

		auto emplace_result_noise_counts = noise_counts.at(sample_idx).at(*noise_indices_it).emplace(noise_split_count_sample, 0);
		emplace_result_noise_counts.first->second++;

		remaining_noise_count -= noise_split_count_sample;
		noise_indices_it++;
	}

	auto emplace_result_noise_counts = noise_counts.at(sample_idx).at(*noise_indices_it).emplace(remaining_noise_count, 0);
	emplace_result_noise_counts.first->second++;

	noise_indices_it++;
	assert(noise_indices_it == noise_indices.end());
}


void VariantClusterGenotyper::sampleHaplotypeFrequencies() {

  	assert((haplotype_frequency_distribution->sumHaplotypeCount() + haplotype_frequency_distribution->sumMissingCount()) == (2 * samples.size()));
  	haplotype_frequency_distribution->sampleFrequencies();
}
