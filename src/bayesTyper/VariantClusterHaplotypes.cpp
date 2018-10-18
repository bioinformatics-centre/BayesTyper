
/*
VariantClusterHaplotypes.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <random>

#include "boost/functional/hash.hpp"

#include "VariantClusterHaplotypes.hpp"
#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "KmerStats.hpp"
#include "Sample.hpp"


uchar VariantClusterHaplotypes::getDiplotypeKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype) {

	uchar diplotype_kmer_multiplicity = 0;

	if (diplotype.first != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_kmer_multiplicities(kmer_idx, diplotype.first);
	} 

	if (diplotype.second != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_kmer_multiplicities(kmer_idx, diplotype.second);
	}

	return diplotype_kmer_multiplicity;
}

uchar VariantClusterHaplotypes::getUniqueKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const Utils::Gender gender) {

	auto kmer_counts = kmers.at(kmer_idx).counts;
	auto multiplicity = getDiplotypeKmerMultiplicity(kmer_idx, diplotype);

	if (kmer_counts) {

		assert(!kmer_counts->hasMulticlusterOccurrence());
		multiplicity += kmer_counts->getInterclusterMultiplicity(gender);
	}

	return multiplicity;
}

uchar VariantClusterHaplotypes::getMulticlusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx, const Utils::Gender gender) {

	auto kmer_counts = kmers.at(kmer_idx).counts;

	assert(kmer_counts);
	assert(kmer_counts->hasMulticlusterOccurrence());

	assert(getDiplotypeKmerMultiplicity(kmer_idx, prev_diplotype) <= kmer_counts->getSampleMultiplicity(sample_idx));

	if (kmer_counts->getSampleCount(sample_idx) == 0) {

		return (getDiplotypeKmerMultiplicity(kmer_idx, diplotype) + kmer_counts->getInterclusterMultiplicity(gender));

	} else {

		return (kmer_counts->getSampleMultiplicity(sample_idx) - getDiplotypeKmerMultiplicity(kmer_idx, prev_diplotype) + getDiplotypeKmerMultiplicity(kmer_idx, diplotype) + kmer_counts->getInterclusterMultiplicity(gender));	
	}
}

uchar VariantClusterHaplotypes::getPreviousMulticlusterKmerMultiplicity(const uint kmer_subset_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx, const Utils::Gender gender) {

	const uint kmer_idx = multicluster_kmer_subset_indices.at(kmer_subset_idx);

	auto kmer_counts = kmers.at(kmer_idx).counts;

	assert(kmer_counts);
	assert(kmer_counts->hasMulticlusterOccurrence());

	assert(kmer_counts->getSampleCount(sample_idx) > 0);
	assert(getDiplotypeKmerMultiplicity(kmer_idx, prev_diplotype) <= sample_multicluster_kmer_multiplicities(kmer_subset_idx, sample_idx));

	return (sample_multicluster_kmer_multiplicities(kmer_subset_idx, sample_idx) - getDiplotypeKmerMultiplicity(kmer_idx, prev_diplotype) + getDiplotypeKmerMultiplicity(kmer_idx, diplotype) + kmer_counts->getInterclusterMultiplicity(gender));
}

void VariantClusterHaplotypes::sampleKmerSubset(mt19937 * prng, const float kmer_subsampling_rate, const uint max_haplotype_variant_kmers, const ushort num_samples) {
	
	assert(unique_kmer_subset_indices.size() <= unique_kmer_indices.size());
	assert(multicluster_kmer_subset_indices.size() <= multicluster_kmer_indices.size());

	unique_kmer_subset_indices.clear();
	multicluster_kmer_subset_indices.clear();

	bernoulli_distribution bernoulli_dist(kmer_subsampling_rate);

	vector<vector<uint> > num_haplotype_variant_subset_kmers(haplotypes.size(), vector<uint>(haplotypes.front().variant_allele_indices.size(), 0));

	shuffle(unique_kmer_indices.begin(), unique_kmer_indices.end(), *prng);

	for (auto & kmer_idx: unique_kmer_indices) {

		if (bernoulli_dist(*prng)) {

			if (!isMaxHaplotypeVariantKmer(&num_haplotype_variant_subset_kmers, max_haplotype_variant_kmers, kmers.at(kmer_idx).variant_haplotype_indices)) {

				unique_kmer_subset_indices.push_back(kmer_idx);	
			}
		}
	}

	shuffle(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end(), *prng);

	for (auto & kmer_idx: multicluster_kmer_indices) {

		if (bernoulli_dist(*prng)) {

			if (!isMaxHaplotypeVariantKmer(&num_haplotype_variant_subset_kmers, max_haplotype_variant_kmers, kmers.at(kmer_idx).variant_haplotype_indices)) {

				multicluster_kmer_subset_indices.push_back(kmer_idx);	
			}
		}
	}

	assert(unique_kmer_subset_indices.size() <= unique_kmer_indices.size());
	assert(multicluster_kmer_subset_indices.size() <= multicluster_kmer_indices.size());

    sample_multicluster_kmer_multiplicities = Utils::MatrixXuchar::Zero(multicluster_kmer_subset_indices.size(), num_samples);

	for (ushort sample_idx = 0; sample_idx < kmer_stats_cache.size(); sample_idx++) {

		kmer_stats_cache.at(sample_idx).update = true;
	}
}

bool VariantClusterHaplotypes::isMaxHaplotypeVariantKmer(vector<vector<uint> > * num_haplotype_variant_subset_kmers, const uint max_haplotype_variant_kmers, const vector<pair<ushort, vector<bool> > > & variant_haplotype_indices) {

	bool is_max_haplotype_variant_kmer = true;

	for (auto & variant_haplotype_idx: variant_haplotype_indices) {

		assert(variant_haplotype_idx.first != Utils::ushort_overflow);

		for (ushort haplotype_idx = 0; haplotype_idx < variant_haplotype_idx.second.size(); haplotype_idx++) {

			if (variant_haplotype_idx.second.at(haplotype_idx) and (num_haplotype_variant_subset_kmers->at(haplotype_idx).at(variant_haplotype_idx.first) < max_haplotype_variant_kmers)) {

				num_haplotype_variant_subset_kmers->at(haplotype_idx).at(variant_haplotype_idx.first)++;
				is_max_haplotype_variant_kmer = false;
			}
		}
	}

	return is_max_haplotype_variant_kmer;
}

bool VariantClusterHaplotypes::isMulticlusterKmerUpdated(const uint kmer_subset_idx, const ushort sample_idx) {

	const uint kmer_idx = multicluster_kmer_subset_indices.at(kmer_subset_idx);

	auto kmer_counts = kmers.at(kmer_idx).counts;

	assert(kmer_counts);
	assert(kmer_counts->hasMulticlusterOccurrence());

	if ((kmer_counts->getSampleCount(sample_idx) > 0) and (kmer_counts->getSampleMultiplicity(sample_idx) != sample_multicluster_kmer_multiplicities(kmer_subset_idx, sample_idx))) {

		return true;	
	}
	
	return false;	
}

void VariantClusterHaplotypes::updateMulticlusterKmerMultiplicities(const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx) {

	if (diplotype != prev_diplotype) {

		kmer_stats_cache.at(sample_idx).update = true;

		for (auto & kmer_idx: multicluster_kmer_indices) {

			uchar cur_multicluster_multiplicity = getDiplotypeKmerMultiplicity(kmer_idx, diplotype);
			uchar prev_multicluster_mulitplicity = getDiplotypeKmerMultiplicity(kmer_idx, prev_diplotype);

			if (cur_multicluster_multiplicity != prev_multicluster_mulitplicity) {

				auto kmer_counts = kmers.at(kmer_idx).counts;
				assert(kmer_counts);

				kmer_counts->reduceSampleMultiplicity(sample_idx, prev_multicluster_mulitplicity);
				kmer_counts->addSampleMultiplicity(sample_idx, cur_multicluster_multiplicity);
			}
		}
	}

	for (uint kmer_subset_idx = 0; kmer_subset_idx < multicluster_kmer_subset_indices.size(); kmer_subset_idx++) {

		const uint kmer_idx = multicluster_kmer_subset_indices.at(kmer_subset_idx);

		auto kmer_counts = kmers.at(kmer_idx).counts;
		assert(kmer_counts);

		if ((getDiplotypeKmerMultiplicity(kmer_idx, diplotype) > 0) and (kmer_counts->getSampleCount(sample_idx) > 0) and (kmer_counts->getSampleMultiplicity(sample_idx) != sample_multicluster_kmer_multiplicities(kmer_subset_idx, sample_idx))) {

			kmer_stats_cache.at(sample_idx).update = true;
		}

		sample_multicluster_kmer_multiplicities(kmer_subset_idx, sample_idx) = kmer_counts->getSampleMultiplicity(sample_idx);
	}
}

void VariantClusterHaplotypes::updateAlleleKmerStats(vector<vector<AlleleKmerStats> > * allele_kmer_stats, const vector<Sample> & samples, const vector<VariantInfo> & variant_cluster_info, const vector<NestedVariantClusterInfo> & nested_variant_cluster_info, const vector<pair<ushort,ushort> > & diplotypes) {	

	assert(kmer_stats_cache.size() == samples.size());

	for (ushort sample_idx = 0; sample_idx < samples.size(); sample_idx++) {

		auto diplotype = diplotypes.at(sample_idx);

		if (kmer_stats_cache.at(sample_idx).update) {

			kmer_stats_cache.at(sample_idx).update = false;

			assert(kmer_stats_cache.at(sample_idx).haplotype_1.size() == kmer_stats_cache.at(sample_idx).haplotype_2.size());

			for (ushort variant_idx = 0; variant_idx < kmer_stats_cache.at(sample_idx).haplotype_1.size(); variant_idx++) {

				kmer_stats_cache.at(sample_idx).haplotype_1.at(variant_idx).reset();
				kmer_stats_cache.at(sample_idx).haplotype_2.at(variant_idx).reset();
			}

			if (diplotype.first != Utils::ushort_overflow) {

				for (auto & kmer_idx: unique_kmer_subset_indices) {

					if (getDiplotypeKmerMultiplicity(kmer_idx, diplotype) > 0) {

						updateKmerStatsCache(kmers.at(kmer_idx), diplotype, sample_idx, getUniqueKmerMultiplicity(kmer_idx, diplotype, samples.at(sample_idx).gender));
					}
				}

				for (auto & kmer_idx: multicluster_kmer_subset_indices) {

					if (getDiplotypeKmerMultiplicity(kmer_idx, diplotype) > 0) {

						updateKmerStatsCache(kmers.at(kmer_idx), diplotype, sample_idx, getMulticlusterKmerMultiplicity(kmer_idx, diplotype, diplotype, sample_idx, samples.at(sample_idx).gender));
					}
				}
			
			} else {

				assert(diplotype.second == Utils::ushort_overflow);		
			}
		} 

		assert(nested_variant_cluster_info.at(sample_idx).nested_kmer_stats.size() <= 2);

		if (diplotype.first != Utils::ushort_overflow) {
	        
	        assert(nested_variant_cluster_info.at(sample_idx).nested_kmer_stats.size() <= 1);
			addHaplotypeKmerStats(allele_kmer_stats, kmer_stats_cache.at(sample_idx).haplotype_1, variant_cluster_info, sample_idx, haplotypes.at(diplotype.first).variant_allele_indices);
		} 

		if (diplotype.second != Utils::ushort_overflow) {

			assert(diplotype.first != Utils::ushort_overflow);
			assert(nested_variant_cluster_info.at(sample_idx).nested_kmer_stats.empty());

			addHaplotypeKmerStats(allele_kmer_stats, kmer_stats_cache.at(sample_idx).haplotype_2, variant_cluster_info, sample_idx, haplotypes.at(diplotype.second).variant_allele_indices);
		}

		for (auto & kmer_stats: nested_variant_cluster_info.at(sample_idx).nested_kmer_stats) {

			addNestedHaplotypeKmerStats(allele_kmer_stats, kmer_stats, variant_cluster_info, sample_idx);			
		}
	}
}

void VariantClusterHaplotypes::updateKmerStatsCache(KmerInfo & kmer_info, const pair<ushort, ushort> & diplotype, const ushort sample_idx, const uchar kmer_multiplicity) {

	assert(diplotype.first != Utils::ushort_overflow);
	assert(kmer_multiplicity > 0);

	double kmer_count = 0;

	if (kmer_info.counts) {

		kmer_count = kmer_info.counts->getSampleCount(sample_idx) / static_cast<double>(kmer_multiplicity);
	}
    
	for (auto & variant_haplotype_idx: kmer_info.variant_haplotype_indices) {
		
		assert(diplotype.first < variant_haplotype_idx.second.size());

		if (variant_haplotype_idx.second.at(diplotype.first)) {

			kmer_stats_cache.at(sample_idx).haplotype_1.at(variant_haplotype_idx.first).addValue(make_pair(kmer_count, true));
		} 

		if (diplotype.second != Utils::ushort_overflow) {

			assert(diplotype.second < variant_haplotype_idx.second.size());

			if (variant_haplotype_idx.second.at(diplotype.second)) {

				kmer_stats_cache.at(sample_idx).haplotype_2.at(variant_haplotype_idx.first).addValue(make_pair(kmer_count, true));
			}
		}
	}
}	

void VariantClusterHaplotypes::addHaplotypeKmerStats(vector<vector<AlleleKmerStats> > * allele_kmer_stats, const vector<KmerStats> & haplotype_kmer_stats_cache, const vector<VariantInfo> & variant_cluster_info, const ushort sample_idx, const vector<ushort> & haplotype_variant_allele_indices) {

	assert(haplotype_kmer_stats_cache.size() == allele_kmer_stats->size());
	assert(haplotype_variant_allele_indices.size() == allele_kmer_stats->size());
	assert(variant_cluster_info.size() == allele_kmer_stats->size());

	ushort last_non_missing_variant_idx = Utils::ushort_overflow;

	for (ushort variant_idx = 0; variant_idx < haplotype_kmer_stats_cache.size(); variant_idx++) {

		auto allele_idx = haplotype_variant_allele_indices.at(variant_idx);
		assert(allele_idx < variant_cluster_info.at(variant_idx).numberOfAlleles());

		if (variant_cluster_info.at(variant_idx).isMissing(allele_idx)) {

			assert(last_non_missing_variant_idx != Utils::ushort_overflow);
			assert(haplotype_kmer_stats_cache.at(variant_idx).getCount() == 0);

			allele_kmer_stats->at(variant_idx).at(sample_idx).addKmerStats(haplotype_kmer_stats_cache.at(last_non_missing_variant_idx), allele_idx);

		} else {

			allele_kmer_stats->at(variant_idx).at(sample_idx).addKmerStats(haplotype_kmer_stats_cache.at(variant_idx), allele_idx);
			last_non_missing_variant_idx = variant_idx;
		}
	}
}      

void VariantClusterHaplotypes::addNestedHaplotypeKmerStats(vector<vector<AlleleKmerStats> > * allele_kmer_stats, const KmerStats & kmer_stats, const vector<VariantInfo> & variant_cluster_info, const ushort sample_idx) {

	assert(variant_cluster_info.size() == allele_kmer_stats->size());

	for (ushort variant_idx = 0; variant_idx < allele_kmer_stats->size(); variant_idx++) {

		assert(variant_cluster_info.at(variant_idx).has_dependency);
		allele_kmer_stats->at(variant_idx).at(sample_idx).addKmerStats(kmer_stats, variant_cluster_info.at(variant_idx).numberOfAlleles() - 1);
	}
}      

