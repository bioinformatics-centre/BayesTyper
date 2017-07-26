
/*
VariantClusterHaplotypes.cpp - This file is part of BayesTyper (v0.9)


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


bool VariantClusterHaplotypes::empty() {

	return (haplotypes.empty() and (haplotype_unique_kmer_multiplicities.cols() == 0) and (haplotype_unique_kmer_multiplicities.rows() == 0) and (haplotype_multicluster_kmer_multiplicities.cols() == 0) and (haplotype_multicluster_kmer_multiplicities.rows() == 0) and unique_kmers.empty() and multicluster_kmers.empty() and unique_kmer_subset_indices.empty() and multicluster_kmer_subset_indices.empty() and (sample_multicluster_kmer_multiplicities.cols() == 0) and (sample_multicluster_kmer_multiplicities.rows() == 0) and variant_kmer_stats_caches.empty());
}

uchar VariantClusterHaplotypes::getDiplotypeUniqueKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype) {

	uchar diplotype_kmer_multiplicity = 0;

	if (diplotype.first != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_unique_kmer_multiplicities(kmer_idx, diplotype.first);
	} 

	if (diplotype.second != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_unique_kmer_multiplicities(kmer_idx, diplotype.second);
	}

	return diplotype_kmer_multiplicity;
}

uchar VariantClusterHaplotypes::getDiplotypeMulticlusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype) {

	uchar diplotype_kmer_multiplicity = 0;

	if (diplotype.first != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_multicluster_kmer_multiplicities(kmer_idx, diplotype.first);
	} 

	if (diplotype.second != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += haplotype_multicluster_kmer_multiplicities(kmer_idx, diplotype.second);
	}

	return diplotype_kmer_multiplicity;
}

uchar VariantClusterHaplotypes::getUniqueKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const Utils::Gender gender) {

	auto multiplicity = getDiplotypeUniqueKmerMultiplicity(kmer_idx, diplotype);

	if (unique_kmers.at(kmer_idx).counts) {

		multiplicity += unique_kmers.at(kmer_idx).counts->getInterclusterMultiplicity(gender);
	} 

	return multiplicity;
}

uchar VariantClusterHaplotypes::getMulticlusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx, const Utils::Gender gender) {

	assert(getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, prev_diplotype) <= multicluster_kmers.at(kmer_idx).counts->getSampleMultiplicity(sample_idx));

	if (multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx) == 0) {

		return (getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, diplotype) + multicluster_kmers.at(kmer_idx).counts->getInterclusterMultiplicity(gender));

	} else {

		return (multicluster_kmers.at(kmer_idx).counts->getSampleMultiplicity(sample_idx) - getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, prev_diplotype) + getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, diplotype) + multicluster_kmers.at(kmer_idx).counts->getInterclusterMultiplicity(gender));
	}
}

uchar VariantClusterHaplotypes::getPreviousMulticlusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx, const Utils::Gender gender) {

	assert(multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx) > 0);
	assert(getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, prev_diplotype) <= sample_multicluster_kmer_multiplicities(kmer_idx, sample_idx));

	return (sample_multicluster_kmer_multiplicities(kmer_idx, sample_idx) - getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, prev_diplotype) + getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, diplotype) + multicluster_kmers.at(kmer_idx).counts->getInterclusterMultiplicity(gender));
}

bool VariantClusterHaplotypes::sampleKmerSubset(mt19937 * prng, const uint max_haplotype_variant_kmers) {
	
	assert(unique_kmer_subset_indices.size() <= unique_kmers.size());
	assert(multicluster_kmer_subset_indices.size() <= multicluster_kmers.size());
	assert(!(haplotypes.empty()));

	if ((unique_kmer_subset_indices.size() == unique_kmers.size()) and (multicluster_kmer_subset_indices.size() == multicluster_kmers.size())) {

		return false;
	}

	unique_kmer_subset_indices.clear();
	multicluster_kmer_subset_indices.clear();

	vector<vector<uint> > num_haplotype_variant_subset_kmers(haplotypes.size(), vector<uint>(haplotypes.front().variant_allele_indices.size(), 0));

	vector<uint> unique_kmer_indices(unique_kmers.size());
	iota(unique_kmer_indices.begin(), unique_kmer_indices.end(), 0);

	shuffle(unique_kmer_indices.begin(), unique_kmer_indices.end(), *prng);

	for (auto & kmer_idx: unique_kmer_indices) {

		if (isKmerInSubset(&num_haplotype_variant_subset_kmers, max_haplotype_variant_kmers, unique_kmers.at(kmer_idx).variant_haplotype_indices)) {

			unique_kmer_subset_indices.push_back(kmer_idx);	
		}
	}

	vector<uint> multicluster_kmer_indices(multicluster_kmers.size());
	iota(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end(), 0);

	shuffle(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end(), *prng);

	for (auto & kmer_idx: multicluster_kmer_indices) {

		if (isKmerInSubset(&num_haplotype_variant_subset_kmers, max_haplotype_variant_kmers, multicluster_kmers.at(kmer_idx).variant_haplotype_indices)) {

			multicluster_kmer_subset_indices.push_back(kmer_idx);	
		}
	}

	assert(unique_kmer_subset_indices.size() <= unique_kmers.size());
	assert(multicluster_kmer_subset_indices.size() <= multicluster_kmers.size());

	for (ushort sample_idx = 0; sample_idx < variant_kmer_stats_caches.size(); sample_idx++) {

		variant_kmer_stats_caches.at(sample_idx).update = true;
	}

	return true;
}

bool VariantClusterHaplotypes::isKmerInSubset(vector<vector<uint> > * num_haplotype_variant_subset_kmers, const uint max_haplotype_variant_kmers, const vector<pair<ushort, vector<bool> > > & variant_haplotype_indices) {

	bool add_kmer = false;

	for (auto & variant_haplotype_idx: variant_haplotype_indices) {

		assert(variant_haplotype_idx.first != Utils::ushort_overflow);

		for (ushort haplotype_idx = 0; haplotype_idx < variant_haplotype_idx.second.size(); haplotype_idx++) {

			if (variant_haplotype_idx.second.at(haplotype_idx) and (num_haplotype_variant_subset_kmers->at(haplotype_idx).at(variant_haplotype_idx.first) < max_haplotype_variant_kmers)) {

				num_haplotype_variant_subset_kmers->at(haplotype_idx).at(variant_haplotype_idx.first)++;
				add_kmer = true;
			}
		}
	}

	return add_kmer;
}

bool VariantClusterHaplotypes::isMulticlusterKmerUpdated(const uint kmer_idx, const ushort sample_idx) {

	if ((multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx) > 0) and (multicluster_kmers.at(kmer_idx).counts->getSampleMultiplicity(sample_idx) != sample_multicluster_kmer_multiplicities(kmer_idx, sample_idx))) {

		return true;	
	}
	
	return false;	
}

void VariantClusterHaplotypes::updateMulticlusterKmerMultiplicities(const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & prev_diplotype, const ushort sample_idx) {

	if (diplotype != prev_diplotype) {

		variant_kmer_stats_caches.at(sample_idx).update = true;

		for (uint kmer_idx = 0; kmer_idx < multicluster_kmers.size(); kmer_idx++) {

			uchar cur_multicluster_multiplicity = getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, diplotype);
			uchar prev_multicluster_mulitplicity = getDiplotypeMulticlusterKmerMultiplicity(kmer_idx, prev_diplotype);

			if (cur_multicluster_multiplicity != prev_multicluster_mulitplicity) {

				multicluster_kmers.at(kmer_idx).counts->reduceSampleMultiplicity(sample_idx, prev_multicluster_mulitplicity);
				multicluster_kmers.at(kmer_idx).counts->addSampleMultiplicity(sample_idx, cur_multicluster_multiplicity);
			}
		}
	}
	
	for (auto & kmer_idx: multicluster_kmer_subset_indices) {

		if ((multicluster_kmers.at(kmer_idx).counts->getSampleCount(sample_idx) > 0) and (multicluster_kmers.at(kmer_idx).counts->getSampleMultiplicity(sample_idx) != sample_multicluster_kmer_multiplicities(kmer_idx, sample_idx))) {

			variant_kmer_stats_caches.at(sample_idx).update = true;
		}

		sample_multicluster_kmer_multiplicities(kmer_idx, sample_idx) = multicluster_kmers.at(kmer_idx).counts->getSampleMultiplicity(sample_idx);
	}
}

void VariantClusterHaplotypes::updateVariantKmerStats(vector<VariantKmerStats> * variant_kmer_stats, const pair<ushort,ushort> & diplotype, const ushort sample_idx, const Utils::Gender gender) {	

	if (variant_kmer_stats_caches.at(sample_idx).update) {

		updateVariantKmerStatsCaches(diplotype, sample_idx, gender);
		variant_kmer_stats_caches.at(sample_idx).update = false;
	} 

	assert(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.size() == variant_kmer_stats->size());
	assert(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_2.size() == variant_kmer_stats->size());

	if ((diplotype.first != Utils::ushort_overflow) and (diplotype.second != Utils::ushort_overflow)) {

		for (ushort variant_idx = 0; variant_idx < variant_kmer_stats->size(); variant_idx++) {

			variant_kmer_stats->at(variant_idx).addAlleleKmerStats(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.at(variant_idx), haplotypes.at(diplotype.first).variant_allele_indices.at(variant_idx), sample_idx);
			variant_kmer_stats->at(variant_idx).addAlleleKmerStats(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_2.at(variant_idx), haplotypes.at(diplotype.second).variant_allele_indices.at(variant_idx), sample_idx);
		}
	
	} else if (diplotype.first != Utils::ushort_overflow) {

		assert(diplotype.second == Utils::ushort_overflow);

		for (ushort variant_idx = 0; variant_idx < variant_kmer_stats->size(); variant_idx++) {

			variant_kmer_stats->at(variant_idx).addAlleleKmerStats(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.at(variant_idx), haplotypes.at(diplotype.first).variant_allele_indices.at(variant_idx), sample_idx);
		}		

	} else {

		assert(diplotype.second == Utils::ushort_overflow);
	}
}

void VariantClusterHaplotypes::updateVariantKmerStatsCaches(const pair<ushort,ushort> & diplotype, const ushort sample_idx, const Utils::Gender gender) {	

	assert(variant_kmer_stats_caches.at(sample_idx).update);
	assert(variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.size() == variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_2.size());

	for (ushort variant_idx = 0; variant_idx < variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.size(); variant_idx++) {

		variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.at(variant_idx).reset();
		variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_2.at(variant_idx).reset();
	}

	if (diplotype.first != Utils::ushort_overflow) {

		for (auto & kmer_idx: unique_kmer_subset_indices) {

			auto kmer_multiplicity = getUniqueKmerMultiplicity(kmer_idx, diplotype, gender);
			updateVariantKmerStatsCache(unique_kmers.at(kmer_idx), diplotype, sample_idx, kmer_multiplicity);
		}

		for (auto & kmer_idx: multicluster_kmer_subset_indices) {

			auto kmer_multiplicity = getMulticlusterKmerMultiplicity(kmer_idx, diplotype, diplotype, sample_idx, gender);
			updateVariantKmerStatsCache( multicluster_kmers.at(kmer_idx), diplotype, sample_idx, kmer_multiplicity);
		}
	
	} else {

		assert(diplotype.second == Utils::ushort_overflow);		
	}
}


void VariantClusterHaplotypes::updateVariantKmerStatsCache(KmerInfo & kmer_info, const pair<ushort, ushort> & diplotype, const ushort sample_idx, const float kmer_multiplicity) {

	assert(diplotype.first != Utils::ushort_overflow);

	for (auto & variant_haplotype_idx: kmer_info.variant_haplotype_indices) {
		
		assert(diplotype.first < variant_haplotype_idx.second.size());

		if (variant_haplotype_idx.second.at(diplotype.first)) {

			variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_1.at(variant_haplotype_idx.first).addKmer(kmer_info.counts, sample_idx, kmer_multiplicity);
		} 

		if (diplotype.second != Utils::ushort_overflow) {

			assert(diplotype.second < variant_haplotype_idx.second.size());

			if (variant_haplotype_idx.second.at(diplotype.second)) {

				variant_kmer_stats_caches.at(sample_idx).allele_kmer_stats_2.at(variant_haplotype_idx.first).addKmer(kmer_info.counts, sample_idx, kmer_multiplicity);
			}
		}
	}
}	


