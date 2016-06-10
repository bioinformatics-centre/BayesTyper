
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

#include "../Eigen/Dense"

#include "boost/functional/hash.hpp"

#include "VariantClusterHaplotypes.hpp"
#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "KmerCoverage.hpp"
#include "Sample.hpp"


VariantClusterHaplotypes::~VariantClusterHaplotypes() {

	for (auto & kmer_count: kmers) {

		if (kmer_count->isEmpty()) {

			delete kmer_count;
		}
	}
}


bool VariantClusterHaplotypes::empty() {

	return (variants.empty() and nested_variant_cluster_indices.empty() and (kmer_haplotype_multiplicities.cols() == 0) and (kmer_haplotype_multiplicities.rows() == 0) and kmers.empty() and unique_kmer_indices.empty() and multicluster_kmer_indices.empty() and haplotype_variant_kmer_indices.empty() and multicluster_kmer_indices_subset.empty() and haplotype_multicluster_kmer_indices.empty() and redundant_multicluster_haplotypes.empty());
}


uchar VariantClusterHaplotypes::getDiplotypeKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype) {

	uchar diplotype_kmer_multiplicity = 0;

	if (diplotype.first != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += kmer_haplotype_multiplicities(kmer_idx, diplotype.first);
	} 

	if (diplotype.second != Utils::ushort_overflow) {

		diplotype_kmer_multiplicity += kmer_haplotype_multiplicities(kmer_idx, diplotype.second);
	}

	return diplotype_kmer_multiplicity;
}


uchar VariantClusterHaplotypes::getDiplotypeInterclusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const Utils::Sex sex) {

	return (getDiplotypeKmerMultiplicity(kmer_idx, diplotype) + kmers.at(kmer_idx)->getInterclusterMultiplicity(sex));
}


uchar VariantClusterHaplotypes::getDiplotypeMultilusterKmerMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & last_diplotype, const ushort sample_idx, const Utils::Sex sex) {

	assert(kmers.at(kmer_idx)->isMulti());

	uchar diplotype_multicluster_kmer_multiplicity = getDiplotypeInterclusterKmerMultiplicity(kmer_idx, diplotype, sex);

	diplotype_multicluster_kmer_multiplicity += static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->getMulticlusterMultiplicity(sample_idx);
	diplotype_multicluster_kmer_multiplicity -= getDiplotypeKmerMultiplicity(kmer_idx, last_diplotype);

	return diplotype_multicluster_kmer_multiplicity;
}


void VariantClusterHaplotypes::resetMulticlusterKmers(const uint max_multicluster_kmers, mt19937 * prng, const ushort num_samples) {

	if ((max_multicluster_kmers > 0) and (!(multicluster_kmer_indices.empty()))) {

		if (multicluster_kmer_indices_subset.empty()) {

			multicluster_kmer_indices_subset.reserve(max_multicluster_kmers);
		}

		if (max_multicluster_kmers < multicluster_kmer_indices.size()) {
			
			multicluster_kmer_indices_subset.clear();
			unordered_set<uint> sampled_multicluster_kmer_indices;
			
			uint max_haplotype_multicluster_kmers = floor(static_cast<float>(max_multicluster_kmers)/haplotype_multicluster_kmer_indices.size());

			for (auto & multicluster_kmer_indices: haplotype_multicluster_kmer_indices) {

				shuffle(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end(), *prng);
				uint num_cur_haplotype_multicluster_kmers = 0;

				auto multicluster_kmer_indices_it = multicluster_kmer_indices.begin();

				while (multicluster_kmer_indices_it != multicluster_kmer_indices.end()) {

					if (num_cur_haplotype_multicluster_kmers == max_haplotype_multicluster_kmers) {

						break;
					}

					auto sampled_multicluster_kmer_indices_emplace = sampled_multicluster_kmer_indices.emplace(*multicluster_kmer_indices_it);
					
					if (sampled_multicluster_kmer_indices_emplace.second) {

						multicluster_kmer_indices_subset.push_back(*sampled_multicluster_kmer_indices_emplace.first);
						num_cur_haplotype_multicluster_kmers++;
					}

					multicluster_kmer_indices_it++;
				}

				sort(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end());
			}

			shuffle(multicluster_kmer_indices.begin(), multicluster_kmer_indices.end(), *prng);
			
			auto multicluster_kmer_indices_it = multicluster_kmer_indices.begin();
			assert(multicluster_kmer_indices_it != multicluster_kmer_indices.end());

			while (multicluster_kmer_indices_subset.size() < max_multicluster_kmers) {

				if (sampled_multicluster_kmer_indices.count(*multicluster_kmer_indices_it) < 1) {

					multicluster_kmer_indices_subset.push_back(*multicluster_kmer_indices_it);
				}

				multicluster_kmer_indices_it++;
				assert(multicluster_kmer_indices_it != multicluster_kmer_indices.end());
			}
		
		} else if (multicluster_kmer_indices_subset.empty()) {

			multicluster_kmer_indices_subset = multicluster_kmer_indices;
		}
	}

	for (auto & kmer_idx: multicluster_kmer_indices_subset) {
		
		assert(kmers.at(kmer_idx)->isMulti());

		for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {

			static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->setIndex(sample_idx, Utils::uint_overflow);
		}
	}
}


bool VariantClusterHaplotypes::hasUpdatedMulticlusterMultiplicity(const ushort sample_idx) {

	for (auto & kmer_idx: multicluster_kmer_indices_subset) {

		assert(kmers.at(kmer_idx)->isMulti());

		if (static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->getIndex(sample_idx) != Utils::uint_overflow) {

			return true;
		}
	}

	return false;
}


void VariantClusterHaplotypes::updateMulticlusterMultiplicityIndices(const uint variant_cluster_index, const ushort sample_idx) {

	for (auto & kmer_idx: multicluster_kmer_indices) {

		assert(kmers.at(kmer_idx)->isMulti());

		if (static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->getIndex(sample_idx) == variant_cluster_index) {

			static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->setIndex(sample_idx, Utils::uint_overflow);
		}
	}
}


bool VariantClusterHaplotypes::isMulticlusterMultiplicityConstant(const pair<ushort,ushort> & diplotype_1, const pair<ushort,ushort> & diplotype_2) {

	if (diplotype_1 == diplotype_2) {

		return true;
	
	} else {

		if (diplotype_1.first == Utils::ushort_overflow) {

			assert(diplotype_1.second == Utils::ushort_overflow);
			assert(diplotype_2.first != Utils::ushort_overflow);

			if (diplotype_2.second != Utils::ushort_overflow) {

				return (haplotype_multicluster_kmer_indices.at(diplotype_2.first).empty() and haplotype_multicluster_kmer_indices.at(diplotype_2.second).empty());

			} else {

				return haplotype_multicluster_kmer_indices.at(diplotype_2.first).empty();
			}

		} else if (diplotype_2.first == Utils::ushort_overflow) {

			assert(diplotype_2.second == Utils::ushort_overflow);
			assert(diplotype_1.first != Utils::ushort_overflow);

			if (diplotype_1.second != Utils::ushort_overflow) {

				return (haplotype_multicluster_kmer_indices.at(diplotype_1.first).empty() and haplotype_multicluster_kmer_indices.at(diplotype_1.second).empty());

			} else {

				return haplotype_multicluster_kmer_indices.at(diplotype_1.first).empty();
			}

		} else if ((diplotype_1.second == Utils::ushort_overflow) and (diplotype_2.second == Utils::ushort_overflow)) {

			if (redundant_multicluster_haplotypes.at(diplotype_1.first).at(diplotype_2.first)) {

				return true;
			}

		} else if ((diplotype_1.second != Utils::ushort_overflow) and (diplotype_2.second != Utils::ushort_overflow)) {

			if (redundant_multicluster_haplotypes.at(diplotype_1.first).at(diplotype_2.first) and redundant_multicluster_haplotypes.at(diplotype_1.second).at(diplotype_2.second)) {

				return true;
			}
			
			if (redundant_multicluster_haplotypes.at(diplotype_1.first).at(diplotype_2.second) and redundant_multicluster_haplotypes.at(diplotype_1.second).at(diplotype_2.first)) {

				return true;
			}
		}
	}

	return false;
}


void VariantClusterHaplotypes::removeDiplotypeMulticlusterMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const ushort sample_idx, const uint variant_cluster_index) {

	assert(kmers.at(kmer_idx)->isMulti());

	static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->setIndex(sample_idx, variant_cluster_index);
	static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->reduceMulticlusterMultiplicity(sample_idx, getDiplotypeKmerMultiplicity(kmer_idx, diplotype));
}

void VariantClusterHaplotypes::removeDiplotypeMulticlusterMultiplicities(const pair<ushort,ushort> & diplotype, const ushort sample_idx, const uint variant_cluster_index) {

	vector<uint> first_dummy_indices;
	vector<uint> second_dummy_indices;

	auto first_multicluster_kmer_indices_it = first_dummy_indices.begin();
	auto first_multicluster_kmer_indices_end = first_dummy_indices.end();

	auto second_multicluster_kmer_indices_it = second_dummy_indices.begin();
	auto second_multicluster_kmer_indices_end = second_dummy_indices.end();

	if (diplotype.first != Utils::ushort_overflow) {

		first_multicluster_kmer_indices_it = haplotype_multicluster_kmer_indices.at(diplotype.first).begin();
		first_multicluster_kmer_indices_end = haplotype_multicluster_kmer_indices.at(diplotype.first).end();
	}

	if (diplotype.second != Utils::ushort_overflow) {

		assert(diplotype.first != Utils::ushort_overflow);

		second_multicluster_kmer_indices_it = haplotype_multicluster_kmer_indices.at(diplotype.second).begin();
		second_multicluster_kmer_indices_end = haplotype_multicluster_kmer_indices.at(diplotype.second).end();			
	}

	while ((first_multicluster_kmer_indices_it != first_multicluster_kmer_indices_end) or (second_multicluster_kmer_indices_it != second_multicluster_kmer_indices_end)) {

		if (first_multicluster_kmer_indices_it == first_multicluster_kmer_indices_end) {

			removeDiplotypeMulticlusterMultiplicity(*second_multicluster_kmer_indices_it, diplotype, sample_idx, variant_cluster_index);
			second_multicluster_kmer_indices_it++;
		
		} else if (second_multicluster_kmer_indices_it == second_multicluster_kmer_indices_end) {

			removeDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, sample_idx, variant_cluster_index);
			first_multicluster_kmer_indices_it++;
		
		} else if (*first_multicluster_kmer_indices_it < *second_multicluster_kmer_indices_it) {

			removeDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, sample_idx, variant_cluster_index);
			first_multicluster_kmer_indices_it++;
		
		} else if (*second_multicluster_kmer_indices_it < *first_multicluster_kmer_indices_it) {

			removeDiplotypeMulticlusterMultiplicity(*second_multicluster_kmer_indices_it, diplotype, sample_idx, variant_cluster_index);
			second_multicluster_kmer_indices_it++;
		
		} else {

			assert(*first_multicluster_kmer_indices_it == *second_multicluster_kmer_indices_it);
			removeDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, sample_idx, variant_cluster_index);

			first_multicluster_kmer_indices_it++;
			second_multicluster_kmer_indices_it++;
		}
	}

	assert(first_multicluster_kmer_indices_it == first_multicluster_kmer_indices_end);
	assert(second_multicluster_kmer_indices_it == second_multicluster_kmer_indices_end);
}


void VariantClusterHaplotypes::addDiplotypeMulticlusterMultiplicity(const uint kmer_idx, const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & last_diplotype, const ushort sample_idx, const uint variant_cluster_index) {

	assert(kmers.at(kmer_idx)->isMulti());

	uchar last_multicluster_multplicity = getDiplotypeKmerMultiplicity(kmer_idx, last_diplotype);
	uchar cur_multicluster_multplicity = getDiplotypeKmerMultiplicity(kmer_idx, diplotype);

	if (last_multicluster_multplicity == cur_multicluster_multplicity) {

		assert(static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->getIndex(sample_idx) == variant_cluster_index);
		static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->setIndex(sample_idx, Utils::uint_overflow); 
	
	} else {

		static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->setIndex(sample_idx, variant_cluster_index);									
	}
	
	static_cast<MultiClusterKmerCounts *>(kmers.at(kmer_idx))->addMulticlusterMultiplicity(sample_idx, cur_multicluster_multplicity);
}

void VariantClusterHaplotypes::addDiplotypeMulticlusterMultiplicities(const pair<ushort,ushort> & diplotype, const pair<ushort,ushort> & last_diplotype, const ushort sample_idx, const uint variant_cluster_index) {

	vector<uint> first_dummy_indices;
	vector<uint> second_dummy_indices;

	auto first_multicluster_kmer_indices_it = first_dummy_indices.begin();
	auto first_multicluster_kmer_indices_end = first_dummy_indices.end();

	auto second_multicluster_kmer_indices_it = second_dummy_indices.begin();
	auto second_multicluster_kmer_indices_end = second_dummy_indices.end();

	if (diplotype.first != Utils::ushort_overflow) {

		first_multicluster_kmer_indices_it = haplotype_multicluster_kmer_indices.at(diplotype.first).begin();
		first_multicluster_kmer_indices_end = haplotype_multicluster_kmer_indices.at(diplotype.first).end();
	}

	if (diplotype.second != Utils::ushort_overflow) {

		assert(diplotype.first != Utils::ushort_overflow);

		second_multicluster_kmer_indices_it = haplotype_multicluster_kmer_indices.at(diplotype.second).begin();
		second_multicluster_kmer_indices_end = haplotype_multicluster_kmer_indices.at(diplotype.second).end();			
	}

	while ((first_multicluster_kmer_indices_it != first_multicluster_kmer_indices_end) or (second_multicluster_kmer_indices_it != second_multicluster_kmer_indices_end)) {

		if (first_multicluster_kmer_indices_it == first_multicluster_kmer_indices_end) {

			addDiplotypeMulticlusterMultiplicity(*second_multicluster_kmer_indices_it, diplotype, last_diplotype, sample_idx, variant_cluster_index);
			second_multicluster_kmer_indices_it++;
		
		} else if (second_multicluster_kmer_indices_it == second_multicluster_kmer_indices_end) {

			addDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, last_diplotype, sample_idx, variant_cluster_index);
			first_multicluster_kmer_indices_it++;
		
		} else if (*first_multicluster_kmer_indices_it < *second_multicluster_kmer_indices_it) {

			addDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, last_diplotype, sample_idx, variant_cluster_index);
			first_multicluster_kmer_indices_it++;
		
		} else if (*second_multicluster_kmer_indices_it < *first_multicluster_kmer_indices_it) {

			addDiplotypeMulticlusterMultiplicity(*second_multicluster_kmer_indices_it, diplotype, last_diplotype, sample_idx, variant_cluster_index);
			second_multicluster_kmer_indices_it++;
		
		} else {

			assert(*first_multicluster_kmer_indices_it == *second_multicluster_kmer_indices_it);
			addDiplotypeMulticlusterMultiplicity(*first_multicluster_kmer_indices_it, diplotype, last_diplotype, sample_idx, variant_cluster_index);

			first_multicluster_kmer_indices_it++;
			second_multicluster_kmer_indices_it++;
		}
	}

	assert(first_multicluster_kmer_indices_it == first_multicluster_kmer_indices_end);
	assert(second_multicluster_kmer_indices_it == second_multicluster_kmer_indices_end);
}


KmerCoverage VariantClusterHaplotypes::getAlleleKmerCoverage(const pair<ushort, ushort> & allele_id, const ushort haplotype, const vector<pair<ushort, ushort> > & haplotype_allele_cover, const ushort num_samples) {

	KmerCoverage allele_kmer_coverage(num_samples);
	bool is_first_coverage = true;

	for (auto & haplotype_allele: haplotype_allele_cover) {

		if (haplotype == haplotype_allele.first) {

			continue;
		}

		if (allele_id.second != haplotype_allele.second) {

			allele_kmer_coverage.mergeInCoverage(getAllelePairKmerCoverage(allele_id.first, haplotype, haplotype_allele.first, num_samples), is_first_coverage);
			is_first_coverage = false;
		}
	}

	if (is_first_coverage) {

		for (auto &kmer_idx: haplotype_variant_kmer_indices.at(haplotype).at(allele_id.first)) {	

			updateAlleleKmerCoverage(&allele_kmer_coverage, kmers.at(kmer_idx), num_samples);
		}
	}

	return allele_kmer_coverage;
}


KmerCoverage VariantClusterHaplotypes::getAllelePairKmerCoverage(const ushort variant_idx, const ushort cur_haplotype, const ushort alternative_haplotype, const ushort num_samples) {	

	assert(cur_haplotype != alternative_haplotype);

	KmerCoverage allele_kmer_coverage(num_samples);

	auto cur_haplotype_kmer_indices_it = haplotype_variant_kmer_indices.at(cur_haplotype).at(variant_idx).begin();
	auto cur_haplotype_kmer_indices_end = haplotype_variant_kmer_indices.at(cur_haplotype).at(variant_idx).end();

	auto alternative_cur_haplotype_kmer_indices_it = haplotype_variant_kmer_indices.at(alternative_haplotype).at(variant_idx).begin();
	auto alternative_cur_haplotype_kmer_indices_end = haplotype_variant_kmer_indices.at(alternative_haplotype).at(variant_idx).end();

	while ((cur_haplotype_kmer_indices_it != cur_haplotype_kmer_indices_end) or (alternative_cur_haplotype_kmer_indices_it != alternative_cur_haplotype_kmer_indices_end)) {

		if (cur_haplotype_kmer_indices_it == cur_haplotype_kmer_indices_end) {

			alternative_cur_haplotype_kmer_indices_it++;
		
		} else if (alternative_cur_haplotype_kmer_indices_it == alternative_cur_haplotype_kmer_indices_end) {

			updateAlleleKmerCoverage(&allele_kmer_coverage, kmers.at(*cur_haplotype_kmer_indices_it), num_samples);
			cur_haplotype_kmer_indices_it++;
		
		} else if (*cur_haplotype_kmer_indices_it < *alternative_cur_haplotype_kmer_indices_it) {

			updateAlleleKmerCoverage(&allele_kmer_coverage, kmers.at(*cur_haplotype_kmer_indices_it), num_samples);
			cur_haplotype_kmer_indices_it++;
		
		} else if (*alternative_cur_haplotype_kmer_indices_it < *cur_haplotype_kmer_indices_it) {

			alternative_cur_haplotype_kmer_indices_it++;
		
		} else {

			assert(*cur_haplotype_kmer_indices_it == *alternative_cur_haplotype_kmer_indices_it);

			if (kmer_haplotype_multiplicities(*cur_haplotype_kmer_indices_it, cur_haplotype) != kmer_haplotype_multiplicities(*cur_haplotype_kmer_indices_it, alternative_haplotype)) {

				updateAlleleKmerCoverage(&allele_kmer_coverage, kmers.at(*cur_haplotype_kmer_indices_it), num_samples);
			}

			cur_haplotype_kmer_indices_it++;
			alternative_cur_haplotype_kmer_indices_it++;
		}
	}

	assert(cur_haplotype_kmer_indices_it == cur_haplotype_kmer_indices_end);
	assert(alternative_cur_haplotype_kmer_indices_it == alternative_cur_haplotype_kmer_indices_end);

	return allele_kmer_coverage;
}


void VariantClusterHaplotypes::updateAlleleKmerCoverage(KmerCoverage * allele_kmer_coverage, KmerCounts * kmer_count, const ushort num_samples) {	

	for (ushort sample_idx = 0; sample_idx < num_samples; sample_idx++) {	

		allele_kmer_coverage->addKmer(sample_idx, kmer_count);
	}
}

