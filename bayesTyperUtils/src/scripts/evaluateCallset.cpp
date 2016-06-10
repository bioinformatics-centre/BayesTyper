
/*
evaluateCallset.cpp - This file is part of BayesTyper (v0.9)


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
#include <iostream>
#include <algorithm>
#include <memory>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"
#include "vcf++/FastaReader.hpp"
#include "vcf++/FastaRecord.hpp"

static const int window_size = 1000;

struct Region {

	uint start;
	uint end;

	Region() : start(0), end(0) {}
	Region(uint start_in, uint end_in) : start(start_in), end(end_in) {}
};

bool operator==(const Region & lhs, const Region & rhs) {
    
    return ((lhs.start == rhs.start) and (lhs.end == rhs.end));
}

bool operator!=(const Region & lhs, const Region & rhs) {
    
    return !(lhs == rhs);
}

struct Haplotype {

	string sequence;
	map<uint, Region> ref_coord_map;
	map<uint, Region>::const_iterator ref_coord_map_hint;
};


struct AlleleInfo {

	uint allele_idx;
	string allele_seq;
	uint ref_len;
	bool ref_call;
	
	AlleleInfo(uint allele_idx_in, string & allele_seq_in, uint ref_len_in, bool ref_call_in) : allele_idx(allele_idx_in), allele_seq(allele_seq_in), ref_len(ref_len_in), ref_call(ref_call_in) {}
};

bool operator==(const AlleleInfo & lhs, const AlleleInfo & rhs) {
    
    return ((lhs.allele_seq == rhs.allele_seq) and (lhs.ref_len == rhs.ref_len) and (lhs.ref_call == rhs.ref_call));
}

bool operator!=(const AlleleInfo & lhs, const AlleleInfo & rhs) {
    
    return !(lhs == rhs);
}

struct VariantInfo {

	const uint position;
	const uint variant_idx;

	vector<AlleleInfo> alleles;

	VariantInfo(const uint position_in, const uint variant_idx_in) : position(position_in), variant_idx(variant_idx_in) {

		alleles.reserve(2);
	}

	uint maxReferenceLength() {

		assert(!(alleles.empty()));

		uint max_ref_len = alleles.front().ref_len;

		for (uint i = 1; i < alleles.size(); i++) {

			max_ref_len = max(max_ref_len, alleles.at(i).ref_len);
		}

		return max_ref_len;
	}

	uint minReferenceLength() {

		assert(!(alleles.empty()));

		uint min_ref_len = alleles.front().ref_len;

		for (uint i = 1; i < alleles.size(); i++) {

			min_ref_len = min(min_ref_len, alleles.at(i).ref_len);
		}

		return min_ref_len;
	}
};


void updateHaplotype(Haplotype * haplotype, const uint position, const AlleleInfo & allele, const string & chromosome_sequence) {

	assert(allele.ref_len > 0);
	
	if (haplotype->ref_coord_map.empty()) {

		assert(allele.ref_len > 0);
		assert(!(allele.allele_seq.empty()));
		assert(haplotype->sequence.empty());

		Region haplotype_region(haplotype->sequence.size(), haplotype->sequence.size() + allele.allele_seq.size() - 1); 

		haplotype->sequence.append(allele.allele_seq);

		haplotype->ref_coord_map.emplace(position, haplotype_region);
		haplotype->ref_coord_map.emplace(position + allele.ref_len - 1, haplotype_region);

	} else {

		if (haplotype->ref_coord_map.rbegin()->first < position) {

			assert(allele.ref_len > 0);
			assert(!(allele.allele_seq.empty()));

			haplotype->sequence.append(chromosome_sequence.substr(haplotype->ref_coord_map.rbegin()->first + 1, position - haplotype->ref_coord_map.rbegin()->first - 1));

			Region haplotype_region(haplotype->sequence.size(), haplotype->sequence.size() + allele.allele_seq.size() - 1); 

			haplotype->sequence.append(allele.allele_seq);

			haplotype->ref_coord_map.emplace(position, haplotype_region);
			haplotype->ref_coord_map.emplace(position + allele.ref_len - 1, haplotype_region);
		}
	}
}

vector<VariantInfo> getCalledVariantInfo(const vector<Variant *> & variants, const string & sample_id, pair<Haplotype, Haplotype> * haplotypes, const string & chromosome_sequence, const uint start_position, const uint end_position) {

	assert(start_position <= end_position);

	if (haplotypes) {

		haplotypes->first.sequence.append(chromosome_sequence.substr(start_position, 1));
		haplotypes->second.sequence.append(chromosome_sequence.substr(start_position, 1));

		haplotypes->first.ref_coord_map.emplace(start_position, Region(0,0));
		haplotypes->second.ref_coord_map.emplace(start_position, Region(0,0));
	}

	vector<VariantInfo> called_variants;
	called_variants.reserve(variants.size());

	uint variant_idx = 0;

	for (auto & variant: variants) {	

		if (!(called_variants.empty())) {

			assert((variant->pos() - 1) > called_variants.back().position);
		}

		auto genotype = variant->getSample(sample_id).genotypeEstimate();

		if (!(genotype.empty())) {

			called_variants.emplace_back(variant->pos() - 1, variant_idx);

			for (auto allele_idx: genotype) {

				Allele reference = variant->ref();
				Allele allele = variant->allele(allele_idx);

				Auxiliaries::rightTrimAllelePair(&reference, &allele);

				assert(!(reference.seq().empty()));
				assert(!(allele.seq().empty()));

				if (allele.isMissing()) {

					reference.seq() = reference.seq().substr(0,1);
					allele.seq().clear();
				}

				called_variants.back().alleles.emplace_back(AlleleInfo(allele_idx, allele.seq(), reference.seq().size(), reference == allele));
			}

			assert(start_position < called_variants.back().position);
			assert((called_variants.back().position + called_variants.back().maxReferenceLength() - 1) < end_position);

			if (haplotypes) {

				assert(variant->getSample(sample_id).isPhased());
				assert(called_variants.back().alleles.size() == 2);

				if (called_variants.back().position <= haplotypes->first.ref_coord_map.rbegin()->first) {

					called_variants.back().alleles.front().allele_seq.clear();
					called_variants.back().alleles.front().ref_len = 1;
				}

				if (called_variants.back().position <= haplotypes->second.ref_coord_map.rbegin()->first) {

					called_variants.back().alleles.back().allele_seq.clear();
					called_variants.back().alleles.back().ref_len = 1;
				}

				updateHaplotype(&(haplotypes->first), called_variants.back().position, called_variants.back().alleles.front(), chromosome_sequence);
				updateHaplotype(&(haplotypes->second), called_variants.back().position, called_variants.back().alleles.back(), chromosome_sequence);
			} 
		}

		variant_idx++;
	}

	if (haplotypes) {

		const uint var_end_position_1 = haplotypes->first.ref_coord_map.rbegin()->first;
		const uint var_end_position_2 = haplotypes->second.ref_coord_map.rbegin()->first;

		assert(var_end_position_1 < end_position);
		assert(var_end_position_2 < end_position);
	
		haplotypes->first.sequence.append(chromosome_sequence.substr(var_end_position_1, end_position - var_end_position_1));
		haplotypes->second.sequence.append(chromosome_sequence.substr(var_end_position_2, end_position - var_end_position_2));

		haplotypes->first.ref_coord_map.emplace(end_position, Region(haplotypes->first.sequence.size() - 1, haplotypes->first.sequence.size() - 1));
		haplotypes->second.ref_coord_map.emplace(end_position, Region(haplotypes->second.sequence.size() - 1, haplotypes->second.sequence.size() - 1));
	}

	return called_variants;
}

pair<bool, bool> isSequenceInHaplotype(const string & sequence, const Haplotype & haplotype, const Region haplotype_region) {

	if (sequence.empty()) {

		return make_pair(true, false);
	}

	assert(haplotype_region.start <= haplotype_region.end);
		
	if (haplotype.sequence.compare(haplotype_region.start, sequence.size(), sequence) == 0) {

		if (sequence.size() == (haplotype_region.end - haplotype_region.start + 1)) {

			return make_pair(true, true);
		
		} else {

			return make_pair(true, false);			
		}
	}

	return make_pair(false, false);
}


pair<typename map<uint, Region>::const_iterator, typename map<uint, Region>::const_iterator> findCoordinates(const uint reference_position, const map<uint, Region> & ref_coord_map, const map<uint, Region>::const_iterator ref_coord_map_hint) {

	assert(ref_coord_map_hint != ref_coord_map.end());
	assert(ref_coord_map_hint->first <= reference_position);

	auto coordinates = make_pair(ref_coord_map.end(), ref_coord_map.end());

	auto ref_coord_map_it = ref_coord_map.find(reference_position);

	if (ref_coord_map_it != ref_coord_map.end()) {

		assert(ref_coord_map_it != ref_coord_map.begin()); 

		auto prev_ref_coord_map_it = ref_coord_map_it;
		prev_ref_coord_map_it--;

		auto next_ref_coord_map_it = ref_coord_map_it;
		next_ref_coord_map_it++;

		assert(next_ref_coord_map_it != ref_coord_map.end());

		if (prev_ref_coord_map_it->second == ref_coord_map_it->second) {

			assert(ref_coord_map_it->second != next_ref_coord_map_it->second);

			coordinates.first = prev_ref_coord_map_it;
			coordinates.second = ref_coord_map_it;

		} else {

			coordinates.first = ref_coord_map_it;
			coordinates.second = next_ref_coord_map_it;
		}
	
	} else {

		ref_coord_map_it = ref_coord_map_hint;
		auto prev_ref_coord_map_it = ref_coord_map_it;

		while (ref_coord_map_it != ref_coord_map.end()) {

			assert(ref_coord_map_it->first != reference_position);

			if (reference_position < ref_coord_map_it->first) {

				coordinates.first = prev_ref_coord_map_it;
				coordinates.second = ref_coord_map_it;

				break;
			}

			prev_ref_coord_map_it = ref_coord_map_it;
			ref_coord_map_it++;
		}
	}

	assert(coordinates.first != coordinates.second);
	assert(coordinates.first != ref_coord_map.end());

	assert(coordinates.first->first <= reference_position);
	assert(reference_position <= coordinates.second->first);

	assert(coordinates.first->second.start <= coordinates.first->second.end);	
	assert(coordinates.second->second.start <= coordinates.second->second.end);	

	return coordinates;
}

Region getHaplotypeRegion(const uint reference_position, const Haplotype & haplotype) {
		
	auto coordinates = findCoordinates(reference_position, haplotype.ref_coord_map, haplotype.ref_coord_map_hint);
	assert(coordinates.first != coordinates.second);

	if (coordinates.first->second == coordinates.second->second) {

		return Region({coordinates.first->second.start, coordinates.first->second.end});	
	}  

	assert(coordinates.second->first != reference_position);

	if (coordinates.first->first == reference_position) {

		return Region({coordinates.first->second.start, coordinates.first->second.end});		
	}

	assert((coordinates.second->first - coordinates.first->first) == (coordinates.second->second.start - coordinates.first->second.end));

	const uint haplotype_position = coordinates.first->second.end + (reference_position - coordinates.first->first);
	return Region({haplotype_position, haplotype_position});		
}

Region getReferenceRegion(const uint reference_position, const Haplotype & haplotype) {

	auto coordinates = findCoordinates(reference_position, haplotype.ref_coord_map, haplotype.ref_coord_map_hint);
	assert(coordinates.first != coordinates.second);

	if (coordinates.first->second == coordinates.second->second) {

		return Region({coordinates.first->first, coordinates.second->first});
	
	} else {

		assert((coordinates.second->first - coordinates.first->first) == (coordinates.second->second.start - coordinates.first->second.end));
		return Region({reference_position, reference_position});		
	}
}

void updateRegionHint(const uint reference_position, Haplotype * haplotype) {

	auto coordinates = findCoordinates(reference_position, haplotype->ref_coord_map, haplotype->ref_coord_map_hint);
	assert(coordinates.first != coordinates.second);

	haplotype->ref_coord_map_hint = coordinates.first;
}

bool isGroundTruthHaplotypeInCallset(vector<VariantInfo>::iterator cs_variant_it, const vector<VariantInfo>::iterator cs_variant_end, Haplotype * cs_haplotype, const Haplotype & gt_haplotype, const Region & gt_reference_region, Region gt_haplotype_region, const string chromosome_sequence, const bool is_first) {

	assert(cs_haplotype->sequence.empty() == cs_haplotype->ref_coord_map.empty());		
	assert(gt_haplotype.sequence.empty() == gt_haplotype.ref_coord_map.empty());

	uint cs_allele_end = 0;

	if (!(cs_haplotype->ref_coord_map.empty())) {

		cs_allele_end = cs_haplotype->ref_coord_map.rbegin()->first;

		if (cs_variant_it != cs_variant_end) {

			if (cs_variant_it->position < gt_reference_region.start) {

				gt_haplotype_region.start = min(gt_haplotype_region.start, getHaplotypeRegion(cs_variant_it->position, gt_haplotype).start);
			} 
		}

		if (gt_reference_region.end < cs_allele_end) {

			gt_haplotype_region.end = max(gt_haplotype_region.end, getHaplotypeRegion(cs_allele_end, gt_haplotype).end);
		} 
	}  

	auto is_cs_sequence_found = isSequenceInHaplotype(cs_haplotype->sequence, gt_haplotype, gt_haplotype_region);
	bool is_cs_sequence_complete = ((is_cs_sequence_found.second) and (gt_reference_region.end <= cs_allele_end));

	if (!is_cs_sequence_found.first) {

		assert(!is_cs_sequence_found.second);
		assert(!is_cs_sequence_complete);

		return is_cs_sequence_complete;
	}

	if (is_cs_sequence_complete) {

		assert(is_cs_sequence_found.first);
		assert(!(cs_haplotype->ref_coord_map.empty()));
			
		return is_cs_sequence_complete;
	}

	auto next_cs_variant_it = cs_variant_it;

	if (!is_first) {

		next_cs_variant_it++;
	}

	assert(!is_cs_sequence_complete);

	while (next_cs_variant_it != cs_variant_end) {

		if (next_cs_variant_it->position <= cs_allele_end) {

			next_cs_variant_it++;	
			continue;	
		}

		if ((next_cs_variant_it->position + next_cs_variant_it->maxReferenceLength() - 1 + window_size) < gt_reference_region.start) {

			next_cs_variant_it++;	
			continue;
		}

		break;	
	}

	if (next_cs_variant_it == cs_variant_end) {

		if (cs_haplotype->ref_coord_map.empty()) {

			cs_haplotype->sequence = chromosome_sequence.at(gt_reference_region.start);
			cs_haplotype->ref_coord_map.emplace(gt_reference_region.start, Region(0,0));	

			is_cs_sequence_complete = isGroundTruthHaplotypeInCallset(next_cs_variant_it, cs_variant_end, cs_haplotype, gt_haplotype, gt_reference_region, gt_haplotype_region, chromosome_sequence, true);
		}
	
		return is_cs_sequence_complete;
	}

	if ((gt_reference_region.end + window_size) < next_cs_variant_it->position) {

		return is_cs_sequence_complete;
	}

	Haplotype cur_cs_haplotype = *cs_haplotype;
	
	if (cur_cs_haplotype.ref_coord_map.empty()) {

		if ((next_cs_variant_it->position + next_cs_variant_it->minReferenceLength() - 1) < gt_reference_region.start) { 	

		 	is_cs_sequence_complete = isGroundTruthHaplotypeInCallset(next_cs_variant_it, cs_variant_end, cs_haplotype, gt_haplotype, gt_reference_region, gt_haplotype_region, chromosome_sequence, false);

			if (is_cs_sequence_complete) {

				return is_cs_sequence_complete;
			}

		} else if (gt_reference_region.start < next_cs_variant_it->position) {

			cs_haplotype->sequence = chromosome_sequence.at(gt_reference_region.start);
			cs_haplotype->ref_coord_map.emplace(gt_reference_region.start, Region(0,0));	

			return isGroundTruthHaplotypeInCallset(next_cs_variant_it, cs_variant_end, cs_haplotype, gt_haplotype, gt_reference_region, gt_haplotype_region, chromosome_sequence, true);
		} 
	}

	auto allele_it = next_cs_variant_it->alleles.begin();

	while (allele_it != next_cs_variant_it->alleles.end()) {

		auto prev_allele_it = next_cs_variant_it->alleles.begin();

		while (prev_allele_it != allele_it) {

			if (*prev_allele_it == *allele_it) {

				break;
			}

			prev_allele_it++;
		}

		if (prev_allele_it != allele_it) {

			allele_it++;
			continue;
		}

		if ((next_cs_variant_it->position + allele_it->ref_len - 1 + window_size) < gt_reference_region.start) {

			allele_it++;
			continue;
		}

		if (cur_cs_haplotype.ref_coord_map.empty() and allele_it->allele_seq.empty()) {

			allele_it++;
			continue;
		}	

		if (!(allele_it->allele_seq.empty() and (next_cs_variant_it->position > cs_allele_end))) {

			*cs_haplotype = cur_cs_haplotype;

			updateHaplotype(cs_haplotype, next_cs_variant_it->position, *allele_it, chromosome_sequence);

			auto is_cs_sequence_complete = isGroundTruthHaplotypeInCallset(next_cs_variant_it, cs_variant_end, cs_haplotype, gt_haplotype, gt_reference_region, gt_haplotype_region, chromosome_sequence, false);

			if (is_cs_sequence_complete) {

				return is_cs_sequence_complete;
			}
		}

		allele_it++;
	}

	return is_cs_sequence_complete;
}

bool evaluateAllele(const uint position, const AlleleInfo & allele, Haplotype * gt_haplotype, vector<VariantInfo>::iterator * cs_variants_start_it, vector<VariantInfo> * cs_variants, const string & chromosome_sequence) {

	assert(allele.ref_len > 0);

	const Region gt_reference_region(getReferenceRegion(position, *gt_haplotype).start, getReferenceRegion(position + allele.ref_len - 1, *gt_haplotype).end);

	assert(gt_reference_region.start <= gt_reference_region.end);
	assert(gt_reference_region.start <= position);
	assert((position + allele.ref_len - 1) <= gt_reference_region.end);

	auto cs_variants_it = *cs_variants_start_it;

	while (cs_variants_it != cs_variants->end()) {

		if (gt_reference_region.start <= (cs_variants_it->position + cs_variants_it->maxReferenceLength() - 1 + window_size)) {

			Haplotype cs_haplotype;
			Region gt_haplotype_region;

			updateRegionHint(min(gt_reference_region.start, cs_variants_it->position), gt_haplotype);

			gt_haplotype_region.start = getHaplotypeRegion(gt_reference_region.start, *gt_haplotype).start;
			gt_haplotype_region.end = getHaplotypeRegion(gt_reference_region.end, *gt_haplotype).end;

			auto is_gt_haplotype_in_cs = isGroundTruthHaplotypeInCallset(cs_variants_it, cs_variants->end(), &cs_haplotype, *gt_haplotype, gt_reference_region, gt_haplotype_region, chromosome_sequence, true); 

			if (!(allele.allele_seq.empty()) and (cs_haplotype.ref_coord_map.count(position) < 1)) {

				is_gt_haplotype_in_cs = false;
			}

			return is_gt_haplotype_in_cs;
		}

		if ((gt_reference_region.end + window_size) < cs_variants_it->position) {

			break;
		}

		*cs_variants_start_it = cs_variants_it;
		cs_variants_it++;
	}

	if (allele.ref_call) {

		return true;
	}

	return false;
}

pair<int, int> calculateFlankEditDistances(const uint position, const Haplotype & haplotype, const string & chromosome_sequence, const uint read_length) {

	auto coordinates = findCoordinates(position, haplotype.ref_coord_map, haplotype.ref_coord_map_hint);

	assert(coordinates.first->first <= position);
	assert(position <= coordinates.second->first);

	if ((coordinates.first->first != position) and (coordinates.first->second == coordinates.second->second)) {

		return make_pair(-1, -1);
	}

	Region haplotype_region = Region({coordinates.first->second.start, coordinates.first->second.end});

	if (coordinates.first->first != position) {

		assert((coordinates.second->first - coordinates.first->first) == (coordinates.second->second.start - coordinates.first->second.end));
		const uint haplotype_position = coordinates.first->second.end + (position - coordinates.first->first);

		haplotype_region.start = haplotype_position; 
		haplotype_region.end = haplotype_position;
	}

	uint edit_distance_snp = 0;
	uint edit_distance_sv = 0;

	auto upstream_ref_coord_map_it = coordinates.first;

	if (coordinates.first->first == position) {

		assert(upstream_ref_coord_map_it != haplotype.ref_coord_map.begin());
		upstream_ref_coord_map_it--;
	}

	while (upstream_ref_coord_map_it != haplotype.ref_coord_map.begin()) {

		assert(upstream_ref_coord_map_it->second.end < haplotype_region.start);
		assert(upstream_ref_coord_map_it->second.start <= upstream_ref_coord_map_it->second.end);

		if ((haplotype_region.start - upstream_ref_coord_map_it->second.end) >= read_length) {

			break;
		}

		auto further_upstream_ref_coord_map_it = upstream_ref_coord_map_it;
		further_upstream_ref_coord_map_it--;

		if (upstream_ref_coord_map_it->second == further_upstream_ref_coord_map_it->second) {

			assert(further_upstream_ref_coord_map_it->first < upstream_ref_coord_map_it->first);
			edit_distance_sv++;

			upstream_ref_coord_map_it--;
		
		} else if (upstream_ref_coord_map_it->second.start < upstream_ref_coord_map_it->second.end) {

			edit_distance_sv++;

		} else {

			assert(upstream_ref_coord_map_it->second.start == upstream_ref_coord_map_it->second.end);

			if (haplotype.sequence.at(upstream_ref_coord_map_it->second.start) != chromosome_sequence.at(upstream_ref_coord_map_it->first)) {

				edit_distance_snp++;
			}
		}

		upstream_ref_coord_map_it--;
	}


	auto downstream_ref_coord_map_it = coordinates.second;

	auto ref_coord_map_last_it = haplotype.ref_coord_map.end();
	ref_coord_map_last_it--;

	if (coordinates.first->second == coordinates.second->second) {

		assert(downstream_ref_coord_map_it != ref_coord_map_last_it);
		downstream_ref_coord_map_it++;
	}

	while (downstream_ref_coord_map_it != ref_coord_map_last_it) {

		assert(haplotype_region.end < downstream_ref_coord_map_it->second.start);
		assert(downstream_ref_coord_map_it->second.start <= downstream_ref_coord_map_it->second.end);

		if ((downstream_ref_coord_map_it->second.start - haplotype_region.end) >= read_length) {

			break;
		}

		auto further_downstream_ref_coord_map_it = downstream_ref_coord_map_it;
		further_downstream_ref_coord_map_it++;

		if (downstream_ref_coord_map_it->second == further_downstream_ref_coord_map_it->second) {

			assert(downstream_ref_coord_map_it->first < further_downstream_ref_coord_map_it->first);
			edit_distance_sv++;

			downstream_ref_coord_map_it++;
	
		} else if (downstream_ref_coord_map_it->second.start < downstream_ref_coord_map_it->second.end) {

			edit_distance_sv++;

		} else {

			assert(downstream_ref_coord_map_it->second.start == downstream_ref_coord_map_it->second.end);

			if (haplotype.sequence.at(downstream_ref_coord_map_it->second.start) != chromosome_sequence.at(downstream_ref_coord_map_it->first)) {

				edit_distance_snp++;
			}
		}

		downstream_ref_coord_map_it++;
	}

	assert((edit_distance_snp + edit_distance_sv) <= (2 * (read_length - 1)));
	return make_pair(edit_distance_snp, edit_distance_sv);
}

void addAlleleStats(unordered_map<string, pair<uint, uint> > * allele_stats, const string & chromosome_id, const string & sample_id, const Auxiliaries::AlleleAttributes & allele_attributes, const pair<int, int> & flank_edit_distances, const bool is_called) {

	JoiningString allele_stats_line_elements('\t');

	allele_stats_line_elements.join(chromosome_id);
	allele_stats_line_elements.join(sample_id);

	allele_stats_line_elements.join(allele_attributes.typeStr());
	allele_stats_line_elements.join(to_string(allele_attributes.length));
	allele_stats_line_elements.join(to_string(allele_attributes.num_ambiguous));
	allele_stats_line_elements.join(to_string(allele_attributes.sv_length));
	
	allele_stats_line_elements.join(to_string(flank_edit_distances.first));
	allele_stats_line_elements.join(to_string(flank_edit_distances.second));

	auto allele_stats_emplace = allele_stats->emplace(allele_stats_line_elements.str(), make_pair(0,0));
	
	allele_stats_emplace.first->second.first++;
	allele_stats_emplace.first->second.second += is_called;
}

void addGenotypeStats(unordered_map<string, pair<uint, uint> > * genotype_stats, const string & chromosome_id, const string & sample_id, const Auxiliaries::AlleleAttributes & allele_attributes, const pair<int, int> & flank_edit_distances_1, const pair<int, int> & flank_edit_distances_2, const pair<uint, uint> & genotype_call_status, const string & genotype_type, const bool is_correct) {

	JoiningString genotype_stats_line_elements('\t');

	genotype_stats_line_elements.join(chromosome_id);
	genotype_stats_line_elements.join(sample_id);

	genotype_stats_line_elements.join(allele_attributes.typeStr());
	genotype_stats_line_elements.join(to_string(allele_attributes.length));
	genotype_stats_line_elements.join(to_string(allele_attributes.num_ambiguous));
	genotype_stats_line_elements.join(to_string(allele_attributes.sv_length));
	
	genotype_stats_line_elements.join(to_string(flank_edit_distances_1.first));
	genotype_stats_line_elements.join(to_string(flank_edit_distances_1.second));
	
	genotype_stats_line_elements.join(to_string(flank_edit_distances_2.first));
	genotype_stats_line_elements.join(to_string(flank_edit_distances_2.second));

	genotype_stats_line_elements.join(to_string(genotype_call_status.first));
	genotype_stats_line_elements.join(to_string(genotype_call_status.second));
	genotype_stats_line_elements.join(genotype_type);

	auto genotype_stats_emplace = genotype_stats->emplace(genotype_stats_line_elements.str(), make_pair(0,0));	

	genotype_stats_emplace.first->second.first++;
	genotype_stats_emplace.first->second.second += is_correct;
}

void evaluateVariants(const vector<Variant *> & gt_variants, const vector<Variant *> & cs_variants, FastaRecord * chromosome_sequence, const uint start_position, const uint end_position, const uint read_length, unordered_map<string, pair<uint, uint> > * allele_stats, unordered_map<string, pair<uint, uint> > * genotype_stats) {

	assert(start_position <= end_position);

	vector<string> sample_ids;

	if (!(gt_variants.empty())) {

		sample_ids = gt_variants.front()->sampleIds();
	
	} else if (!(cs_variants.empty())) {

		sample_ids = cs_variants.front()->sampleIds();
	}

	for (auto & sample_id: sample_ids) {

		pair<Haplotype, Haplotype> gt_haplotypes;

		auto gt_called_variants = getCalledVariantInfo(gt_variants, sample_id, &gt_haplotypes, chromosome_sequence->seq(), start_position, end_position);
		auto cs_called_variants = getCalledVariantInfo(cs_variants, sample_id, nullptr, chromosome_sequence->seq(), start_position, end_position);

		assert(gt_called_variants.size() == gt_variants.size());
		assert(cs_called_variants.size() <= cs_variants.size());

		assert(!(gt_haplotypes.first.sequence.empty()));
		assert(!(gt_haplotypes.second.sequence.empty()));

		assert(!(gt_haplotypes.first.ref_coord_map.empty()));
		assert(!(gt_haplotypes.second.ref_coord_map.empty()));

		gt_haplotypes.first.ref_coord_map_hint = gt_haplotypes.first.ref_coord_map.begin();
		gt_haplotypes.second.ref_coord_map_hint = gt_haplotypes.second.ref_coord_map.begin();

		auto gt_called_variants_it = gt_called_variants.begin();
		
		auto cs_called_variants_start_it_1 = cs_called_variants.begin();
		auto cs_called_variants_start_it_2 = cs_called_variants.begin();

		while (gt_called_variants_it != gt_called_variants.end()) {
	
			assert(gt_called_variants_it->alleles.size() == 2);

			if (!(gt_called_variants_it->alleles.front().allele_seq.empty())) {

				bool is_called = evaluateAllele(gt_called_variants_it->position, gt_called_variants_it->alleles.front(), &(gt_haplotypes.first), &cs_called_variants_start_it_1, &cs_called_variants, chromosome_sequence->seq());

				auto allele_attributes = Auxiliaries::alleleAttributes(gt_variants.at(gt_called_variants_it->variant_idx)->allele(gt_called_variants_it->alleles.front().allele_idx), gt_variants.at(gt_called_variants_it->variant_idx)->ref());
				
				auto flank_edit_distances = calculateFlankEditDistances(gt_called_variants_it->position, gt_haplotypes.first, chromosome_sequence->seq(), read_length);

				addAlleleStats(allele_stats, chromosome_sequence->id(), sample_id, allele_attributes, flank_edit_distances, is_called);
			}

			if (!(gt_called_variants_it->alleles.back().allele_seq.empty()) and (gt_called_variants_it->alleles.front() != gt_called_variants_it->alleles.back())) { 

				bool is_called = evaluateAllele(gt_called_variants_it->position, gt_called_variants_it->alleles.back(), &(gt_haplotypes.second), &cs_called_variants_start_it_2, &cs_called_variants, chromosome_sequence->seq());

				auto allele_attributes = Auxiliaries::alleleAttributes(gt_variants.at(gt_called_variants_it->variant_idx)->allele(gt_called_variants_it->alleles.back().allele_idx), gt_variants.at(gt_called_variants_it->variant_idx)->ref());

				auto flank_edit_distances = calculateFlankEditDistances(gt_called_variants_it->position, gt_haplotypes.second, chromosome_sequence->seq(), read_length);

				addAlleleStats(allele_stats, chromosome_sequence->id(), sample_id, allele_attributes, flank_edit_distances, is_called);
			}

			gt_called_variants_it++;
		}

		gt_haplotypes.first.ref_coord_map_hint = gt_haplotypes.first.ref_coord_map.begin();
		gt_haplotypes.second.ref_coord_map_hint = gt_haplotypes.second.ref_coord_map.begin();

		auto cs_called_variants_it = cs_called_variants.begin();
		
		cs_called_variants_start_it_1 = cs_called_variants.begin();
		cs_called_variants_start_it_2 = cs_called_variants.begin();

		while (cs_called_variants_it != cs_called_variants.end()) {

			auto flank_edit_distances_1 = calculateFlankEditDistances(cs_called_variants_it->position, gt_haplotypes.first, chromosome_sequence->seq(), read_length);
			auto flank_edit_distances_2 = calculateFlankEditDistances(cs_called_variants_it->position, gt_haplotypes.second, chromosome_sequence->seq(), read_length);

			bool is_correct_1 = false;
			bool is_correct_2 = false;

			const vector<AlleleInfo> cur_cs_alleles = cs_called_variants_it->alleles;

			assert(!(cur_cs_alleles.empty()));
			assert(cur_cs_alleles.size() <= 2);

			cs_called_variants_it->alleles.clear();
			cs_called_variants_it->alleles.push_back(cur_cs_alleles.front());

			is_correct_1 = evaluateAllele(cs_called_variants_it->position, cs_called_variants_it->alleles.front(), &(gt_haplotypes.first), &cs_called_variants_start_it_1, &cs_called_variants, chromosome_sequence->seq());

			auto allele_attributes_1 = Auxiliaries::alleleAttributes(cs_variants.at(cs_called_variants_it->variant_idx)->allele(cur_cs_alleles.front().allele_idx), cs_variants.at(cs_called_variants_it->variant_idx)->ref());

			if (cur_cs_alleles.size() == 1) {

				if (!is_correct_1) {

					is_correct_1 = evaluateAllele(cs_called_variants_it->position, cs_called_variants_it->alleles.front(), &(gt_haplotypes.second), &cs_called_variants_start_it_2, &cs_called_variants, chromosome_sequence->seq());
				}

				addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_1, flank_edit_distances_1, flank_edit_distances_2, make_pair(1, static_cast<uint>(is_correct_1)), allele_attributes_1.typeStr(), is_correct_1);
			
			} else {

				auto allele_attributes_2 = Auxiliaries::alleleAttributes(cs_variants.at(cs_called_variants_it->variant_idx)->allele(cur_cs_alleles.back().allele_idx), cs_variants.at(cs_called_variants_it->variant_idx)->ref());

				auto genotype_type = allele_attributes_1.typeStr();

				if (allele_attributes_1.type != allele_attributes_2.type) {

					genotype_type = "Multi";
				}

				cs_called_variants_it->alleles.clear();
				cs_called_variants_it->alleles.push_back(cur_cs_alleles.back());

				is_correct_2 = evaluateAllele(cs_called_variants_it->position, cs_called_variants_it->alleles.front(), &(gt_haplotypes.second), &cs_called_variants_start_it_2, &cs_called_variants, chromosome_sequence->seq());

				const uint genotype_correct_call_status_1 = static_cast<uint>(is_correct_1) + static_cast<uint>(is_correct_2);

				if ((is_correct_1 and is_correct_2) or (cur_cs_alleles.front() == cur_cs_alleles.back())) {

					addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_1, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_1), genotype_type, is_correct_1);
					addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_2, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_1), genotype_type, is_correct_2);

				} else {

					cs_called_variants_it->alleles.clear();
					cs_called_variants_it->alleles.push_back(cur_cs_alleles.back());

					bool swap_is_correct_1 = evaluateAllele(cs_called_variants_it->position, cs_called_variants_it->alleles.front(), &(gt_haplotypes.first), &cs_called_variants_start_it_1, &cs_called_variants, chromosome_sequence->seq());

					cs_called_variants_it->alleles.clear();
					cs_called_variants_it->alleles.push_back(cur_cs_alleles.front());

					bool swap_is_correct_2 = evaluateAllele(cs_called_variants_it->position, cs_called_variants_it->alleles.front(), &(gt_haplotypes.second), &cs_called_variants_start_it_2, &cs_called_variants, chromosome_sequence->seq());

					const uint genotype_correct_call_status_2 = static_cast<uint>(swap_is_correct_1) + static_cast<uint>(swap_is_correct_2);
					
					if (genotype_correct_call_status_1 >= genotype_correct_call_status_2) {

						addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_1, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_1), genotype_type, is_correct_1);
						addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_2, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_1), genotype_type, is_correct_2);

					} else {

						addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_2, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_2), genotype_type, swap_is_correct_1);
						addGenotypeStats(genotype_stats, chromosome_sequence->id(), sample_id, allele_attributes_1, flank_edit_distances_1, flank_edit_distances_2, make_pair(2, genotype_correct_call_status_2), genotype_type, swap_is_correct_2);
					}
				}
			}
			
			cs_called_variants_it->alleles = cur_cs_alleles;
			cs_called_variants_it++;
		}	
	}

	for (auto & var: gt_variants) {

		delete var;
	}

	for (auto & var: cs_variants) {

		delete var;
	}
}


int main(int argc, char const *argv[]) {

    if (argc != 6) {

        std::cout << "USAGE: evaluateCallset <ground_truth> <callset> <genome> <read_length> <output_prefix>" << std::endl;
        return 1;
    }

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") evaluateCallset script ..." << endl;

	GenotypedVcfFileReader vcf_reader_gt(argv[1], true);
	GenotypedVcfFileReader vcf_reader_cs(argv[2], true);

	vcf_reader_gt.metaData().infoDescriptors().clear();
	vcf_reader_cs.metaData().infoDescriptors().clear();

	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader_gt.metaData()), {"GT"});
	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader_cs.metaData()), {"GT"});

    unordered_map<string, FastaRecord*> genome_seqs;
    FastaReader genome_reader(argv[3]);
    FastaRecord * cur_fasta_rec;

    while (genome_reader.getNextRecord(&cur_fasta_rec)) {

        assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
    }

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << genome_seqs.size() << " chromosomes." << endl;

	uint read_length = stoi(argv[4]);
	
	vector<string> contig_ids_gt;
	contig_ids_gt.reserve(vcf_reader_gt.metaData().contigs().size());

	vector<string> contig_ids_cs;
	contig_ids_cs.reserve(vcf_reader_cs.metaData().contigs().size());

	for (auto & contig: vcf_reader_gt.metaData().contigs()) {

		contig_ids_gt.push_back(contig.id());
	}

	for (auto & contig: vcf_reader_cs.metaData().contigs()) {

		contig_ids_cs.push_back(contig.id());
	}

	auto contig_ids_gt_it = contig_ids_gt.begin();

	while (contig_ids_gt_it != contig_ids_gt.end()) {

		if (find(contig_ids_cs.begin(), contig_ids_cs.end(), *contig_ids_gt_it) == contig_ids_cs.end()) {

			contig_ids_gt_it = contig_ids_gt.erase(contig_ids_gt_it);
		
		} else {

			contig_ids_gt_it++;
		}
	}

	auto contig_ids_cs_it = contig_ids_cs.begin();

	while (contig_ids_cs_it != contig_ids_cs.end()) {

		if (find(contig_ids_gt.begin(), contig_ids_gt.end(), *contig_ids_cs_it) == contig_ids_gt.end()) {

			contig_ids_cs_it = contig_ids_cs.erase(contig_ids_cs_it);
		
		} else {

			contig_ids_cs_it++;
		}
	}

	assert(contig_ids_gt == contig_ids_cs);

	for (auto & contig: contig_ids_gt) {

		assert(genome_seqs.count(contig) > 0);
	}

	auto sample_ids_gt = vcf_reader_gt.metaData().sampleIds();

	for (auto &sample_id: sample_ids_gt) {

		if (!(vcf_reader_cs.metaData().hasSampleId(sample_id))) {

			vcf_reader_gt.metaData().rmSample(sample_id);
		}
	}

	auto sample_ids_cs = vcf_reader_cs.metaData().sampleIds();

	for (auto &sample_id: sample_ids_cs) {

		if (!(vcf_reader_gt.metaData().hasSampleId(sample_id))) {

			vcf_reader_cs.metaData().rmSample(sample_id);
		}
	}

	assert(vcf_reader_gt.metaData().numSamples() == vcf_reader_cs.metaData().numSamples());

	cout << "[" << Utils::getLocalTime() << "] Evaluating variants on " << contig_ids_gt.size() << " chromosomes and " << vcf_reader_gt.metaData().numSamples() << " samples (intersection bewteen input sets).\n" << endl;

	unordered_map<string, pair<uint, uint> > allele_stats;
	unordered_map<string, pair<uint, uint> > genotype_stats;

	Variant * cur_var_gt;
	Variant * cur_var_cs;

	bool cur_var_gt_parsed = vcf_reader_gt.getNextVariant(&cur_var_gt);
	bool cur_var_cs_parsed = vcf_reader_cs.getNextVariant(&cur_var_cs);

	vector<Variant *> variants_gt;
	vector<Variant *> variants_cs;

	uint num_variants_cs = 0;
	uint num_variants_gt = 0;

	uint cluster_start_position = 0;
	uint cluster_end_position = 0;

	for (auto &contig_id: contig_ids_gt) {

		cluster_start_position = 0;
		cluster_end_position = 0;

		auto chromosome_sequence_it = genome_seqs.find(contig_id);
		assert(chromosome_sequence_it != genome_seqs.end());

		assert(vcf_reader_gt.metaData().getContig(contig_id).typeStr() == vcf_reader_cs.metaData().getContig(contig_id).typeStr());

		while (cur_var_gt_parsed) {

			assert(!(cur_var_gt->ref().seq().empty()));

			if (cur_var_gt->chrom() != contig_id) {

				break;
			}

			if (!(variants_gt.empty())) {

				while (cur_var_cs_parsed) {

					assert(!(cur_var_cs->ref().seq().empty()));

					if (cur_var_cs->chrom() != contig_id) { 

						break;
					}

					if (cluster_end_position >= cur_var_cs->pos()) {

						assert(cluster_start_position <= cur_var_cs->pos());

						variants_cs.push_back(cur_var_cs);
						cluster_end_position = max(cluster_end_position, cur_var_cs->pos() + static_cast<uint>(cur_var_cs->ref().seq().size()) - 1 + window_size);

						num_variants_cs++;
						cur_var_cs_parsed = vcf_reader_cs.getNextVariant(&cur_var_cs);

					} else {

						break;
					}
				}
			}

			if (variants_gt.empty()) {

				assert(cluster_end_position < cur_var_gt->pos());
				variants_gt.push_back(cur_var_gt);

				cluster_start_position = cluster_end_position + 1;
				cluster_end_position = cur_var_gt->pos() + static_cast<uint>(cur_var_gt->ref().seq().size()) - 1 + window_size;
			
			} else if (cluster_end_position >= cur_var_gt->pos()) {

				variants_gt.push_back(cur_var_gt);
				cluster_end_position = max(cluster_end_position, cur_var_gt->pos() + static_cast<uint>(cur_var_gt->ref().seq().size()) - 1 + window_size);

			} else {

				assert(0 < cluster_start_position);
				evaluateVariants(variants_gt, variants_cs, chromosome_sequence_it->second, cluster_start_position - 1, cur_var_gt->pos() - 2, read_length, &allele_stats, &genotype_stats);

				variants_gt.clear();
				variants_cs.clear();

				assert(cluster_end_position < cur_var_gt->pos());
				assert(window_size <= cluster_end_position);

				variants_gt.push_back(cur_var_gt);

				cluster_start_position = cluster_end_position + 1 - window_size;
				cluster_end_position = cur_var_gt->pos() + static_cast<uint>(cur_var_gt->ref().seq().size()) - 1 + window_size;
			}

			num_variants_gt++;
			cur_var_gt_parsed = vcf_reader_gt.getNextVariant(&cur_var_gt);
		}

		while (cur_var_cs_parsed) {

			assert(!(cur_var_cs->ref().seq().empty()));

			if (cluster_start_position == 0) {

				assert(!(variants_gt.empty()));
				cluster_start_position = cur_var_cs->pos() - 1;
			}

			if (cur_var_cs->chrom() != contig_id) {

				break;
			}

			variants_cs.push_back(cur_var_cs);

			num_variants_cs++;
			cur_var_cs_parsed = vcf_reader_cs.getNextVariant(&cur_var_cs);
		}

		if ((!(variants_gt.empty())) or (!(variants_cs.empty()))) {

			assert(0 < cluster_start_position);
			evaluateVariants(variants_gt, variants_cs, chromosome_sequence_it->second, cluster_start_position - 1, chromosome_sequence_it->second->seq().size() - 1, read_length, &allele_stats, &genotype_stats);
		}

		variants_gt.clear();
		variants_cs.clear();

		cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig_id << endl;
	}

	assert(!cur_var_cs_parsed);
	assert(!(vcf_reader_cs.getNextVariant(&cur_var_cs)));

	assert(!cur_var_gt_parsed);
	assert(!(vcf_reader_gt.getNextVariant(&cur_var_gt)));

    for (auto & genome_seq: genome_seqs) {

        delete genome_seq.second;
    }

	cout << "\n[" << Utils::getLocalTime() << "] Evaluated " << num_variants_cs << " callset variants using " <<  num_variants_gt << " ground truth variants." << endl;


	ofstream allele_stats_writer(string(argv[5]) + "_allele_stats.txt");
	allele_stats_writer << "NumAlleles\tNumCalledAlleles\tChromosome\tSampleId\tAlleleType\tAlleleLength\tNumAmbiguous\tAlleleSVLength\tFlankEditDistanceSNP\tFlankEditDistanceSV\n";

	for (auto &stats: allele_stats) {

		allele_stats_writer << stats.second.first << "\t" << stats.second.second << "\t" << stats.first << "\n";
	}

	allele_stats_writer.close();

	ofstream genotype_stats_writer(string(argv[5]) + "_genotype_stats.txt");
	genotype_stats_writer << "NumGenotypedAlleles\tNumCorrectGenotypedAlleles\tChromosome\tSampleId\tAlleleType\tAlleleLength\tNumAmbiguous\tAlleleSVLength\tFlankEditDistanceSNP1\tFlankEditDistanceSV1\tFlankEditDistanceSNP2\tFlankEditDistanceSV2\tGenotypeCallStatus\tGenotypeCorrectCallStatus\tGenotypeType\n";

	for (auto &stats: genotype_stats) {

		genotype_stats_writer << stats.second.first << "\t" << stats.second.second << "\t" << stats.first << "\n";
	}

	genotype_stats_writer.close();


	cout << "[" << Utils::getLocalTime() << "] Finished BayesTyperUtils evaluateCallset\n" << endl;

	return 0;
}




