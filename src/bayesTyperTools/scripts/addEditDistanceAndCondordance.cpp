
/*
addEditDistanceAndCondordance.cpp - This file is part of BayesTyper (v1.1)


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

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"
#include "FastaReader.hpp"
#include "FastaRecord.hpp"


pair<string, string> getAllelePair(Variant & variant, const uint allele_idx) {

	Allele ref_allele = variant.ref();
	Allele allele = variant.allele(allele_idx);

    Auxiliaries::rightTrimAllelePair(&ref_allele, &allele);

    assert(!(ref_allele.seq().empty()));
    assert(!(allele.seq().empty()));

    return make_pair(ref_allele.seq(), allele.seq());
}

vector<vector<bool> > getNestedIndex(vector<Variant *> & variants, const string sample_id) {

	if (variants.empty()) {

		return vector<vector<bool> >();
	}

	const uint ploidy = variants.front()->getSample(sample_id).genotypeEstimate().size();

	vector<vector<bool> > nested_index(ploidy, vector<bool>(variants.size(), false));
	vector<uint> max_allele_end_pos(ploidy, 0);

	for (uint i = 0; i < variants.size(); i++) {

    	Sample & sample = variants.at(i)->getSample(sample_id);

		assert(sample.ploidy() != Sample::Ploidy::Polyploid);
		assert(sample.isPhased());
		assert(sample.genotypeEstimate().size() == ploidy);

        for (uint genotyped_allele_idx = 0; genotyped_allele_idx < sample.genotypeEstimate().size(); genotyped_allele_idx++) {

			if (variants.at(i)->pos() <= max_allele_end_pos.at(genotyped_allele_idx)) {

				nested_index.at(genotyped_allele_idx).at(i) = true;

			} else {

	    		auto allele_pair = getAllelePair(*(variants.at(i)), sample.genotypeEstimate().at(genotyped_allele_idx));
	            
	            max_allele_end_pos.at(genotyped_allele_idx) = max(max_allele_end_pos.at(genotyped_allele_idx), static_cast<uint>(variants.at(i)->pos() + allele_pair.first.size() - 1));
			}
        }
	}

	return nested_index;
}

pair<uint, bool> calculateFlankEditDistance(typename vector<Variant *>::iterator variants_it, vector<Variant *> & variants, const vector<vector<bool> > & nested_index, const string sample_id, const uint genotyped_allele_idx, const uint read_length) {

	assert(variants_it != variants.end());
	assert(genotyped_allele_idx < nested_index.size());

	if (nested_index.at(genotyped_allele_idx).at(variants_it - variants.begin())) {

		return make_pair(0, false);
	}

	uint flank_edit_distance = 0;

	Sample & sample = (*variants_it)->getSample(sample_id);
	assert(sample.genotypeEstimate().size() == nested_index.size());

	auto main_allele_pair = getAllelePair(**variants_it, sample.genotypeEstimate().at(genotyped_allele_idx));

	auto variants_rit = variants.rbegin();

	while (variants_rit != variants.rend()) {

		if (*variants_rit == *variants_it) {

			break;
		}

		variants_rit++;		
	}

	assert(variants_rit != variants.rend());

	uint downstream_start_position = (*variants_rit)->pos();
	uint left_flank_length = main_allele_pair.second.size() + 1;

	variants_rit++;

	while (variants_rit != variants.rend()) {

		if (!(nested_index.at(genotyped_allele_idx).at(variants.rend() - variants_rit - 1))) {

			Sample & cur_sample = (*variants_rit)->getSample(sample_id);
			assert(cur_sample.genotypeEstimate().size() == nested_index.size());

			if (cur_sample.genotypeEstimate().at(genotyped_allele_idx) > 0) {

		    	auto cur_allele_pair = getAllelePair(*(*variants_rit), cur_sample.genotypeEstimate().at(genotyped_allele_idx));

		    	assert(((*variants_rit)->pos() + cur_allele_pair.first.size() - 1) < downstream_start_position);

				left_flank_length += downstream_start_position - ((*variants_rit)->pos() + cur_allele_pair.first.size() - 1) - 1;    	

				if (left_flank_length <= read_length) {

					flank_edit_distance++;
				}

				left_flank_length += cur_allele_pair.second.size();
				downstream_start_position = (*variants_rit)->pos();	
			}
		}

		if (read_length < left_flank_length) {

			break;
		}

		variants_rit++;
	}


	uint upstream_end_position = (*variants_it)->pos() + main_allele_pair.first.size() - 1;
	uint right_flank_length = main_allele_pair.second.size() + 1;

	variants_it++;

	while (variants_it != variants.end()) {

		if (!(nested_index.at(genotyped_allele_idx).at(variants_it - variants.begin()))) {

	    	assert(upstream_end_position < (*variants_it)->pos());

			Sample & cur_sample = (*variants_it)->getSample(sample_id);
			assert(cur_sample.genotypeEstimate().size() == nested_index.size());

			if (cur_sample.genotypeEstimate().at(genotyped_allele_idx) > 0) {

		    	auto cur_allele_pair = getAllelePair(*(*variants_it), cur_sample.genotypeEstimate().at(genotyped_allele_idx));

				right_flank_length += (*variants_it)->pos() - upstream_end_position - 1;

				if (right_flank_length <= read_length) {

					flank_edit_distance++;
				}

				right_flank_length += cur_allele_pair.second.size();
				downstream_start_position = (*variants_it)->pos() + cur_allele_pair.first.size() - 1;
			}
		}

		if (read_length < right_flank_length) {

			break;
		}

		variants_it++;
	}	

	return make_pair(flank_edit_distance, true);
}

void addMeanFlankEditDistance(Variant * cs_variant, typename vector<Variant *>::iterator gt_variants_it, vector<Variant *> * gt_variants, const vector<vector<bool> > & gt_nested_index, const string sample_id, const uint read_length) {

	float mean_flank_edit_distance = -1;

	if (gt_variants_it != gt_variants->end()) {

		if (!(gt_nested_index.empty())) {

	    	pair<uint, bool> sum_flank_edit_distance(0, true);

            for (uint genotyped_allele_idx = 0; genotyped_allele_idx < gt_nested_index.size(); genotyped_allele_idx++) {

            	auto flank_edit_distance = calculateFlankEditDistance(gt_variants_it, *gt_variants, gt_nested_index, sample_id, genotyped_allele_idx, read_length);

            	sum_flank_edit_distance.first += flank_edit_distance.first;
            	sum_flank_edit_distance.second = flank_edit_distance.second;

            	if (!(sum_flank_edit_distance.second)) {

            		break;
            	}
            }

            if (sum_flank_edit_distance.second) {

            	mean_flank_edit_distance = sum_flank_edit_distance.first/static_cast<float>(gt_nested_index.size());
            }
        }

	} 

	if (cs_variant and (gt_variants_it != gt_variants->end())) {

		if (cs_variant->pos() == (*gt_variants_it)->pos()) {

			cs_variant->getSample(sample_id).info().setValue<float>("MFED", mean_flank_edit_distance);

		} else {

			cs_variant->getSample(sample_id).info().setValue<float>("MFED", -1);
		}

		(*gt_variants_it)->getSample(sample_id).info().setValue<float>("MFED", mean_flank_edit_distance);

	} else if (cs_variant) {

		assert(Utils::floatCompare(mean_flank_edit_distance, -1));
		cs_variant->getSample(sample_id).info().setValue<float>("MFED", mean_flank_edit_distance);
	
	} else {

		assert(gt_variants_it != gt_variants->end());
		(*gt_variants_it)->getSample(sample_id).info().setValue<float>("MFED", mean_flank_edit_distance);
	} 
}

void addGenotypeConcordance(Variant * gt_variant, Variant * cs_variant, const string sample_id) {

	assert(gt_variant or cs_variant);

	if (gt_variant and cs_variant) {

		assert(gt_variant->pos() == cs_variant->pos());

    	Sample * gt_sample = &(gt_variant->getSample(sample_id));
    	Sample * cs_sample = &(cs_variant->getSample(sample_id));

		assert(gt_sample->ploidy() != Sample::Ploidy::Polyploid);
		assert(cs_sample->ploidy() != Sample::Ploidy::Polyploid);
    
    	assert(gt_sample->genotypeEstimate().size() <= 2);
    	assert(cs_sample->genotypeEstimate().size() <= 2);

    	string genotype_concordance = "I";

    	if (gt_sample->genotypeEstimate().size() == cs_sample->genotypeEstimate().size()) {

    		if (gt_sample->genotypeEstimate().size() == 2) {

    			auto gt_allele_pair_1 = getAllelePair(*gt_variant, gt_sample->genotypeEstimate().front());
    			auto cs_allele_pair_1 = getAllelePair(*cs_variant, cs_sample->genotypeEstimate().front());

    			auto gt_allele_pair_2 = getAllelePair(*gt_variant, gt_sample->genotypeEstimate().back());
    			auto cs_allele_pair_2 = getAllelePair(*cs_variant, cs_sample->genotypeEstimate().back());

    			if (((gt_allele_pair_1 == cs_allele_pair_1) and (gt_allele_pair_2 == cs_allele_pair_2)) or ((gt_allele_pair_1 == cs_allele_pair_2) and (gt_allele_pair_2 == cs_allele_pair_1)))  {

    				genotype_concordance = "T";

    			} else if (((gt_allele_pair_1 == cs_allele_pair_1) or (gt_allele_pair_2 == cs_allele_pair_2)) or ((gt_allele_pair_1 == cs_allele_pair_2) or (gt_allele_pair_2 == cs_allele_pair_1))) {

    				genotype_concordance = "P";

    			} else {

    				genotype_concordance = "F";
    			}

    		} else if (gt_sample->genotypeEstimate().size() == 1) {

    			auto gt_allele_pair_1 = getAllelePair(*gt_variant, gt_sample->genotypeEstimate().front());
    			auto cs_allele_pair_1 = getAllelePair(*cs_variant, cs_sample->genotypeEstimate().front());

    			if (gt_allele_pair_1 ==  cs_allele_pair_1) {

    				genotype_concordance = "T";

    			} else {

    				genotype_concordance = "F";
    			}
    		} 

    	} else {
    		
    		genotype_concordance = "F";
    	}

		gt_sample->info().setValue<string>("GTCO", genotype_concordance);
		cs_sample->info().setValue<string>("GTCO", genotype_concordance);

	} else {

		Variant * variant = nullptr;

		if (gt_variant) {

			variant = gt_variant;
		
		} else {

			assert(cs_variant);
			variant = cs_variant;			
		}

    	Sample * sample = &(variant->getSample(sample_id));

		assert(sample->ploidy() != Sample::Ploidy::Polyploid);    
    	assert(sample->genotypeEstimate().size() <= 2);

    	if (sample->genotypeEstimate().size() > 0) {

	    	if (*(max_element(sample->genotypeEstimate().begin(), sample->genotypeEstimate().end())) == 0) {

				sample->info().setValue<string>("GTCO", "T");
	    	
	    	} else if (*(min_element(sample->genotypeEstimate().begin(), sample->genotypeEstimate().end())) == 0) {

				sample->info().setValue<string>("GTCO", "P");	    		
	    	
	    	} else {

				sample->info().setValue<string>("GTCO", "F");	    		
	    	}

    	} else {

			sample->info().setValue<string>("GTCO", "I");    		
    	}
    }
}

void addEditDistanceAndCondordance(vector<Variant *> * gt_variants, vector<Variant *> * cs_variants, const vector<string> sample_ids, const uint read_length) {

	for (auto & sample_id: sample_ids) {

		auto gt_nested_index = getNestedIndex(*gt_variants, sample_id);

		auto gt_variants_it = gt_variants->begin();
		auto cs_variants_it = cs_variants->begin();

		uint gt_last_position = 0;
		uint cs_last_position = 0;

		while ((gt_variants_it != gt_variants->end()) or (cs_variants_it != cs_variants->end())) {

			if ((gt_variants_it != gt_variants->end()) and (cs_variants_it != cs_variants->end())) {

				if ((*gt_variants_it)->pos() == (*cs_variants_it)->pos()) {

					assert(gt_last_position < (*gt_variants_it)->pos());
					assert(cs_last_position < (*cs_variants_it)->pos());

					gt_last_position = (*gt_variants_it)->pos();
					cs_last_position = (*cs_variants_it)->pos();

					addMeanFlankEditDistance(*cs_variants_it, gt_variants_it, gt_variants, gt_nested_index, sample_id, read_length);
					addGenotypeConcordance(*gt_variants_it, *cs_variants_it, sample_id);

					gt_variants_it++;
					cs_variants_it++;
				
				} else if ((*gt_variants_it)->pos() < (*cs_variants_it)->pos()) {

					assert(gt_last_position < (*gt_variants_it)->pos());
					gt_last_position = (*gt_variants_it)->pos();

					addMeanFlankEditDistance(nullptr, gt_variants_it, gt_variants, gt_nested_index, sample_id, read_length);
					addGenotypeConcordance(*gt_variants_it, nullptr, sample_id);

					gt_variants_it++;

				} else {

					assert(cs_last_position < (*cs_variants_it)->pos());
					cs_last_position = (*cs_variants_it)->pos();

					addMeanFlankEditDistance(*cs_variants_it, gt_variants->end(), gt_variants, gt_nested_index, sample_id, read_length);
					addGenotypeConcordance(nullptr, *cs_variants_it, sample_id);

					assert(((*cs_variants_it)->pos() < (*gt_variants_it)->pos()));
					cs_variants_it++;
				}

			} else if (gt_variants_it != gt_variants->end()) {

				assert(gt_last_position < (*gt_variants_it)->pos());
				gt_last_position = (*gt_variants_it)->pos();

				addMeanFlankEditDistance(nullptr, gt_variants_it, gt_variants, gt_nested_index, sample_id, read_length);
				addGenotypeConcordance(*gt_variants_it, nullptr, sample_id);
				
				gt_variants_it++;

			} else {

				assert(cs_last_position < (*cs_variants_it)->pos());
				cs_last_position = (*cs_variants_it)->pos();

				addMeanFlankEditDistance(*cs_variants_it, gt_variants_it, gt_variants, gt_nested_index, sample_id, read_length);
				addGenotypeConcordance(nullptr, *cs_variants_it, sample_id);

				cs_variants_it++;
			}
		}
	}
}

int main(int argc, char const *argv[]) {

    if (argc != 6) {

        std::cout << "USAGE: addEditDistanceAndCondordance <ground_truth> <callset> <ground_truth_output_prefix> <callset_output_prefix> <read_length>" << std::endl;
        return 1;
    }

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") addEditDistanceAndCondordance script ...\n" << endl;

	GenotypedVcfFileReader gt_vcf_reader(argv[1], true);
	GenotypedVcfFileReader cs_vcf_reader(argv[2], true);
	
	assert(gt_vcf_reader.metaData().contigs() == cs_vcf_reader.metaData().contigs());
	assert(gt_vcf_reader.metaData().sampleIds() == cs_vcf_reader.metaData().sampleIds());

	cout << "[" << Utils::getLocalTime() << "] Adding mean flank edit distance (MFED) and genotype concordance (GTCO) to " << gt_vcf_reader.metaData().sampleIds().size() << " samples ...\n" << endl;

	auto gt_output_meta_data = gt_vcf_reader.metaData();
	auto cs_output_meta_data = cs_vcf_reader.metaData();

	gt_output_meta_data.formatDescriptors().emplace("GTCO", Attribute::DetailedDescriptor("GTCO", Attribute::Number::One, Attribute::Type::String, "Genotype concordance"));
	gt_output_meta_data.formatDescriptors().emplace("MFED", Attribute::DetailedDescriptor("MFED", Attribute::Number::One, Attribute::Type::Float, "Mean flank edit distance"));

	cs_output_meta_data.formatDescriptors().emplace("GTCO", Attribute::DetailedDescriptor("GTCO", Attribute::Number::One, Attribute::Type::String, "Genotype concordance"));
	cs_output_meta_data.formatDescriptors().emplace("MFED", Attribute::DetailedDescriptor("MFED", Attribute::Number::One, Attribute::Type::Float, "Mean flank edit distance"));

	VcfFileWriter gt_vcf_writer(string(argv[3]) + ".vcf", gt_output_meta_data, true);
	VcfFileWriter cs_vcf_writer(string(argv[4]) + ".vcf", cs_output_meta_data, true);

	const uint read_length = stoi(argv[5]);
	const vector<string> sample_ids = gt_vcf_reader.metaData().sampleIds();

	Variant * gt_cur_var;
	Variant * cs_cur_var;

	gt_vcf_reader.getNextVariant(&gt_cur_var);
	cs_vcf_reader.getNextVariant(&cs_cur_var);

	vector<Variant *> gt_variants;
	vector<Variant *> cs_variants;

	uint num_cs_variants = 0;
	uint num_gt_variants = 0;

	uint cluster_end_position = 0;

	for (auto &contig_id: gt_vcf_reader.metaData().contigs()) {

		cluster_end_position = 0;

		while (gt_cur_var) {

			assert(!(gt_cur_var->ref().seq().empty()));

			if (gt_cur_var->chrom() != contig_id.id()) {

				break;
			}

			while (cs_cur_var) {

				assert(!(cs_cur_var->ref().seq().empty()));

				if (cs_cur_var->chrom() != contig_id.id()) { 

					break;
				}

				if ((cs_variants.empty() and gt_variants.empty() and (cs_cur_var->pos() < gt_cur_var->pos())) or (cs_cur_var->pos() <= cluster_end_position)) {

					cs_variants.push_back(cs_cur_var);
					cluster_end_position = max(cluster_end_position, cs_cur_var->pos() + static_cast<uint>(cs_cur_var->ref().seq().size()) - 1 + read_length);

					num_cs_variants++;
					cs_vcf_reader.getNextVariant(&cs_cur_var);

				} else {

					break;
				}
			}

			if ((cs_variants.empty() and gt_variants.empty()) or (gt_cur_var->pos() <= cluster_end_position)) {

				gt_variants.push_back(gt_cur_var);
				cluster_end_position = max(cluster_end_position, gt_cur_var->pos() + static_cast<uint>(gt_cur_var->ref().seq().size()) - 1 + read_length);

				num_gt_variants++;
				gt_vcf_reader.getNextVariant(&gt_cur_var);

			} else {

				addEditDistanceAndCondordance(&gt_variants, &cs_variants, sample_ids, read_length);

				for (auto & gt_variant: gt_variants) {

					gt_vcf_writer.write(gt_variant);
					delete gt_variant;
				}

				for (auto & cs_variant: cs_variants) {

					cs_vcf_writer.write(cs_variant);
					delete cs_variant;
				}

				gt_variants.clear();
				cs_variants.clear();
			}
		}

		while (cs_cur_var) {

			assert(!(cs_cur_var->ref().seq().empty()));

			if (cs_cur_var->chrom() != contig_id.id()) {

				break;
			}

			cs_variants.push_back(cs_cur_var);

			num_cs_variants++;
			cs_vcf_reader.getNextVariant(&cs_cur_var);
		}

		addEditDistanceAndCondordance(&gt_variants, &cs_variants, sample_ids, read_length);

		for (auto & gt_variant: gt_variants) {

			gt_vcf_writer.write(gt_variant);
			delete gt_variant;
		}

		for (auto & cs_variant: cs_variants) {

			cs_vcf_writer.write(cs_variant);
			delete cs_variant;
		}

		gt_variants.clear();
		cs_variants.clear();

		cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig_id.id() << endl;
	}

	assert(!cs_cur_var);
	assert(!gt_cur_var);

	assert(!(cs_vcf_reader.getNextVariant(&cs_cur_var)));
	assert(!(gt_vcf_reader.getNextVariant(&gt_cur_var)));

	cout << "\n[" << Utils::getLocalTime() << "] Added MFED and GTCO to genotypes on " << num_gt_variants << " ground truth variants." << endl;
	cout << "[" << Utils::getLocalTime() << "] Added MFED and GTCO to genotypes on " << num_cs_variants << " callset variants." << endl;
	cout << endl;

	return 0;
}




