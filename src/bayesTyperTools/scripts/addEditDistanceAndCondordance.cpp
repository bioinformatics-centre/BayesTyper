
/*
addEditDistanceAndCondordance.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <iterator>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"
#include "FastaReader.hpp"
#include "FastaRecord.hpp"


struct AltAlleleInfo {

	const uint pos;
	const uint ref_len;
	const uint alt_len;

	AltAlleleInfo(const uint pos_in, const uint ref_len_in, const uint alt_len_in) : pos(pos_in), ref_len(ref_len_in), alt_len(alt_len_in) {}
};


pair<string, string> getAllelePair(Variant & variant, const uint allele_idx) {

	Allele ref_allele = variant.ref();
	Allele allele = variant.allele(allele_idx);

    Auxiliaries::rightTrimAllelePair(&ref_allele, &allele);

    assert(!ref_allele.seq().empty());
    assert(!allele.seq().empty());

    return make_pair(ref_allele.seq(), allele.seq());
}

vector<vector<AltAlleleInfo> > getHaplotypesAlleleInfo(vector<Variant *> & variants, const string sample_id) {

	if (variants.empty()) {

		return vector<vector<AltAlleleInfo> >(1, vector<AltAlleleInfo>());
	}

	const uint ploidy = variants.front()->getSample(sample_id).genotypeEstimate().size();

	vector<vector<AltAlleleInfo> > haplotype_allele_info(ploidy, vector<AltAlleleInfo>());

	for (auto & variant: variants) {

    	Sample & sample = variant->getSample(sample_id);

		assert(sample.ploidy() != Sample::Ploidy::Polyploid);
		assert(sample.isPhased());
		assert(sample.genotypeEstimate().size() == ploidy);

        for (uint genotyped_allele_idx = 0; genotyped_allele_idx < sample.genotypeEstimate().size(); genotyped_allele_idx++) {

        	if ((sample.genotypeEstimate().at(genotyped_allele_idx) > 0) and !variant->allele(sample.genotypeEstimate().at(genotyped_allele_idx)).isMissing()) {

				Allele ref_allele = variant->ref();
				Allele alt_allele = variant->allele(sample.genotypeEstimate().at(genotyped_allele_idx));

			    auto pos_shift = Auxiliaries::fullTrimAllelePair(&ref_allele, &alt_allele);

			    assert(!ref_allele.seq().empty() or !alt_allele.seq().empty());

	        	if (!haplotype_allele_info.at(genotyped_allele_idx).empty()) {

	        		assert((haplotype_allele_info.at(genotyped_allele_idx).back().pos + haplotype_allele_info.at(genotyped_allele_idx).back().ref_len) <= (variant->pos() + pos_shift.first));
	        	}

	        	haplotype_allele_info.at(genotyped_allele_idx).emplace_back(variant->pos() + pos_shift.first, ref_allele.seq().size(), alt_allele.seq().size());
	        }
        }
	}

	return haplotype_allele_info;
}


pair<uint, bool> calculateFlankEditDistance(Variant * variant, const vector<AltAlleleInfo> & haplotype_allele_info, const uint read_length) {

	assert(variant);

	uint flank_edit_distance = 0;

	if (haplotype_allele_info.empty()) {

		return make_pair(flank_edit_distance, true);
	}

	vector<AltAlleleInfo>::const_iterator haplotype_allele_info_it = lower_bound(haplotype_allele_info.begin(), haplotype_allele_info.end(), variant->pos(), [](const AltAlleleInfo & elem, const uint & val) { return elem.pos < val;});
	reverse_iterator<vector<AltAlleleInfo>::const_iterator> haplotype_allele_info_rit(haplotype_allele_info_it);

	if (haplotype_allele_info_rit != haplotype_allele_info.crend()) {

		assert(haplotype_allele_info_rit->pos < variant->pos());

		if (variant->pos() < (haplotype_allele_info_rit->pos + haplotype_allele_info_rit->ref_len)) {

			return make_pair(flank_edit_distance, false);
		}
	}

	uint upstream_end_position = variant->pos();
	uint right_flank_length = 0;

	while (haplotype_allele_info_it != haplotype_allele_info.cend()) {

    	assert(upstream_end_position <= haplotype_allele_info_it->pos);

		right_flank_length += haplotype_allele_info_it->pos - upstream_end_position;    	

		if (right_flank_length < read_length) {

			flank_edit_distance++;
		
		} else {

			break;
		}

		right_flank_length += haplotype_allele_info_it->alt_len;
		upstream_end_position = haplotype_allele_info_it->pos + haplotype_allele_info_it->ref_len;	

		haplotype_allele_info_it++;
	}	

	uint downstream_start_position = variant->pos();
	uint left_flank_length = 1;

	while (haplotype_allele_info_rit != haplotype_allele_info.crend()) {

    	assert((haplotype_allele_info_rit->pos + haplotype_allele_info_rit->ref_len) <= downstream_start_position);

		left_flank_length += downstream_start_position - (haplotype_allele_info_rit->pos + haplotype_allele_info_rit->ref_len);    	

		if (left_flank_length < read_length) {

			flank_edit_distance++;
		
		} else {

			break;
		}

		left_flank_length += haplotype_allele_info_rit->alt_len;
		downstream_start_position = haplotype_allele_info_rit->pos;	

		haplotype_allele_info_rit++;
	}

	return make_pair(flank_edit_distance, true);
}

void addMaxFlankEditDistance(Variant * variant, const vector<vector<AltAlleleInfo> > & haplotypes_allele_info, const string sample_id, const uint read_length) {

	assert(variant);

	pair<uint, bool> max_flank_edit_distance(0, false);

    for (auto & haplotype_allele_info: haplotypes_allele_info) {

    	auto flank_edit_distance = calculateFlankEditDistance(variant, haplotype_allele_info, read_length);

    	if (flank_edit_distance.second) {

    		max_flank_edit_distance.first = max(max_flank_edit_distance.first, flank_edit_distance.first);
    		max_flank_edit_distance.second = true;
    	}
    }

    if (max_flank_edit_distance.second) {

		variant->getSample(sample_id).info().setValue<float>("MED", max_flank_edit_distance.first);    		
    		    
    } else {

    	variant->getSample(sample_id).info().setValue<float>("MED", -1);
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

    			if (gt_allele_pair_1 == cs_allele_pair_1) {

    				genotype_concordance = "T";

    			} else {

    				genotype_concordance = "F";
    			}
    		} 

    	} else {
    		
    		genotype_concordance = "F";
    	}

    	auto gt_sample_gtco = gt_sample->info().getValue<string>("GTCO");

    	if (gt_sample_gtco.second) {

    		if (genotype_concordance == "T") {

	    		gt_sample->info().setValue<string>("GTCO", genotype_concordance);

    		} else if ((genotype_concordance == "P") and ((gt_sample_gtco.first == "F") | (gt_sample_gtco.first == "I"))) {

	    		gt_sample->info().setValue<string>("GTCO", genotype_concordance);
    		
    		} else if (gt_sample_gtco.first == "I") {

	    		gt_sample->info().setValue<string>("GTCO", genotype_concordance);
    		}
    	
    	} else {

    		gt_sample->info().setValue<string>("GTCO", genotype_concordance);
    	}

		assert(!cs_sample->info().getValue<string>("GTCO").second);
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

		// auto gt_haplotypes_allele_info = getHaplotypesAlleleInfo(*gt_variants, sample_id);

		auto gt_variants_it = gt_variants->begin();
		auto cs_variants_it = cs_variants->begin();

		uint gt_last_position = 0;
		uint cs_last_position = 0;

		while ((gt_variants_it != gt_variants->end()) or (cs_variants_it != cs_variants->end())) {

			if ((gt_variants_it != gt_variants->end()) and (cs_variants_it != cs_variants->end())) {

				if ((*gt_variants_it)->pos() == (*cs_variants_it)->pos()) {

					// addMaxFlankEditDistance(*gt_variants_it, gt_haplotypes_allele_info, sample_id, read_length);
					// addMaxFlankEditDistance(*cs_variants_it, gt_haplotypes_allele_info, sample_id, read_length);

					addGenotypeConcordance(*gt_variants_it, *cs_variants_it, sample_id);

					assert(cs_last_position <= (*cs_variants_it)->pos());
					cs_last_position = (*cs_variants_it)->pos();
					
					cs_variants_it++;

					assert(gt_last_position < (*gt_variants_it)->pos());

					if (cs_variants_it != cs_variants->end()) {

						assert(cs_last_position <= (*cs_variants_it)->pos());

						if (cs_last_position < (*cs_variants_it)->pos()) {

							gt_last_position = (*gt_variants_it)->pos();
							gt_variants_it++;
						}
					
					} else {

						gt_last_position = (*gt_variants_it)->pos();
						gt_variants_it++;
					}

				} else if ((*gt_variants_it)->pos() < (*cs_variants_it)->pos()) {

					// addMaxFlankEditDistance(*gt_variants_it, gt_haplotypes_allele_info, sample_id, read_length);
					addGenotypeConcordance(*gt_variants_it, nullptr, sample_id);

					assert(gt_last_position < (*gt_variants_it)->pos());
					gt_last_position = (*gt_variants_it)->pos();

					gt_variants_it++;

				} else {

					assert(((*cs_variants_it)->pos() < (*gt_variants_it)->pos()));

					// addMaxFlankEditDistance(*cs_variants_it, gt_haplotypes_allele_info, sample_id, read_length);
					addGenotypeConcordance(nullptr, *cs_variants_it, sample_id);

					assert(cs_last_position <= (*cs_variants_it)->pos());
					cs_last_position = (*cs_variants_it)->pos();
					
					cs_variants_it++;
				}

			} else if (gt_variants_it != gt_variants->end()) {

				// addMaxFlankEditDistance(*gt_variants_it, gt_haplotypes_allele_info, sample_id, read_length);
				addGenotypeConcordance(*gt_variants_it, nullptr, sample_id);

				assert(gt_last_position < (*gt_variants_it)->pos());
				gt_last_position = (*gt_variants_it)->pos();				

				gt_variants_it++;

			} else {


				// addMaxFlankEditDistance(*cs_variants_it, gt_haplotypes_allele_info, sample_id, read_length);
				addGenotypeConcordance(nullptr, *cs_variants_it, sample_id);

				assert(cs_last_position <= (*cs_variants_it)->pos());
				cs_last_position = (*cs_variants_it)->pos();

				cs_variants_it++;
			}
		}
	}
}

int main(int argc, char const *argv[]) {

    if (argc != 6) {

        std::cout << "USAGE: addEditDistanceAndCondordance <truth_variant_file> <callset_variant_file> <truth_output_prefix> <callset_output_prefix> <read_length>" << std::endl;
        return 1;
    }

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") addEditDistanceAndCondordance script ...\n" << endl;

    cout << "\nWARNING: edit distance calculation is currently disabled\n" << endl;

	GenotypedVcfFileReader gt_vcf_reader(argv[1], true);
	GenotypedVcfFileReader cs_vcf_reader(argv[2], true);
	
	assert(gt_vcf_reader.metaData().sampleIds() == cs_vcf_reader.metaData().sampleIds());

	cout << "[" << Utils::getLocalTime() << "] Adding maximum edit distance (MED) and genotype concordance (GTCO) to " << gt_vcf_reader.metaData().sampleIds().size() << " samples ...\n" << endl;

	auto gt_output_meta_data = gt_vcf_reader.metaData();
	auto cs_output_meta_data = cs_vcf_reader.metaData();

	const uint read_length = stoi(argv[5]);

	// gt_output_meta_data.formatDescriptors().emplace("MED", Attribute::DetailedDescriptor("MED", Attribute::Number::One, Attribute::Type::Float, "Maximum edit distance to reference within a window of size " + to_string(read_length * 2 - 1) + " centered at position"));
	gt_output_meta_data.formatDescriptors().emplace("GTCO", Attribute::DetailedDescriptor("GTCO", Attribute::Number::One, Attribute::Type::String, "Genotype concordance"));

	// cs_output_meta_data.formatDescriptors().emplace("MED", Attribute::DetailedDescriptor("MED", Attribute::Number::One, Attribute::Type::Float, "Maximum edit distance to reference within a window of size " + to_string(read_length * 2 - 1) + " centered at position"));
	cs_output_meta_data.formatDescriptors().emplace("GTCO", Attribute::DetailedDescriptor("GTCO", Attribute::Number::One, Attribute::Type::String, "Genotype concordance"));

	VcfFileWriter gt_vcf_writer(string(argv[3]) + ".vcf.gz", gt_output_meta_data, true);
	VcfFileWriter cs_vcf_writer(string(argv[4]) + ".vcf.gz", cs_output_meta_data, true);

	const vector<Contig> contigs_merged = Auxiliaries::mergeContigs(gt_vcf_reader.metaData().contigs(), cs_vcf_reader.metaData().contigs());
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

	for (auto & contig: contigs_merged) {

		cluster_end_position = 0;

		while (gt_cur_var) {

			assert(!gt_cur_var->ref().seq().empty());

			if (gt_cur_var->chrom() != contig.id()) {

				break;
			}

			while (cs_cur_var) {

				assert(!cs_cur_var->ref().seq().empty());

				if (cs_cur_var->chrom() != contig.id()) { 

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

			assert(!cs_cur_var->ref().seq().empty());

			if (cs_cur_var->chrom() != contig.id()) {

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

		cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
	}

	assert(!cs_cur_var);
	assert(!gt_cur_var);

	assert(!cs_vcf_reader.getNextVariant(&cs_cur_var));
	assert(!gt_vcf_reader.getNextVariant(&gt_cur_var));

	cout << "\n[" << Utils::getLocalTime() << "] Added MED and GTCO to genotypes on " << num_gt_variants << " ground truth variants." << endl;
	cout << "[" << Utils::getLocalTime() << "] Added MED and GTCO to genotypes on " << num_cs_variants << " callset variants." << endl;
	cout << endl;

	return 0;
}




