
/*
annotate.cpp - This file is part of BayesTyper (v1.1)


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
#include <unordered_map>
#include <set>
#include <algorithm>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"

#include "edlib/edlib.h"

namespace Annotate {

	struct AnnotatedAllele {

		uint trimmed_pos;
		string ref_seq;
		string alt_seq;
		vector<string> var_ids;
	};


	float calcMatchScore(uint ref1_len, uint alt1_len, uint ref2_len, uint alt2_len, uint ref_edit, uint alt_edit) {

		float score = 1 - (ref_edit + alt_edit)/static_cast<float>(max(ref1_len, ref2_len) + max(alt1_len, alt2_len));
		assert(score >= 0);
		assert(score <= 1);
		return score;
	}

	int getNCount(const string & seq) {

		return static_cast<int>(count(seq.begin(), seq.end(), 'N'));
	}

	int edlibAlignSafe(const string & string1, const string & string2) {

		if (string1.empty() and string2.empty()) {

			return 0;

		} else if (string1.empty()) {

			return string2.size();

		} else if (string2.empty()) {

			return string1.size();

		} else {

			int n_diff = abs(getNCount(string1) - getNCount(string2));
			auto ed_dist_temp = edlibAlign(string1.c_str(), string1.size(), string2.c_str(), string2.size(), edlibDefaultAlignConfig()).editDistance;
			assert(ed_dist_temp >= n_diff);

			return ed_dist_temp - n_diff;
		}
	}

	void annotate(const string & vcf_filename, const string & annotation_filename, const string & output_prefix, const float match_threshold, const float window_size_scale, const bool overwrite_prev_anno) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") annotate ...\n" << endl;

		bool self_anno_mode = (vcf_filename == annotation_filename);

		assert(window_size_scale >= 1);

		cout << "\n[" << Utils::getLocalTime() << "] Match threshold: " << match_threshold << endl;
		cout << "\n[" << Utils::getLocalTime() << "] Window size scale: " << window_size_scale << endl;
		cout << "\n[" << Utils::getLocalTime() << "] Overwrite previous annotation: " << overwrite_prev_anno << "\n" << endl;
		cout << "\n[" << Utils::getLocalTime() << "] Self annotate mode: " << self_anno_mode << "\n" << endl;

		GenotypedVcfFileReader callset_vcf_reader(vcf_filename, true);
		VcfFileReader annotation_vcf_reader(annotation_filename, true);

		auto annotation_contigs = annotation_vcf_reader.metaData().contigs();
		assert(annotation_contigs == callset_vcf_reader.metaData().contigs());

		auto output_meta_data = callset_vcf_reader.metaData();

		if (!output_meta_data.infoDescriptors().count("AAI")) {

			vector<pair<string, string> > aai_descriptor_elems({make_pair("ID","AAI"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele annotation")});
			output_meta_data.infoDescriptors().emplace("AAI", Attribute::DetailedDescriptor(aai_descriptor_elems));
		}

		VcfFileWriter output_vcf(output_prefix + ".vcf", output_meta_data, true);

		Variant * cur_annotation_var;
		assert(annotation_vcf_reader.getNextVariant(&cur_annotation_var));

		Variant * cur_callset_var;
		assert(callset_vcf_reader.getNextVariant(&cur_callset_var));

		uint num_annotation_variants = 0;
		uint num_annotation_alleles = 0;
		uint num_callset_variants = 0;
		uint num_callset_alleles = 0;
		uint num_annotated_callset_alleles = 0;

		for (auto & annotation_contig : annotation_contigs) {

			// Build contig var index for donor
			map<uint,vector<AnnotatedAllele> > annotation_allele_map;
			cout << "\n[" << Utils::getLocalTime() << "] Indexing annotated variants on " << annotation_contig.id() << "\n" << endl;

			do {

				if (!cur_annotation_var) {

					break;
				}

				if (cur_annotation_var->chrom() != annotation_contig.id()) {

					break;
				}

				++num_annotation_variants;

				for (uint alt_idx = 0; alt_idx < cur_annotation_var->numAlts(); ++alt_idx) {

					if (!cur_annotation_var->alt(alt_idx).isMissing()) {

						assert(!(cur_annotation_var->alt(alt_idx).isID()));

						++num_annotation_alleles;

						auto ref_trimmed = cur_annotation_var->ref();
						auto alt_trimmed = cur_annotation_var->alt(alt_idx);

						uint trimmed_pos = Auxiliaries::fullTrimAllelePair(&ref_trimmed, &alt_trimmed).first + cur_annotation_var->pos();
						assert(trimmed_pos <= annotation_contig.length());

						assert(!cur_annotation_var->ids().empty());

						AnnotatedAllele ann_all = AnnotatedAllele{trimmed_pos, ref_trimmed.seq(), alt_trimmed.seq(), cur_annotation_var->ids()};
						auto empl_res_start = annotation_allele_map.emplace(trimmed_pos, vector<AnnotatedAllele> { ann_all });

						if (!empl_res_start.second) {

							empl_res_start.first->second.emplace_back(ann_all);
						}

						if (!ref_trimmed.seq().empty()) {

							auto empl_res_end = annotation_allele_map.emplace(trimmed_pos + ref_trimmed.seq().size() - 1, vector<AnnotatedAllele> { ann_all });

							if (!empl_res_end.second) {

								empl_res_end.first->second.emplace_back(ann_all);
							}
						}
					}
				}

				delete cur_annotation_var;

			} while (annotation_vcf_reader.getNextVariant(&cur_annotation_var));

			cout << "\n[" << Utils::getLocalTime() << "] Searching for annotated variants in callset ...\n" << endl;

			do {

				if (!cur_callset_var) {

					break;
				}

				if (cur_callset_var->chrom() != annotation_contig.id()) {

					break;
				}

				num_callset_variants++;

				if (overwrite_prev_anno) {

					cur_callset_var->setIds(vector<string>{});
				}

				for (uint alt_idx = 0; alt_idx < cur_callset_var->numAlts(); ++alt_idx) {

					if (overwrite_prev_anno) {

						cur_callset_var->alt(alt_idx).info().rm("AAI");
					}

					++num_callset_alleles;

					set<string> all_ann_ids;

					if (!cur_callset_var->alt(alt_idx).isMissing()) {

						assert(!(cur_callset_var->alt(alt_idx).isID()));

						auto ref_trimmed = cur_callset_var->ref();
						auto alt_trimmed = cur_callset_var->alt(alt_idx);

						uint trimmed_pos = Auxiliaries::fullTrimAllelePair(&ref_trimmed, &alt_trimmed).first + cur_callset_var->pos();
						assert(trimmed_pos <= annotation_contig.length());

						string ref_trimmed_seq = ref_trimmed.seq();
						string alt_trimmed_seq = alt_trimmed.seq();

						uint window_size = ceil(window_size_scale * static_cast<float>(max(ref_trimmed_seq.size(), alt_trimmed_seq.size())));

						auto anno_all_iter = annotation_allele_map.lower_bound(static_cast<int>(trimmed_pos) - static_cast<int>(window_size));

						while (anno_all_iter != annotation_allele_map.end()) {

							if (anno_all_iter->first >= (trimmed_pos + ref_trimmed_seq.size() + window_size)) {

								break;
							}

							for (auto & anno_all : anno_all_iter->second) {

								if (self_anno_mode and trimmed_pos == anno_all.trimmed_pos and ref_trimmed_seq == anno_all.ref_seq and alt_trimmed_seq == anno_all.alt_seq) {

									continue;
								}

								if (ref_trimmed_seq.size() == 1 and alt_trimmed_seq.size() == 1) { //SNV case

									if (trimmed_pos == anno_all_iter->first and ref_trimmed_seq == anno_all.ref_seq and alt_trimmed_seq == anno_all.alt_seq) {

										for (auto & var_id : anno_all.var_ids) {

											all_ann_ids.insert(var_id);
										}
									}

								} else { // Other variants

									// Upper bounds on score
									int ref_edit_dis = abs(static_cast<int>(ref_trimmed_seq.size()) - static_cast<int>(anno_all.ref_seq.size()));
									int alt_edit_dis = abs(static_cast<int>(alt_trimmed_seq.size()) - static_cast<int>(anno_all.alt_seq.size()));

									if (calcMatchScore(ref_trimmed_seq.size(), alt_trimmed_seq.size(), anno_all.ref_seq.size(), anno_all.alt_seq.size(), ref_edit_dis, alt_edit_dis) < match_threshold) {

										continue;
									}

									ref_edit_dis = edlibAlignSafe(ref_trimmed_seq, anno_all.ref_seq);
									assert(ref_edit_dis >= 0);

									if (calcMatchScore(ref_trimmed_seq.size(), alt_trimmed_seq.size(), anno_all.ref_seq.size(), anno_all.alt_seq.size(), ref_edit_dis, alt_edit_dis) < match_threshold) {

										continue;
									}

									alt_edit_dis = edlibAlignSafe(alt_trimmed_seq, anno_all.alt_seq);
									assert(alt_edit_dis >= 0);

									if (calcMatchScore(ref_trimmed_seq.size(), alt_trimmed_seq.size(), anno_all.ref_seq.size(), anno_all.alt_seq.size(), ref_edit_dis, alt_edit_dis) >= match_threshold) {

										for (auto & var_id : anno_all.var_ids) {

											all_ann_ids.insert(var_id);
										}
									}
								}
							}

							++anno_all_iter;
						}
					}

					if (!all_ann_ids.empty()) {

						++num_annotated_callset_alleles;
					}

					for (auto & all_ann_id : all_ann_ids) {

						cur_callset_var->addId(all_ann_id);
					}

					auto cur_aai = cur_callset_var->alt(alt_idx).info().getValue<string>("AAI");

					if (cur_aai.second and cur_aai.first != ".") {

						auto cur_aai_split = Utils::splitString(cur_aai.first, ':');
						for (auto & ann : cur_aai_split) {

							assert(ann != ".");
							all_ann_ids.insert(ann);
						}
					}

					JoiningString aai_js(':');
					for (auto & all_ann_id : all_ann_ids) {

						aai_js.join(all_ann_id);
					}

					string aai_string = aai_js.str();

					if (aai_string.empty()) {

						aai_string = ".";
					}

					cur_callset_var->alt(alt_idx).info().setValue<string>("AAI", aai_string);
				}

				output_vcf.write(cur_callset_var);
				delete cur_callset_var;

				if (num_callset_variants % 100000 == 0) {

					cout << "\n[" << Utils::getLocalTime() << "] Completed annotating " << num_callset_variants << " callset variants ...\n" << endl;
				}

			} while (callset_vcf_reader.getNextVariant(&cur_callset_var));
		}

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_annotation_variants << " annotation variant(s) and " << num_callset_variants << " callset variant(s)" <<  endl;
		cout << "[" << Utils::getLocalTime() << "] Parsed " << num_annotation_alleles << " annotation allele(s) and " << num_callset_alleles << " callset allele(s)" <<  endl;
		cout << "\n[" << Utils::getLocalTime() << "] " << num_annotated_callset_alleles << " callset alleles(s) were annotated " <<  endl;
		cout << endl;
	}
}
