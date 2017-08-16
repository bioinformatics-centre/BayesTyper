
/*
transferAlleleFreqFields.cpp - This file is part of BayesTyper (v1.1)


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

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"

int main(int argc, char const *argv[]) {

	if (argc != 4) {

		std::cout << "USAGE: transferAlleleFreqFields <target_vcf> <donor_vcf> <out_prefix>" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") transferAlleleFreqFields script ...\n" << endl;

	GenotypedVcfFileReader target_vcf_reader(argv[1], true);
	VcfFileReader donor_vcf_reader(argv[2], true);
	auto target_contigs_iter = target_vcf_reader.metaData().contigs().begin();
	auto donor_contigs_iter = donor_vcf_reader.metaData().contigs().begin();

	while (target_contigs_iter != target_vcf_reader.metaData().contigs().end() and donor_contigs_iter != donor_vcf_reader.metaData().contigs().end()) {

		assert(*target_contigs_iter == *donor_contigs_iter);
		++target_contigs_iter;
		++donor_contigs_iter;
	}

	// Auxiliaries::removeNonRelevantInfoDescriptors(&(donor_vcf_reader.metaData()), info_fields);

	auto donor_contigs = donor_vcf_reader.metaData().contigs();

	auto output_meta_data = target_vcf_reader.metaData();
	vector<pair<string, string> > gdk_af_descriptor_elems({make_pair("ID","GDK_AF"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele frequency in gdk callset")});
	vector<pair<string, string> > gdk_ac_descriptor_elems({make_pair("ID","GDK_AC"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele count in gdk callset")});
	vector<pair<string, string> > gdk_an_descriptor_elems({make_pair("ID","GDK_AN"), make_pair("Number","A"), make_pair("Type","String"), make_pair("Description","Allele number in gdk callset")});

	output_meta_data.infoDescriptors().emplace("GDK_AF", Attribute::DetailedDescriptor(gdk_af_descriptor_elems));
	output_meta_data.infoDescriptors().emplace("GDK_AC", Attribute::DetailedDescriptor(gdk_ac_descriptor_elems));
	output_meta_data.infoDescriptors().emplace("GDK_AN", Attribute::DetailedDescriptor(gdk_an_descriptor_elems));

	VcfFileWriter output_vcf(string(argv[3]) + ".vcf", output_meta_data, true);

	Variant * cur_donor_var;
	assert(donor_vcf_reader.getNextVariant(&cur_donor_var));

	Variant * cur_target_var;
	assert(target_vcf_reader.getNextVariant(&cur_target_var));

	uint num_donor_variants = 0;
	uint num_target_variants = 0;
	uint num_target_variants_w_donor = 0;
	uint num_target_alleles = 0;
	uint num_matched_target_alleles = 0;

	for (auto & donor_contig : donor_contigs) {

		// Build contig var index for donor
		map<uint,Variant *> donor_var_map;

		do {

			if (!cur_donor_var) {

				break;
			}

			if (cur_donor_var->chrom() != donor_contig.id()) {

				break;
			}

			num_donor_variants++;
			donor_var_map.emplace(cur_donor_var->pos(), cur_donor_var); // Assumes only one variant per position

		} while (donor_vcf_reader.getNextVariant(&cur_donor_var));

		do {

			if (!cur_target_var) {

				break;
			}

			if (cur_target_var->chrom() != donor_contig.id()) {

				break;
			}

			num_target_variants++;
			auto donor_var_fi = donor_var_map.find(cur_target_var->pos());

			if (donor_var_fi != donor_var_map.end()) {

				auto donor_var_match = donor_var_fi->second;
				num_target_variants_w_donor++;
				auto donor_an = donor_var_match->info().getValue<int>("AN");
				assert(donor_an.second);

				for (uint target_alt_idx = 0; target_alt_idx < cur_target_var->numAlts(); ++target_alt_idx) {

					++num_target_alleles;

					Allele cur_target_ref = cur_target_var->ref();
					Allele cur_target_alt = cur_target_var->alt(target_alt_idx);

					Auxiliaries::rightTrimAllelePair(&cur_target_ref, &cur_target_alt);

					bool found_match = false;

					for (uint donor_alt_idx = 0; donor_alt_idx < donor_var_match->numAlts(); ++donor_alt_idx) {

						Allele cur_donor_ref = donor_var_match->ref();
						Allele cur_donor_alt = donor_var_match->alt(donor_alt_idx);

						Auxiliaries::rightTrimAllelePair(&cur_donor_ref, &cur_donor_alt);

						if (cur_target_ref.seq() == cur_donor_ref.seq() and cur_target_alt.seq() == cur_donor_alt.seq()) {

							found_match = true;
							++num_matched_target_alleles;

							auto donor_af = donor_var_match->alt(donor_alt_idx).info().getValue<float>("AF");
							auto donor_ac = donor_var_match->alt(donor_alt_idx).info().getValue<int>("AC");
							assert(donor_af.second);
							assert(donor_ac.second);

							assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AF", to_string(donor_af.first)));
							assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AC", to_string(donor_ac.first)));
							assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AN", to_string(donor_an.first)));
							break;
						}
					}

					if (!found_match) {

						// cout << "no match " << endl;
						// cout << cur_target_var->vcf(target_vcf_reader.metaData()) << endl;
						// cout << donor_var_match->vcf(donor_vcf_reader.metaData()) << endl;
						assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AF", string("NA")));
						assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AC", string("NA")));
						assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AN", string("NA")));
					}
				}

			} else {

				// assert(cur_target_var->info().setValue("GDK_AN", string("NA")));

				for (uint target_alt_idx = 0; target_alt_idx < cur_target_var->numAlts(); ++target_alt_idx) {

					++num_target_alleles;
					assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AF", string("NA")));
					assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AC", string("NA")));
					assert(cur_target_var->alt(target_alt_idx).info().setValue("GDK_AN", string("NA")));
				}
			}

			output_vcf.write(cur_target_var);
			delete cur_target_var;

		} while (target_vcf_reader.getNextVariant(&cur_target_var));

		for (auto & donor_var : donor_var_map) {

			delete donor_var.second;
		}
	}

	cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_donor_variants << " donor variant(s) and " << num_target_variants << " target variant(s)" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] \t" << num_target_variants_w_donor << " target variant(s) matched a donor variant" <<  endl;
	cout << "\n[" << Utils::getLocalTime() << "] \t" << num_matched_target_alleles << " target alleles(s) out of " << num_target_alleles << " matched a donor allele" <<  endl;

	return 0;
}
