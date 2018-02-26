
/*
transferAlleleCallsetOrigin.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <vector>
#include <numeric>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Stats.hpp"
#include "Auxiliaries.hpp"


int main(int argc, char const *argv[]) {

	if (argc != 4) {

		std::cout << "USAGE: transferAlleleCallsetOrigin <input> <combine_output> <output_prefix>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") transferAlleleCallsetOrigin script ...\n" << endl;

	GenotypedVcfFileReader vcf_reader(argv[1], true);
	auto output_meta_data = vcf_reader.metaData();

	VcfFileReader vcf_reader_aco(argv[2], true);

	assert(vcf_reader.metaData().contigs() == vcf_reader_aco.metaData().contigs());

	output_meta_data.infoDescriptors().emplace("ACO", vcf_reader_aco.metaData().infoDescriptors().at("ACO"));

	VcfFileWriter vcf_writer(string(argv[3]) + ".vcf", output_meta_data, true);
	
	Variant * cur_var_aco;
	vcf_reader_aco.getNextVariant(&cur_var_aco);

	Variant * cur_var;
	vcf_reader.getNextVariant(&cur_var);

	uint num_variants = 0;

	for (auto & contig: vcf_reader.metaData().contigs()) {

		unordered_map<string, vector<pair<string, string> > > aco_index;

		while (cur_var_aco and cur_var) {

			if (cur_var_aco->chrom() != contig.id()) {

				break;
			}

			if (cur_var_aco->chrom() == cur_var->chrom()) {

				for (uint alt_allele_idx = 0; alt_allele_idx < cur_var_aco->numAlts(); alt_allele_idx++) {

					Allele ref_allele = cur_var_aco->ref();
					Allele alt_allele = cur_var_aco->alt(alt_allele_idx);

				    auto pos_shift = Auxiliaries::fullTrimAllelePair(&ref_allele, &alt_allele);

				    assert(!(ref_allele.seq().empty()) or !(alt_allele.seq().empty()));

				    auto aco_value = alt_allele.info().getValue<string>("ACO");
				    assert(aco_value.second);

				    string aco_index_key = cur_var_aco->chrom() + "_" + to_string(cur_var_aco->pos() + pos_shift.first) + "_" + to_string(ref_allele.seq().size());

				    auto aco_index_it = aco_index.emplace(aco_index_key, vector<pair<string, string> >());
				    aco_index_it.first->second.emplace_back(alt_allele.seq(), aco_value.first);
				}
			}

			delete cur_var_aco;
			vcf_reader_aco.getNextVariant(&cur_var_aco);
		}

		while (cur_var) {

			if (cur_var->chrom() != contig.id()) {

				break;
			}

			num_variants++;

			for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

				Allele ref_allele = cur_var->ref();
				Allele alt_allele = cur_var->alt(alt_allele_idx);

			    auto pos_shift = Auxiliaries::fullTrimAllelePair(&ref_allele, &alt_allele);

			    assert(!(ref_allele.seq().empty()) or !(alt_allele.seq().empty()));

				cur_var->alt(alt_allele_idx).info().setValue<string>("ACO", ".");

			    string aco_index_key = cur_var->chrom() + "_" + to_string(cur_var->pos() + pos_shift.first) + "_" + to_string(ref_allele.seq().size());

			    auto aco_index_it = aco_index.find(aco_index_key);

			    if (aco_index_it != aco_index.end()) {

			    	int match_idx = -1;

			    	for (uint aco_seq_idx = 0; aco_seq_idx < aco_index_it->second.size(); aco_seq_idx++) {

			    		if (alt_allele.seq() == aco_index_it->second.at(aco_seq_idx).first) {

			    			assert(match_idx == -1);
			    			match_idx = aco_seq_idx;
			    		}
			    	}

			    	assert(match_idx < static_cast<int>(aco_index_it->second.size()));

			    	if (match_idx != -1) {

				    	cur_var->alt(alt_allele_idx).info().setValue<string>("ACO", aco_index_it->second.at(match_idx).second);
			    	}
			    } 
			}

			vcf_writer.write(cur_var);

			delete cur_var;
			vcf_reader.getNextVariant(&cur_var);
		}

		cout << "[" << Utils::getLocalTime() << "] Finished chromosome " << contig.id() << endl;
	}

	if (cur_var_aco) {

		delete cur_var_aco;
	}

	if (cur_var) {

		delete cur_var;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Added allele call-set origin (ACO) to " << num_variants << "." << endl;
	cout << endl;

	return 0;
}
