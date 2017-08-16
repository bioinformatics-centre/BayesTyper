
/*
convertSeqToAlleleId.cpp - This file is part of BayesTyper (v1.1)


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
#include "FastaRecord.hpp"
#include "FastaReader.hpp"
#include "Stats.hpp"
#include "Auxiliaries.hpp"


int main(int argc, char const *argv[]) {

	if (argc != 5) {

		std::cout << "USAGE: convertSeqToAlleleId <input> <output_prefix> <genome> <min_sv_length>" << std::endl;
		return 1;
	}

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") convertSeqToAlleleId script ...\n" << endl;

	VcfFileReader vcf_reader(string(argv[1]), true);
	auto output_meta_data = vcf_reader.metaData();

	output_meta_data.infoDescriptors().emplace("END", Attribute::DetailedDescriptor("END", Attribute::Number::A, Attribute::Type::Int, "End position of the variant described in this record"));
	output_meta_data.infoDescriptors().emplace("SVTYPE", Attribute::DetailedDescriptor("SVTYPE", Attribute::Number::A, Attribute::Type::String, "Type of structural variant"));

	VcfFileWriter vcf_writer(string(argv[2]) + ".vcf", output_meta_data, true);

    cout << "[" << Utils::getLocalTime() << "] Parsing reference genome fasta ..." << endl;

    unordered_map<string, FastaRecord*> genome_seqs;
    FastaReader genome_reader(argv[3]);
    FastaRecord * cur_fasta_rec;

	while (genome_reader.getNextRecord(&cur_fasta_rec)) {

        assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
        genome_seqs.at(cur_fasta_rec->id())->convertToUppercase();
    }

    cout << "[" << Utils::getLocalTime() << "] Parsed " << genome_seqs.size() << " chromosome(s)\n" << endl;

	const uint min_sv_length = stoi(argv[4]);

    Variant * cur_var;
    string cur_chromosome = "";

    auto genome_seqs_it = genome_seqs.find(cur_chromosome);
    assert(genome_seqs_it == genome_seqs.end());

	uint num_variants = 0;

	uint num_converted_deletions = 0;
	uint num_converted_duplications = 0;
	uint num_converted_inversions = 0;

    while (vcf_reader.getNextVariant(&cur_var)) {

        num_variants++;

        if (cur_var->chrom() != cur_chromosome) {

            genome_seqs_it = genome_seqs.find(cur_var->chrom());
            assert(genome_seqs_it != genome_seqs.end());

            cur_chromosome = cur_var->chrom();
        }

        vector<uint> excluded_alt_indices;
        excluded_alt_indices.reserve(cur_var->numAlts());

		for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

			assert(!(cur_var->alt(alt_allele_idx).isID()));

			auto allele_attributes = Auxiliaries::alleleAttributes(cur_var->alt(alt_allele_idx), cur_var->ref());

			if (allele_attributes.type == Auxiliaries::Type::Deletion) {

				Allele ref_allele = cur_var->ref();
	        	Allele alt_allele = cur_var->alt(alt_allele_idx);

	            Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);

	            assert(ref_allele.seq().size() > 1);
	            assert(alt_allele.seq().size() == 1);

	            if (min_sv_length <= (ref_allele.seq().size() - 1)) {

	            	num_converted_deletions++;

	            	cur_var->alt(alt_allele_idx).seq() = "<DEL>";
	            	cur_var->alt(alt_allele_idx).info().setValue<int>("END", cur_var->pos() + ref_allele.seq().size() - 1);
	            	cur_var->alt(alt_allele_idx).info().setValue<string>("SVTYPE", "DEL");

	            } else {

	            	excluded_alt_indices.push_back(alt_allele_idx);
	            }

			} else if (allele_attributes.type == Auxiliaries::Type::Insertion) {

				Allele ref_allele = cur_var->ref();
	        	Allele alt_allele = cur_var->alt(alt_allele_idx);

	            Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);

	            assert(ref_allele.seq().size() == 1);
	            assert(alt_allele.seq().size() > 1);

	            if ((min_sv_length <= (alt_allele.seq().size() - 1)) and (genome_seqs_it->second->seq().substr(cur_var->pos() - 1, alt_allele.seq().size()) == alt_allele.seq())) {

	            	num_converted_duplications++;

	            	cur_var->alt(alt_allele_idx).seq() = "<DUP>";
	            	cur_var->alt(alt_allele_idx).info().setValue<int>("END", cur_var->pos() + alt_allele.seq().size() - 1);
	            	cur_var->alt(alt_allele_idx).info().setValue<string>("SVTYPE", "DUP");

	            } else {

	            	excluded_alt_indices.push_back(alt_allele_idx);
	            }

			} else if (allele_attributes.type == Auxiliaries::Type::Inversion) {

				Allele ref_allele = cur_var->ref();
	        	Allele alt_allele = cur_var->alt(alt_allele_idx);

	            Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);

	            assert(ref_allele.seq().size() == alt_allele.seq().size());

	            if ((min_sv_length <= (ref_allele.seq().size() - 1)) and (ref_allele.seq().substr(1, string::npos) == Auxiliaries::reverseComplementSequence(alt_allele.seq().substr(1, string::npos)))) {

	            	num_converted_inversions++;

	            	cur_var->alt(alt_allele_idx).seq() = "<INV>";
	            	cur_var->alt(alt_allele_idx).info().setValue<int>("END", cur_var->pos() + ref_allele.seq().size() - 1);
	            	cur_var->alt(alt_allele_idx).info().setValue<string>("SVTYPE", "INV");

	            } else {

	            	excluded_alt_indices.push_back(alt_allele_idx);
	            }

			} else {

            	excluded_alt_indices.push_back(alt_allele_idx);				
			}
		}         

        assert(cur_var->numAlts() > 0);
        assert(cur_var->numAlts() >= excluded_alt_indices.size());

        if (!(excluded_alt_indices.empty())) {

            cur_var->removeAlts(excluded_alt_indices);
        }

        if (cur_var->numAlts() > 0) {

        	cur_var->ref().seq() = cur_var->ref().seq().front();
            vcf_writer.write(cur_var);
        }

		delete cur_var;

		if ((num_variants % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
		}
	}

    cout << "\n[" << Utils::getLocalTime() << "] Converted " << num_converted_deletions + num_converted_duplications + num_converted_inversions << " alternative allele(s):\n" << endl;

    cout << "\t- <DEL>: " << num_converted_deletions << endl;
    cout << "\t- <DUP>: " << num_converted_duplications << endl;
    cout << "\t- <INV>: " << num_converted_inversions << endl;

	cout << endl;

	return 0;
}
