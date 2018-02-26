
/*
convertSeqToAlleleId.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

    output_meta_data.infoDescriptors().clear();
	output_meta_data.infoDescriptors().emplace("END", Attribute::DetailedDescriptor("END", Attribute::Number::One, Attribute::Type::Int, "End position of the variant described in this record"));
	output_meta_data.infoDescriptors().emplace("SVTYPE", Attribute::DetailedDescriptor("SVTYPE", Attribute::Number::One, Attribute::Type::String, "Type of structural variant"));
	output_meta_data.infoDescriptors().emplace("CIPOS", Attribute::DetailedDescriptor("CIPOS", Attribute::Number::Two, Attribute::Type::Int, "Confidence interval around POS for imprecise variants"));
	output_meta_data.infoDescriptors().emplace("CIEND", Attribute::DetailedDescriptor("CIEND", Attribute::Number::Two, Attribute::Type::Int, "Confidence interval around END for imprecise variants"));

	VcfFileWriter vcf_writer(string(argv[2]) + ".vcf", output_meta_data, true);

    cout << "[" << Utils::getLocalTime() << "] Parsing reference genome fasta ..." << endl;

    unordered_map<string, FastaRecord*> genome_seqs;
    FastaReader genome_reader(argv[3]);
    FastaRecord * cur_fasta_rec;

	while (genome_reader.getNextRecord(&cur_fasta_rec)) {

        string cur_fasta_id = Utils::splitString(cur_fasta_rec->id(), '\t').front();

        assert(genome_seqs.emplace(cur_fasta_id, cur_fasta_rec).second);
        genome_seqs.at(cur_fasta_id)->convertToUppercase();
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

        Auxiliaries::rightTrimVariant(cur_var);

        multimap<uint, Variant> converted_variants;

        for (uint alt_allele_idx = 0; alt_allele_idx < cur_var->numAlts(); alt_allele_idx++) {

            bool alt_allele_excluded = true;

            uint cur_pos = cur_var->pos();

    		assert(!(cur_var->alt(alt_allele_idx).isID()));

            Allele cur_ref_allele = cur_var->ref();
            Allele cur_alt_allele = cur_var->alt(alt_allele_idx);

            Auxiliaries::rightTrimAllelePair(&cur_ref_allele, &cur_alt_allele);

            assert(!(cur_ref_allele.seq().empty()));
            assert(!(cur_alt_allele.seq().empty()));

            if ((cur_ref_allele.seq().size() == cur_alt_allele.seq().size()) and (cur_ref_allele.seq().size() == 1)) {

                continue;
            }

            AttributeSet cur_info;

            if ((cur_ref_allele.seq().size() > 1) and (cur_alt_allele.seq().size() == 1) and (cur_ref_allele.seq().front() == cur_alt_allele.seq().front())) {

    	        if (min_sv_length <= (cur_ref_allele.seq().size() - 1)) {

    	        	alt_allele_excluded = false;
                    num_converted_deletions++;

    	        	cur_info.setValue<int>("END", cur_pos + cur_ref_allele.seq().size() - 1);
    	        	cur_info.setValue<string>("SVTYPE", "DEL");

                    cur_ref_allele.seq() = cur_ref_allele.seq().front();
                    cur_alt_allele = Allele("<DEL>");
    	        } 
            
            } else if ((cur_ref_allele.seq().size() == 1) and (cur_alt_allele.seq().size() > 1) and (cur_ref_allele.seq().front() == cur_alt_allele.seq().front())) {

                if ((min_sv_length <= (cur_alt_allele.seq().size() - 1)) and (genome_seqs_it->second->seq().substr(cur_pos - 1, cur_alt_allele.seq().size()) == cur_alt_allele.seq())) {

            		alt_allele_excluded = false;
                    num_converted_duplications++;

                	cur_info.setValue<int>("END", cur_pos + cur_alt_allele.seq().size() - 1);
                	cur_info.setValue<string>("SVTYPE", "DUP");

                    cur_ref_allele.seq() = cur_ref_allele.seq().front();
                    cur_alt_allele = Allele("<DUP>");
                }

            } else if (cur_ref_allele.seq().size() == cur_alt_allele.seq().size()) {

                if ((cur_ref_allele.seq().front() == cur_alt_allele.seq().front()) and (min_sv_length <= (cur_ref_allele.seq().size() - 1)) and (cur_ref_allele.seq().substr(1, string::npos) == Auxiliaries::reverseComplementSequence(cur_alt_allele.seq().substr(1, string::npos)))) {

                    alt_allele_excluded = false;
                    num_converted_inversions++;

                    cur_alt_allele = Allele("<INV>");
                    cur_info.setValue<int>("END", cur_pos + cur_ref_allele.seq().size() - 1);
                    cur_info.setValue<string>("SVTYPE", "INV");

                    cur_ref_allele.seq() = cur_ref_allele.seq().front();
                
                } else if ((cur_ref_allele.seq().front() != cur_alt_allele.seq().front()) and (min_sv_length <= cur_ref_allele.seq().size()) and (cur_ref_allele.seq() == Auxiliaries::reverseComplementSequence(cur_alt_allele.seq()))) {

                    alt_allele_excluded = false;
                    num_converted_inversions++;

                    cur_pos--;                

                    cur_alt_allele = Allele("<INV>");
                    cur_info.setValue<int>("END", cur_pos + cur_ref_allele.seq().size());
                    cur_info.setValue<string>("SVTYPE", "INV");

                    cur_ref_allele.seq() = genome_seqs_it->second->seq().at(cur_pos - 1);                
                }
    		}         

            if (!alt_allele_excluded) {

                assert(cur_alt_allele.isID());

                cur_info.setValue<string>("CIPOS", "0,0");
                cur_info.setValue<string>("CIEND", "0,0");

                converted_variants.emplace(cur_pos, Variant(cur_chromosome, cur_pos, cur_ref_allele, vector<Allele>(1, cur_alt_allele), cur_info));
            }
        }

        for (auto & converted_variant: converted_variants) {

            vcf_writer.write(&(converted_variant.second));
        }

		delete cur_var;

		if ((num_variants % 1000000) == 0) {

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
