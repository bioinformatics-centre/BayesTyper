
/*
generateDiplotypes.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <random>

#include "VcfFile.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"
#include "Auxiliaries.hpp"
#include "FastaRecord.hpp"
#include "FastaReader.hpp"

const vector<char> canonical_bases = {'A','C','G','T'};
const uint seed = 12345678;


struct Diplotype {

    pair<uint, uint> max_allele_end_pos;
    pair<FastaRecord *, FastaRecord *> haplotype_seqs;

    ofstream * genome_outfile;
};


int main(int argc, char const *argv[]) {

	if ((argc != 4) and (argc != 5)) {

		std::cout << "USAGE: generateDiplotypes <variant_file> <genome_file> <output_prefix> (<sample_ambiguous>)" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") generateDiplotypes script ...\n" << endl;

    mt19937 prng(seed);

    unordered_map<string,FastaRecord*> ref_genome_seqs;
    FastaReader genome_reader(argv[2]);
    FastaRecord * cur_fasta_rec;

    cout << "[" << Utils::getLocalTime() << "] Parsing reference genome ..." << endl;

    while (genome_reader.getNextRecord(&cur_fasta_rec)) {

    	cur_fasta_rec->convertToUppercase();

    	if (argc == 5) {

	    	cur_fasta_rec->sampleAmbiguousBases(canonical_bases, &prng);
    	}

        assert(ref_genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
    }

    cout << "[" << Utils::getLocalTime() << "] Parsed " << ref_genome_seqs.size() << " chromosome(s)" << endl;

	GenotypedVcfFileReader vcf_reader(string(argv[1]), true);

	vcf_reader.metaData().infoDescriptors().clear();
	vcf_reader.metaData().filterDescriptors().clear();
	vcf_reader.metaData().miscMeta().clear();

	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader.metaData()), {"GT"});

    unordered_map<string, Diplotype> sample_diplotypes;

    for (auto & sample_id: vcf_reader.metaData().sampleIds()) {

        auto sample_diplotypes_it = sample_diplotypes.emplace(sample_id, Diplotype());
        assert(sample_diplotypes_it.second);

        sample_diplotypes_it.first->second.max_allele_end_pos = make_pair(0,0);
        sample_diplotypes_it.first->second.genome_outfile = new ofstream(string(argv[3]) + "_" + sample_id + ".fa");

	    if (!sample_diplotypes_it.first->second.genome_outfile->is_open()) {

	        cerr << "\nERROR: Unable to write file " << string(argv[3]) + "_" + sample_id + ".fa" << "\n" << endl;
	        exit(1);
	    }
    }

    assert(vcf_reader.metaData().sampleIds().size() == sample_diplotypes.size());

	cout << "\n[" << Utils::getLocalTime() << "] Generating diplotypes for " << sample_diplotypes.size() << " sample(s) ...\n" << endl;

    auto contigs = vcf_reader.metaData().contigs();

	for (auto & contig: contigs) {

		assert(ref_genome_seqs.count(contig.id()) > 0);
	} 

	uint num_variants = 0;

	Variant * cur_var;
	vcf_reader.getNextVariant(&cur_var);

    uniform_int_distribution<int> canonical_bases_sampler(0, canonical_bases.size() - 1);

	for (auto & contig: contigs) {		

		auto ref_genome_seqs_it = ref_genome_seqs.find(contig.id());
		assert(ref_genome_seqs_it != ref_genome_seqs.end());

		for (auto & diplotype: sample_diplotypes) {

			diplotype.second.max_allele_end_pos.first = 0;
			diplotype.second.max_allele_end_pos.second = 0;

			diplotype.second.haplotype_seqs.first = new FastaRecord(contig.id() + "_1", ref_genome_seqs_it->second->seq().size() * 1.1);
			diplotype.second.haplotype_seqs.second = new FastaRecord(contig.id() + "_2",  ref_genome_seqs_it->second->seq().size() * 1.1);
		}

		uint prev_position = 0; 

		while (cur_var) {

			if (cur_var->chrom() != contig.id()) {

				break;
			}

			num_variants++;

			assert(prev_position < cur_var->pos());
			prev_position = cur_var->pos();

    		assert(cur_var->pos() <= ref_genome_seqs_it->second->seq().size());

    		if (cur_var->ref().seq().find_first_of("N") == string::npos) {

				assert(cur_var->ref().seq() == ref_genome_seqs_it->second->seq().substr(cur_var->pos() - 1, cur_var->ref().seq().size()));
    		}

    		assert(cur_var->numAlts() > 0);

			vector<pair<uint, string> > trimmed_alt_alleles;
			trimmed_alt_alleles.reserve(cur_var->numAlts());

			for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

            	Allele ref_allele = cur_var->ref();
            	Allele alt_allele = cur_var->alt(alt_idx);

            	assert(!ref_allele.isMissing());
            	assert(!alt_allele.isMissing());

                Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);                
                trimmed_alt_alleles.emplace_back(ref_allele.seq().size(), alt_allele.seq());

			    auto trimmed_alt_alleles_seq_it = trimmed_alt_alleles.back().second.begin();
			    assert(trimmed_alt_alleles_seq_it != trimmed_alt_alleles.back().second.end());

    			if (argc == 5) {

				    while (trimmed_alt_alleles_seq_it != trimmed_alt_alleles.back().second.end()) {

				        if (*trimmed_alt_alleles_seq_it == 'N') {

				    		*trimmed_alt_alleles_seq_it = canonical_bases.at(canonical_bases_sampler(prng));
				        }

				        trimmed_alt_alleles_seq_it++;
				    }
				}
			}

			for (auto & diplotype: sample_diplotypes) {

				Sample & cur_sample = cur_var->getSample(diplotype.first);

				assert(cur_sample.callStatus() == Sample::CallStatus::Complete);
				assert(cur_sample.ploidy() == Sample::Ploidy::Diploid);
				assert(cur_sample.isPhased());

				auto cur_genotype = cur_sample.genotypeEstimate();
				assert(cur_genotype.size() == 2);

				if ((cur_genotype.front() > 0) and (diplotype.second.max_allele_end_pos.first < cur_var->pos())) {

					diplotype.second.haplotype_seqs.first->appendSeq(string(ref_genome_seqs_it->second->seq().begin() + diplotype.second.max_allele_end_pos.first, ref_genome_seqs_it->second->seq().begin() + cur_var->pos() - 1));
					diplotype.second.haplotype_seqs.first->appendSeq(trimmed_alt_alleles.at(cur_genotype.front() - 1).second);

  					diplotype.second.max_allele_end_pos.first = cur_var->pos() + trimmed_alt_alleles.at(cur_genotype.front() - 1).first - 1;
    			}

				if ((cur_genotype.back() > 0) and (diplotype.second.max_allele_end_pos.second < cur_var->pos())) {

					diplotype.second.haplotype_seqs.second->appendSeq(string(ref_genome_seqs_it->second->seq().begin() + diplotype.second.max_allele_end_pos.second, ref_genome_seqs_it->second->seq().begin() + cur_var->pos() - 1));
					diplotype.second.haplotype_seqs.second->appendSeq(trimmed_alt_alleles.at(cur_genotype.back() - 1).second);

  					diplotype.second.max_allele_end_pos.second = cur_var->pos() + trimmed_alt_alleles.at(cur_genotype.back() - 1).first - 1;
    			}
			}

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}

			delete cur_var;
			vcf_reader.getNextVariant(&cur_var);
		}

		if (prev_position > 0) {

			for (auto & diplotype: sample_diplotypes) {
		
	    		assert(diplotype.second.max_allele_end_pos.first <= ref_genome_seqs_it->second->seq().size());
				assert(diplotype.second.max_allele_end_pos.second <= ref_genome_seqs_it->second->seq().size());

				diplotype.second.haplotype_seqs.first->appendSeq(string(ref_genome_seqs_it->second->seq().begin() + diplotype.second.max_allele_end_pos.first, ref_genome_seqs_it->second->seq().end()));
				diplotype.second.haplotype_seqs.second->appendSeq(string(ref_genome_seqs_it->second->seq().begin() + diplotype.second.max_allele_end_pos.second, ref_genome_seqs_it->second->seq().end()));

		    	if (argc == 5) {

					assert(diplotype.second.haplotype_seqs.first->seq().find_first_of("N") == string::npos);
					assert(diplotype.second.haplotype_seqs.second->seq().find_first_of("N") == string::npos);
				}

				*diplotype.second.genome_outfile << diplotype.second.haplotype_seqs.first->str();
				*diplotype.second.genome_outfile << diplotype.second.haplotype_seqs.second->str();
			}
		}

		for (auto & diplotype: sample_diplotypes) {

			delete diplotype.second.haplotype_seqs.first;
			delete diplotype.second.haplotype_seqs.second;
		}
	}

	for (auto &chr_seq: ref_genome_seqs) {

        delete chr_seq.second;
    }

	for (auto & diplotype: sample_diplotypes) {

		diplotype.second.genome_outfile->close();
		delete diplotype.second.genome_outfile;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Wrote diplotypes using " << num_variants << " variants" << endl;
    cout << endl;

	return 0;
}
