
/*
prepareSimulationFemaleSamples.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"
#include "vcf++/FastaRecord.hpp"
#include "vcf++/FastaReader.hpp"

const vector<char> canonical_bases = {'A','C','G','T'};
const uint seed = 12345678;

int main(int argc, char const *argv[]) {

	if (argc != 5) {

		std::cout << "USAGE: prepareSimulationFemaleSamples <vcf> <genome> <samples> <output_prefix>" << std::endl;
		return 1;
	}

	cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") prepareSimulationFemaleSamples script ...\n" << endl;

    mt19937 prng(seed);

    unordered_map<string,FastaRecord*> ref_genome_seqs;
    FastaReader genome_reader(argv[2]);
    FastaRecord * cur_fasta_rec;

    while (genome_reader.getNextRecord(&cur_fasta_rec)) {

    	cur_fasta_rec->convertToUppercase();
    	cur_fasta_rec->sampleAmbiguousBases(canonical_bases, &prng);

        assert(ref_genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);

        cout << "[" << Utils::getLocalTime() << "] Parsed sequence: " << cur_fasta_rec->id() << endl;
    }

    unordered_set<string> samples;

    ifstream samples_file(argv[3]);
    string line;

    while (getline(samples_file, line)) {

    	if (!(line.empty())) {

	        assert(samples.insert(line).second);
    	}
    }

    samples_file.close();

	cout << "\n[" << Utils::getLocalTime() << "] Preparing " << samples.size() << " sample(s) for simulation ..." << endl;

	string vcf_filename(argv[1]);
	GenotypedVcfFileReader vcf_reader(vcf_filename, true);

	vcf_reader.metaData().infoDescriptors().clear();
	vcf_reader.metaData().filterDescriptors().clear();
	vcf_reader.metaData().miscMeta().clear();

	assert(vcf_reader.metaData().filterDescriptors().emplace("REFCALL", Attribute::Descriptor("REFCALL", "Reference call")).second);

	Auxiliaries::removeNonRelevantFormatDescriptors(&(vcf_reader.metaData()), {"GT"});

    auto sample_ids = vcf_reader.metaData().sampleIds();

    unordered_map<string, pair<uint, uint> > sample_last_allele_end_position;
    unordered_map<string, vector<pair<FastaRecord *, FastaRecord *> > > sample_genome_seqs;
    unordered_map<string, ofstream *> sample_genome_writers;

    uint num_rm_samples = 0;

    for (auto &sample_id: sample_ids) {

        if (samples.count(sample_id) < 1) {

            vcf_reader.metaData().rmSample(sample_id);
            num_rm_samples++;
        
        } else {

        	assert(sample_last_allele_end_position.emplace(sample_id, make_pair(0,0)).second);
			assert(sample_genome_seqs.emplace(sample_id, vector<pair<FastaRecord *, FastaRecord *> >()).second);
			assert(sample_genome_writers.emplace(sample_id, new ofstream(string(argv[4]) + "_" + sample_id + "_genome.fa")).second);
        }
    }

    sample_ids = vcf_reader.metaData().sampleIds();

    cout << "[" << Utils::getLocalTime() << "] Removed " << num_rm_samples << " sample(s)\n" << endl;

    assert(sample_ids.size() == samples.size());
    assert(sample_ids.size() == sample_last_allele_end_position.size());
    assert(sample_ids.size() == sample_genome_seqs.size());
    assert(sample_ids.size() == sample_genome_writers.size());

    uniform_int_distribution<int> canonical_bases_sampler(0, canonical_bases.size() - 1);

    VcfFileWriter vcf_writer(string(argv[4]) + ".vcf", vcf_reader.metaData(), true);

	Variant * cur_var;

	uint num_variants = 0;
	uint num_alternative_calls = 0;

	string cur_chr = "";
    uint prev_position = 0;

    unordered_set<string> parsed_chr;
	FastaRecord * ref_chr_seq = nullptr;

	while (vcf_reader.getNextVariant(&cur_var)) {

		if ((vcf_reader.metaData().getContig(cur_var->chrom()).type() == Contig::Type::Autosomal) or (vcf_reader.metaData().getContig(cur_var->chrom()).type() == Contig::Type::ChrX)) {

			assert(!(cur_var->chrom().empty()));

			if (cur_chr != cur_var->chrom()) {

				if (ref_chr_seq) {

					assert(!(cur_chr.empty()));

					for (auto &sample: sample_genome_seqs) {

						auto sample_last_allele_end_position_it = sample_last_allele_end_position.find(sample.first);
						assert(sample_last_allele_end_position_it != sample_last_allele_end_position.end());

                		assert(sample_last_allele_end_position_it->second.first <= ref_chr_seq->seq().size());
                		assert(sample_last_allele_end_position_it->second.second <= ref_chr_seq->seq().size());

                		assert(sample.second.size() == parsed_chr.size());

						sample.second.back().first->appendSeq(string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.first, ref_chr_seq->seq().end()));
						sample.second.back().second->appendSeq(string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.second, ref_chr_seq->seq().end()));

						assert(sample.second.back().first->seq().find_first_of("N") == string::npos);
						assert(sample.second.back().second->seq().find_first_of("N") == string::npos);

						auto sample_genome_writers_it = sample_genome_writers.find(sample.first);
						assert(sample_genome_writers_it != sample_genome_writers.end());

						*(sample_genome_writers_it->second) << sample.second.back().first->str();
						*(sample_genome_writers_it->second) << sample.second.back().second->str();

				    	delete sample.second.back().first;
				   		delete sample.second.back().second;
					}					
				
				} else {

					assert(cur_chr.empty());
				}

				for (auto &end_positions: sample_last_allele_end_position) {

					end_positions.second.first = 0;
					end_positions.second.second = 0;
				}

				cur_chr = cur_var->chrom();
				prev_position = 0;
				
				assert(parsed_chr.insert(cur_var->chrom()).second);

				auto ref_genome_seqs_it = ref_genome_seqs.find(cur_var->chrom());
				assert(ref_genome_seqs_it != ref_genome_seqs.end());

				ref_chr_seq = ref_genome_seqs_it->second;

				for (auto &sample: sample_genome_seqs) {

					sample.second.emplace_back(new FastaRecord(cur_var->chrom() + "_1", ref_chr_seq->seq().size() * 1.05), new FastaRecord(cur_var->chrom() + "_2", ref_chr_seq->seq().size() * 1.05));
				}
			}

			num_variants++;

			assert(prev_position < cur_var->pos());
			prev_position = cur_var->pos();

    		assert(cur_var->pos() <= ref_chr_seq->seq().size());

    		if (cur_var->ref().seq().find_first_of("N") == string::npos) {

				assert(cur_var->ref().seq() == ref_chr_seq->seq().substr(cur_var->pos() - 1, cur_var->ref().seq().size()));
    		}

			vector<pair<uint, string> > alt_alleles;
			alt_alleles.reserve(cur_var->numAlts());

			for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

            	Allele ref_allele = cur_var->ref();
            	Allele alt_allele = cur_var->alt(alt_idx);

            	assert(!(ref_allele.isMissing()));
            	assert(!(alt_allele.isMissing()));

                Auxiliaries::rightTrimAllelePair(&ref_allele, &alt_allele);
                alt_alleles.emplace_back(ref_allele.seq().size(), alt_allele.seq());

			    auto alt_allele_seq_it = alt_alleles.back().second.begin();
			    assert(alt_allele_seq_it != alt_alleles.back().second.end());

			    if (*alt_allele_seq_it == 'N') {

			    	if (alt_alleles.back().first != alt_alleles.back().second.size()) {

			    		assert((alt_alleles.back().first == 1) or (alt_alleles.back().second.size() == 1));
			    		*alt_allele_seq_it = ref_chr_seq->seq().at(cur_var->pos() - 1);

			    	} else {	

			    		*alt_allele_seq_it = canonical_bases.at(canonical_bases_sampler(prng));
			    	}
			    } 

			    alt_allele_seq_it++;

			    while (alt_allele_seq_it != alt_alleles.back().second.end()) {

			        if (*alt_allele_seq_it == 'N') {

			    		*alt_allele_seq_it = canonical_bases.at(canonical_bases_sampler(prng));			            
			        }

			        alt_allele_seq_it++;
			    }
			} 

			bool has_alternative_call = false;
			
			for (auto &sample_id: sample_ids) {

				Sample * cur_sample = &(cur_var->getSample(sample_id));

				assert(cur_sample->ploidy() == Sample::Ploidy::Diploid);
				assert(cur_sample->callStatus() == Sample::CallStatus::Complete);
				assert(cur_sample->isInformative());

				vector<ushort> cur_genotype = cur_sample->genotypeEstimate();
				assert(cur_genotype.size() == 2);

				vector<ushort> new_genotype = cur_genotype;

				auto sample_genome_seqs_it = sample_genome_seqs.find(sample_id);
				assert(sample_genome_seqs_it != sample_genome_seqs.end());

				auto sample_last_allele_end_position_it = sample_last_allele_end_position.find(sample_id);
				assert(sample_last_allele_end_position_it != sample_last_allele_end_position.end());

				if (cur_var->pos() <= sample_last_allele_end_position_it->second.first) {

					Allele missing_allele = Allele();
					cur_var->addAlt(missing_allele);

					new_genotype.front() = cur_var->numAlts();
				
				} else if (cur_genotype.front() > 0) {

					has_alternative_call = true;

					string haplotype_sequence = string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.first, ref_chr_seq->seq().begin() + cur_var->pos() - 1); 
					haplotype_sequence += alt_alleles.at(cur_genotype.front() - 1).second;
					sample_genome_seqs_it->second.back().first->appendSeq(haplotype_sequence);

  					sample_last_allele_end_position_it->second.first = cur_var->pos() + alt_alleles.at(cur_genotype.front() - 1).first - 1;
    			}

				if (cur_var->pos() <= sample_last_allele_end_position_it->second.second) {

					Allele missing_allele = Allele();
					cur_var->addAlt(missing_allele);

					new_genotype.back() = cur_var->numAlts();
				
				} else if (cur_genotype.back() > 0) {

					has_alternative_call = true;

					string haplotype_sequence = string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.second, ref_chr_seq->seq().begin() + cur_var->pos() - 1); 
					haplotype_sequence += alt_alleles.at(cur_genotype.back() - 1).second;
					sample_genome_seqs_it->second.back().second->appendSeq(haplotype_sequence);

  					sample_last_allele_end_position_it->second.second = cur_var->pos() + alt_alleles.at(cur_genotype.back() - 1).first - 1;
    			}

    			cur_sample->newGenotypeEstimate(new_genotype);
			}

			cur_var->setIds({});
			cur_var->setQual({0, false});

			if (has_alternative_call) {

				num_alternative_calls++;
				cur_var->setFilters({"PASS"});
			
			} else {

				cur_var->setFilters({"REFCALL"});
			}

			if ((num_variants % 100000) == 0) {

				std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
			}

			vcf_writer.write(cur_var);
		}

		delete cur_var;
	}

 	assert(ref_chr_seq);
	assert(!(cur_chr.empty()));

	for (auto &sample: sample_genome_seqs) {

		auto sample_last_allele_end_position_it = sample_last_allele_end_position.find(sample.first);
		assert(sample_last_allele_end_position_it != sample_last_allele_end_position.end());

		assert(sample_last_allele_end_position_it->second.first <= ref_chr_seq->seq().size());
		assert(sample_last_allele_end_position_it->second.second <= ref_chr_seq->seq().size());

		assert(sample.second.size() == parsed_chr.size());

		sample.second.back().first->appendSeq(string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.first, ref_chr_seq->seq().end()));
		sample.second.back().second->appendSeq(string(ref_chr_seq->seq().begin() + sample_last_allele_end_position_it->second.second, ref_chr_seq->seq().end()));

		assert(sample.second.back().first->seq().find_first_of("N") == string::npos);
		assert(sample.second.back().second->seq().find_first_of("N") == string::npos);

		auto sample_genome_writers_it = sample_genome_writers.find(sample.first);
		assert(sample_genome_writers_it != sample_genome_writers.end());

		*(sample_genome_writers_it->second) << sample.second.back().first->str();
		*(sample_genome_writers_it->second) << sample.second.back().second->str();

    	delete sample.second.back().first;
   		delete sample.second.back().second;

       	sample_genome_writers_it->second->close();
       	delete sample_genome_writers_it->second;
	}	

	for (auto &chr_seq: ref_genome_seqs) {

        delete chr_seq.second;
    }

	cout << "\n[" << Utils::getLocalTime() << "] Wrote " << sample_ids.size() << " diploid genomes" << endl;

	cout << "\n[" << Utils::getLocalTime() << "] Prepared " << num_variants << " variants for simulation:\n" << endl;
    cout << "\t- Reference calls: " << num_variants - num_alternative_calls << endl;
    cout << "\t- Alternative allele calls: " << num_alternative_calls << endl;

    cout << endl;
	return 0;
}
