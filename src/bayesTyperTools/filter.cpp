
/*
filter.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <algorithm>

#include "VcfFile.hpp"
#include "Stats.hpp"
#include "Auxiliaries.hpp"
#include "FastaRecord.hpp"
#include "FastaReader.hpp"
#include "JoiningString.hpp"

#include "filter.hpp"

static const float fak_beta_value = 0.275;
static const uint min_inbreeding_samples = 10;
static const uint min_homozygote_samples = 10;

namespace Filter {

    pair<float, bool> calcInbreedingCoefficient(Variant & cur_var, const Contig::Type chrom_type, const ushort num_indepedent_samples, const regex & indepedent_sample_regex) {

        if (chrom_type == Contig::Type::Autosomal) {

            auto inbreeding_stats = Stats::calcInbreedingStats(cur_var, indepedent_sample_regex);
            assert(inbreeding_stats.num_samples == num_indepedent_samples);

            if (!(inbreeding_stats.is_fixed)) {

                return make_pair(abs(inbreeding_stats.inbreeding_coef), true);
            }
        }

        return make_pair(0, false);
    }

    vector<bool> getFilteredHomopolymerAlleles(Variant & cur_var, const string & chrom_seq, const uint max_homopolymer_length) {

        auto homopolymer_info = Auxiliaries::getHomopolymerInfo(cur_var.pos(), chrom_seq);   

        if (homopolymer_info.first > max_homopolymer_length) {

            return Auxiliaries::getHomopolymerAlleles(cur_var, homopolymer_info.second, 1);

        } else {

            return vector<bool>(cur_var.numAlls(), false);
        }
    }

    void filter(const string & vcf_filename, const string & genome_filename, const string & output_prefix, const string & kmer_coverage_filename, string indepedent_samples_regex_str, FilterValues filter_values) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") filter ...\n" << endl;
        cout << "[" << Utils::getLocalTime() << "] Parsing reference genome fasta ..." << endl;

        unordered_map<string, FastaRecord*> genome_seqs;
        FastaReader genome_reader(genome_filename);
        FastaRecord * cur_fasta_rec;

        while (genome_reader.getNextRecord(&cur_fasta_rec)) {

            assert(genome_seqs.emplace(cur_fasta_rec->id(), cur_fasta_rec).second);
            genome_seqs.at(cur_fasta_rec->id())->convertToUppercase();
        }

        cout << "[" << Utils::getLocalTime() << "] Parsed " << genome_seqs.size() << " chromosome(s)\n" << endl;

        GenotypedVcfFileReader vcf_reader(vcf_filename, true);
        
        auto sample_ids = vcf_reader.metaData().sampleIds();
        assert(!(sample_ids.empty()));

        uint num_indepedent_samples = 0;

        regex indepedent_sample_regex(indepedent_samples_regex_str);

        for (auto & sample_id: sample_ids) {

            if (regex_match(sample_id, indepedent_sample_regex)) {

                num_indepedent_samples++;
            }
        }

        if (num_indepedent_samples >= min_inbreeding_samples) {

            assert(filter_values.max_inbreeding_coef >= 0);
            assert(filter_values.max_inbreeding_coef <= 1);

            cout << "[" << Utils::getLocalTime() << "] Calculating inbreeding coefficient using " << num_indepedent_samples << " independent samples\n" << endl;
            cout << "[" << Utils::getLocalTime() << "] Adding absolute inbreeding coefficient filter (> " << filter_values.max_inbreeding_coef << ")" << endl;
        }

        if (sample_ids.size() >= min_homozygote_samples) {

            assert(filter_values.min_homozygote >= 0);

            cout << "[" << Utils::getLocalTime() << "] Adding number of homozygote unfiltered genotypes filter (< " << filter_values.min_homozygote << ")" << endl;
        }

        if (filter_values.max_homopolymer_length != -1) {

            assert(filter_values.max_homopolymer_length >= 0);   

            cout << "[" << Utils::getLocalTime() << "] Adding homopolymer length filter (> " << filter_values.max_homopolymer_length << ")" << endl;
        }

        assert(filter_values.min_nak_value >= 0);   

        cout << "[" << Utils::getLocalTime() << "] Adding number of allele kmer (NAK) filter (< " << filter_values.min_nak_value << ")" << endl;

        filter_values.min_sample_fak_values = vector<float>(sample_ids.size(), 0);

        if (!(kmer_coverage_filename.empty())) {

            ifstream kmer_coverage_file(kmer_coverage_filename);
            assert(kmer_coverage_file.is_open());

            string line;

            while (kmer_coverage_file.good()) {

                getline(kmer_coverage_file, line);

                if (line.empty()) {

                    continue;
                }

                auto line_split = Utils::splitString(line, '\t');
                assert(line_split.size() == 3); 

                if (line_split.front() == "Sample") {

                    continue;
                }

                auto sample_ids_it = find(sample_ids.begin(), sample_ids.end(), line_split.front());

                if (sample_ids_it != sample_ids.end()) {

                    assert(floatCompare(filter_values.min_sample_fak_values.at(sample_ids_it - sample_ids.begin()), 0));
                    filter_values.min_sample_fak_values.at(sample_ids_it - sample_ids.begin()) = 1 - exp(-(fak_beta_value * stof(line_split.at(1))));

                    assert(filter_values.min_sample_fak_values.at(sample_ids_it - sample_ids.begin()) >= 0);
                    assert(filter_values.min_sample_fak_values.at(sample_ids_it - sample_ids.begin()) <= 1);
                }
            }
    
            kmer_coverage_file.close();
        }

        cout << "[" << Utils::getLocalTime() << "] Adding sample specific observed kmer fraction (FAK) filter (" << sample_ids.front() << " < " << filter_values.min_sample_fak_values.front();
        
        for (ushort sample_idx = 1; sample_idx < sample_ids.size(); sample_idx++) {

            assert(filter_values.min_sample_fak_values.at(sample_idx) >= 0);
            assert(filter_values.min_sample_fak_values.at(sample_idx) <= 1);

            cout << ", " << sample_ids.at(sample_idx) << " < " << filter_values.min_sample_fak_values.at(sample_idx);
        }

        cout << ")" << endl;

        assert(filter_values.min_gpp_value >= 0);
        assert(filter_values.min_gpp_value <= 1);

        cout << "[" << Utils::getLocalTime() << "] Adding genotype posterior probability (GPP) filter (< " << filter_values.min_gpp_value << ")" << endl;
        cout << endl;

        auto output_meta_data = vcf_reader.metaData();
        Auxiliaries::removeNonRelevantFilterDescriptors(&output_meta_data, {"UV"});

        string all_samples_filtered_label("ASF");
        output_meta_data.filterDescriptors().emplace(all_samples_filtered_label, Attribute::Descriptor(all_samples_filtered_label, "All samples filtered"));

        output_meta_data.formatDescriptors().emplace("SAF", Attribute::DetailedDescriptor({make_pair("ID","SAF"), make_pair("Number","R"), make_pair("Type","String"), make_pair("Description","Sample allele filter ('P': Pass, 'F': Filtered)")}));

        VcfFileWriter vcf_writer(output_prefix + ".vcf", output_meta_data, true);

        ulong num_variants = 0;
        map<string, ulong> variant_filter_labels;

        ulong num_genotypes = 0;
        ulong num_genotypes_unfiltered = 0;
        ulong num_genotypes_filtered = 0;

        vector<uint> max_filtered_allele_end_pos(sample_ids.size(), 0);

        Variant * cur_var;

        string cur_chromosome = "";
        Contig::Type cur_chromosome_type = Contig::Type::Unknown;

        auto genome_seqs_it = genome_seqs.find(cur_chromosome);
        assert(genome_seqs_it == genome_seqs.end());

        while (vcf_reader.getNextVariant(&cur_var)) {

            num_variants++;

            if (find(cur_var->filters().begin(), cur_var->filters().end(), "UV") != cur_var->filters().end()) {
                
                cur_var->setFilters({"UV"});

                auto variant_filter_labels_it = variant_filter_labels.emplace("UV", 0);
                variant_filter_labels_it.first->second++;

                vcf_writer.write(cur_var);
                delete cur_var;

                continue;            
            }

            Auxiliaries::resetFilters(cur_var);

            if (cur_var->chrom() != cur_chromosome) {

                max_filtered_allele_end_pos = vector<uint>(sample_ids.size(), 0);
                genome_seqs_it = genome_seqs.find(cur_var->chrom());

                cur_chromosome = cur_var->chrom();
                cur_chromosome_type = vcf_reader.metaData().getContig(cur_var->chrom()).type();
            }

            assert(genome_seqs_it != genome_seqs.end());

            vector<bool> filtered_alleles(cur_var->numAlls(), false);

            if (filter_values.max_homopolymer_length != -1) {

                assert(filter_values.max_homopolymer_length >= 0);
                filtered_alleles = getFilteredHomopolymerAlleles(*cur_var, genome_seqs_it->second->seq(), filter_values.max_homopolymer_length);
            }

            if (num_indepedent_samples >= min_inbreeding_samples) {

                auto inbreeding_coefficient = calcInbreedingCoefficient(*cur_var, cur_chromosome_type, num_indepedent_samples, indepedent_sample_regex);

                if (inbreeding_coefficient.second and (inbreeding_coefficient.first > filter_values.max_inbreeding_coef)) {

                    filtered_alleles = vector<bool>(cur_var->numAlls(), true);
                }
            }
            
            if (sample_ids.size() >= min_homozygote_samples) {

                uint num_homozygote_genotypes = 0;

                for (ushort sample_idx = 0; sample_idx < sample_ids.size(); sample_idx++) {

                    Sample * cur_sample = &(cur_var->getSample(sample_ids.at(sample_idx)));
                    
                    assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);
                    assert((cur_sample->callStatus() == Sample::CallStatus::Complete) or (cur_sample->callStatus() == Sample::CallStatus::Missing));
                    
                    if (cur_sample->ploidy() == Sample::Ploidy::Zeroploid) {

                        num_homozygote_genotypes++;                 
                        continue;
                    }

                    if (cur_sample->callStatus() == Sample::CallStatus::Complete) {

                        assert(!(cur_sample->genotypeEstimate().empty()));
                        assert(cur_sample->genotypeEstimate().size() <= 2);

                        if (cur_sample->genotypeEstimate().front() == cur_sample->genotypeEstimate().back()) {

                            num_homozygote_genotypes++;
                        }
                    }   
                }

                if (num_homozygote_genotypes < filter_values.min_homozygote) {

                    filtered_alleles = vector<bool>(cur_var->numAlls(), true);
                }
            }

            assert(filtered_alleles.size() == cur_var->numAlls());

            for (ushort sample_idx = 0; sample_idx < sample_ids.size(); sample_idx++) {

                Sample * cur_sample = &(cur_var->getSample(sample_ids.at(sample_idx)));
                
                assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);
                assert((cur_sample->callStatus() == Sample::CallStatus::Complete) or (cur_sample->callStatus() == Sample::CallStatus::Missing));
                
                if (cur_sample->ploidy() == Sample::Ploidy::Zeroploid) {

                    continue;
                }

                num_genotypes++;

                assert(cur_sample->alleleInfo().size() == filtered_alleles.size());
                assert(cur_sample->alleleInfo().size() == cur_var->numAlls());

                for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

                    auto saf_value = cur_sample->alleleInfo().at(allele_idx).getValue<string>("SAF");
                    assert(saf_value.second);
                    assert(saf_value.first == "P");

                    auto app_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("APP");
                    assert(app_value.second);

                    bool is_filtered = filtered_alleles.at(allele_idx);

                    if (!is_filtered) {

                        if (cur_var->pos() <= max_filtered_allele_end_pos.at(sample_idx)) {

                            assert(Auxiliaries::hasMissing(*cur_var));
                            is_filtered = true;
                        
                        } else if (!(cur_var->allele(allele_idx).isMissing())) {

                            if (!(Utils::floatCompare(app_value.first, 0))) {

                                assert(app_value.first > 0);

                                auto nak_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("NAK");                  
                                assert(nak_value.second);

                                auto fak_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("FAK");                  
                                assert(fak_value.second);

                                assert(nak_value.first >= 0);
               
                                if (nak_value.first < filter_values.min_nak_value) {

                                    is_filtered = true;
                                
                                } else if ((fak_value.first >= 0) and (fak_value.first < filter_values.min_sample_fak_values.at(sample_idx))) {
                                    
                                    is_filtered = true;
                                }
                            }
                        }
                    }

                    if (is_filtered) {

                        assert(!(cur_sample->alleleInfo().at(allele_idx).setValue<string>("SAF", "F")));
                    }
                }

                for (auto allele_idx: cur_sample->genotypeEstimate()) {

                    auto saf_value = cur_sample->alleleInfo().at(allele_idx).getValue<string>("SAF");
                    assert(saf_value.second);
                    assert((saf_value.first == "P") or (saf_value.first == "F"));

                    if (!(cur_var->allele(allele_idx).isMissing()) and (saf_value.first == "F")) {

                        Allele cur_ref_allele = cur_var->ref();
                        Allele cur_alt_allele = cur_var->allele(allele_idx);

                        Auxiliaries::rightTrimAllelePair(&cur_ref_allele, &cur_alt_allele);

                        assert(!(cur_ref_allele.seq().empty()));
                        assert(!(cur_alt_allele.seq().empty()));

                        max_filtered_allele_end_pos.at(sample_idx) = max(max_filtered_allele_end_pos.at(sample_idx), static_cast<uint>(cur_var->pos() + cur_ref_allele.seq().size() - 1));
                    }
                }

                auto max_gpp = Auxiliaries::getMaxGenotypePosterior(*cur_sample);
                assert(max_gpp.second or (cur_sample->callStatus() == Sample::CallStatus::Missing));

                if (max_gpp.first < filter_values.min_gpp_value) {

                    cur_sample->clearGenotypeEstimate();                    
                }   

                for (auto allele_idx: cur_sample->genotypeEstimate()) {

                    auto saf_value = cur_sample->alleleInfo().at(allele_idx).getValue<string>("SAF");
                    assert(saf_value.second);

                    if (saf_value.first == "F") {

                        cur_sample->clearGenotypeEstimate();
                        break;
                    } 
                }

                if (cur_sample->callStatus() == Sample::CallStatus::Missing) {

                    num_genotypes_filtered++;

                } else {

                    assert(cur_sample->callStatus() == Sample::CallStatus::Complete);
                    num_genotypes_unfiltered++;
                }
            }

            Auxiliaries::updateAlleleStatsAndCallProb(cur_var);

            auto an_value = cur_var->info().getValue<int>("AN");
            assert(an_value.second);

            if (an_value.first == 0) {

                cur_var->addFilter(all_samples_filtered_label);
            
            } else if (cur_var->filters().empty()) {

                cur_var->setFilters({"PASS"});                
            }

            JoiningString filter_string(';');
            filter_string.join(cur_var->filters());

            auto variant_filter_labels_it = variant_filter_labels.emplace(filter_string.str(), 0);
            variant_filter_labels_it.first->second++;

            vcf_writer.write(cur_var);
            delete cur_var;

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }
        }

        for (auto & genome_seq : genome_seqs) {

            delete genome_seq.second;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Completed filtering of " << num_variants << " variants:\n\n";

        ulong variant_filter_label_sum = 0;

        for (auto variant_filter_label: variant_filter_labels) {

            cout << "\t- " << variant_filter_label.second << " variants were labeled with " << variant_filter_label.first << endl;
            variant_filter_label_sum += variant_filter_label.second;
        }

        assert(num_variants == variant_filter_label_sum);

        cout << "\n[" << Utils::getLocalTime() << "] And out of " << num_genotypes << " genotypes " << num_genotypes_filtered << " were filtered";
        cout << endl;

        assert(num_genotypes == (num_genotypes_unfiltered + num_genotypes_filtered));
    }
}
