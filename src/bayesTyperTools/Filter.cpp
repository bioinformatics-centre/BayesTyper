
/*
Filter.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "Filter.hpp"


static const float fak_beta_value = 0.275;

namespace Filter {

    void filter(const string & variant_file, const string & outfile, const string & kmer_coverage_filename, FilterValues filter_values, const uint min_filter_samples) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperTools (" << BT_VERSION << ") filter ...\n" << endl;

        GenotypedVcfFileReader vcf_reader(variant_file, true);
        
        auto sample_ids = vcf_reader.metaData().sampleIds();
        assert(!(sample_ids.empty()));

        if (sample_ids.size() >= min_filter_samples) {

            assert(filter_values.min_homozygote >= 0);

            cout << "[" << Utils::getLocalTime() << "] Adding number of homozygote genotypes filter (< " << filter_values.min_homozygote << ")" << endl;
        }

        assert(filter_values.min_nak_value >= 0);   

        cout << "[" << Utils::getLocalTime() << "] Adding number of allele kmer (NAK) filter (< " << filter_values.min_nak_value << ")" << endl;

        assert(filter_values.min_gpp_value >= 0);
        assert(filter_values.min_gpp_value <= 1);

        cout << "[" << Utils::getLocalTime() << "] Adding genotype posterior probability (GPP) filter (< " << filter_values.min_gpp_value << ")" << endl;

        filter_values.min_sample_fak_values = vector<float>(sample_ids.size(), 0);

        if (!(kmer_coverage_filename.empty())) {

            ifstream kmer_coverage_infile(kmer_coverage_filename);
            assert(kmer_coverage_infile.is_open());

            string line;

            while (kmer_coverage_infile.good()) {

                getline(kmer_coverage_infile, line);

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
    
            kmer_coverage_infile.close();
        }

        cout << "[" << Utils::getLocalTime() << "] Adding sample specific observed kmer fraction (FAK) filter (" << sample_ids.front() << " < " << filter_values.min_sample_fak_values.front();
        
        for (ushort sample_idx = 1; sample_idx < sample_ids.size(); sample_idx++) {

            assert(filter_values.min_sample_fak_values.at(sample_idx) >= 0);
            assert(filter_values.min_sample_fak_values.at(sample_idx) <= 1);

            cout << ", " << sample_ids.at(sample_idx) << " < " << filter_values.min_sample_fak_values.at(sample_idx);
        }

        cout << ")\n" << endl;

        auto output_meta_data = vcf_reader.metaData();

        output_meta_data.filterDescriptors().erase("HOM");

        if (sample_ids.size() >= min_filter_samples) {

            assert(output_meta_data.filterDescriptors().emplace("HOM", Attribute::Descriptor("HOM", "Less than " + to_string(filter_values.min_homozygote) + " homozygote genotypes (calculated before other filters)")).second);
        }
        
        assert(output_meta_data.filterDescriptors().count("AN0") > 0);
        assert(output_meta_data.formatDescriptors().count("SAF") > 0);

        VcfFileWriter vcf_writer(outfile, output_meta_data, true);

        ulong num_variants = 0;
        map<string, ulong> variant_filter_labels;

        ulong num_genotypes = 0;
        ulong num_genotypes_unfiltered = 0;
        ulong num_genotypes_filtered = 0;

        Variant * cur_var;

        while (vcf_reader.getNextVariant(&cur_var)) {

            num_variants++;
            Auxiliaries::resetFilters(cur_var);
            
            if (sample_ids.size() >= min_filter_samples) {

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

                    cur_var->addFilter("HOM");
                }
            }

            for (ushort sample_idx = 0; sample_idx < sample_ids.size(); sample_idx++) {

                Sample * cur_sample = &(cur_var->getSample(sample_ids.at(sample_idx)));
                
                assert(cur_sample->ploidy() != Sample::Ploidy::Polyploid);
                assert((cur_sample->callStatus() == Sample::CallStatus::Complete) or (cur_sample->callStatus() == Sample::CallStatus::Missing));
                
                if (cur_sample->ploidy() == Sample::Ploidy::Zeroploid) {

                    continue;
                }

                num_genotypes++;

                assert(cur_sample->alleleInfo().size() == cur_var->numAlls());

                for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

                    auto saf_value = cur_sample->alleleInfo().at(allele_idx).getValue<string>("SAF");
                    
                    if (!(saf_value.second)) {

                        continue;
                    }

                    assert(saf_value.first == "P");

                    auto app_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("APP");
                    assert(app_value.second);

                    if (!(Utils::floatCompare(app_value.first, 0))) {

                        assert(app_value.first > 0);

                        auto nak_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("NAK");                  
                        assert(nak_value.second);

                        auto fak_value = cur_sample->alleleInfo().at(allele_idx).getValue<float>("FAK");                  
                        assert(fak_value.second);

                        assert(nak_value.first >= 0);
       
                        if (Utils::floatLess(nak_value.first, filter_values.min_nak_value)) {

                            assert(!(cur_sample->alleleInfo().at(allele_idx).setValue<string>("SAF", "F")));
                        
                        } else if ((fak_value.first >= 0) and Utils::floatLess(fak_value.first, filter_values.min_sample_fak_values.at(sample_idx))) {
                            
                            assert(!(cur_sample->alleleInfo().at(allele_idx).setValue<string>("SAF", "F")));
                        }
                    }
                }

                auto max_gpp = Auxiliaries::getMaxGenotypePosterior(*cur_sample);
                assert(max_gpp.second or (cur_sample->callStatus() == Sample::CallStatus::Missing));

                if (Utils::floatLess(max_gpp.first, filter_values.min_gpp_value)) {

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

                cur_var->addFilter("AN0");
            
            } else if (cur_var->filters().empty()) {

                cur_var->setFilters({"PASS"});                
            }

            auto cur_filters = cur_var->filters();
            sort(cur_filters.begin(), cur_filters.end());

            JoiningString filter_string(';');
            filter_string.join(cur_filters);

            auto variant_filter_labels_it = variant_filter_labels.emplace(filter_string.str(), 0);
            variant_filter_labels_it.first->second++;

            vcf_writer.write(cur_var);
            delete cur_var;

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }
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
