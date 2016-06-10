
/*
simplifyDoublePolymorphisms.cpp - This file is part of BayesTyper (v0.9)


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


#include <sstream>
#include <math.h>
#include <thread>
#include <map>

#include "vcf++/VcfFile.hpp"
#include "vcf++/Variant.hpp"
#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

void writeVariant(VcfFileWriter * output_vcf, Variant * cur_var) {

    cur_var->removeRedundantAlts();
    Auxiliaries::rightTrimVariant(cur_var);

    assert(cur_var->numAlts() > 0);
        
    cur_var->setIds({});
    cur_var->setQual({0, false});
    cur_var->setFilters({});

    output_vcf->write(cur_var);
}

void appendSeqToAlleles(Variant * cur_var, const string & new_seq) {

    for (uint allele_idx = 0; allele_idx < cur_var->numAlls(); allele_idx++) {

        assert(!(cur_var->allele(allele_idx).isMissing()));
        cur_var->allele(allele_idx).seq().append(new_seq);
    }
}

void addVariantToBuffer(map<uint, Variant *> * variant_buffer, Variant * cur_var) {

    auto variant_buffer_it = variant_buffer->find(cur_var->pos());

    if (variant_buffer_it != variant_buffer->end()) {

        uint min_reference_length = min(variant_buffer_it->second->ref().seq().size(), cur_var->ref().seq().size());
        assert(variant_buffer_it->second->ref().seq().substr(0, min_reference_length) == cur_var->ref().seq().substr(0, min_reference_length));

        if (variant_buffer_it->second->ref().seq().size() < cur_var->ref().seq().size()) {

            appendSeqToAlleles(variant_buffer_it->second, cur_var->ref().seq().substr(min_reference_length));

        } else if (variant_buffer_it->second->ref().seq().size() > cur_var->ref().seq().size()) {

            appendSeqToAlleles(cur_var, variant_buffer_it->second->ref().seq().substr(min_reference_length));
        }

        assert(variant_buffer_it->second->ref().seq() == cur_var->ref().seq());

        for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

            variant_buffer_it->second->addAlt(cur_var->alt(alt_idx));
        } 

        delete cur_var;
    
    } else {

        assert(variant_buffer->emplace(cur_var->pos(), cur_var).second);
    }
}

Variant * simplifyNucleotidePolymorphism(Variant * cur_var, const uint alt_idx, Allele & trimmed_ref_allele, Allele & trimmed_alt_allele, const uint trimmed_left_bases) {

    if (trimmed_ref_allele.seq().substr(1, trimmed_ref_allele.seq().size() - 2) == trimmed_alt_allele.seq().substr(1, trimmed_alt_allele.seq().size() - 2)) {

        cur_var->alt(alt_idx).seq().at(trimmed_left_bases + trimmed_alt_allele.seq().size() - 1) = cur_var->ref().seq().at(trimmed_left_bases + trimmed_ref_allele.seq().size() - 1);

        return new Variant(cur_var->chrom(), cur_var->pos() + trimmed_left_bases + trimmed_ref_allele.seq().size() - 1, Allele(trimmed_ref_allele.seq().substr(trimmed_ref_allele.seq().size() - 1)), vector<Allele>(1, Allele(trimmed_alt_allele.seq().substr(trimmed_alt_allele.seq().size() - 1))), cur_var->info());
    
    } else {

        return nullptr;
    }
}

Variant * simplifyMixedPolymorphism(Variant * cur_var, const uint alt_idx, Allele & trimmed_ref_allele, Allele & trimmed_alt_allele, const uint trimmed_left_bases) {

    uint min_allele_length = min(trimmed_ref_allele.seq().size(), trimmed_alt_allele.seq().size());

    if (trimmed_ref_allele.seq().substr(1, min_allele_length - 1) == trimmed_alt_allele.seq().substr(1, min_allele_length - 1)) {

        if (trimmed_ref_allele.seq().size() > trimmed_alt_allele.seq().size()) {

            cur_var->alt(alt_idx).seq().append(cur_var->ref().seq().substr(cur_var->alt(alt_idx).seq().size(), string::npos));

        } else {

            assert(trimmed_ref_allele.seq().size() < trimmed_alt_allele.seq().size());
            cur_var->alt(alt_idx).seq().erase(cur_var->ref().seq().size(), string::npos);
        }

        return new Variant(cur_var->chrom(), cur_var->pos() + trimmed_left_bases + min_allele_length - 1, Allele(trimmed_ref_allele.seq().substr(min_allele_length - 1, string::npos)), vector<Allele>(1, Allele(trimmed_alt_allele.seq().substr(min_allele_length - 1, string::npos))), cur_var->info());
    
    } else {

        return nullptr;
    } 
}

int main(int argc, char const *argv[]) {

    if (argc != 3) {

        std::cout << "USAGE: simplifyDoublePolymorphisms <variants> <output_prefix>" << std::endl;
        return 1;
    }

    cout << "\n[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") simplifyDoublePolymorphisms script...\n" << endl;

    cout << "\n[" << Utils::getLocalTime() << "] WARNING: Not all types of mixed double polymorphisms are at the moment supported.\n" << endl;

    VcfFileReader vcf_reader(argv[1], true);

    vcf_reader.metaData().filterDescriptors().clear();
    vcf_reader.metaData().infoDescriptors().clear();
    vcf_reader.metaData().formatDescriptors().clear();

    VcfFileWriter output_vcf(string(argv[2]) + ".vcf", vcf_reader.metaData(), true);
    
    Variant * cur_var;
    map<uint, Variant *> variant_buffer;

    string cur_chrom = "";

    uint num_variants = 0;
    uint num_alt_alleles = 0;
    uint num_np_alt_alleles = 0;
    uint num_mp_alt_alleles = 0;

    while (vcf_reader.getNextVariant(&cur_var)) {

        if (cur_chrom != cur_var->chrom()) {

            auto variant_buffer_it = variant_buffer.begin();

            while (!(variant_buffer.empty())) {

                writeVariant(&output_vcf, variant_buffer_it->second);
                delete variant_buffer_it->second;

                variant_buffer_it = variant_buffer.erase(variant_buffer_it);
            }

            assert(variant_buffer_it == variant_buffer.end());
        }

        cur_chrom = cur_var->chrom();
        uint cur_pos = cur_var->pos();

        num_variants++;
        num_alt_alleles += cur_var->numAlts();

        Auxiliaries::rightTrimVariant(cur_var);

        if (cur_var->ref().seq().size() > 1) {

            for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

                if (cur_var->alt(alt_idx).seq().size() > 1) {

                    Allele ref_allele = cur_var->ref().seq();
                    Allele alt_allele = cur_var->alt(alt_idx).seq();
                    assert(!(alt_allele.isMissing()));

                    auto trimmed_bases = Auxiliaries::partialTrimAllelePair(&ref_allele, &alt_allele);

                    if ((ref_allele.seq().size() > 1) and (alt_allele.seq().size() > 1)) {
                        
                        assert(ref_allele.seq().front() != alt_allele.seq().front());
                        assert(ref_allele.seq().back() != alt_allele.seq().back());

                        if (ref_allele.seq().size() == alt_allele.seq().size()) {

                            auto new_var = simplifyNucleotidePolymorphism(cur_var, alt_idx, ref_allele, alt_allele, trimmed_bases.first);

                            if (new_var) {

                                num_np_alt_alleles++;                                
                                addVariantToBuffer(&variant_buffer, new_var);
                            }

                        } else {

                            cout << endl;
                            cout << ref_allele.seq() << " " << alt_allele.seq() << endl;

                            auto new_var = simplifyMixedPolymorphism(cur_var, alt_idx, ref_allele, alt_allele, trimmed_bases.first);

                            if (new_var) {

                                num_mp_alt_alleles++;
                                addVariantToBuffer(&variant_buffer, new_var);
                            }
                        } 
                    }
                }
            }
        }   

        addVariantToBuffer(&variant_buffer, cur_var);

        auto variant_buffer_it = variant_buffer.begin();

        while (!(variant_buffer.empty())) {

            if (variant_buffer_it->second->pos() > cur_pos) {

                break;
            }

            writeVariant(&output_vcf, variant_buffer_it->second);
            delete variant_buffer_it->second;

            variant_buffer_it = variant_buffer.erase(variant_buffer_it);
        }

        if ((num_variants % 100000) == 0) {

            std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
        }
    }

    auto variant_buffer_it = variant_buffer.begin();

    while (!(variant_buffer.empty())) {

        writeVariant(&output_vcf, variant_buffer_it->second);
        delete variant_buffer_it->second;

        variant_buffer_it = variant_buffer.erase(variant_buffer_it);
    }

    cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants and " << num_alt_alleles << " alternative alleles:" <<  endl;
    cout << "\n\t- Number of simplified double nucleotide polymorphisms (e.g. AGA/GGT -> A/G and A/T): " << num_np_alt_alleles << endl;
    cout << "\t- Number of simplified mixed double polymorphisms (e.g. AGAA/GGA -> A/G and G/GA): " << num_mp_alt_alleles << endl;
    cout << endl;

    return 0;
}
