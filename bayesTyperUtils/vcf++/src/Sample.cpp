
/*
Sample.cpp - This file is part of BayesTyper (v0.9)


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


#include <regex>
#include <algorithm>

#include "Utils.hpp"
#include "Sample.hpp"
#include "JoiningString.hpp"

Sample::Sample() {

    call_status_ = CallStatus::Missing;
    ploidy_ = Ploidy::Zeroploid;
    phased_ = false;

    num_alleles = 0;

    var_has_missing_allele = false;
}

Sample::Sample(const string & sample_str, const string & format_str, const map<string, Attribute::DetailedDescriptor> & format_dscrs, ushort num_alleles_in, bool var_has_missing_allele_in) {

    num_alleles = num_alleles_in;
    phased_ = false;
    var_has_missing_allele = var_has_missing_allele_in;

    auto sample_split_str = Utils::splitString(sample_str, ':');
    auto format_str_split =  Utils::splitString(format_str, ':');

    assert(sample_split_str.size() <= format_str_split.size());
    
    assert(!(format_str_split.empty()));
    assert(format_str_split.front() == "GT");

    if (sample_split_str.empty()) {

        assert(format_str_split.size() == 1);

        ploidy_ = Ploidy::Zeroploid;
        call_status_ = CallStatus::Complete;

    } else if (sample_split_str.front().size() == 0) {

        ploidy_ = Ploidy::Zeroploid;
        call_status_ = CallStatus::Complete;

    } else {

        vector<string> gt_str_split;

        if (sample_split_str.front().find("/") != string::npos) {

            phased_ = false;
            gt_str_split = Utils::splitString(sample_split_str.front(), '/');

        } else {

            phased_ = true;
            gt_str_split = Utils::splitString(sample_split_str.front(), '|');
        }

        for (auto & all_idx : gt_str_split) {

            if (all_idx != ".") {

                genotype_estimate.push_back(stoi(all_idx));
                assert(genotype_estimate.back() < Utils::ushort_overflow);
            }
        }

        uchar gt_ploidy = gt_str_split.size();
        assert(gt_ploidy > 0);

        if (gt_ploidy == genotype_estimate.size()) {

            call_status_ = CallStatus::Complete;

        } else if (genotype_estimate.size() == 0) {

            call_status_ = CallStatus::Missing;

        } else {

            call_status_ = CallStatus::Partial;
        }

        switch (gt_ploidy) {

            case 1 : {

                ploidy_ = Ploidy::Haploid;

                assert((call_status_ == CallStatus::Complete) or (call_status_ == CallStatus::Missing));

                allele_info = vector<AttributeSet>(num_alleles);
                genotype_info = vector<AttributeSet>(num_alleles);

                break;
            }

            case 2 : {

                ploidy_ = Ploidy::Diploid;

                allele_info = vector<AttributeSet>(num_alleles);
                genotype_info = vector<AttributeSet>(numPossibleGenotypes());

                break;
            }

            default : {

                assert(gt_ploidy > 2);
                ploidy_ = Ploidy::Polyploid;
                cout << "WARNING: Ploidy greater than two detected ..." << endl;
            }
        }
    }

    for (uint sample_split_str_idx = 1; sample_split_str_idx < sample_split_str.size(); sample_split_str_idx++) {

        if (sample_split_str.at(sample_split_str_idx) == ".") {

            continue;
        }

        auto fmt_str_dscr_find_res = format_dscrs.find(format_str_split.at(sample_split_str_idx));

        if (fmt_str_dscr_find_res != format_dscrs.end()) {

            auto current_format_str_dscr = fmt_str_dscr_find_res->second;

            assert(current_format_str_dscr.type() != Attribute::Type::Flag);

            switch (current_format_str_dscr.number()) {

                case Attribute::Number::G : {

                    auto format_val_str_split = Utils::splitString(sample_split_str.at(sample_split_str_idx), ',');
                    assert(genotype_info.size() == format_val_str_split.size());

                    for (uint gt_idx = 0; gt_idx < genotype_info.size(); gt_idx++) {

                        genotype_info.at(gt_idx).add(current_format_str_dscr.id(), format_val_str_split.at(gt_idx), current_format_str_dscr.type());
                    }

                    break;
                }

                case Attribute::Number::R : {

                    auto format_val_str_split = Utils::splitString(sample_split_str.at(sample_split_str_idx), ',');
                    assert(allele_info.size() == format_val_str_split.size());

                    for (uint allele_idx = 0; allele_idx < allele_info.size(); allele_idx++) {

                        allele_info.at(allele_idx).add(current_format_str_dscr.id(), format_val_str_split.at(allele_idx), current_format_str_dscr.type());
                    }

                    break;
                }

                default : {

                    assert((current_format_str_dscr.number() == Attribute::Number::One) or (current_format_str_dscr.number() == Attribute::Number::Dot));
                    info_.add(current_format_str_dscr.id(), sample_split_str.at(sample_split_str_idx), current_format_str_dscr.type());
                }
            }
        }
    }
}

Sample::Ploidy Sample::ploidy() const {

    return ploidy_;
}

Sample::CallStatus Sample::callStatus() const {

    return call_status_;
}

const vector<ushort> & Sample::genotypeEstimate() const {

    // assert(!is_missing);
    return genotype_estimate;
}

bool Sample::isPhased() const {

    return phased_;
}

bool Sample::isInformative() const {

    if (call_status_ == CallStatus::Missing or ploidy_ == Ploidy::Zeroploid or ploidy_ == Ploidy::Polyploid) {

        return false;

    } else {

        return true;
    }
}

bool Sample::varHasMissingAllele() const {

    return var_has_missing_allele;
}

bool Sample::alleleIsMissing(ushort all_idx) const {

    if (varHasMissingAllele()) {

        return (all_idx + 1 == num_alleles);

    } else {

        return false;
    }
}

AttributeSet & Sample::info() {

    return info_;
}

vector<AttributeSet> & Sample::alleleInfo() {

    return allele_info;
}

vector<AttributeSet> & Sample::genotypeInfo() {

    return genotype_info;
}

string Sample::vcf(const vector<Attribute::DetailedDescriptor> & format_descriptors) {

    JoiningString sample_element(':');
    assert(format_descriptors.front().id() == "GT");

    JoiningString sample_gt_element;

    if (phased_) {

        sample_gt_element.setDelim('|');

    } else {

        sample_gt_element.setDelim('/');
    }

    for (auto & allele_idx : genotype_estimate) {

        sample_gt_element.join(to_string(allele_idx));
    }

    if ((call_status_ == CallStatus::Partial and ploidy_ == Ploidy::Diploid) or (call_status_ == CallStatus::Missing and ploidy_ == Ploidy::Haploid)) {

        sample_gt_element.join(".");

    } else if (call_status_ == CallStatus::Missing and ploidy_ == Ploidy::Diploid) {

        sample_gt_element.join(".");
        sample_gt_element.join(".");
    }

    sample_element.join(sample_gt_element.str());

    for (uint format_descriptor_idx = 1; format_descriptor_idx < format_descriptors.size(); format_descriptor_idx++) {

        auto current_format_descriptor = format_descriptors.at(format_descriptor_idx);

        switch (current_format_descriptor.number()) {

            case Attribute::Number::G : {

                sample_element.join(joinVectorAttribute(genotype_info, current_format_descriptor.id()));
                break;
            }

            case Attribute::Number::R : {

                sample_element.join(joinVectorAttribute(allele_info, current_format_descriptor.id()));
                break;
            }

            default : {

                assert(current_format_descriptor.number() == Attribute::Number::One);
                auto info_get = info_.getValue(current_format_descriptor.id());

                if (info_get.second) {

                    sample_element.join(info_get.first.str());

                } else {

                    sample_element.join(".");
                }
            }
        }
    }

    return sample_element.str();
}

string Sample::joinVectorAttribute(const vector<AttributeSet> & vector_attribute, const string & id) {

    if (vector_attribute.size() > 0) {

        auto first_get = vector_attribute.front().getValue(id);

        if (first_get.second) {

            JoiningString sample_gt_element(',');

            for (auto & att_set : vector_attribute) {

                auto att_get = att_set.getValue(id);
                assert(att_get.second);
                sample_gt_element.join(att_get.first.str());
            }

            return sample_gt_element.str();

        } else {

            return ".";
        }

    } else {

        return ".";
    }
}

uint Sample::numPossibleGenotypes() {

    switch (ploidy_) {

        case Ploidy::Haploid : {

            return num_alleles;
        }

        case Ploidy::Diploid : {

            return (num_alleles * num_alleles + num_alleles)/2;
        }

        default : {

            assert(false);
        }
    }
}

void Sample::newGenotypeEstimate(const vector<ushort> & new_genotype_estimate) {

    genotype_estimate = new_genotype_estimate;

    assert(ploidy_ != Ploidy::Polyploid);

    if (ploidy_ == Ploidy::Zeroploid) {

        assert(genotype_estimate.empty());
        assert(call_status_ == CallStatus::Complete);

    } else if (ploidy_ == Ploidy::Haploid) {

        if (genotype_estimate.empty()) {

            call_status_ = CallStatus::Missing;
        
        } else {

            assert(genotype_estimate.size() == 1);
            call_status_ = CallStatus::Complete;
        }

    } else {

        assert(ploidy_ == Ploidy::Diploid);

        if (genotype_estimate.empty()) {

            call_status_ = CallStatus::Missing;

        } else if (genotype_estimate.size() == 1) {

            call_status_ = CallStatus::Partial;
        
        } else {

            assert(genotype_estimate.size() == 2);
            call_status_ = CallStatus::Complete;
        }    
    }
}

void Sample::addAllele(const bool allele_is_missing) {

    num_alleles++;

    allele_info.clear();
    genotype_info.clear();   
    
    if (allele_is_missing) {

        assert(!var_has_missing_allele);
        var_has_missing_allele = true;
    } 
}

void Sample::removeAllele(const ushort removed_allele, const bool allele_is_missing) { 

    num_alleles--; 

    allele_info.clear();
    genotype_info.clear(); 

    if (allele_is_missing) {

        var_has_missing_allele = false;
    }

    vector<ushort> new_genotype_estimate;
    new_genotype_estimate.reserve(genotype_estimate.size());

    for (auto &allele: genotype_estimate) {

        if (allele > removed_allele) {

            new_genotype_estimate.push_back(allele - 1);

        } else if (allele < removed_allele) {

            new_genotype_estimate.push_back(allele);
        }
    }

    newGenotypeEstimate(new_genotype_estimate);
}

void Sample::clearGenotypeEstimate() {

    genotype_estimate.clear();
    call_status_ = CallStatus::Missing;
}

vector<ushort> Sample::oneToTwoDimIdx(uint linear_idx) const {

    ushort i = floor((2 * num_alleles + 1 - sqrt(pow(2 * num_alleles + 1, 2) - 8 * linear_idx)) / 2);
    ushort j = linear_idx - num_alleles * i + i * (i + 1) / 2;

    assert(i <= j);
    return {i, j};
}

uint Sample::twoToOneDimIdx(const vector<ushort> & two_d_idx) const {

    assert(two_d_idx.size() == 2);
    assert(two_d_idx.front() <= two_d_idx.back());

    return (two_d_idx.back() + num_alleles * two_d_idx.front() - two_d_idx.front() * (two_d_idx.front() + 1) / 2);
}
