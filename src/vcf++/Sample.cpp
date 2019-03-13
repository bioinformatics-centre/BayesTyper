
/*
Sample.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

    call_status_ = CallStatus::Complete;
    ploidy_ = Ploidy::Zeroploid;

    phased_ = false;

    num_alleles = 0;
    var_has_missing_allele = false;
}

Sample::Sample(const string & sample_str, const string & format_str, const map<string, Attribute::DetailedDescriptor> & format_dscrs, ushort num_alleles_in, bool var_has_missing_allele_in) {

    call_status_ = CallStatus::Complete;
    ploidy_ = Ploidy::Zeroploid;
    
    phased_ = false;

    num_alleles = num_alleles_in;
    var_has_missing_allele = var_has_missing_allele_in;

    auto sample_split_str = Utils::splitString(sample_str, ':');
    auto format_str_split =  Utils::splitString(format_str, ':');

    assert(sample_split_str.size() <= format_str_split.size());

    assert(!format_str_split.empty());
    assert(format_str_split.front() == "GT");

    if (sample_split_str.empty()) {

        assert(format_str_split.size() == 1);

    } else if (sample_split_str.front().size() > 0) {

        vector<string> gt_str_split;

        if (sample_split_str.front().find("/") != string::npos) {

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

        const uchar gt_ploidy = gt_str_split.size();
        
        assert(gt_ploidy > 0);
        assert(genotype_estimate.size() <= gt_ploidy);
        
        if (genotype_estimate.size() == 0) {

            call_status_ = CallStatus::Missing;

        } else if (genotype_estimate.size() < gt_ploidy) {

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
            assert(current_format_str_dscr.id() == format_str_split.at(sample_split_str_idx));            
            
            assert(current_format_str_dscr.type() != Attribute::Type::Flag);
            assert(current_format_str_dscr.number() != Attribute::Number::Zero);

            auto format_str_val_split = Utils::splitString(sample_split_str.at(sample_split_str_idx), ',');

            switch (current_format_str_dscr.number()) {

                case Attribute::Number::G : {

                    parseFormatAttributeVector(&genotype_info, current_format_str_dscr, format_str_val_split, 0);
                    break;
                }

                case Attribute::Number::R : {

                    parseFormatAttributeVector(&allele_info, current_format_str_dscr, format_str_val_split, 0);
                    break;
                }

                case Attribute::Number::A : {

                    parseFormatAttributeVector(&allele_info, current_format_str_dscr, format_str_val_split, 1);
                    break;
                }

                default : {

                    if (format_str_val_split.size() > 1) {

                        if (find_if_not(format_str_val_split.begin(), format_str_val_split.end(), [](string str){return str == ".";}) != format_str_val_split.end()) {

                            assert(info_.add(current_format_str_dscr.id(), sample_split_str.at(sample_split_str_idx), Attribute::Type::String));
                        }

                    } else {

                        if (sample_split_str.at(sample_split_str_idx) != ".") {

                            assert(info_.add(current_format_str_dscr.id(), sample_split_str.at(sample_split_str_idx), current_format_str_dscr.type()));
                        }
                    }
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

    assert(!format_descriptors.empty());
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

    JoiningString sample_element(':');
    sample_element.join(sample_gt_element.str());

    for (uint format_descriptor_idx = 1; format_descriptor_idx < format_descriptors.size(); format_descriptor_idx++) {

        auto current_format_descriptor = format_descriptors.at(format_descriptor_idx);

        assert(current_format_descriptor.type() != Attribute::Type::Flag);
        assert(current_format_descriptor.number() != Attribute::Number::Zero);

        switch (current_format_descriptor.number()) {

            case Attribute::Number::G : {

                sample_element.join(writeFormatAttributeVector(genotype_info, current_format_descriptor.id(), 0));
                break;
            }

            case Attribute::Number::R : {

                sample_element.join(writeFormatAttributeVector(allele_info, current_format_descriptor.id(), 0));
                break;
            }

            case Attribute::Number::A : {

                sample_element.join(writeFormatAttributeVector(allele_info, current_format_descriptor.id(), 1));
                break;
            }

            default : {

                auto info_value = info_.getValue(current_format_descriptor.id());

                if (info_value.second) {

                    sample_element.join(info_value.first.str());

                } else {

                    sample_element.join(".");
                }
            }
        }
    }

    if (sample_element.empty()) {

        return ".";
    
    } else {

        return sample_element.str();
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
    genotype_info.clear();

    if (ploidy_ != Ploidy::Zeroploid) {
        
        assert(ploidy_ != Ploidy::Polyploid);
        allele_info.emplace_back(AttributeSet());
    }

    if (allele_is_missing) {

        assert(!var_has_missing_allele);
        var_has_missing_allele = true;
    }
}

void Sample::removeAllele(const ushort removed_allele, const bool allele_is_missing, const bool convert_genotype_to_ref) {

    num_alleles--;
    genotype_info.clear();

    if (ploidy_ != Ploidy::Zeroploid) {
        
        assert(ploidy_ != Ploidy::Polyploid);
        allele_info.erase(allele_info.begin() + removed_allele);
    }

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

        } else if (convert_genotype_to_ref) {

            new_genotype_estimate.push_back(0);
        }
    }

    newGenotypeEstimate(new_genotype_estimate);
}

void Sample::clearGenotypeEstimate() {

    genotype_estimate.clear();
    call_status_ = CallStatus::Missing;
}

vector<ushort> Sample::oneToTwoDimIdx(const uint one_d_idx) const {

    uint idx_sum = 0;
    ushort second_idx = 0;

    while ((idx_sum + second_idx) < one_d_idx) {

        second_idx++;
        idx_sum += second_idx;
    }

    ushort first_idx = one_d_idx - idx_sum;

    assert(first_idx <= second_idx);
    return {first_idx, second_idx};
}

uint Sample::twoToOneDimIdx(const vector<ushort> & two_d_idx) const {

    assert(two_d_idx.size() == 2);
    assert(two_d_idx.front() <= two_d_idx.back());

    return ((two_d_idx.back() * (two_d_idx.back() + 1)) / 2 + two_d_idx.front());
}

void Sample::parseFormatAttributeVector(vector<AttributeSet> * attributes, const Attribute::DetailedDescriptor & attribute_dscr, const vector<string> & attribute_str, const uint first_attribute_idx) {

    assert(first_attribute_idx < attributes->size());
    assert(attribute_str.size() == (attributes->size() - first_attribute_idx));
    
    auto attribute_str_it = attribute_str.begin();

    for (uint attribute_idx = first_attribute_idx; attribute_idx < attributes->size(); attribute_idx++) {

        if (*attribute_str_it != ".") {

            assert(attributes->at(attribute_idx).add(attribute_dscr.id(), *attribute_str_it, attribute_dscr.type()));
        }

        attribute_str_it++;
    }

    assert(attribute_str_it == attribute_str.end());
}

string Sample::writeFormatAttributeVector(const vector<AttributeSet> & attributes, const string & id, const uint first_attribute_idx) {

    if (attributes.size() > 0) {

        assert(first_attribute_idx < attributes.size());

        bool found_format_id = false;

        JoiningString format_element_value(',');

        for (uint attribute_idx = first_attribute_idx; attribute_idx < attributes.size(); attribute_idx++) {

            auto format_value = attributes.at(attribute_idx).getValue(id);
            
            if (format_value.second) {
    
                found_format_id = true;
                format_element_value.join(format_value.first.str());

            } else {

                format_element_value.join(".");
            }
        }

        if (found_format_id) {

            return format_element_value.str();            

        } else {

            return ".";
        }

    } else {

        return ".";
    }
}

