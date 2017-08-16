
/*
Variant.cpp - This file is part of BayesTyper (v1.1)


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


#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>

#include "Utils.hpp"
#include "Variant.hpp"
#include "JoiningString.hpp"

Variant::Variant() {}

Variant::Variant(const string & chrom_in, const uint pos_in, const Allele & ref_in, const vector<Allele> & alts_in, const AttributeSet & info_in) {

    _chrom = chrom_in;
    _pos = pos_in;
    _alleles.push_back(ref_in);

    for (auto & alt : alts_in) {

        _alleles.push_back(alt);
    }

    _qual = make_pair(0, false);
    _info = info_in;

    samples_parsed = false;
}

Variant::Variant(const vector<string> & var_line, VcfMetaData & vcf_meta_data) {

    assert(var_line.size() >= 8);

    _chrom = var_line.at(0);
    _pos = stoi(var_line.at(1));
    _ids = parseField(var_line.at(2));
    _alleles.push_back(Allele(var_line.at(3)));

    assert(!(_alleles.front().isID()));
    assert(!(_alleles.front().isMissing()));

    auto alt_strs = Utils::splitString(var_line.at(4), ',');

    for (uint i = 0; i < alt_strs.size(); i++) {

        assert(alt_strs.at(i).size() > 0);
        _alleles.push_back(Allele(alt_strs.at(i)));
    }

    assert(_alleles.size() > 1);
    assert(_alleles.size() < Utils::ushort_overflow);
    assert(!(var_line.at(5).empty()));

    if (var_line.at(5) == ".") {

        _qual = make_pair(0, false);

    } else {

        _qual = make_pair(stof(var_line.at(5)), true);
    }

    _filters = parseField(var_line.at(6));

    auto info_strs = Utils::splitString(var_line.at(7), ';');

    for (auto & info_str : info_strs) {

        auto info_str_split = Utils::splitString(info_str, '=');
        assert(info_str_split.size() <= 2);

        if (info_str_split.back() == ".") {

            continue;
        }

        auto info_descriptor_get = vcf_meta_data.getInfoDescriptor(info_str_split.front());

        if (info_descriptor_get.second) {

            auto info_descriptor = info_descriptor_get.first;
            assert(info_descriptor.number() != Attribute::Number::G);

            auto info_str_val_split = Utils::splitString(info_str_split.back(), ',');

            switch (info_descriptor.number()) {

                case Attribute::Number::R :

                    assert(info_descriptor.type() != Attribute::Type::Flag);
                    assert(info_str_split.size() == 2);

                    assert(info_str_val_split.size() == numAlls());

                    for (uint allele_idx = 0; allele_idx < numAlls(); allele_idx++) {

                        assert(allele(allele_idx).info().add(info_str_split.front(), info_str_val_split.at(allele_idx), info_descriptor.type()));
                    }

                    break;

                case Attribute::Number::A : {

                    assert(info_descriptor.type() != Attribute::Type::Flag);
                    assert(info_str_split.size() == 2);

                    assert(info_str_val_split.size() == numAlts());

                    for (uint alt_allele_idx = 0; alt_allele_idx < numAlts(); alt_allele_idx++) {

                        assert(alt(alt_allele_idx).info().add(info_str_split.front(), info_str_val_split.at(alt_allele_idx), info_descriptor.type()));
                    }

                    break;
                }

                case Attribute::Number::Zero : {

                    assert(info_descriptor.type() == Attribute::Type::Flag);

                    assert(info_str_split.size() == 1);
                    assert(info_str_val_split.size() == 1);

                    assert(_info.addFlag(info_str_split.front()));

                    break;
                }

                default : {

                    assert(info_descriptor.type() != Attribute::Type::Flag);
                    assert(info_str_split.size() == 2);

                    if (info_str_val_split.size() > 1) {

                        assert(_info.add(info_str_split.front(), info_str_split.back(), Attribute::Type::String));

                    } else {

                        assert(_info.add(info_str_split.front(), info_str_split.back(), info_descriptor.type()));
                    }
                }
            }
        }
    }

    if (var_line.size() > 8) {

        samples_parsed = false;

        assert(vcf_meta_data.numSamples() <= var_line.size() - 9);
        assert(sample_info.sample_ids.empty());

        for (uint var_line_idx = 9; var_line_idx < var_line.size(); var_line_idx++) {

            auto sample_id = vcf_meta_data.sampleIdxToId(var_line_idx - 9);

            if (sample_id.second) {

                sample_info.sample_ids.push_back(sample_id.first);
                sample_info.sample_strs.push_back(var_line.at(var_line_idx));
            }
        }

        assert(sample_info.sample_ids.size() == sample_info.sample_strs.size());

        sample_info.fmt_str = var_line.at(8);
        sample_info.fmt_dscrs = vcf_meta_data.formatDescriptors();

    } else {

        samples_parsed = true;
    }
}

bool Variant::isGenotyped() const {

    if (samples_parsed) {

        return !_samples.empty();

    } else {

        return !sample_info.sample_ids.empty();
    }
}

const string & Variant::chrom() const {

    return _chrom;
}

uint Variant::pos() const {

    return _pos;
}

const vector<string> & Variant::ids() const {

    return _ids;
}

Allele & Variant::ref() {

    return _alleles.front();
}

Allele & Variant::allele(uint allele_idx) {

    return _alleles.at(allele_idx);
}

Allele & Variant::alt(uint alt_idx) {

    return _alleles.at(alt_idx + 1);
}

uint Variant::numAlls() const {

    return _alleles.size();
}

uint Variant::numAlts() const {

    return _alleles.size() - 1;
}

pair<float, bool> Variant::qual() const {

    return _qual;
}

const vector<string> & Variant::filters() const {

    return _filters;
}

void Variant::setChrom(string chrom) {

    _chrom = chrom;
}

void Variant::setPos(uint pos) {

    _pos = pos;
}

void Variant::addId(const string & id) {

    addElementToField(id, &_ids);
}

void Variant::setIds(const vector<string> & ids) {

    _ids = ids;
}

void Variant::setRef(const Allele & ref_in) {

    ref() = ref_in;
}

void Variant::setAlt(const Allele & alt_in, uint alt_idx) {

    alt(alt_idx) = alt_in;
}

pair<Allele &, bool> Variant::addAlt(const Allele & new_alt) {

    parseSamplesIfNotParsed();

    auto new_alt_find = find(_alleles.begin() + 1, _alleles.end(), new_alt);

    if (new_alt_find == _alleles.end()) {

        _alleles.push_back(new_alt);

        for (auto &sample: _samples) {

            sample.second.addAllele(new_alt.isMissing());
        }

        return pair<Allele &, bool>(_alleles.back(), true);

    } else {

        return pair<Allele &, bool>(*new_alt_find, false);
    }
}

void Variant::removeAlts(vector<uint> alt_indices, const bool convert_genotype_to_ref) {

    parseSamplesIfNotParsed();

    sort(alt_indices.begin(), alt_indices.end());
    auto ait = alt_indices.rbegin();

    while (ait != alt_indices.rend()) {

        for (auto &sample: _samples) {

            sample.second.removeAllele(*ait + 1, _alleles.at(*ait + 1).isMissing(), convert_genotype_to_ref);
        }

        _alleles.erase(_alleles.begin() + *ait + 1);
        ait++;
    }

    assert(!(_alleles.empty()));
}

void Variant::removeRedundantAlts() {

    vector<uint> redundant_alt_indices;
    redundant_alt_indices.reserve(numAlts());

    for (uint first_alt_idx = 1; first_alt_idx < numAlts(); first_alt_idx++) {

        for (uint second_alt_idx = 0; second_alt_idx < first_alt_idx; second_alt_idx++) {

            if (alt(first_alt_idx) == alt(second_alt_idx)) {

                redundant_alt_indices.push_back(first_alt_idx);
                break;
            }

        }
    }

    removeAlts(redundant_alt_indices);
}

void Variant::setQual(const pair<float, bool> & qual) {

    _qual = qual;
}

void Variant::addFilter(const string & filter) {

    addElementToField(filter, &_filters);
}

void Variant::setFilters(const vector<string> & filters) {

    _filters = filters;
}

AttributeSet & Variant::info() {

    return _info;
}

vector<string> Variant::parseField(const string & field_str) {

    vector<string> field = Utils::splitString(field_str, ';');
    assert(!(field.empty()));

    auto field_it = field.begin();
    assert(field_it != field.end());

    while (field_it != field.end()) {

        if (*field_it == ".") {

            field_it = field.erase(field_it);

        } else {

            field_it++;
        }
    }

    return field;
}

void Variant::addElementToField(const string & element, vector<string> * field) {

    assert(!(element.empty()));
    assert(element != ".");

    if (find(field->begin(), field->end(), element) == field->end()) {

        field->push_back(element);
    }
}

string Variant::writeField(const vector<string> & field) {

    if (!(field.empty())) {

        JoiningString field_elements(';');

        unordered_set<string> elements;

        for (auto &element: field) {

            assert(!(element.empty()));
            assert(element != ".");

            if (elements.insert(element).second) {

                field_elements.join(element);
            }
        }

        return field_elements.str();

    } else {

        return ".";
    }
}

string Variant::vcf(VcfMetaData & vcf_meta_data) {

    JoiningString vcf_line('\t');
    vcf_line.join(_chrom);
    vcf_line.join(to_string(_pos));
    vcf_line.join(writeField(_ids));

    assert(!(ref().isID()));

    vcf_line.join(ref().seq());

    assert(numAlts() > 0);
    JoiningString alts_element(',');

    for (uint alt_idx = 0; alt_idx < numAlts(); alt_idx++) {

        alts_element.join(alt(alt_idx).seq());
    }

    vcf_line.join(alts_element.str());

    if (_qual.second) {

        vcf_line.join(to_string(_qual.first));

    } else {

        vcf_line.join(".");
    }

    sort(_filters.begin(), _filters.end());
    vcf_line.join(writeField(_filters));

    JoiningString info_element(';');

    for (auto & info_id_descriptor : vcf_meta_data.infoDescriptors()) {

        JoiningString info_element_value(',');
        bool found_info_id = false;

        auto info_descriptor = info_id_descriptor.second;
        assert(info_descriptor.number() != Attribute::Number::G);

        switch (info_descriptor.number()) {

            case Attribute::Number::R : {

                assert(info_descriptor.type() != Attribute::Type::Flag);

                for (uint allele_idx = 0; allele_idx < numAlls(); allele_idx++) {

                    auto info_value = allele(allele_idx).info().getValue(info_descriptor.id());
                    
                    if (info_value.second) {

                        found_info_id = true;
                        info_element_value.join(info_value.first.str());

                    } else {

                        info_element_value.join(".");                            
                    }
                }

                break;
            }

            case Attribute::Number::A : {

                assert(info_descriptor.type() != Attribute::Type::Flag);

                for (uint alt_allele_idx = 0; alt_allele_idx < numAlts(); alt_allele_idx++) {

                    auto info_value = alt(alt_allele_idx).info().getValue(info_descriptor.id());

                    if (info_value.second) {

                        found_info_id = true;
                        info_element_value.join(info_value.first.str());

                    } else {

                        info_element_value.join(".");                            
                    }
                }

                break;
            }

            case Attribute::Number::Zero : {

                assert(info_descriptor.type() == Attribute::Type::Flag);

                auto info_value = _info.getValue(info_descriptor.id());

                if (info_value.second) {

                    found_info_id = true;
                    assert(info_value.first.str().empty());
                }

                break;
            }

            default : {

                assert(info_descriptor.type() != Attribute::Type::Flag);

                auto info_value = _info.getValue(info_descriptor.id());

                if (info_value.second) {

                    found_info_id = true;
                    info_element_value.join(info_value.first.str());
                }
            }
        }

        if (found_info_id) {

            if (info_element_value.empty()) {

                info_element.join(info_descriptor.id());

            } else {

                auto info_element_key_value = info_descriptor.id() + "=" + info_element_value.str();
                info_element.join(info_element_key_value);
            }
        }
    }

    if (info_element.empty()) {

        vcf_line.join(".");

    } else {

        vcf_line.join(info_element.str());
    }

    if (isGenotyped()) {

        if (!samples_parsed) {

            // If meta data format descriptors has changed, samples need to be parsed for consistency
            if (vcf_meta_data.formatDescriptors() != sample_info.fmt_dscrs) {

                parseSamplesIfNotParsed();
            }
        }

        if (samples_parsed) {

            // Assemble format str
            set<string> fmt_str_ids;
            for (auto & sample : _samples) {

                for (auto & fmt_dscr : vcf_meta_data.formatDescriptors()) {

                    if (sample.second.info().hasValue(fmt_dscr.first)) {

                        fmt_str_ids.insert(fmt_dscr.first);
                    }

                    if (!(sample.second.alleleInfo().empty())) {

                        assert(!(sample.second.genotypeInfo().empty()));

                        if (sample.second.alleleInfo().front().hasValue(fmt_dscr.first) or sample.second.genotypeInfo().front().hasValue(fmt_dscr.first)) {

                            fmt_str_ids.insert(fmt_dscr.first);
                        }

                    } else {

                        assert(sample.second.genotypeInfo().empty());
                    }
                }

                if (fmt_str_ids.size() + 1 == vcf_meta_data.formatDescriptors().size()) {

                    break;
                }
            }

            assert(fmt_str_ids.find("GT") == fmt_str_ids.end());

            JoiningString fmt_str(':');
            fmt_str.join("GT");

            vector<Attribute::DetailedDescriptor> cur_fmt_dscrs;
            cur_fmt_dscrs.reserve(fmt_str_ids.size());

            auto fmt_dscr_get_gt = vcf_meta_data.getFormatDescriptor("GT");
            assert(fmt_dscr_get_gt.second);
            cur_fmt_dscrs.push_back(fmt_dscr_get_gt.first);

            for (auto & fmt_str_id : fmt_str_ids) {

                fmt_str.join(fmt_str_id);
                auto fmt_dscr_get = vcf_meta_data.getFormatDescriptor(fmt_str_id);
                assert(fmt_dscr_get.second);
                cur_fmt_dscrs.push_back(fmt_dscr_get.first);
            }

            vcf_line.join(fmt_str.str());

            for (auto & sample_id : vcf_meta_data.sampleIds()) {

                vcf_line.join(getSample(sample_id).vcf(cur_fmt_dscrs));
            }

        } else {

            vcf_line.join(sample_info.fmt_str);

            for (auto & sample_id : vcf_meta_data.sampleIds()) {

                auto sample_idx = vcf_meta_data.sampleIdToIdx(sample_id);
                assert(sample_idx.second);
                assert(sample_info.sample_ids.at(sample_idx.first) == sample_id);

                vcf_line.join(sample_info.sample_strs.at(sample_idx.first));
            }
        }
    }

    return vcf_line.str();
}

void Variant::parseSamplesIfNotParsed() {

    if (!samples_parsed) {

        for (uint sample_idx = 0; sample_idx < sample_info.sample_ids.size(); sample_idx++) {

            assert(_samples.emplace(sample_info.sample_ids.at(sample_idx), Sample(sample_info.sample_strs.at(sample_idx), sample_info.fmt_str, sample_info.fmt_dscrs, numAlls(), _alleles.back().isMissing())).second);
        }

        samples_parsed = true;

        sample_info.sample_ids.clear();
        sample_info.sample_strs.clear();
    }
}

uint Variant::numSamples() {

    if (samples_parsed) {

        assert(sample_info.sample_ids.empty());
        assert(sample_info.sample_strs.empty());

        return _samples.size();

    } else {

        assert(_samples.empty());
        return sample_info.sample_ids.size();
    }
}

vector<string> Variant::sampleIds() {

    if (samples_parsed) {

        vector<string> sample_ids;
        sample_ids.reserve(numSamples());

        for (auto & sample : _samples) {

            sample_ids.push_back(sample.first);
        }

        return sample_ids;

    } else {

        return sample_info.sample_ids;
    }
}

Sample & Variant::getSample(const string & sample_id) {

    parseSamplesIfNotParsed();
    auto smpl_find = _samples.find(sample_id);
    assert(smpl_find != _samples.end());
    return smpl_find->second;
}

bool Variant::addSample(const string & sample_id, const Sample & sample) {

    parseSamplesIfNotParsed();
    return(_samples.emplace(sample_id, sample).second);
}

bool Variant::hasSample(const string & sample_id) {

    parseSamplesIfNotParsed();
    return (_samples.find(sample_id) != _samples.end());
}

void Variant::clearSamples() {

    if (samples_parsed) {

        _samples.clear();
        assert(sample_info.sample_ids.empty());
        assert(sample_info.sample_strs.empty());

    } else {

        assert(_samples.empty());
        sample_info.sample_ids.clear();
        sample_info.sample_strs.clear();
    }
}
