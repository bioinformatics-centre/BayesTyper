
/*
VcfMetaData.cpp - This file is part of BayesTyper (v1.1)


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

#include "Utils.hpp"
#include "VcfMetaData.hpp"
#include "JoiningString.hpp"
#include "Contig.hpp"

VcfMetaData::VcfMetaData() {

    header_added = false;
    add_sample_allowed = true;
    num_samples = 0;
    num_columns = 0;
}

uint VcfMetaData::numSamples() const {

    assert(header_added);
    return num_samples;
}

uint VcfMetaData::numColumns() const {

    assert(header_added);
    return num_columns;
}

vector<string> VcfMetaData::splitStringSkipQuotes(const string & str, char delim) {

    vector<string> elems;
    elems.reserve(100);

    string buf;
    bool in_quote = false;

    for (uint i = 0; i < str.size(); i++) {

        if (str.at(i) == '\"') {

            if (in_quote) {

                in_quote = false;

            } else {

                in_quote = true;
            }
        }

        if (str.at(i) == delim and !buf.empty() and !in_quote) {

            elems.push_back(buf);
            buf.clear();

        } else {

            buf += str.at(i);
        }
    }

    if (!buf.empty()) {

        elems.push_back(buf);
    }

    assert(!in_quote);

    return elems;
}

string VcfMetaData::trimString(const string & str, const string & trim_str) {

    auto start = str.find_first_not_of(trim_str);
    auto end = str.find_last_not_of(trim_str);

    if (start == string::npos and end == string::npos) {

        return "";

    } else {

        return str.substr(start, end - start + 1);
    }
}

vector<pair<string,string> > VcfMetaData::metaDataLineToTokens(const string & meta_data_line) {

    auto meta_data_line_trim_split = splitStringSkipQuotes(trimString(meta_data_line, "><"), ',');

    vector<pair<string,string> > tokens;

    for (auto & token : meta_data_line_trim_split) {

        auto token_split = splitStringSkipQuotes(token, '=');

        token_split.back() = trimString(token_split.back(), "\"");

        assert(token_split.size() == 2);
        tokens.emplace_back(token_split.front(), token_split.back());
    }

    return tokens;
}

bool VcfMetaData::addLine(const string & meta_data_line) {

    if (meta_data_line.substr(0,12) == "##fileformat") {

        auto filefmt_str_split = Utils::splitString(meta_data_line, '=');
        assert(filefmt_str_split.size() == 2);

        _format = filefmt_str_split.back();

    } else if (meta_data_line.substr(0,6) == "##INFO") {

        Attribute::DetailedDescriptor info_descriptor(metaDataLineToTokens(meta_data_line.substr(7)));
        assert(info_descriptors.emplace(info_descriptor.id(), info_descriptor).second);

    } else if (meta_data_line.substr(0,8) == "##FILTER") {

        Attribute::Descriptor filter_descriptor(metaDataLineToTokens(meta_data_line.substr(9)));
        assert(filter_descriptors.emplace(filter_descriptor.id(), filter_descriptor).second);

    } else if (meta_data_line.substr(0,8) == "##FORMAT") {

        Attribute::DetailedDescriptor format_descriptor(metaDataLineToTokens(meta_data_line.substr(9)));
        assert(format_descriptors.emplace(format_descriptor.id(), format_descriptor).second);

    } else if (meta_data_line.substr(0,6) == "#CHROM") {

        assert(sample_id_to_idx.size() == 0);

        header_added = true;
        auto header_split = Utils::splitString(meta_data_line, '\t');
        assert(header_split.size() >= 8);

        num_columns = header_split.size();

        for (uint header_idx = 9; header_idx < header_split.size(); header_idx++) {

            addSample(header_split.at(header_idx));
        }

    } else if (meta_data_line.substr(0,8) == "##contig") {

        auto contig_line_trim_split = splitStringSkipQuotes(trimString(meta_data_line.substr(9), "><"), ',');

        _contigs.emplace_back(contig_line_trim_split);
        assert(contig_id_to_idx.emplace(_contigs.back().id(), _contigs.size() - 1).second);

    } else {

        auto misc_meta_str_split = Utils::splitString(meta_data_line, '=', 2);
        assert(misc_meta_str_split.size() == 2);
        misc_meta.emplace(misc_meta_str_split.front().substr(2),misc_meta_str_split.back());
    }

    return header_added;
}

void VcfMetaData::clearSamples() {

    num_samples = 0;
    sample_id_to_idx.clear();
    sample_idx_to_id.clear();

    add_sample_allowed = true;
}

void VcfMetaData::addSample(const string & sample_id) {

    assert(add_sample_allowed);

    num_samples++;

    assert(sample_id_to_idx.emplace(sample_id, num_samples - 1).second);
    sample_idx_to_id.emplace(num_samples - 1, sample_id);

    assert(sample_id_to_idx.size() == num_samples);
    assert(sample_idx_to_id.size() == num_samples);
}

void VcfMetaData::rmSample(const string & sample_id) {

    auto sample_idx = sampleIdToIdx(sample_id);

    assert(sample_idx.second);
    assert(sample_id_to_idx.erase(sample_id) == 1);
    assert(sample_idx_to_id.erase(sample_idx.first) == 1);

    num_samples--;

    assert(sample_id_to_idx.size() == num_samples);
    assert(sample_idx_to_id.size() == num_samples);

    add_sample_allowed = false;
}

bool VcfMetaData::hasSampleId(const string & sample_id) const {

    return (sample_id_to_idx.find(sample_id) != sample_id_to_idx.end());
}

vector<string> VcfMetaData::sampleIds() const {

    auto sample_id_iter = sample_id_to_idx.begin();

    vector<string> sample_ids(numSamples());
    for (uint sample_idx = 0; sample_idx < numSamples(); sample_idx++) {

        sample_ids.at(sample_idx) = sample_id_iter->first;
        sample_id_iter++;
    }

    assert(sample_id_iter == sample_id_to_idx.end());

    return sample_ids;
}

pair<string,bool> VcfMetaData::sampleIdxToId(const uint sample_idx) const {

    auto sample_id_to_idx_find_res = sample_idx_to_id.find(sample_idx);

    if (sample_id_to_idx_find_res == sample_idx_to_id.end()) {

        return make_pair("", false);

    } else {

        return make_pair(sample_id_to_idx_find_res->second, true);
    }
}

pair<uint,bool> VcfMetaData::sampleIdToIdx(const string & sampleId) const {

    auto sample_id_to_idx_find_res = sample_id_to_idx.find(sampleId);

    if (sample_id_to_idx_find_res == sample_id_to_idx.end()) {

        return make_pair(0, false);

    } else {

        return make_pair(sample_id_to_idx_find_res->second, true);
    }
}

pair<Attribute::DetailedDescriptor,bool> VcfMetaData::getInfoDescriptor(const string & id) const {

    auto descriptor_find_res = info_descriptors.find(id);

    if (descriptor_find_res == info_descriptors.end()) {

        return make_pair(Attribute::DetailedDescriptor(), false);

    } else {

        return make_pair(descriptor_find_res->second, true);
    }
}

pair<Attribute::Descriptor,bool> VcfMetaData::getFilterDescriptor(const string & id) const {

    auto descriptor_find_res = filter_descriptors.find(id);
    if (descriptor_find_res == filter_descriptors.end()) {

        return make_pair(Attribute::Descriptor(), false);

    } else {

        return make_pair(descriptor_find_res->second, true);
    }
}

pair<Attribute::DetailedDescriptor,bool> VcfMetaData::getFormatDescriptor(const string & id) const {

    auto descriptor_find_res = format_descriptors.find(id);
    if (descriptor_find_res == format_descriptors.end()) {

        return make_pair(Attribute::DetailedDescriptor(), false);

    } else {

        return make_pair(descriptor_find_res->second, true);
    }
}

const vector<Contig> & VcfMetaData::contigs() {

    return _contigs;
}

uint VcfMetaData::getContigIndex(const string & id) {

    auto contig_id_to_idx_it = contig_id_to_idx.find(id);
    assert(contig_id_to_idx_it != contig_id_to_idx.end());
    assert(id == _contigs.at(contig_id_to_idx_it->second).id());

    return contig_id_to_idx_it->second;
}

const Contig & VcfMetaData::getContig(const string & id) {

    return _contigs.at(getContigIndex(id));
}

multimap<string,string> & VcfMetaData::miscMeta() {

    return misc_meta;
}

map<string, Attribute::DetailedDescriptor> & VcfMetaData::infoDescriptors() {

    return info_descriptors;
}

map<string, Attribute::Descriptor> & VcfMetaData::filterDescriptors() {

    return filter_descriptors;
}

map<string, Attribute::DetailedDescriptor> & VcfMetaData::formatDescriptors() {

    return format_descriptors;
}

void VcfMetaData::setFormat(const string & format_in) {

    _format = format_in;
}

string VcfMetaData::vcf() const {

    JoiningString meta_str('\n');

    meta_str.join("##fileformat=" + _format);

    for (auto & misc_meta_line : misc_meta) {

        meta_str.join("##" + misc_meta_line.first + "=" + misc_meta_line.second);
    }

    for (auto & contig : _contigs) {

        meta_str.join("##contig=<ID=" + contig.id() + ",length=" + to_string(contig.length()) + ">");
    }

    for (auto & filter_descriptor : filter_descriptors) {

        meta_str.join("##FILTER=" + filter_descriptor.second.str());
    }

    for (auto & info_descriptor : info_descriptors) {

        meta_str.join("##INFO=" + info_descriptor.second.str());
    }

    if (numSamples() > 0) {

        for (auto & format_descriptor : format_descriptors) {

            meta_str.join("##FORMAT=" + format_descriptor.second.str());
        }
    }

    JoiningString header_str('\t');
    header_str.join("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

    if (numSamples() > 0) {

        header_str.join("FORMAT");
        for (auto & sample_id : sampleIds()) {

            header_str.join(sample_id);
        }
    }

    meta_str.join(header_str.str());
    return meta_str.str();
}
