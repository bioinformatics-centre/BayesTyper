
/*
VcfMetaData.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef VCFMETADATA
#define VCFMETADATA

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "assert.h"

#include "Utils.hpp"
#include "Attribute.hpp"
#include "Contig.hpp"

using namespace std;
using namespace Utils;

class VcfMetaData {

  public:

    VcfMetaData();

    uint numSamples() const;
    uint numColumns() const;

    bool addLine(const string &);
    bool hasSampleId(const string &) const;
    void addSample(const string &);
    void rmSample(const string &);
    void clearSamples();

    pair<uint,bool> sampleIdToIdx(const string &) const;
    pair<string,bool> sampleIdxToId(const uint) const;
    vector<string> sampleIds() const;

    pair<Attribute::DetailedDescriptor, bool> getInfoDescriptor(const string &) const;
    pair<Attribute::Descriptor, bool> getFilterDescriptor(const string &) const;
    pair<Attribute::DetailedDescriptor, bool> getFormatDescriptor(const string &) const;

    const vector<Contig> & contigs();
    uint getContigIndex(const string &);
    const Contig & getContig(const string &);

    multimap<string,string> & miscMeta();

    map<string, Attribute::DetailedDescriptor> & infoDescriptors();
    map<string, Attribute::Descriptor> & filterDescriptors();
    map<string, Attribute::DetailedDescriptor> & formatDescriptors();

    void setFormat(const string &);

    string vcf() const;

  private:

    vector<string> splitStringSkipQuotes(const string &, char);
    string trimString(const string &, const string &);
    vector<pair<string,string> > getDescriptorTokens(const string &);
    vector<pair<string,string> > metaDataLineToTokens(const string &);
    
    string _format;

    vector<Contig> _contigs;
    unordered_map<string, uint> contig_id_to_idx;

    multimap<string,string> misc_meta;

    map<string, Attribute::Descriptor> filter_descriptors;
    map<string, Attribute::DetailedDescriptor> info_descriptors;
    map<string, Attribute::DetailedDescriptor> format_descriptors;

    map<uint, string> sample_idx_to_id;
    map<string, uint> sample_id_to_idx;

    uint num_samples;
    uint num_columns;

    bool header_added;
    bool add_sample_allowed;
};

#endif
