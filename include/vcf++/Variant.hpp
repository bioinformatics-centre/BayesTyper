
/*
Variant.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef VARIANT
#define VARIANT

#include <string>
#include <vector>
#include <unordered_map>
#include "assert.h"

#include "Allele.hpp"
#include "Sample.hpp"
#include "VcfMetaData.hpp"
#include "AttributeSet.hpp"

using namespace std;

typedef unsigned int uint;

class Variant {

public:

    Variant();
    Variant(const string &, const uint, const Allele &, const vector<Allele> &, const AttributeSet &);
    Variant(const vector<string> &, VcfMetaData &);

    const string & chrom() const;
    uint pos() const;
    const vector<string> & ids() const;
    Allele & ref();
    Allele & alt(uint);
    Allele & allele(uint);
    uint numAlls() const;
    uint numAlts() const;
    pair<float, bool> qual() const;
    const vector<string> & filters() const;

    void setChrom(string);
    void setPos(uint);

    void addId(const string &);
    void setIds(const vector<string> &);

    void setRef(const Allele &);

    void setAlt(const Allele &, uint);
    pair<Allele &, bool> addAlt(const Allele &);
    void removeAlts(vector<uint>, const bool convert_genotype_to_ref = false);
    void removeRedundantAlts();

    void setQual(const pair<float, bool> &);

    void addFilter(const string &);
    void setFilters(const vector<string> &);

    AttributeSet & info();

    vector<string> sampleIds();
    bool hasSample(const string &);
    Sample & getSample(const string &);
    bool addSample(const string &, const Sample &);
    void clearSamples();

    bool isGenotyped() const;
    uint numSamples();

    string vcf(VcfMetaData &);

protected:

    struct DelayedSampleParsingInfo {

        vector<string> sample_ids;
        vector<string> sample_strs;
        string fmt_str;
        map<string, Attribute::DetailedDescriptor> fmt_dscrs;
    };

    void parseSamplesIfNotParsed();

    string _chrom;
    uint _pos;
    vector<string> _ids;
    vector<Allele> _alleles;
    pair<float, bool> _qual;
    vector<string> _filters;

    vector<string> parseField(const string &);
    void addElementToField(const string &, vector<string> * field);
    string writeField(const vector<string> &);

    AttributeSet _info;

    bool is_genotyped;
    bool samples_parsed;
    DelayedSampleParsingInfo sample_info;
    unordered_map<string, Sample> _samples;
};

#endif
