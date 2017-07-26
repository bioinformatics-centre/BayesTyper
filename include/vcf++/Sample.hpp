
/*
Sample.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef SAMPLE
#define SAMPLE

#include <unordered_set>
#include <vector>
#include <map>

#include "Utils.hpp"
#include "Attribute.hpp"
#include "AttributeSet.hpp"

class Sample {

public:

    enum class CallStatus {Missing, Partial, Complete};
    enum class Ploidy {Zeroploid, Haploid, Diploid, Polyploid};

    Sample();
    Sample(const string &, const string &, const map<string, Attribute::DetailedDescriptor> &, ushort, bool);

    vector<ushort> oneToTwoDimIdx(uint) const;
    uint twoToOneDimIdx(const vector<ushort> &) const;
    uint numPossibleGenotypes();

    Ploidy ploidy() const;
    CallStatus callStatus() const;
    const vector<ushort> & genotypeEstimate() const;

    AttributeSet & info();
    vector<AttributeSet> & alleleInfo();
    vector<AttributeSet> & genotypeInfo();

    void addAllele(const bool);
    void removeAllele(const ushort, const bool, const bool);

    void newGenotypeEstimate(const vector<ushort> &);
    void clearGenotypeEstimate();

    bool isPhased() const;
    bool isInformative() const;
    bool varHasMissingAllele() const;
    bool alleleIsMissing(ushort) const;

    string vcf(const vector<Attribute::DetailedDescriptor> &);

private:

    string joinVectorAttribute(const vector<AttributeSet> &, const string &, const bool);

    CallStatus call_status_;
    Ploidy ploidy_;
    bool phased_;

    ushort num_alleles;

    vector<ushort> genotype_estimate;

    AttributeSet info_;
    vector<AttributeSet> allele_info;
    vector<AttributeSet> genotype_info;

    bool var_has_missing_allele;
};

#endif
