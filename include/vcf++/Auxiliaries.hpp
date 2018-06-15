
/*
Auxiliaries.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __vcf__Auxiliaries_hpp
#define __vcf__Auxiliaries_hpp

#include <iostream>

#include "Utils.hpp"
#include "Allele.hpp"
#include "Variant.hpp"

namespace Auxiliaries {

    enum class Type {Reference, SNP, Insertion, Deletion, Inversion, Complex, Missing};

    class AlleleAttributes {

        public:

            AlleleAttributes(const Type type_in, const uint length_in, const uint num_ambiguous_in, const int sv_length_in) : type(type_in), length(length_in), num_ambiguous(num_ambiguous_in), sv_length(sv_length_in) {}
            string typeStr() const;

        	const Type type;
        	const uint length;
            const uint num_ambiguous;
        	const int sv_length;
    };

    uint leftTrimAllelePair(Allele *, Allele *, const bool);
    uint rightTrimAllelePair(Allele *, Allele *);

    pair<uint, uint> partialTrimAllelePair(Allele *, Allele *);
    pair<uint, uint> fullTrimAllelePair(Allele *, Allele *);

    AlleleAttributes alleleAttributes(Allele &, Allele &);
    string variantType(Variant &);

    pair<string, bool> variantOrigins(Variant &);

    bool isInversion(Allele &, Allele &, const float, const uint);
    string reverseComplementSequence(const string &);

    uint rightTrimVariant(Variant *);

    bool hasMissing(Variant &);

    bool isAnnotated(Variant &);
    bool isAnnotated(Allele &);

    bool hasRepeat(Variant &);
    bool hasRepeat(Allele &);

    uint repeatLength(Allele &, const string &);

    bool hasAmbiguous(Variant &);
    bool hasAmbiguous(Allele &);

    bool isAlleleCalled(Allele &, const float);

    pair<float, bool> getMaxGenotypePosterior(Sample &);
    pair<int, bool> getMaxGenotypePosteriorIndex(Sample &);

    void resetFilters(Variant * variant);
    void updateAlleleStatsAndCallProb(Variant * variant);

    vector<uint> getCalledAlleleIdxsSorted(Variant &, const float);
    vector<uint> getCalledAlleleIdxsSorted(Sample &, const float);

    vector<uint> getNonZeroProbAlleleIdxsSorted(Variant &);
    vector<uint> getNonZeroProbAlleleIdxsSorted(Sample &);

    void removeNonRelevantFilterDescriptors(VcfMetaData *, const unordered_set<string> &);
    void removeNonRelevantInfoDescriptors(VcfMetaData *, const unordered_set<string> &);
    void removeNonRelevantFormatDescriptors(VcfMetaData *, const unordered_set<string> &);

    pair<uint, string> getHomopolymerInfo(uint, const string &);
    vector<bool> getHomopolymerAlleles(Variant &, const string &, const uint);

    vector<Contig> mergeContigs(const vector<Contig> &, const vector<Contig> &);
}

ostream& operator<< (ostream &, Auxiliaries::Type);
bool operator==(const Auxiliaries::AlleleAttributes &, const Auxiliaries::AlleleAttributes &);
bool operator!=(const Auxiliaries::AlleleAttributes &, const Auxiliaries::AlleleAttributes &);


#endif
