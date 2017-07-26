
/*
Trio.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef TRIO
#define TRIO

#include <vector>
#include <iostream>
#include <string>

#include "Utils.hpp"
#include "Variant.hpp"
#include "Sample.hpp"
#include "VcfMetaData.hpp"

class Trio {

    public:

        struct TrioInfo {

            string id;

            string father;
            string mother;
            string child;
        };

        Trio(Variant &, const TrioInfo &);

        Sample father;
        Sample mother;
        Sample child;

        bool is_filtered;
        bool is_diploid;
        bool is_informative;
        bool is_concordant;
        bool is_exclusively_child_heterozygote;
        bool is_parents_bi_allelelic_heterozygote;
        bool is_reference_call;
        bool has_called_missing;

        bool isFiltered();
        bool isDiploid();
        bool isInformative();
        bool isConcordant();
        bool isExclusivelyChildHeterozygote();
        bool isParentsBiAllelicHeterozygote();
        bool isReferenceCall();
        bool hasCalledMissing();

        struct DeNovoEvent {

            string child_id;

            uint ancestral_allele_idx;
            uint de_novo_allele_idx;
        };

        pair<bool, DeNovoEvent> de_novo_event;

        bool isDeNovo();
        DeNovoEvent deNovoEvent();

        static vector<TrioInfo> parsePedigree(const VcfMetaData &, const string &);
        static vector<TrioInfo> parseGenomeDKPedigree(const VcfMetaData &);
};


#endif
