
/*
FastaRecord.cpp - This file is part of BayesTyper (v0.9)


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


#include <random>
#include <algorithm>
#include <string>

#include "Utils.hpp"

#include "FastaRecord.hpp"

FastaRecord::FastaRecord(const string & id_in, uint expected_seq_len) {

    id_ = Utils::splitString(id_in, ' ').front();
    seq_.reserve(expected_seq_len);
}

FastaRecord::FastaRecord(const string & id_in, const string & seq_in) {

    id_ = Utils::splitString(id_in, ' ').front();
    seq_ = seq_in;
}

const string & FastaRecord::id() {

    return id_;
}

const string & FastaRecord::seq() {

    return seq_;
}

void FastaRecord::appendSeq(const string & seq_in) {

    seq_ += seq_in;
}

void FastaRecord::shrinkSeqToFit() {

    seq_.shrink_to_fit();
}

string FastaRecord::stripUCSCId(const string & id) const {

    if (id.substr(0,3) == "chr") {

        return id.substr(3);

    } else {

        return id;
    }
}

string FastaRecord::str() const {

    return ">" + id_ + "\n" + seq_ + "\n";
}

void FastaRecord::sampleAmbiguousBases(const vector<char> & canonical_bases, mt19937 * prng) {

    uniform_int_distribution<int> canonical_bases_sampler(0, canonical_bases.size() - 1);

    auto seq_it = seq_.begin();

    while (seq_it != seq_.end()) {

        if ((*seq_it == 'N') or (*seq_it == 'n')) {

            *seq_it = canonical_bases.at(canonical_bases_sampler(*prng));
        }

        seq_it++;
    }
}

void FastaRecord::convertToUppercase() {

    std::transform(seq_.begin(), seq_.end(), seq_.begin(), ::toupper);
}
