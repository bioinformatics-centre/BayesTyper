
/*
FastaReader.cpp - This file is part of BayesTyper (v0.9)


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


#include <assert.h>

#include "FastaReader.hpp"

FastaReader::FastaReader(const string & fasta_filename) {

    fasta_file.open(fasta_filename);
    assert(fasta_file.is_open());
    last_line_read = !getline(fasta_file, cur_line);
}

FastaReader::~FastaReader() {

    fasta_file.close();
}

bool FastaReader::getNextRecord(FastaRecord ** fasta_rec) {

    if (last_line_read) {

        return false;

    } else {

        assert(cur_line.at(0) == '>');
        *fasta_rec = new FastaRecord(cur_line.substr(1), 300000000);

        while (last_line_read = !getline(fasta_file, cur_line), !last_line_read) {

            assert(!cur_line.empty());
            if (cur_line.at(0) == '>') {

                assert(!(*fasta_rec)->seq().empty());
                break;

            } else {

                (*fasta_rec)->appendSeq(cur_line);
            }
        }

        (*fasta_rec)->shrinkSeqToFit();

        return true;
    }
}
