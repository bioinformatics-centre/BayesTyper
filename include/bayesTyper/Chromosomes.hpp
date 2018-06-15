
/*
Chromosomes.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__Chromosomes_hpp
#define __bayesTyper__Chromosomes_hpp

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Utils.hpp"


typedef vector<pair<string, string> >::const_iterator SequenceIterator;

class Chromosomes {

  public:

    Chromosomes(const string &, const bool);

    void addFasta(const string &, const bool);
    void addSequence(const pair<string, string> &, const bool);

    SequenceIterator cbegin() const;
    SequenceIterator cend() const;
    SequenceIterator find(const string &) const;

    bool isDecoy(const string &) const;

    uint getTotalCount() const;
    uint getDecoyCount() const;

    ulong getTotalLength() const;
    ulong getDecoyLength() const;

    Utils::ChromClass getClass(string) const; 

    void convertToUpper();    

  protected:

    void parseFasta(const string &, const bool);

  	vector<pair<string, string> > _chromosomes;

    unordered_map<string, uint> order;
    unordered_set<string> decoys;

    ulong total_length;
    ulong decoy_length;
};


#endif
