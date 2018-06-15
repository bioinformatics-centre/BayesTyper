
/*
AttributeSet.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __vcf__AttributeSet_hpp
#define __vcf__AttributeSet_hpp

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "assert.h"

#include "Attribute.hpp"

using namespace std;

class AttributeSet {

  public:

    bool add(const string &, const string &, Attribute::Type);
    template<typename DataType>
    bool setValue(const string &, DataType);

    bool addFlag(const string & key);

    bool rm(const string &);
    bool hasValue(const string &) const;
    pair<Attribute::Value,bool> getValue(const string &) const;

    template<typename DataType>
    pair<DataType,bool> getValue(const string & key) const;

  private:

    unordered_map<string, Attribute::Value> attributes;
};

#endif
