
/*
AttributeFilter.cpp - This file is part of BayesTyper (v1.1)


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


#include "AttributeFilter.hpp"

AttributeFilter::AttributeFilter(const string & id_in) {

    id_ = id_in;
}

AttributeFilter::~AttributeFilter() {}


const string & AttributeFilter::id() const {

    return id_;
}

FlagAttributeFilter::FlagAttributeFilter(const string & id_in) : AttributeFilter(id_in) {}

FlagAttributeFilter::~FlagAttributeFilter(){};

bool FlagAttributeFilter::pass(const AttributeSet & att_set) {

    return att_set.hasValue(id());
}

ValueAttributeFilter::ValueAttributeFilter(const string & id_in, Attribute::Value value_in, Attribute::CmpOp * pass_operator_in) : AttributeFilter(id_in) {

    value_ = value_in;
    pass_operator = pass_operator_in;
}

ValueAttributeFilter::~ValueAttributeFilter(){}

const Attribute::Value & ValueAttributeFilter::value() const {

    return value_;
}

void ValueAttributeFilter::setValue(Attribute::Value value_in) {

    value_ = value_in;
}

bool ValueAttributeFilter::pass(const AttributeSet & att_set) {

    auto get_result = att_set.getValue(id_);

    if (get_result.second) {

        return (*pass_operator)(get_result.first, value_);

    } else {

        return false;
    }
}
