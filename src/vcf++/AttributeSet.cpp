
/*
AttributeSet.cpp - This file is part of BayesTyper (v1.1)


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


#include <iostream>
#include "AttributeSet.hpp"

bool AttributeSet::add(const string & key, const string & value_str, Attribute::Type type) {

    Attribute::Value ins_val(type, value_str);
    auto ins_res = attributes.insert(make_pair(key, ins_val));

    if (ins_res.second) {

        return true;

    } else {

        ins_res.first->second = ins_val;
        return false;
    }
}

template<typename DataType>
bool AttributeSet::setValue(const string & key, DataType value) {

    auto ins_res = attributes.insert(make_pair(key, Attribute::Value(value)));

    if (ins_res.second) {

        return true;

    } else {

        ins_res.first->second = value;
        return false;
    }
}

bool AttributeSet::addFlag(const string & key) {

    return attributes.insert(make_pair(key, Attribute::Value())).second;
}

bool AttributeSet::rm(const string & key) {

    return (attributes.erase(key) == 1);
}

bool AttributeSet::hasValue(const string & key) const {

    return attributes.find(key) != attributes.end();
}

pair<Attribute::Value,bool> AttributeSet::getValue(const string & key) const {

    auto find_res = attributes.find(key);

    if (find_res == attributes.end()) {

        return make_pair(Attribute::Value(), false);

    } else {

        return make_pair(find_res->second, true);
    }
}

template<typename DataType>
pair<DataType,bool> AttributeSet::getValue(const string & key) const {

    auto find_res = attributes.find(key);

    if (find_res == attributes.end()) {

        return make_pair(DataType(), false);

    } else {

        return make_pair(boost::get<DataType>(find_res->second.data()), true);
    }
}

template bool AttributeSet::setValue<float>(const string & key, float value);
template bool AttributeSet::setValue<int>(const string & key, int value);
template bool AttributeSet::setValue<string>(const string & key, string value);
template bool AttributeSet::setValue<char>(const string & key, char value);

template pair<float,bool> AttributeSet::getValue<float>(const string & key) const;
template pair<int,bool> AttributeSet::getValue<int>(const string & key) const;
template pair<string,bool> AttributeSet::getValue<string>(const string & key) const;
template pair<char,bool> AttributeSet::getValue<char>(const string & key) const;
