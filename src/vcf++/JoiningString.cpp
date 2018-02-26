
/*
JoiningString.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "JoiningString.hpp"

JoiningString::JoiningString() {

    add_delim = false;
}

JoiningString::JoiningString(char delim) {

    add_delim = false;
    delim_ = delim;
    str_.reserve(1000);
}

void JoiningString::join(const string & join_str) {

    if (add_delim) {

        str_ += delim_;

    } else {

        add_delim = true;
    }

    str_ += join_str;
}

void JoiningString::join(const vector<string> & join_strs) {

    for (auto & join_str : join_strs) {

        join(join_str);
    }
}

const string & JoiningString::str() const {

    return str_;
}

bool JoiningString::empty() const {

    return str_.empty();
}

char JoiningString::delim() const {

    return delim_;
}

void JoiningString::setDelim(char delim) {

    delim_ = delim;
}
