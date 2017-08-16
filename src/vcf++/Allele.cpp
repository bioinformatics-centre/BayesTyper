
/*
Allele.cpp - This file is part of BayesTyper (v1.1)


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


#include "Allele.hpp"

#include <iostream>
#include <algorithm>

Allele::Allele(const string & seq) {

  _seq = seq;
  assert(!(_seq.empty()));

  if (_seq.front() == '<') {

    assert(_seq.back() == '>');

    _is_id = true;
    _is_missing = false;

  } else if ((_seq.front() == '[') or (_seq.front() == ']') or (_seq.back() == '[') or (_seq.back() == ']')) {

    _is_id = true;
    _is_missing = false;

  } else if (_seq == "*") {

    _is_id = false;
    _is_missing = true;
  
  } else {

    _is_id = false;
    _is_missing = false;    

    transform(_seq.begin(), _seq.end(), _seq.begin(), ::toupper);
    assert(_seq.find_first_not_of("ACGTRYSWKMBDHVN") == string::npos);
  }
}


string & Allele::seq() {

  return _seq;
}

bool Allele::isID() const {

  return _is_id;
}

bool Allele::isMissing() const {

  return _is_missing;
}

AttributeSet & Allele::info() {

  return _info;
}

bool operator==(const Allele & lhs, const Allele & rhs) {

    return ((lhs._seq == rhs._seq) and (lhs._is_id == rhs._is_id) and (lhs._is_missing == rhs._is_missing));
}

bool operator!=(const Allele & lhs, const Allele & rhs) {

    return !(lhs == rhs);
}
