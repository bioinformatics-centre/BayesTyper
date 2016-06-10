
/*
Contig.cpp - This file is part of BayesTyper (v0.9)


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


#include <string>
#include <vector>
#include <algorithm>

#include "Contig.hpp"

Contig::Contig(const vector<string> & contig_info_line) {

    assert(contig_info_line.size() >= 2);

    _id = "";
    _length = 0;

    for (auto &contig_info: contig_info_line) {

        auto contig_info_split = Utils::splitString(contig_info, '=');
        assert(contig_info_split.size() == 2);

        if (contig_info_split.front() == "ID") {

            _id = contig_info_split.back();

        } else if (contig_info_split.front() == "length") {

            _length = stoi(contig_info_split.back());
        }
    }

    assert(!(_id.empty()));
    assert(_length > 0);

    string id_copy = _id; 
    transform(id_copy.begin(), id_copy.end(), id_copy.begin(), ::tolower);

    if (id_copy.size() > 2) {

        if (id_copy.substr(0,3) == "chr") {

            id_copy = id_copy.substr(3);
        }
    }

    assert(!(id_copy.empty()));

    bool is_pure_digit = true;

    for (auto &symbol: id_copy) {

        if (!(isdigit(symbol))) {

            is_pure_digit = false;
        }
    }

    if (is_pure_digit) {

        _type = Type::Autosomal;
    
    } else if (id_copy == "x") {

        _type = Type::ChrX;

    } else if (id_copy == "y") {

        _type = Type::ChrY;

    } else if ((id_copy == "m") or (id_copy == "mt")) {

        _type = Type::Mitochondrial;

    } else {

        _type = Type::Unknown;
    }
}

string Contig::id() const {

    return _id;
}

uint Contig::length() const {

    return _length;
}

Contig::Type Contig::type() const {

    return _type;
}

string Contig::typeStr() const {

    stringstream type_ss;
    type_ss << _type;
    return type_ss.str();
}

bool operator==(const Contig & lhs, const Contig & rhs) {

    return ((lhs.id() == rhs.id()) and (lhs.length() == rhs.length()) and (lhs.type() == rhs.type()));
}

bool operator!=(const Contig & lhs, const Contig & rhs) {

    return !(lhs == rhs);
}

ostream& operator<< (ostream & os, Contig::Type type) {

    switch (type) {

        case Contig::Type::Autosomal:

            return os << "Autosomal";

        case Contig::Type::ChrX:

            return os << "ChrX";

        case Contig::Type::ChrY:

            return os << "ChrY";

        case Contig::Type::Mitochondrial:

            return os << "Mitochondrial";

        case Contig::Type::Unknown:

            return os << "Unknown";

        default:

            assert(false);
    }

    return os;
}
