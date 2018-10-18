
/*
Attribute.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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
#include <sstream>

#include "JoiningString.hpp"
#include "CompareOperators.hpp"

#include "Attribute.hpp"


vector<string> Attribute::Descriptor::type_labels = {"Float", "Integer", "String", "Character", "Flag"};
vector<string> Attribute::Descriptor::number_labels = {"0", "1", "2", "3", "4", "R", "A", "G", "."};

Attribute::Descriptor::Descriptor() {

    _id = "";
    _description = "";
}

Attribute::Descriptor::Descriptor(const string & id_in, const string & description_in) {

    _id = id_in;
    _description = description_in;
}

Attribute::Descriptor::Descriptor(const vector<pair<string,string> > & descriptor_elems) {

    for (auto & descriptor_elem : descriptor_elems) {

        if (descriptor_elem.first == "ID") {

            _id = descriptor_elem.second;

        } else if (descriptor_elem.first == "Description") {

            _description = descriptor_elem.second;
        }
    }

    assert(!_id.empty());
}

string Attribute::Descriptor::id() const {

    assert(!_id.empty());
    return _id;
}

string Attribute::Descriptor::description() const {

    return _description;
}


string Attribute::Descriptor::str() const {

    assert(!_id.empty());

    JoiningString descriptor_str(',');
    descriptor_str.join("ID=" + _id);

    descriptor_str.join("Description=\"" + _description + "\"");

    return "<" + descriptor_str.str() + ">";
}

bool Attribute::Descriptor::operator==(const Attribute::Descriptor & rhs) const {

    if ((this->_id == rhs._id) and (this->_description == rhs._description)) {

        return true;

    } else {

        return false;
    }
}

bool Attribute::Descriptor::operator!=(const Attribute::Descriptor & rhs) const {

    return !(*this == rhs);
}

Attribute::DetailedDescriptor::DetailedDescriptor() {

    has_number = false;
    has_type = false;
}

Attribute::DetailedDescriptor::DetailedDescriptor(const string & id, Number number, Type type, const string & description) {

    has_number = true;
    has_type = true;

    _id = id;
    _number = number;
    _type = type;
    _description = description;
}

Attribute::DetailedDescriptor::DetailedDescriptor(const vector<pair<string,string> > & descriptor_elems) : Descriptor(descriptor_elems) {

    has_number = false;
    has_type = false;

    for (auto & descriptor_elem : descriptor_elems) {

        if (descriptor_elem.first == "Number") {

            has_number = true;
            if (descriptor_elem.second == "0") {

                _number = Attribute::Number::Zero;

            } else if (descriptor_elem.second == "1") {

                _number = Attribute::Number::One;

            } else if (descriptor_elem.second == "2") {

                _number = Attribute::Number::Two;

            } else if (descriptor_elem.second == "3") {

                _number = Attribute::Number::Three;              

            } else if (descriptor_elem.second == "4") {

                _number = Attribute::Number::Four;

            } else if (descriptor_elem.second == "R") {

                _number = Attribute::Number::R;

            } else if (descriptor_elem.second == "A") {

                _number = Attribute::Number::A;

            } else if (descriptor_elem.second == "G") {

                _number = Attribute::Number::G;

            } else if (descriptor_elem.second == ".") {

                _number = Attribute::Number::Dot;

            } else {

                assert(false);
            }

        } else if (descriptor_elem.first == "Type") {

            has_type = true;
            if (descriptor_elem.second == "Float") {

                _type = Attribute::Type::Float;

            } else if (descriptor_elem.second == "Integer") {

                _type = Attribute::Type::Int;

            } else if (descriptor_elem.second == "String") {

                _type = Attribute::Type::String;

            } else if (descriptor_elem.second == "Character") {

                _type = Attribute::Type::Char;

            } else if (descriptor_elem.second == "Flag") {

                _type = Attribute::Type::Flag;

            } else {

                assert(false);
            }
        }
    }

    assert(has_number);
    assert(has_type);
}

Attribute::Number Attribute::DetailedDescriptor::number() const {

    assert(has_number);
    return _number;
}

Attribute::Type Attribute::DetailedDescriptor::type() const {

    assert(has_type);
    return _type;
}

string Attribute::DetailedDescriptor::str() const {

    assert(!_id.empty());
    assert(has_number);
    assert(has_type);

    JoiningString descriptor_str(',');
    descriptor_str.join("ID=" + _id);
    descriptor_str.join("Number=" + number_labels.at(char(_number)));
    descriptor_str.join("Type=" + type_labels.at(char(_type)));
    descriptor_str.join("Description=\"" + _description + "\"");

    return "<" + descriptor_str.str() + ">";
}

bool Attribute::DetailedDescriptor::operator==(const Attribute::DetailedDescriptor & rhs) const {

    const Attribute::Descriptor & base_rhs = rhs;

    if ((dynamic_cast<const Attribute::Descriptor &>(*this) == base_rhs) and (this->_number == rhs._number) and (this->_type == rhs._type)) {

        return true;

    } else {

        return false;
    }
}

bool Attribute::DetailedDescriptor::operator!=(const Attribute::DetailedDescriptor & rhs) const {

    return !(*this == rhs);
}

Attribute::Value::Value() {

    type_ = Type::Flag;
    data_ = string("");
}

Attribute::Value::Value(Type type_in, const string & data_str) {

    assert(!data_str.empty());

    type_ = type_in;

    switch (type_) {

        case Type::Float :

        data_ = stof(data_str);
        break;

        case Type::Int :

        data_ = stoi(data_str);
        break;

        case Type::String :

        data_ = data_str;
        break;

        case Type::Char :

        assert(data_str.size() == 1);
        data_ = char(data_str.at(0));
        break;

        default :

        assert(false);
    }
}

Attribute::Value::Value(float val) {

    type_ = Attribute::Type::Float;
    data_ = val;
}

Attribute::Value::Value(int val) {

    type_ = Attribute::Type::Int;
    data_ = val;
}

Attribute::Value::Value(const string & val) {

    assert(!val.empty());

    type_ = Attribute::Type::String;
    data_ = val;
}

Attribute::Value::Value(char val) {

    type_ = Attribute::Type::Char;
    data_ = val;
}

Attribute::Type Attribute::Value::type() const {

    return type_;
}

Attribute::Value::Data Attribute::Value::data() const {

    return data_;
}

string Attribute::Value::str() const {

    stringstream ss;
    ss << data_;
    return ss.str();
}

bool Attribute::operator==(const Attribute::Value & lhs, const Attribute::Value & rhs) {

    if (lhs.type() == rhs.type()) {

        Attribute::IsEqCmpOp isEq;
        return isEq(lhs, rhs);

    } else {

        return false;
    }
}


bool Attribute::operator!=(const Attribute::Value & lhs, const Attribute::Value & rhs) {

    return !(lhs == rhs);
}
