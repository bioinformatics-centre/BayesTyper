
/*
Attribute.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef ATTRIBUTE
#define ATTRIBUTE

#include <vector>
#include <string>
#include "assert.h"

#include "boost/variant.hpp"

using namespace std;

namespace Attribute {

    enum class Type {Float, Int, String, Char, Flag};
    enum class Number {Zero, One, Two, Three, Four, R, A, G, Dot};

    class Descriptor {

    public:

        Descriptor();
        Descriptor(const string &, const string &);
        Descriptor(const vector<pair<string,string> > &);
        string id() const;
        string description() const;
        virtual string str() const;

        bool operator==(const Descriptor &) const;
        bool operator!=(const Descriptor &) const;

    protected:

        string _id;
        string _description;

        static vector<string> type_labels;
        static vector<string> number_labels;
    };

    class DetailedDescriptor : public Descriptor {

      public:

        DetailedDescriptor();
        DetailedDescriptor(const vector<pair<string,string> > &);
        DetailedDescriptor(const string &, Number, Type, const string &);

        Number number() const;
        Type type() const;
        string str() const;

        bool operator==(const DetailedDescriptor &) const;
        bool operator!=(const DetailedDescriptor &) const;

      private:

        Number _number;
        Type _type;

        bool has_number;
        bool has_type;
    };

    class Value {

        public:

            typedef boost::variant<float, int, std::string, char> Data;

            Value();
            Value(Type, const string &);
            Value(float);
            Value(int);
            Value(const string &);
            Value(char);

            Type type() const;
            Data data() const;
            string str() const;

        private:

            Type type_;
            Data data_;
    };

    bool operator==(const Attribute::Value &, const Attribute::Value &);
    bool operator!=(const Attribute::Value &, const Attribute::Value &);

};


#endif
