
/*
AttributeFilter.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef ATTRIBUTE_FILTER
#define ATTRIBUTE_FILTER

#include <string>

#include "CompareOperators.hpp"
#include "Attribute.hpp"
#include "AttributeSet.hpp"

class AttributeFilter {

    public:

        AttributeFilter(const string &);
        virtual ~AttributeFilter();
        const string & id() const;
        virtual bool pass(const AttributeSet &) = 0;

    protected:

        string id_;
};

class FlagAttributeFilter : public AttributeFilter {

    public:

        FlagAttributeFilter(const string &);
        ~FlagAttributeFilter();

        bool pass(const AttributeSet &);
};

class ValueAttributeFilter : public AttributeFilter {

    public:

        ValueAttributeFilter(const string &, Attribute::Value, Attribute::CmpOp *);
        ~ValueAttributeFilter();

        const string & id() const;
        const Attribute::Value & value() const;
        void setValue(Attribute::Value);
        bool pass(const AttributeSet &);

    private:

        Attribute::Value value_;
        Attribute::CmpOp * pass_operator;
};

#endif
