
/*
OptionsContainer.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__OptionsContainer_hpp
#define __bayesTyper__OptionsContainer_hpp

#include <string>
#include <vector>
#include <time.h>
#include <unordered_map>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;


class OptionsContainer {

    private:

        struct OptionValueBase {

            virtual ~OptionValueBase() {};
        };

        template<typename ValueType>
        struct OptionValue : public OptionValueBase {

            OptionValue(ValueType value_in) : value(value_in) {}
            ValueType value;
        };

        std::unordered_map<std::string, std::pair<OptionValueBase*, std::string> > options;

        const string version;
        const string start_time;

    public:

        inline OptionsContainer(const std::string &, const std::string &);
        inline ~OptionsContainer();

        template<typename ValueType>
        void parseValue(std::string, ValueType);

        template<typename ValueType>
        void parseValuePair(std::string, std::string);

        inline void parseFiles(std::string, std::string);
        inline void parseRegions(std::string, std::string);

        template<typename ValueType>
        const ValueType & getValue(std::string) const;

        template<typename ValueType>
        void updateValue(std::string, ValueType);

        inline std::string writeHeader();
};

#include "OptionsContainer.tpp"

#endif