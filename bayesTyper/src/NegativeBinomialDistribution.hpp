
/*
NegativeBinomialDistribution.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__NegativeBinomialDistribution_hpp
#define __bayesTyper__NegativeBinomialDistribution_hpp

#include <unordered_map>

typedef unsigned int uint;
typedef unsigned long int ulong;

class NegativeBinomialDistribution {


    public:

        static std::pair<float, float> methodOfMomentsEst(const std::unordered_map<uint,ulong> &);

        NegativeBinomialDistribution();
        NegativeBinomialDistribution(std::pair<float, float>);

        float p() const;
        float size() const;
        float mean() const;
        float var() const;

        void p(float);
        void size(float);
        float pmf(uint) const;
        float pmf(uint, uint) const;
        float logPmf(uint) const;
        float logPmf(uint, uint) const;

        uint quantile(float) const;

    private:

        float logBinomialCoefTerm(uint, uint) const;

        float p_;
        float size_;
};

#endif
