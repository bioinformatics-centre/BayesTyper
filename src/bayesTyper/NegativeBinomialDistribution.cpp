
/*
NegativeBinomialDistribution.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <cmath>
#include <assert.h>
#include <iostream>
#include <limits>

#include "Utils.hpp"
#include "NegativeBinomialDistribution.hpp"

static const double max_p = 0.99;


std::pair<double, double> NegativeBinomialDistribution::methodOfMomentsEst(const std::vector<ulong> & samples) {

    double mean = 0;
    double var = 0;

    double num_obs = 0;

    for (uint i = 0; i < samples.size(); i++) {

        mean += static_cast<double>(i * samples.at(i));
        num_obs += samples.at(i);
    }

    mean /= num_obs;

    assert(num_obs > 1);

    for (uint i = 0; i < samples.size(); i++) {

        var += std::pow(static_cast<double>(i) - mean, 2) * static_cast<double>(samples.at(i));
    }

    var /= num_obs - 1;

    if (max_p < (mean / var)) {

        var = mean / max_p;
    }

    double p = mean / var;
    double size = std::pow(mean, 2)/(var - mean);

    return std::make_pair(p, size);
}

NegativeBinomialDistribution::NegativeBinomialDistribution() {

    p_ = max_p;
    size_ = p_/(1 - p_);
}

NegativeBinomialDistribution::NegativeBinomialDistribution(std::pair<double, double> parameters) {

    assert(parameters.first > 0);
    assert(parameters.first < 1);
    assert(parameters.second > 0);

    p_ = parameters.first;
    size_ = parameters.second;
}

double NegativeBinomialDistribution::p() const {

    return p_;
}

double NegativeBinomialDistribution::size() const {

    return size_;
}

double NegativeBinomialDistribution::mean() const {

    return size_ * (1 - p_)/p_;
}

double NegativeBinomialDistribution::var() const {

    return size_ * (1 - p_)/std::pow(p_,2);
}

void NegativeBinomialDistribution::p(double p) {

    assert(p > 0);
    assert(p < 1);

    p_ = p;
}

void NegativeBinomialDistribution::size(double size) {

    assert(size > 0);

    size_ = size;
}

double NegativeBinomialDistribution::logBinomialCoefTerm(uint obs, uint size_scale) const {

    return std::lgamma(obs + size_ * size_scale) - std::lgamma(size_ * size_scale) - std::lgamma(obs + 1);
}

double NegativeBinomialDistribution::pmf(uint obs) const {

    return pmf(obs, 1);
}

double NegativeBinomialDistribution::pmf(uint obs, uint size_scale) const {

    assert(size_scale > 0);
    return std::exp(logBinomialCoefTerm(obs, size_scale)) *  std::pow(p_, size_ * size_scale) * std::pow((1 - p_), obs);
}

double NegativeBinomialDistribution::logPmf(uint obs) const {

    return logPmf(obs, 1);
}

double NegativeBinomialDistribution::logPmf(uint obs, uint size_scale) const {

    assert(size_scale > 0);
    return logBinomialCoefTerm(obs, size_scale) + std::log(p_) * size_ * size_scale + std::log(1 - p_) * obs;
}

uint NegativeBinomialDistribution::quantile(double quantile) const {

    assert(quantile > 0);
    assert(quantile <= 1);

    double log_quantile = log(quantile);
    assert(std::isfinite(log_quantile));

    uint count = 0;
    double cum_log_prob = logPmf(count);

    while (cum_log_prob < log_quantile) {

        count++;
        cum_log_prob = Utils::logAddition(cum_log_prob, logPmf(count));
        assert(std::isfinite(cum_log_prob));
    }

    return count;
}
