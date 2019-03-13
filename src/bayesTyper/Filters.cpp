
/*
Filters.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "Filters.hpp"


static const float observed_kmer_beta_value = 0.275;

Filters::Filters(const OptionsContainer & options_container, const vector<vector<NegativeBinomialDistribution> > & genomic_count_distributions) {

    min_genotype_posterior = options_container.getValue<float>("min-genotype-posterior");
    min_number_of_kmers = options_container.getValue<float>("min-number-of-kmers");

    min_fraction_observed_kmers = vector<float>(genomic_count_distributions.size(), 0);

    if (!options_container.getValue<bool>("disable-observed-kmers")) {

        for (ushort sample_idx = 0; sample_idx < genomic_count_distributions.size(); sample_idx++) {

            assert(genomic_count_distributions.at(sample_idx).size() == 1);

            min_fraction_observed_kmers.at(sample_idx) = 1 - exp(-(observed_kmer_beta_value * genomic_count_distributions.at(sample_idx).front().mean()));

            assert(min_fraction_observed_kmers.at(sample_idx) >= 0);
            assert(min_fraction_observed_kmers.at(sample_idx) <= 1);
        }
    }
}

float Filters::minGenotypePosterior() const {

    return min_genotype_posterior;
}

float Filters::minNumberOfKmers() const {

    return min_number_of_kmers;
}

float Filters::minFractionObservedKmers(const ushort sample_idx) const {

    return min_fraction_observed_kmers.at(sample_idx);
}




