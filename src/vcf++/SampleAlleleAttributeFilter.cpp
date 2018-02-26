
/*
SampleAlleleAttributeFilter.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "SampleAlleleAttributeFilter.hpp"


SampleAlleleAttributeFilter::SampleAlleleAttributeFilter(AttributeFilter * att_flt_in, Attribute::ReductionOp * red_op_in) {

    att_flt = att_flt_in;
    red_op = red_op_in;
}

bool SampleAlleleAttributeFilter::pass(Sample & sample) const {

    if (!sample.isInformative()) {

        return false;
    }

    auto gt_est = sample.genotypeEstimate();

    vector<bool> filter_pass_status;
    for (auto all_idx : gt_est) {

        if (!sample.alleleIsMissing(all_idx)) {

            filter_pass_status.push_back(att_flt->pass(sample.alleleInfo().at(all_idx)));
        }
    }

    if (filter_pass_status.empty()) {

        return true;

    } else if (filter_pass_status.size() == 1) {

        return filter_pass_status.front();

    } else {

        assert(filter_pass_status.size() == 2);
        return (*red_op)(filter_pass_status.front(),filter_pass_status.back());
    }
}
