
/*
SparsityEstimator.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <unordered_set>

#include "SparsityEstimator.hpp"
#include "DiscreteSampler.hpp"


SparsityEstimator::SparsityEstimator(const uint prng_seed) {

    prng = mt19937(prng_seed);
}

unordered_set<ushort> SparsityEstimator::estimateMinimumColumnCover(const Utils::MatrixXuchar & data_matrix, Utils::RowVectorXbool * uncovered_rows) {

    assert(data_matrix.cols() < Utils::ushort_overflow);    
    assert(data_matrix.rows() == uncovered_rows->size());

    unordered_set<ushort> min_column_cover;

    while (uncovered_rows->sum() > 0) {
        
        Utils::RowVectorXuint column_row_cover = (*uncovered_rows).cast<uint>() * data_matrix.cast<uint>();
        assert(column_row_cover.size() == data_matrix.cols());

        uint max_row_cover = column_row_cover.maxCoeff();
        assert(max_row_cover > 0);
    
        DiscreteSampler column_sampler(column_row_cover.size());

        vector<ushort> max_row_cover_column_indices;
        max_row_cover_column_indices.reserve(column_row_cover.size());
        
        for (ushort column_idx = 0; column_idx < column_row_cover.size(); column_idx++) {
            
            if (column_row_cover[column_idx] == max_row_cover) {
                
                column_sampler.addOutcome(1);
                max_row_cover_column_indices.push_back(column_idx);
            }
        }
        	
        const ushort sampled_column_idx = max_row_cover_column_indices.at(column_sampler.sample(&prng)); 
        assert(min_column_cover.insert(sampled_column_idx).second);

        *uncovered_rows = *uncovered_rows - (uncovered_rows->array() * (data_matrix.col(sampled_column_idx).transpose().cast<bool>().array())).matrix();
    }

    return min_column_cover;
}

