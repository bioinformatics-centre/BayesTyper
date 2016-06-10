
/*
SparsityEstimator.cpp - This file is part of BayesTyper (v0.9)


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


SparsityEstimator::SparsityEstimator(const uint prng_seed) {

    prng = mt19937(prng_seed);
}

uint SparsityEstimator::estimateMinimumSetCover(Eigen::MatrixXuchar const & data_element_map, Eigen::RowVectorXbool * non_zero_counts) {
    
    assert(data_element_map.rows() == non_zero_counts->size());
    assert(in_minimum_set.empty());

    in_minimum_set = std::vector<bool>(data_element_map.cols(), false);
    uint minimum_set_size = 0;

    Eigen::MatrixXbool binary_data_element_map((data_element_map.array() > 0).matrix().cast<bool>());

    while (non_zero_counts->sum() > 0) {
        
        Eigen::ColVectorXbool covered_data(*non_zero_counts * binary_data_element_map);

        uint max_val = covered_data.maxCoeff();
        assert(max_val > 0);
        
        vector<uint> max_indices;
        max_indices.reserve(covered_data.size());
        
        for (uint i = 0; i < covered_data.size(); i++) {
            
            if (max_val == covered_data[i]) {
                
                max_indices.push_back(i);
            }
        }
        	
        uniform_int_dist.param(uniform_int_distribution<>::param_type(0, max_indices.size()-1));
        uint max_pos = max_indices.at(uniform_int_dist(prng));
        
        assert(!in_minimum_set.at(max_pos));

        in_minimum_set.at(max_pos) = true;
        minimum_set_size++;

        *non_zero_counts = (non_zero_counts->array() - non_zero_counts->array() * binary_data_element_map.col(max_pos).transpose().array()).matrix();
    }

    return minimum_set_size;
}


vector<bool> & SparsityEstimator::getMinimumSet() {

    return in_minimum_set;
}

