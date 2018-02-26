
/*
FrequencyDistribution.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random> 
#include <algorithm> 

#include "boost/math/special_functions/gamma.hpp"

#include "FrequencyDistribution.hpp"
#include "Utils.hpp"

const double FrequencyDistribution::dirichlet_parameter = 1;

FrequencyDistribution::FrequencyDistribution(const ushort num_elements, const uint prng_seed) {

    prng = mt19937(prng_seed);

	observation_counts = vector<ushort>(num_elements, 0); 
	frequencies = vector<double>(num_elements, double(1)/num_elements); 
    non_zero_frequencies = vector<bool>(num_elements, true); 
}


void FrequencyDistribution::reset() {

    auto num_elements = observation_counts.size();

    observation_counts = vector<ushort>(num_elements, 0); 
    frequencies = vector<double>(num_elements, double(1)/num_elements); 
    non_zero_frequencies = vector<bool>(num_elements, true); 
}


ushort FrequencyDistribution::getNumElements() {

    return frequencies.size();
}


void FrequencyDistribution::incrementObservationCount(const ushort element_idx) {

	observation_counts.at(element_idx)++;
}


pair<bool, double> FrequencyDistribution::getElementFrequency(const ushort element_idx) {

	return pair<bool, double>(non_zero_frequencies.at(element_idx), frequencies.at(element_idx));
}


void FrequencyDistribution::sampleFrequencies(const uint sum_observation_counts) {

	double norm_const = 0;

	for (uint i = 0; i < frequencies.size(); i++) {

		gamma_dist.param(gamma_distribution<>::param_type(observation_counts.at(i) + 1, 1));

        frequencies.at(i) = gamma_dist(prng);		
        norm_const += frequencies.at(i);

        observation_counts.at(i) = 0;
	}

    for (auto &freq: frequencies) {

        freq /= norm_const;
    }
}


SparseFrequencyDistribution::SparseFrequencyDistribution(const ushort num_elements, uint prng_seed, double sparsity_in) : FrequencyDistribution(num_elements, prng_seed), sparsity(sparsity_in) {

    assert(sparsity > 0); 
    assert(sparsity < 1);

    for (uint i = 0; i < frequencies.size(); i++) {

        zero_count_indices.insert(i);
    }
}


void SparseFrequencyDistribution::reset() {

    FrequencyDistribution::reset();

    plus_count_indices.clear();
    zero_count_indices.clear();

    for (uint i = 0; i < frequencies.size(); i++) {

        zero_count_indices.insert(i);
    }
}


void SparseFrequencyDistribution::updateCachedSimplexProbVector(vector<double> * simplex_prob_vector, const uint total_number_of_observations, const ushort count_plus_size) {

    assert(total_number_of_observations > 0);
    assert(count_plus_size > 0);
    assert(count_plus_size <= frequencies.size());

    simplex_prob_vector->reserve(frequencies.size() - count_plus_size + 1);
                
	// Cardinality of equivalence class of size |s+| is zero
	double cardinal_eq_z_log = 0;

	// Calculate probability of member of equivalence class
	double prob_z_log = count_plus_size*log(sparsity) + (frequencies.size()-count_plus_size)*log(1-sparsity);

	// Probability of assignment given the binary vector
	double prob_t_log = boost::math::lgamma(count_plus_size*dirichlet_parameter) - boost::math::lgamma(total_number_of_observations + count_plus_size*dirichlet_parameter);

	// Full probability of the binary vector		
	double prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
	double row_sum = prob_eq_z_log;

	simplex_prob_vector->push_back(row_sum);

	for (ushort j = count_plus_size + 1; j < frequencies.size() + 1; j++) {
	
		// Calculate cardinality of equivalence class
		cardinal_eq_z_log = boost::math::lgamma(frequencies.size()-count_plus_size+1)-(boost::math::lgamma(j-count_plus_size+1)+boost::math::lgamma(frequencies.size() - j + 1));
	
		// Calculate probability of member of equivalence class
		prob_z_log = j*log(sparsity) + (frequencies.size()-j)*log(1-sparsity);
	
		// Probability of assignment given the binary vector
		prob_t_log = boost::math::lgamma(j*dirichlet_parameter) - boost::math::lgamma(total_number_of_observations + j*dirichlet_parameter);
	
		// Full probability of the binary std::vector<char> v;
		prob_eq_z_log = cardinal_eq_z_log + prob_z_log + prob_t_log;
	    
        row_sum += log(1 + exp(prob_eq_z_log - row_sum));
		simplex_prob_vector->push_back(row_sum);
        
        if (Utils::doubleCompare(simplex_prob_vector->back(), *(simplex_prob_vector->rbegin() + 1))) {

            break;     
        }
	}

	// Row-normalise and transform back from log-space
	for (auto &prob: *simplex_prob_vector) {
	
	    prob = exp(prob - row_sum);
	}

    assert(Utils::doubleCompare(simplex_prob_vector->back(), 1));
}


void SparseFrequencyDistribution::incrementObservationCount(const ushort element_idx) {

    if (observation_counts.at(element_idx) == 0) {

        assert(plus_count_indices.insert(element_idx).second);
        assert(zero_count_indices.erase(element_idx));
    } 

    observation_counts.at(element_idx)++;
}


void SparseFrequencyDistribution::sampleFrequencies(const uint sum_observation_counts) {

    if (cached_simplex_prob_vectors.count(sum_observation_counts) > 0) {

        if (cached_simplex_prob_vectors.at(sum_observation_counts).count(plus_count_indices.size()-1) < 1) {

            auto emplaced_prob_vector = cached_simplex_prob_vectors.at(sum_observation_counts).emplace(pair<ushort, vector<double> >(plus_count_indices.size()-1, vector<double>()));
            assert(emplaced_prob_vector.second);

            updateCachedSimplexProbVector(&(emplaced_prob_vector.first->second), sum_observation_counts, plus_count_indices.size());
        }
    
    } else {

        assert(cached_simplex_prob_vectors.emplace(pair<uint, unordered_map<ushort, vector<double> > >(sum_observation_counts, unordered_map<ushort, vector<double> >())).second);
        auto emplaced_prob_vector = cached_simplex_prob_vectors.at(sum_observation_counts).emplace(pair<ushort, vector<double> >(plus_count_indices.size()-1, vector<double>()));
        assert(emplaced_prob_vector.second);

        updateCachedSimplexProbVector(&(emplaced_prob_vector.first->second), sum_observation_counts, plus_count_indices.size());
    }

    auto prob_vector = &cached_simplex_prob_vectors.at(sum_observation_counts).at(plus_count_indices.size()-1);
    ushort simplex_size = ushort(upper_bound(prob_vector->begin(), prob_vector->end(), generate_canonical<double,std::numeric_limits<double>::digits>(prng)) - prob_vector->begin()) + plus_count_indices.size();

    assert(simplex_size > 0);
    
    // Sample gammas observed alleles
    double norm_const = 0;
    
    for (auto &plus_count_idx: plus_count_indices) {
                
        gamma_dist.param(gamma_distribution<>::param_type(observation_counts.at(plus_count_idx) + dirichlet_parameter, 1));
                
        frequencies.at(plus_count_idx) = gamma_dist(prng);
        norm_const += frequencies.at(plus_count_idx);

        non_zero_frequencies.at(plus_count_idx) = true;
    }
    
    gamma_dist.param(gamma_distribution<>::param_type(dirichlet_parameter, 1));

    // Sample gammas for the expanded set of alleles
    while (plus_count_indices.size() < simplex_size) {
           
        uniform_int_dist.param(uniform_int_distribution<>::param_type(0, zero_count_indices.size()-1));            
        uint sampled_position = uniform_int_dist(prng);

        auto zero_count_indices_iter = zero_count_indices.begin();      
        uint zero_count_indices_position = 0; 

        while (zero_count_indices_position < sampled_position) {

            zero_count_indices_iter++;
            zero_count_indices_position++;
        }

        assert(zero_count_indices_position == sampled_position);
        assert(zero_count_indices_iter != zero_count_indices.end()); 
        assert(observation_counts.at(*zero_count_indices_iter) == 0);
                
        frequencies.at(*zero_count_indices_iter) = gamma_dist(prng);
        norm_const += frequencies.at(*zero_count_indices_iter);

        assert(frequencies.at(*zero_count_indices_iter) > 0);

        non_zero_frequencies.at(*zero_count_indices_iter) = true;

        assert(plus_count_indices.insert(*zero_count_indices_iter).second);
        assert(zero_count_indices.erase(*zero_count_indices_iter));   
    }

    assert((zero_count_indices.size() + plus_count_indices.size()) == frequencies.size());

    for (auto &zero_count_idx: zero_count_indices) {

        frequencies.at(zero_count_idx) = 0;
        non_zero_frequencies.at(zero_count_idx) = false;

        observation_counts.at(zero_count_idx) = 0;
    }

    for (auto &plus_count_idx: plus_count_indices) {

        frequencies.at(plus_count_idx) /= norm_const;
        assert(zero_count_indices.insert(plus_count_idx).second);

        assert(non_zero_frequencies.at(plus_count_idx) == true);

        observation_counts.at(plus_count_idx) = 0;
    }

    assert(zero_count_indices.size() == frequencies.size());

    plus_count_indices.clear();
}

