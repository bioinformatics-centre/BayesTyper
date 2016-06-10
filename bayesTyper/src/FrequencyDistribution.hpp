
/*
FrequencyDistribution.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__FrequencyDistribution_hpp
#define __bayesTyper__FrequencyDistribution_hpp

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random> 

#include "Utils.hpp"

class FrequencyDistribution {

	public:

		FrequencyDistribution(const ushort, const uint);

		virtual void reset();

		pair<bool, double> getElementFrequency(const ushort);
		ushort getNumElements();

		virtual ~FrequencyDistribution() {};
		virtual void incrementObservationCount(const ushort);		
		virtual void sampleFrequencies(const uint);

	protected:

		static const double dirichlet_parameter;

		vector<ushort> observation_counts;
		vector<double> frequencies;	
		vector<bool> non_zero_frequencies;

		mt19937 prng;
		gamma_distribution<> gamma_dist;
};


class SparseFrequencyDistribution : public FrequencyDistribution {

	public:
		
		SparseFrequencyDistribution(const ushort, const uint, const double);

		void reset();

		void incrementObservationCount(const ushort);		
		void sampleFrequencies(const uint);

	private:

		void updateCachedSimplexProbVector(vector<double> *, const uint, const ushort);

		const double sparsity;
		
		unordered_map<uint, unordered_map<ushort, vector<double> > > cached_simplex_prob_vectors; 

	    unordered_set<ushort> plus_count_indices;
    	unordered_set<ushort> zero_count_indices; 

		uniform_int_distribution<> uniform_int_dist;
};

#endif	