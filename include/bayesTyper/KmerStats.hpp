
/*
KmerStats.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__KmerStats_hpp
#define __bayesTyper__KmerStats_hpp

#include <vector>

#include "Utils.hpp"
#include "KmerCounts.hpp"
#include "LinearMap.hpp"

using namespace std;

class KmerStats {

	public:

		KmerStats();

		void reset();
		void addValue(const pair<double, bool> &);

		uint getCount() const;		
		pair<double, bool> getFraction() const;
		pair<double, bool> getMean() const;
		pair<double, bool> getVariance() const;

	private:

		uint count;
		double fraction;
		double mean;
		double M2;
};

class AlleleKmerStats {

	public:

		vector<KmerStats> count_stats;
		vector<KmerStats> fraction_stats;
		vector<KmerStats> mean_stats;

		AlleleKmerStats() {};
		AlleleKmerStats(const ushort);

		void addKmerStats(const KmerStats &, const ushort);
};

#endif