
/*
KmerStats.hpp - This file is part of BayesTyper (v0.9)


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
		void addKmer(KmerCounts * const, const ushort, const float);

		float getCount() const;		
		pair<float, bool> getFraction() const;
		pair<float, bool> getMean() const;

	protected:

		float count;
		float fraction;
		float mean;
};

class FixedKmerStats {

	public:

		FixedKmerStats();
		FixedKmerStats(const float, const float, const float);

		float getCount() const;		
		float getFraction() const;
		float getMean() const;

		float count;
		float fraction;
		float mean;
};

class MedianKmerStats {

	public:

		MedianKmerStats();

		void addKmerStats(const KmerStats &);

		pair<float, bool> getMedianCount();		
		pair<float, bool> getMedianFraction();
		pair<float, bool> getMedianMean();

	private:

		struct floatLess {

			bool operator() (const float & first, const float & second) const {

				if (Utils::floatCompare(first, second)) {

					return false;

				} else {

					return first < second;
				}
			}
		};

		pair<float, bool> median(SortedLinearMap<float, uint, floatLess> &, const uint);

		SortedLinearMap<float, uint, floatLess> count_values;
		uint num_count_values;

		SortedLinearMap<float, uint, floatLess> fraction_values;
		uint num_fraction_values;

		SortedLinearMap<float, uint, floatLess> mean_values;
		uint num_mean_values;
};

class VariantKmerStats {

	public:

		vector<vector<MedianKmerStats> > allele_kmer_stats;

		VariantKmerStats();
		VariantKmerStats(const ushort, const ushort);

		void addAlleleKmerStats(const KmerStats &, const ushort, const ushort);
};

#endif