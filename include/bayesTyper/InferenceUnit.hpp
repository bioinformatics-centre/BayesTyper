
/*
VariantClusterGroup.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__InferenceUnit_hpp
#define __bayesTyper__InferenceUnit_hpp

#include <vector>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include "VariantClusterGroup.hpp"


class InferenceUnit {

	public:

		uint index;

		string cluster_options_header;

		uint num_variants;
		uint num_variant_clusters;

		ulong num_path_kmers;

		std::vector<VariantClusterGroup*> variant_cluster_groups;

		InferenceUnit() {}
		InferenceUnit(const uint index_in) : index(index_in), num_variants(0), num_variant_clusters(0), num_path_kmers(0) {}

	private:

		friend class boost::serialization::access;

		template<class Archive>
	    void serialize(Archive & ar, const unsigned int version) {

	    	ar & index;
	    	ar & cluster_options_header;
	    	ar & num_variants;
	    	ar & num_variant_clusters;
	    	ar & num_path_kmers;
	    	ar & variant_cluster_groups;
	    }

};

#endif