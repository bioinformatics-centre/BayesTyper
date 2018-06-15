
/*
VariantInfo.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__VariantInfo_hpp
#define __bayesTyper__VariantInfo_hpp

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include "Utils.hpp"

using namespace std;


class AlleleInfo {

public:

	uint ref_length;
	string sequence;
	string aco_att;

	AlleleInfo() {}
	AlleleInfo(const uint ref_length_in, const string & sequence_in, const string & aco_att_in) : ref_length(ref_length_in), sequence(sequence_in), aco_att(aco_att_in) {}

private:

	friend class boost::serialization::access;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {

    	ar & ref_length;
    	ar & sequence;
    	ar & aco_att;
    }
};


class VariantInfo {

public:

	uint position;
	string id;

	bool has_dependency;
	vector<AlleleInfo> alt_alleles;

	VariantInfo() {}
	VariantInfo(const uint position_in, const string & id_in, const bool has_dependency_in, const vector<AlleleInfo> & alt_alleles_in) : position(position_in), id(id_in), has_dependency(has_dependency_in), alt_alleles(alt_alleles_in) {}

	ushort numberOfAlleles() const {

		return (1 + static_cast<ushort>(has_dependency) + alt_alleles.size());
	} 

	bool isMissing(const ushort allele_idx) const {

		assert(allele_idx < numberOfAlleles());

		if (has_dependency and (allele_idx == numberOfAlleles() - 1)) {

			return true;

		} else {

			return false;
		}
	}

	uint maxReferenceLength() const {

		uint max_ref_length = 0;

		for (auto & alt_allele: alt_alleles) {

			max_ref_length = max(max_ref_length, alt_allele.ref_length);
		}

		return max_ref_length;
	} 	

private:

	friend class boost::serialization::access;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {

    	ar & position;
    	ar & id;
    	ar & has_dependency;
    	ar & alt_alleles;
    }
};


#endif