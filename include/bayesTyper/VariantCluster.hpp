
/*
VariantCluster.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef __bayesTyper__VariantCluster_hpp
#define __bayesTyper__VariantCluster_hpp

#include <map>
#include <set>

#include "Utils.hpp"

using namespace std;

class VariantFileParser;
class VariantClusterGraph;

class VariantCluster {

	friend class VariantFileParser;
	friend class VariantClusterGraph;
	friend class VariantClusterGroup;

	protected:

		struct Variant {

			const uint id;
			const bool has_dependency;

			Utils::VariantType type;			
			uint max_reference_length;
			uint num_redundant_nucleotides;
			
			vector<ushort> excluded_alt_alleles;
			vector<pair<ushort, pair<uint, const string> > > alternative_alleles;

			Variant(const uint id_in, const bool has_dependency_in) : id(id_in), has_dependency(has_dependency_in) {

				type = Utils::VariantType::Unsupported;
				max_reference_length = 0;
				num_redundant_nucleotides = Utils::uint_overflow;
			}
		};

		struct ContainedCluster {

			const uint cluster_idx;
			const uint left_flank;
			const uint right_flank;

			ContainedCluster(uint cluster_idx_in, uint left_flank_in, uint right_flank_in) : cluster_idx(cluster_idx_in), left_flank(left_flank_in), right_flank(right_flank_in) {}
		};

	private:

		struct ContainedClusterCompare {
    		
    		bool operator() (const ContainedCluster & lhs, const ContainedCluster & rhs) const {
        	   	
    			assert(lhs.left_flank <= lhs.right_flank);
    			assert(rhs.left_flank <= rhs.right_flank);

    			assert((lhs.right_flank < rhs.left_flank) or (rhs.right_flank < lhs.left_flank));

        		return lhs.left_flank < rhs.left_flank;
    		}
    	};

    public:

		uint cluster_idx;
		uint left_flank;
		uint right_flank;

		string chromosome_id;
		string chromosome_name;
		Utils::ChromosomeClass chromosome_class;

		map<uint, Variant> variants;
		set<ContainedCluster, ContainedClusterCompare> contained_clusters;

		void mergeVariantClusters(const VariantCluster &);
};

#endif