
/*
VariantClusterGraphPath.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <functional>

#include "VariantClusterGraphPath.hpp"
#include "Nucleotide.hpp"


static const uchar min_observed_kmers = 2;

VariantClusterGraphPath::VariantClusterGraphPath(const ushort num_variants) {

	path.reserve(num_variants * 2);

	kmer_score.first = 0;
	kmer_score.second = 0;
}

void VariantClusterGraphPath::addVertex(const uint vertex_idx, const VariantClusterGraphVertex & vertex, KmerBloom<Utils::kmer_size> * kmer_bloom) {

	if (path.empty()) {

		assert(vertex_idx == 0);
	
	} else {

		assert(path.back().index < vertex_idx);
	}

	path.emplace_back(vertex_idx, &vertex); 

    if (vertex.is_disconnected) {

    	if (vertex.sequence.empty()) {

    		path.back().num_observed_kmers = min_observed_kmers;
    	}

    	kmer_pair.reset();
    } 

    bitset<2> nt_bits;

    assert((vertex.sequence.size() % 2) == 0);
    auto sequence_it = vertex.sequence.begin();

    while (sequence_it != vertex.sequence.end()) {

        nt_bits.set(0, *sequence_it);
        sequence_it++;

        nt_bits.set(1, *sequence_it);
        sequence_it++;

        if (kmer_pair.move(make_pair(nt_bits, true))) {

            updateScore(kmer_bloom->lookup(kmer_pair.getLexicographicalLowestKmer()), (sequence_it - vertex.sequence.begin()) / 2);
        }
    }
}

void VariantClusterGraphPath::updateScore(const bool is_kmer_observed, uint cur_sequence_length) {

	assert(cur_sequence_length > 0);
	assert(!(path.empty()));

	if (is_kmer_observed) {

		kmer_score.first++;

		auto path_rit = path.rbegin();
		assert(path_rit != path.rend());

	    if ((cur_sequence_length > 1) or !(path_rit->vertex->is_first_nucleotides_redundant)) {

	    	if (path_rit->num_observed_kmers < min_observed_kmers) {

	    		path_rit->num_observed_kmers++;
	    	}
	    }
		    
		path_rit++;

		while (path_rit != path.rend()) {

	        if ((Utils::kmer_size <= cur_sequence_length) or (path_rit->num_observed_kmers == min_observed_kmers)) {

	            break;
	        }

	    	if (path_rit->num_observed_kmers < min_observed_kmers) {

	    		path_rit->num_observed_kmers++;
	    	}

		    cur_sequence_length += (path_rit->vertex->sequence.size() / 2);

			path_rit++;
		}
	}

	kmer_score.second++;
}

const vector<typename VariantClusterGraphPath::PathVertexInfo> & VariantClusterGraphPath::getPath() const {

	return path;
}

double VariantClusterGraphPath::getKmerScore() const {

	assert(kmer_score.first <= kmer_score.second);

	if (kmer_score.second > 0) {
		
		return (kmer_score.first / static_cast<double>(kmer_score.second));
	
	} else {

		return 1;
	} 
}

uint VariantClusterGraphPath::getVertexScore(const vector<bool> & covered_observed_vertices, const bool is_complete) const {

	uint observed_vertex_score = 0;

	for (auto & vertex: path) {

		assert(vertex.num_observed_kmers <= min_observed_kmers);

		if ((vertex.num_observed_kmers == min_observed_kmers) and !(covered_observed_vertices.at(vertex.index))) {

			observed_vertex_score++;
		}
	}

	if (!is_complete) {

		auto path_rit = path.rbegin();
		assert(path_rit != path.rend());

		uint cur_sequence_length = 0;

		while (path_rit != path.rend()) {

	        if (((Utils::kmer_size - 1) <= cur_sequence_length) or (path_rit->num_observed_kmers == min_observed_kmers)) {

	            break;
	        }

			if (!(covered_observed_vertices.at(path_rit->index))) {

				observed_vertex_score++;
			}

		    cur_sequence_length += (path_rit->vertex->sequence.size() / 2);
			path_rit++;
		}
	}

	assert(observed_vertex_score <= path.size());
	return observed_vertex_score;
}

void VariantClusterGraphPath::updateObservedCoveredVertices(vector<bool> * covered_observed_vertices, const bool is_complete) const {

	for (auto & vertex: path) {

		assert(vertex.num_observed_kmers <= min_observed_kmers);

		if (vertex.num_observed_kmers == min_observed_kmers) {

			covered_observed_vertices->at(vertex.index) = true;
		}
	}

	if (!is_complete) {

		auto path_rit = path.rbegin();
		assert(path_rit != path.rend());

		uint cur_sequence_length = 0;

		while (path_rit != path.rend()) {

	        if (((Utils::kmer_size - 1) <= cur_sequence_length) or (path_rit->num_observed_kmers == min_observed_kmers)) {

	            break;
	        }

			covered_observed_vertices->at(path_rit->index) = true;

		    cur_sequence_length += (path_rit->vertex->sequence.size() / 2);
			path_rit++;
		}
	}
}

