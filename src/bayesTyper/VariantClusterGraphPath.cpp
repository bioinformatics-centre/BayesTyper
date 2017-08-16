
/*
VariantClusterGraphPath.cpp - This file is part of BayesTyper (v1.1)


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


#include "VariantClusterGraphPath.hpp"


template <uchar kmer_size>
VariantClusterGraphPath<kmer_size>::VariantClusterGraphPath(const ushort num_variants) {

	path.reserve(num_variants * 2);

	path_score.first = 0;
	path_score.second = 0;
}

template <uchar kmer_size>
void VariantClusterGraphPath<kmer_size>::addVertex(const uint vertex_idx, const VariantClusterGraphVertex & vertex) {

	if (path.empty()) {

		assert(vertex_idx == 0);
	
	} else {

		assert(path.back().first < vertex_idx);
	}

	path.emplace_back(vertex_idx, &vertex);
}

template <uchar kmer_size>
const vector<pair<uint, const VariantClusterGraphVertex *> > & VariantClusterGraphPath<kmer_size>::getPath() const {

	return path;
}

template <uchar kmer_size>
void VariantClusterGraphPath<kmer_size>::updateScore(KmerCounts * kmer_counts, const ushort sample_idx) {

	if (kmer_counts) {

		if (kmer_counts->getSampleCount(sample_idx) > 0) {

			path_score.first++;
		}
	}

	path_score.second++;
}

template <uchar kmer_size>
double VariantClusterGraphPath<kmer_size>::getPathKmerScore() const {

	assert(path_score.first <= path_score.second);

	if (path_score.second > 0) {
		
		return (path_score.first / static_cast<double>(path_score.second));
	
	} else {

		return 1;
	} 
}

template <uchar kmer_size>
uint VariantClusterGraphPath<kmer_size>::getPathVertexScore(const unordered_set<uint> & covered_vertices) const {

	uint path_vertex_score = 0;

	for (auto & vertex: path) {

		assert(vertex.second);

		if (covered_vertices.count(vertex.first) < 1) {

			path_vertex_score++;
		}
	}

	return path_vertex_score;
}

template <uchar kmer_size>
void VariantClusterGraphPath<kmer_size>::updateCoveredVertices(unordered_set<uint> * covered_vertices) const {

	for (auto & vertex: path) {

		assert(vertex.second);
		covered_vertices->insert(vertex.first);
	}
}

template <uchar kmer_size>
bool VariantClusterGraphPath<kmer_size>::operator == (const VariantClusterGraphPath<kmer_size> & rhs) const {

	auto path_it_1 = this->path.crbegin();
	auto path_it_2 = rhs.path.crbegin();

	assert(path_it_1 != this->path.crend());
	assert(path_it_2 != rhs.path.crend());

	auto sequence_it_1 = path_it_1->second->sequence.crbegin();
	auto sequence_it_2 = path_it_2->second->sequence.crbegin();

	bool is_disconnected_1 = false;
	bool is_disconnected_2 = false;

	while (true) {

		while (sequence_it_1 == path_it_1->second->sequence.crend()) {
			
			assert(path_it_1 != this->path.crend());

			if (path_it_1->second->is_disconnected) {

				is_disconnected_1 = true;
			}

			path_it_1++;

			if (path_it_1 != this->path.crend()) {

				sequence_it_1 = path_it_1->second->sequence.crbegin();
			
			} else {

				break;
			}
		}

		while (sequence_it_2 == path_it_2->second->sequence.crend()) {

			assert(path_it_2 != rhs.path.crend());

			if (path_it_2->second->is_disconnected) {

				is_disconnected_2 = true;
			}

			path_it_2++;

			if (path_it_2 != rhs.path.crend()) {

				sequence_it_2 = path_it_2->second->sequence.crbegin();
			
			} else {

				break;
			}
		}
		
		if (is_disconnected_1 != is_disconnected_2) {

			return false;
		}	

		is_disconnected_1 = false;
		is_disconnected_2 = false;

		if ((path_it_1 == this->path.crend()) or (path_it_2 == rhs.path.crend())) {

			break;
		}	

		if (sequence_it_1 == sequence_it_2) {

			sequence_it_1 = path_it_1->second->sequence.crend();
			sequence_it_2 = path_it_2->second->sequence.crend();
		}	

		while ((sequence_it_1 != path_it_1->second->sequence.crend()) and (sequence_it_2 != path_it_2->second->sequence.crend())) {

			if (*sequence_it_1 != *sequence_it_2) {

				return false;
			}

			sequence_it_1++;
			sequence_it_2++;
		}
	}

	if ((path_it_1 != this->path.crend()) or (path_it_2 != rhs.path.crend())) {

		return false;
	}

	return true;
}


template class VariantClusterGraphPath<31>;
template class VariantClusterGraphPath<39>;
template class VariantClusterGraphPath<47>;
template class VariantClusterGraphPath<55>;
template class VariantClusterGraphPath<63>;

