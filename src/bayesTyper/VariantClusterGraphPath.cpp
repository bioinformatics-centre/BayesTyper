

#include "VariantClusterGraphPath.hpp"


template <uchar kmer_size>
VariantClusterGraphPath<kmer_size>::VariantClusterGraphPath(const ushort num_variants) {

	path.reserve(num_variants * 2);

	path_score.first = 0;
	path_score.second = 0;
}

template <uchar kmer_size>
bool VariantClusterGraphPath<kmer_size>::operator == (const VariantClusterGraphPath<kmer_size> & rhs) const {

	auto path_it_1 = this->path.crbegin();
	auto path_it_2 = rhs.path.crbegin();

	assert(path_it_1 != this->path.crend());
	assert(path_it_2 != rhs.path.crend());

	auto sequence_it_1 = (*path_it_1)->sequence.crbegin();
	auto sequence_it_2 = (*path_it_2)->sequence.crbegin();

	bool is_disconnected_1 = false;
	bool is_disconnected_2 = false;

	while (true) {

		while (sequence_it_1 == (*path_it_1)->sequence.crend()) {
			
			assert(path_it_1 != this->path.crend());

			if ((*path_it_1)->is_disconnected) {

				is_disconnected_1 = true;
			}

			path_it_1++;

			if (path_it_1 != this->path.crend()) {

				sequence_it_1 = (*path_it_1)->sequence.crbegin();
			
			} else {

				break;
			}
		}

		while (sequence_it_2 == (*path_it_2)->sequence.crend()) {

			assert(path_it_2 != rhs.path.crend());

			if ((*path_it_2)->is_disconnected) {

				is_disconnected_2 = true;
			}

			path_it_2++;

			if (path_it_2 != rhs.path.crend()) {

				sequence_it_2 = (*path_it_2)->sequence.crbegin();
			
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

			sequence_it_1 = (*path_it_1)->sequence.crend();
			sequence_it_2 = (*path_it_2)->sequence.crend();
		}	

		while ((sequence_it_1 != (*path_it_1)->sequence.crend()) and (sequence_it_2 != (*path_it_2)->sequence.crend())) {

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

template <uchar kmer_size>
void VariantClusterGraphPath<kmer_size>::addVertex(const VariantClusterGraphVertex & vertex) {

	path.push_back(&vertex);
}

template <uchar kmer_size>
const vector<const VariantClusterGraphVertex *> & VariantClusterGraphPath<kmer_size>::getPath() const {

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
uint VariantClusterGraphPath<kmer_size>::getPathVertexScore(const unordered_set<const VariantClusterGraphVertex *> & covered_vertices) const {

	uint path_vertex_score = 0;

	for (auto & vertex: path) {

		if (covered_vertices.count(vertex) < 1) {

			path_vertex_score++;
		}
	}

	return path_vertex_score;
}

template <uchar kmer_size>
void VariantClusterGraphPath<kmer_size>::updateCoveredVertices(unordered_set<const VariantClusterGraphVertex *> * covered_vertices) const {

	for (auto & vertex: path) {

		covered_vertices->insert(vertex);
	}
}


template class VariantClusterGraphPath<31>;
template class VariantClusterGraphPath<39>;
template class VariantClusterGraphPath<47>;
template class VariantClusterGraphPath<55>;
template class VariantClusterGraphPath<63>;

