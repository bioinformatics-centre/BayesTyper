
template<uchar kmer_size>
VariantClusterGraphKmerPath<kmer_size>::VariantClusterGraphKmerPath(const uint expected_path_length, const bool is_dummy_in) : is_dummy(is_dummy_in) {

	kmers.reserve(expected_path_length);
	path.reserve(expected_path_length);

	num_nucleotides = 0;
	num_reference_alleles = 0;
}


template<uchar kmer_size>
bool VariantClusterGraphKmerPath<kmer_size>::operator == (const VariantClusterGraphKmerPath<kmer_size> & rhs) const {

	auto lvit = this->path.crbegin();
	auto rvit = rhs.path.crbegin();

	auto lsit = (*lvit)->sequences.crbegin();
	auto rsit = (*rvit)->sequences.crbegin();

	auto lkit = lsit->crbegin();
	auto rkit = rsit->crbegin();

	bool left_kmer_overlap = true;
	bool right_kmer_overlap = true;

	while ((lvit != this->path.crend()) and (rvit != rhs.path.crend())) {

		while (lsit == (*lvit)->sequences.crend()) {

			lvit++;

			if (lvit != this->path.crend()) {

				lsit = (*lvit)->sequences.crbegin();
				lkit = lsit->crbegin();
			
			} else {

				break;
			}
		}

		while (rsit == (*rvit)->sequences.crend()) {

			rvit++;

			if (rvit != rhs.path.crend()) {

				rsit = (*rvit)->sequences.crbegin();
				rkit = rsit->crbegin();
			
			} else {

				break;
			}
		}

		if ((lvit == this->path.crend()) or (rvit == rhs.path.crend())) {

			continue;
		}

		while ((lsit != (*lvit)->sequences.crend()) and (rsit != (*rvit)->sequences.crend())) {

			while (lkit == lsit->crend()) {

				lsit++;

				if (lsit != (*lvit)->sequences.crend()) {

					left_kmer_overlap = false;
					lkit = lsit->crbegin();
				
				} else {

					break;
				}
			}

			while (rkit == rsit->crend()) {

				rsit++;

				if (rsit != (*rvit)->sequences.crend()) {

					right_kmer_overlap = false;
					rkit = rsit->crbegin();
				
				} else {

					break;
				}
			}

			if ((lsit == (*lvit)->sequences.crend()) or (rsit == (*rvit)->sequences.crend())) {

				continue;
			}	

			if (left_kmer_overlap != right_kmer_overlap) {

				return false;
			}		

			if (lkit == rkit) {

				left_kmer_overlap = true;
				right_kmer_overlap = true;
				
				lkit = lsit->crend();
				rkit = rsit->crend();
				continue;
			}	

			while ((lkit != lsit->crend()) and (rkit != rsit->crend())) {

				if (*lkit != *rkit) {

					return false;
				} 

				left_kmer_overlap = true;
				right_kmer_overlap = true;

				lkit++;
				rkit++;
			}
		}
	}

	if ((lvit != this->path.crend()) or (rvit != rhs.path.crend())) {

		return false;
	}

	return true;
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::addVertex(VariantClusterGraphVertex & vertex, const uint variant_size) {

	path.push_back(&vertex);

	for (auto &vit: vertex.variant_ids) {

	    assert(vit.allele_id.first < Utils::ushort_overflow);
	    assert(vit.allele_id.second < Utils::ushort_overflow);
		assert(processed_variants.count(vit.allele_id.first) < 1);

	    if (vit.is_allele_vertex) {

	    	assert(&vit == &(vertex.variant_ids.front()));
	    	assert(variants.emplace(vit.allele_id, &vertex).second);

	    	if (vit.allele_id.second == 0) {

		    	num_reference_alleles++;
	    	}
	    }

	    auto running_variants_it = running_variants.find(vit.allele_id.first);

	    if (running_variants_it != running_variants.end()) {

	    	assert(!vit.is_allele_vertex);
 			assert((running_variants_it->second.last_nucleotide - kmer_size + 1) == num_nucleotides);

			running_variants_it->second.last_nucleotide += variant_size;

	    } else {

	    	assert((vit.allele_id.second == 0) or vit.is_allele_vertex);
			assert(running_variants.emplace(vit.allele_id.first, VariantScore(num_nucleotides + kmer_size + variant_size - 1)).second);
	    }
	}

	addNestedVariantClusterIndices(vertex);
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::addNestedVariantClusterIndices(VariantClusterGraphVertex &) {}


template<uchar kmer_size>
const vector<typename VariantClusterGraphKmerPath<kmer_size>::KmerMultiplicities *> & VariantClusterGraphKmerPath<kmer_size>::getKmerVertexMultiplicities() {

	return kmers;
}


template<uchar kmer_size>
const vector<const VariantClusterGraphVertex *> & VariantClusterGraphKmerPath<kmer_size>::getPath() {

	return path;
}


template<uchar kmer_size>
const unordered_map<pair<ushort, ushort>, const VariantClusterGraphVertex *, boost::hash<pair<ushort, ushort> > > & VariantClusterGraphKmerPath<kmer_size>::getVariants() {

	return variants;
}


template<uchar kmer_size>
uint VariantClusterGraphKmerPath<kmer_size>::getNumberOfVariantKmers(const ushort variant_idx) {

	auto processed_variants_it = processed_variants.find(variant_idx);
	auto running_variants_it = running_variants.find(variant_idx);

	if (processed_variants_it != processed_variants.end()) {

		assert(running_variants_it == running_variants.end());
		return processed_variants_it->second.score_counts.second;
	}

	assert(running_variants_it != running_variants.end());
	return running_variants_it->second.score_counts.second;
}


template<uchar kmer_size>
uint VariantClusterGraphKmerPath<kmer_size>::getNumberOfReferenceAlleles() {

	return num_reference_alleles;
}


template<uchar kmer_size>
typename VariantClusterGraphKmerPath<kmer_size>::KmerMultiplicities * VariantClusterGraphKmerPath<kmer_size>::newKmerVertex(const uint vertex_size) {

	kmers.emplace_back(new KmerMultiplicities());	
	kmers.back()->reserve(vertex_size);

	return kmers.back();
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::addKmerVertex(KmerMultiplicities * kmers_in) {
	
	kmers.push_back(kmers_in);	
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::newKmer(KmerCounts * kmer_counts, const ushort sample_idx, const bool is_shared) {

	updateScore(kmer_counts, sample_idx);

	if (!is_shared) {

	    kmers.back()->emplace_back(kmer_pair.getLexicographicalLowestKmer(), KmerInfo(kmer_counts));
	
	} else {

		assert(kmers.back()->back().first == kmer_pair.getLexicographicalLowestKmer());
	}
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::updateScore(KmerCounts * kmer_counts, const ushort sample_idx) {

	bool has_positive_count = false;

	if (is_dummy) {

		has_positive_count = true;

	} else {

		if (kmer_counts) {

			has_positive_count = (kmer_counts->getCount(sample_idx) > 0);
		}
	}

	auto vit = running_variants.begin();

    while (vit != running_variants.end()) {

		if (vit->second.last_nucleotide <= num_nucleotides) {

			assert(processed_variants.insert(*vit).second);
			vit = running_variants.erase(vit);

			continue;
		}

		vit->second.score_counts.first += has_positive_count;
		vit->second.score_counts.second++;

		vit++;
	}
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::incrementNucleotide() {

	num_nucleotides++;
}


template<uchar kmer_size>
bool VariantClusterGraphKmerPath<kmer_size>::isDummy() {

	return is_dummy;
}


template<uchar kmer_size>
pair<double, uint> VariantClusterGraphKmerPath<kmer_size>::getScore() {

	return calculatePathScore(unordered_set<const VariantClusterGraphVertex *>());
}



template<uchar kmer_size>
pair<double, uint> VariantClusterGraphKmerPath<kmer_size>::getScore(const unordered_set<const VariantClusterGraphVertex *> & covered_vertices) {

	return calculatePathScore(covered_vertices);
}


template<uchar kmer_size>
pair<double, uint> VariantClusterGraphKmerPath<kmer_size>::calculatePathScore(const unordered_set<const VariantClusterGraphVertex *> & covered_vertices) {

	assert(!(processed_variants.empty()) or !(running_variants.empty()));

	double score = 0;

	for (auto &vit: processed_variants) {

		score += calculateVariantScore(vit.second.score_counts);
	}

	for (auto &vit: running_variants) {

		score += calculateVariantScore(vit.second.score_counts);
	}

	uint num_unique_covered = 0;

	for (auto &var: variants) {

		if (covered_vertices.count(var.second) < 1) {

			if (is_dummy) {

				if (var.first.second == 0) {

					num_unique_covered += variants.size();
				
				} else {

					num_unique_covered++;
				}

			} else {

				num_unique_covered++;
			}		
		}
	}

	return make_pair(score/(processed_variants.size() + running_variants.size()), num_unique_covered);
}


template<uchar kmer_size>
double VariantClusterGraphKmerPath<kmer_size>::calculateVariantScore(pair<uint, uint> & score_counts) {

	assert(score_counts.first <= score_counts.second);

	if (score_counts.second > 0) {
		
		return score_counts.first/static_cast<double>(score_counts.second);

	} else {

		return 1;
	}	
}


template<uchar kmer_size>
void VariantClusterGraphKmerPath<kmer_size>::addCoveredVertices(unordered_set<const VariantClusterGraphVertex *> * covered_vertices) {

	for (auto &var: variants) {

		covered_vertices->insert(var.second);
	}
}


template<uchar kmer_size>
VariantClusterGraphFullKmerPath<kmer_size>::VariantClusterGraphFullKmerPath(const uint expected_path_length, const bool is_dummy_in) : VariantClusterGraphKmerPath<kmer_size>(expected_path_length, is_dummy_in) {

	num_multicluster_kmers = 0;
}


template<uchar kmer_size>
void VariantClusterGraphFullKmerPath<kmer_size>::addNestedVariantClusterIndices(VariantClusterGraphVertex & vertex) {

	if (vertex.isComplex()) {

		assert(!(vertex.variant_ids.empty()));

		ComplexVariantClusterGraphVertex * complex_vertex = static_cast<ComplexVariantClusterGraphVertex * >(&vertex);
		assert(!(complex_vertex->nested_variant_cluster_indices.empty()));

		nested_variant_cluster_indices.insert(complex_vertex->nested_variant_cluster_indices.begin(), complex_vertex->nested_variant_cluster_indices.end());
	}
}


template<uchar kmer_size>
const unordered_set<uint> & VariantClusterGraphFullKmerPath<kmer_size>::getNestedVariantClusterIndices() {

	return nested_variant_cluster_indices;
}


template<uchar kmer_size>
uint VariantClusterGraphFullKmerPath<kmer_size>::getNumberOfMultiClusterKmers() {

	return num_multicluster_kmers;
}


template<uchar kmer_size>
void VariantClusterGraphFullKmerPath<kmer_size>::newKmer(KmerCounts * kmer_counts, const ushort sample_idx, const bool is_shared) {

	VariantClusterGraphKmerPath<kmer_size>::updateScore(kmer_counts, sample_idx);

	bool add_kmer = true;

	if (kmer_counts) {

		if (kmer_counts->isExcluded()) {

			add_kmer = false;
		
		} else if (kmer_counts->isMulti()) { 

			num_multicluster_kmers++;
		}
	}

	if (add_kmer) {

		if (!is_shared) {

		    VariantClusterGraphKmerPath<kmer_size>::kmers.back()->emplace_back(VariantClusterGraphKmerPath<kmer_size>::kmer_pair.getLexicographicalLowestKmer(), KmerInfo(kmer_counts));

		    KmerInfo * cur_kmer_info = &(VariantClusterGraphKmerPath<kmer_size>::kmers.back()->back().second);
		    cur_kmer_info->variant_ids.reserve(VariantClusterGraphKmerPath<kmer_size>::running_variants.size());

	    	for (auto &variant: VariantClusterGraphKmerPath<kmer_size>::running_variants) {

				assert(variant.second.last_nucleotide > VariantClusterGraphKmerPath<kmer_size>::num_nucleotides);
				cur_kmer_info->variant_ids.emplace_back(variant.first);
			}

		} else {

			assert(VariantClusterGraphKmerPath<kmer_size>::kmers.back()->back().first == VariantClusterGraphKmerPath<kmer_size>::kmer_pair.getLexicographicalLowestKmer());
		    assert(VariantClusterGraphKmerPath<kmer_size>::kmers.back()->back().second.kmer_counts == kmer_counts);
		}
	}
}
