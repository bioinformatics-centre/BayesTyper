


template<uchar size>
pair<bool, bool> Kmer<size>::move(const char nt) {

	auto nt_bits = DNAtoBits(nt);

	if (!nt_bits.first) {

		this->reset();
		
		return pair<bool, bool>(false, false);
	}

	return pair<bool, bool>(true, this->shift(nt_bits.second));
}


template<uchar size>
pair<bool, bitset<2> > Kmer<size>::DNAtoBits(const char nt) {

	bitset<2> nt_bits;

    if ((nt == 'A') or (nt == 'a')) {

        nt_bits.set(0, 0);
        nt_bits.set(1, 0);

    } else if ((nt == 'C') or (nt == 'c')) {

        nt_bits.set(0, 1);
        nt_bits.set(1, 0);
    
    } else if ((nt == 'G') or (nt == 'g')) {

        nt_bits.set(0, 0);
        nt_bits.set(1, 1);
    
    } else if ((nt == 'T') or (nt == 't')) {

        nt_bits.set(0, 1);
        nt_bits.set(1, 1);
    
    } else {

        return pair<bool, bitset<2> >(false, nt_bits);
    }

    return pair<bool, bitset<2> >(true, nt_bits);
}


template<uchar size>
bitset<size*2> & Kmer<size>::getKmer() {

	return bitmer;
}


template<uchar size>
KmerForward<size>::KmerForward() : end_position(size * 2 - 2) {

	cur_position = 0;
	is_complete = false;
}


template<uchar size>
bool KmerForward<size>::shift(bitset<2> nt_bits) {

	if (is_complete) {

		(this->bitmer)>>=2;
		update(nt_bits);

		return is_complete;

	} else {

		if (cur_position != end_position) {	

			update(nt_bits);			
			cur_position += 2;

			return is_complete;
		
		} else {

			update(nt_bits);			
			is_complete = true;

			return is_complete;		
		} 
	}
}


template<uchar size>
void KmerForward<size>::update(bitset<2> nt_bits) {

    this->bitmer.set(cur_position, nt_bits[0]);
	this->bitmer.set(cur_position + 1, nt_bits[1]);
}


template<uchar size>
void KmerForward<size>::shrink(const uchar new_size) {
	
	if (is_complete) {
		
		assert(cur_position == end_position);
		(this->bitmer)>>=2;		
	}

	assert(cur_position >= new_size * 2);
	uchar shift = cur_position - new_size * 2 + 2;

	(this->bitmer) >>= shift;		
	cur_position -= shift;

	is_complete = false;
}


template<uchar size>
void KmerForward<size>::reset() {

	cur_position = 0;
	is_complete = false;
}


template<uchar size>
KmerReverseComplement<size>::KmerReverseComplement() : end_position(0) {

	cur_position = size * 2 - 2;
	is_complete = false;
}


template<uchar size>
bool KmerReverseComplement<size>::shift(bitset<2> nt_bits) {

	if (is_complete) {

		(this->bitmer)<<=2;
		update(nt_bits);

		return is_complete;

	} else {

		if (cur_position != end_position) {	

			update(nt_bits);			
			cur_position -= 2;

			return is_complete;
		
		} else {

			update(nt_bits);			
			is_complete = true;

			return is_complete;		
		} 
	}
}


template<uchar size>
void KmerReverseComplement<size>::update(bitset<2> nt_bits) {

    this->bitmer.set(cur_position, ~nt_bits[0]);
	this->bitmer.set(cur_position + 1, ~nt_bits[1]);
}


template<uchar size>
void KmerReverseComplement<size>::shrink(const uchar new_size) {

	if (is_complete) {
		
		assert(cur_position == end_position);
		(this->bitmer)<<=2;		
	}

	assert((size * 2 - 2 - cur_position) >= new_size * 2);
	uchar shift = size * 2 - 2 - cur_position - new_size * 2 + 2;

	(this->bitmer) <<= shift;		
	cur_position += shift;

	is_complete = false;
}


template<uchar size>
void KmerReverseComplement<size>::reset() {

	cur_position = size * 2 - 2;
	is_complete = false;
}


template<uchar size>
pair<bool, bool> KmerPair<size>::move(const char nt) {

	return KmerForward<size>::move(nt);
}


template<uchar size>
bool KmerPair<size>::shift(bitset<2> nt_bits) {

	auto is_complete_fw = KmerForward<size>::shift(nt_bits);
	auto is_complete_rc = KmerReverseComplement<size>::shift(nt_bits);

	assert(is_complete_fw == is_complete_rc);

	return is_complete_fw;
}


template<uchar size>
void KmerPair<size>::shrink(const uchar new_size) {

	KmerForward<size>::shrink(new_size);
	KmerReverseComplement<size>::shrink(new_size);
}


template<uchar size>
void KmerPair<size>::reset() {

	KmerForward<size>::reset();
	KmerReverseComplement<size>::reset();
}


template<uchar size>
bitset<size*2> & KmerPair<size>::getForwardKmer() {

	return KmerForward<size>::getKmer();
}


template<uchar size>
bitset<size*2> & KmerPair<size>::getReverseComplementKmer() {

	return KmerReverseComplement<size>::getKmer();
}


template<uchar size>
bitset<size*2> & KmerPair<size>::getLexicographicalLowestKmer() {

	uchar start_position = 0;
	uchar end_position = size * 2 - 2;

	while (start_position <= end_position) {

		if (KmerForward<size>::getKmer()[start_position + 1] < KmerReverseComplement<size>::getKmer()[start_position + 1]) {

			return KmerForward<size>::getKmer();
		
		} else if (KmerForward<size>::getKmer()[start_position + 1] > KmerReverseComplement<size>::getKmer()[start_position + 1]) {

			return KmerReverseComplement<size>::getKmer();
		} 

		if (KmerForward<size>::getKmer()[start_position] < KmerReverseComplement<size>::getKmer()[start_position]) {

			return KmerForward<size>::getKmer();
		
		} else if (KmerForward<size>::getKmer()[start_position] > KmerReverseComplement<size>::getKmer()[start_position]) {

			return KmerReverseComplement<size>::getKmer();
		} 

		start_position += 2;
	}

	return KmerForward<size>::getKmer();
}


template<uchar size>
bool KmerPair<size>::operator == (const KmerPair<size> & rhs) const {

	if (this->KmerForward<size>::cur_position != rhs.KmerForward<size>::cur_position) {

		return false;
	}

	uchar cur_end_position = this->KmerForward<size>::cur_position;

	if (this->KmerForward<size>::is_complete) {

		assert(cur_end_position == this->KmerForward<size>::end_position);
		cur_end_position += 2;
	}

	for (ushort i = 0; i < cur_end_position; i++) {

		if (this->KmerForward<size>::bitmer[i] != rhs.KmerForward<size>::bitmer[i]) {

			return false;
		}
	}

	return true;
}


