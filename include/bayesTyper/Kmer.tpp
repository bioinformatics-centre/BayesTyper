
/*
Kmer.tpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


template<uchar size>
const bitset<size*2> & Kmer<size>::getKmer() {

	return bitmer;
}


template<uchar size>
KmerForward<size>::KmerForward() : end_position(size * 2 - 2) {

	cur_position = 0;
	is_complete = false;
}

template<uchar size>
bool KmerForward<size>::move(bitset<2> nt_bits) {

	if (is_complete) {

		(this->bitmer)>>=2;
		update(nt_bits);

	} else {

		if (cur_position != end_position) {	

			update(nt_bits);			
			cur_position += 2;
		
		} else {

			update(nt_bits);			
			is_complete = true;
		} 
	}

	return is_complete;		
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
bool KmerReverseComplement<size>::move(bitset<2> nt_bits) {

	if (is_complete) {

		(this->bitmer)<<=2;
		update(nt_bits);

	} else {

		if (cur_position != end_position) {	

			update(nt_bits);			
			cur_position -= 2;
		
		} else {

			update(nt_bits);			
			is_complete = true;
		} 
	}

	return is_complete;		
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
bool KmerPair<size>::move(bitset<2> nt_bits) {

	auto is_complete_fw = KmerForward<size>::move(nt_bits);
	auto is_complete_rc = KmerReverseComplement<size>::move(nt_bits);

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
bool KmerPair<size>::isComplete() {

	assert(KmerForward<size>::is_complete == KmerReverseComplement<size>::is_complete);
	return KmerForward<size>::is_complete;
}

template<uchar size>
const bitset<size*2> & KmerPair<size>::getForwardKmer() {

	return KmerForward<size>::getKmer();
}

template<uchar size>
const bitset<size*2> & KmerPair<size>::getReverseComplementKmer() {

	return KmerReverseComplement<size>::getKmer();
}

template<uchar size>
const bitset<size*2> & KmerPair<size>::getLexicographicalLowestKmer() {

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


