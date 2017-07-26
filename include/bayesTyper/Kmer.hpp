
/*
Kmer.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__Kmer_hpp
#define __bayesTyper__Kmer_hpp

#include <bitset>

#include "Utils.hpp"

template<uchar size>
class Kmer {

	public:

		Kmer() {};
		virtual ~Kmer() {}

		const bitset<size*2> & getKmer();

	protected:
		
		bitset<size*2> bitmer;
};
 

template<uchar size>
class KmerForward : public Kmer<size> {

	public:

		KmerForward();
		~KmerForward() {}

		bool move(bitset<2>);
		void shrink(const uchar);
		void reset();

	protected:

		uchar cur_position;
		const uchar end_position;
		bool is_complete;

	private:

		void update(bitset<2>);
};


template<uchar size>
class KmerReverseComplement : public Kmer<size> {
	
	public:

		KmerReverseComplement();
		~KmerReverseComplement() {}

		bool move(bitset<2>);
		void shrink(const uchar);
		void reset();	

	protected:

		uchar cur_position;
		const uchar end_position;
		bool is_complete;

	private:

		void update(bitset<2>);		
};


template<uchar size>
class KmerPair : public KmerForward<size>, public KmerReverseComplement<size> {

	public:

		bool operator == (const KmerPair<size> &) const;

		bool move(bitset<2>);
		void shrink(const uchar);
		void reset();

		const bitset<size*2> & getForwardKmer();
		const bitset<size*2> & getReverseComplementKmer();
		const bitset<size*2> & getLexicographicalLowestKmer();
};


#include "Kmer.tpp"


#endif