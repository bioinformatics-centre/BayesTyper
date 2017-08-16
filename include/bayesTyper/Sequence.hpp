
/*
Sequence.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef __bayesTyper__Sequence_hpp
#define __bayesTyper__Sequence_hpp

#include <bitset>

#include "Utils.hpp"

namespace Sequence {

	template<uchar size = 1>
	pair<bitset<size * 2>, bool> ntToBit(const char nt) {

		bitset<size * 2> bits;

	    if ((nt == 'A') or (nt == 'a')) {

	        bits.set(0, 0);
	        bits.set(1, 0);

	    } else if ((nt == 'C') or (nt == 'c')) {

	        bits.set(0, 1);
	        bits.set(1, 0);
	    
	    } else if ((nt == 'G') or (nt == 'g')) {

	        bits.set(0, 0);
	        bits.set(1, 1);
	    
	    } else if ((nt == 'T') or (nt == 't')) {

	        bits.set(0, 1);
	        bits.set(1, 1);
	    
	    } else {

	        return make_pair(bits, false);
	    }

	    return make_pair(bits, true);
	}

	template<uchar size>
	pair<bitset<size * 2>, bool> ntToBit(const string & nt_seq) {

		assert(nt_seq.size() == size);
		bitset<size * 2> bit_seq;

		for (ushort i = 0; i < (size * 2); i += 2) {

			auto bits = ntToBit<1>(nt_seq.at(i));

			if (bits.second) {

		    	bit_seq.set(i, bits.first[0]);
		    	bit_seq.set(i + 1, bits.first[1]);
			
			} else {

		        return make_pair(bit_seq, false);
			}
		}

	    return make_pair(bit_seq, true);
	}

	template<uchar size>
	string bitToNt(const bitset<size * 2> & bit_seq) {

		string nt_seq;
		nt_seq.reserve(size);

		for (ushort i = 0; i < (size * 2); i += 2) {

			if (bit_seq[i] == bit_seq[i + 1]) {

				if (bit_seq[i] == 0) {

					nt_seq += "A";
				
				} else {

					nt_seq += "T";
				}				
			
			} else {

				if (bit_seq[i] == 0) {

					nt_seq += "G";
				
				} else {

					nt_seq += "C";
				}					
			}
		}

	    return nt_seq;
	}

	template<uchar size>
	uchar gcBiasBin(const bitset<size * 2> & bit_seq, const uchar num_bins) {

		assert(num_bins >= 1);
		assert(num_bins <= size);

		if (num_bins == 1) {

			return 0;
		}

		uchar gc_count = 0;

		for (ushort i = 0; i < (size * 2); i += 2) {

			if (bit_seq[i] != bit_seq[i + 1]) {

				gc_count++;
			}
		}

		assert(gc_count <= size);

		if (gc_count == 0) {

			return 0;

		} else if (gc_count == size) {

			return num_bins - 1;

		} else {

			return ceil((gc_count * num_bins)/static_cast<float>(size)) - 1;
		}
	}
}

#endif