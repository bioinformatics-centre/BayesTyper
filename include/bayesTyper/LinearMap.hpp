
/*
LinearMap.hpp - This file is part of BayesTyper (v1.1)


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


#ifndef __bayesTyper__LinearMap_hpp
#define __bayesTyper__LinearMap_hpp

#include <vector>
#include <algorithm>
#include <functional>
#include <iterator>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

template<typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType> >
class LinearMap {

	public: 

		LinearMap();

		void reserve(const uint);
		bool empty();
		uint size();

		typedef typename std::pair<KeyType,ValueType> MapElement;
		typedef typename std::vector<MapElement>::iterator iterator;
		
		iterator erase(iterator);
		
		iterator begin();
		iterator end();

	protected:

		struct PairComparator {

			KeyComparator keycomparator;
			
			bool operator() (const MapElement & lhs, const KeyType & rhs) const {

				return keycomparator(lhs.first, rhs);				
		  	}	

			bool operator() (const KeyType & lhs, const MapElement & rhs) const {

				return keycomparator(lhs, rhs.first);				
		  	}

			bool operator() (const MapElement & lhs, const MapElement & rhs) const {

				return keycomparator(lhs.first, rhs.first);				
		  	}
		};

		std::vector<std::pair<KeyType,ValueType> > map;
};


template<typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType> >
class SortedLinearMap : public LinearMap<KeyType, ValueType, KeyComparator> {

	public: 

		SortedLinearMap();

		std::pair<typename LinearMap<KeyType, ValueType, KeyComparator>::iterator, bool> insert(std::pair<KeyType,ValueType>);
		typename LinearMap<KeyType, ValueType, KeyComparator>::iterator find(KeyType);
};


template<typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType> >
class PartialSortedLinearMap : public LinearMap<KeyType, ValueType, KeyComparator> {

	public: 

		PartialSortedLinearMap();

		std::pair<typename LinearMap<KeyType, ValueType, KeyComparator>::iterator, bool> insert(std::pair<KeyType,ValueType>, const bool);
		typename LinearMap<KeyType, ValueType, KeyComparator>::iterator find(KeyType);
		void sort();


	private:

		uint sorted_segment_length;
};


#include "LinearMap.tpp"

#endif