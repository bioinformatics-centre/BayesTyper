
/*
Combinator.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "Combinator.hpp"

vector<vector<uint> > Combinator::enumerateCombinations(uint cardinality, uint complexity) {

	vector<uint> start_path;
	vector<vector<uint> > all_paths;

	enumerationRecursion(0, cardinality, complexity, start_path, all_paths);

	return all_paths;
}

vector<vector<bool> > Combinator::enumerateBinaryCombinations(uint cardinality) {

	vector<bool> start_path;
	vector<vector<bool> > all_paths;

	binaryEnumerationRecursion(0, cardinality, start_path, all_paths);

	return all_paths;
}


void Combinator::enumerationRecursion(uint depth, uint cardinality, uint complexity, vector<uint> path, vector<vector<uint> > & all_paths) {

	uint local_depth = depth + 1;

	if (local_depth <= cardinality) {

		for (uint i = 0; i < complexity; i++) {

			vector<uint> local_path = path;
			local_path.push_back(i);
			enumerationRecursion(local_depth, cardinality, complexity, local_path, all_paths);
		}

	} else {

		all_paths.push_back(path);
	}
}

void Combinator::binaryEnumerationRecursion(uint depth, uint cardinality, vector<bool> path, vector<vector<bool> > & all_paths) {

	uint local_depth = depth + 1;

	if (local_depth <= cardinality) {

		for (uint i = 0; i < 2; i++) {

			vector<bool> local_path = path;

			if (i == 0) {

				local_path.push_back(false);

			} else {

				local_path.push_back(true);
			}

			binaryEnumerationRecursion(local_depth, cardinality, local_path, all_paths);
		}

	} else {

		all_paths.push_back(path);
	}
}
