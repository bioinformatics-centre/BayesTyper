
/*
PartitionGraph.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__PartitionGraph_hpp
#define __bayesTyper__PartitionGraph_hpp

#include <vector>

#include "boost/graph/adjacency_list.hpp"

#include "Utils.hpp"
#include "KmerCounts.hpp"

class PartitionGraph {

	private:
	
		typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, KmerCounts *> Graph;
		typedef boost::graph_traits<Graph>::vertex_descriptor VertexType;

		void depthFirstSearch(vector<uint> *, const uint, const VertexType &) const;

	public:

		Graph graph;

		PartitionGraph(const uint);
		
		bool addEdge(const uint, const uint, KmerCounts *);
		uint numberOfVertices() const;

		uint connectedComponents(vector<uint> * component_indices) const;
};

#endif