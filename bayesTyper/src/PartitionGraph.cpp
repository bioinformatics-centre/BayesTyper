
/*
PartitionGraph.cpp - This file is part of BayesTyper (v0.9)


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


#include <vector>
#include <queue>

#include "PartitionGraph.hpp"
#include "Utils.hpp"


PartitionGraph::PartitionGraph(const uint num_vertices) : graph(num_vertices) {}


bool PartitionGraph::addEdge(const uint first_vertex_idx, const uint second_vertex_idx, KmerCounts * kmer_counts) {

    return boost::add_edge(first_vertex_idx, second_vertex_idx, kmer_counts, graph).second;
}


uint PartitionGraph::numberOfVertices() const {

    return num_vertices(graph);
}


uint PartitionGraph::connectedComponents(vector<uint> * component_indices) const {

    assert(component_indices->empty());
    *component_indices = vector<uint>(num_vertices(graph), Utils::uint_overflow);

    uint component_idx = 0;
    auto vit = boost::vertices(graph); 

    while (vit.first != vit.second) {

        if (component_indices->at(*vit.first) == Utils::uint_overflow) {

            depthFirstSearch(component_indices, component_idx, *vit.first);
            component_idx++;
        } 

        vit.first++;
    }

    return component_idx;
}


void PartitionGraph::depthFirstSearch(vector<uint> * component_indices, const uint component_idx, const VertexType & start_vertex) const {

    queue<VertexType> dfs_queue; 
    dfs_queue.push(start_vertex);

    while (!(dfs_queue.empty())) {

        auto cur_vertex = dfs_queue.front();
        dfs_queue.pop();

        if (component_indices->at(cur_vertex) == Utils::uint_overflow) {

            component_indices->at(cur_vertex) = component_idx;

            auto eit_out = boost::out_edges(cur_vertex, graph);
            
            while (eit_out.first != eit_out.second) {

                assert(graph[*eit_out.first]);
                assert(!(graph[*eit_out.first]->isEmpty()));

                if (!(graph[*eit_out.first]->isExcluded())) {

                    dfs_queue.push(boost::target(*eit_out.first, graph));
                }

                eit_out.first++;
            }           

        } else {

            assert(component_indices->at(cur_vertex) == component_idx);
        }
    }
}
