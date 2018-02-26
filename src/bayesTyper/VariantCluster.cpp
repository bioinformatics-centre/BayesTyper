
/*
VariantCluster.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include <map>
#include <unordered_set>
#include <math.h>

#include "VariantCluster.hpp"
#include "Utils.hpp"


void VariantCluster::mergeVariantClusters(const VariantCluster & in_variant_cluster) {

    assert(cluster_idx != in_variant_cluster.cluster_idx);
    assert(chromosome_id == in_variant_cluster.chromosome_id);
    assert(chromosome_name == in_variant_cluster.chromosome_name);
    assert(chromosome_class == in_variant_cluster.chromosome_class);
    assert(contained_clusters.empty() and in_variant_cluster.contained_clusters.empty());

    left_flank = min(left_flank, in_variant_cluster.left_flank);
    right_flank = max(right_flank, in_variant_cluster.right_flank);

    for (auto &vit: in_variant_cluster.variants) {

        assert(variants.insert(vit).second);
    }
}

