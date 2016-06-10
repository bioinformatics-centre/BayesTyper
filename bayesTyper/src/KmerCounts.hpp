
/*
KmerCounts.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef __bayesTyper__KmerCounts_hpp
#define __bayesTyper__KmerCounts_hpp

#include <vector>
#include <mutex>
#include <random>

#include "Utils.hpp"

using namespace std;

class KmerCounts {

    public:

        KmerCounts();
        virtual ~KmerCounts() {};

        bool hasMaxCount();
        bool hasClusterOccurrence();
        bool hasMulticlusterOccurrence();
        bool hasIncompleteClusterOccurrence();

        bool hasDecoyOccurrence();
        bool hasMaxMultiplicity();
        bool hasMulticlusterGroupOccurrence();
        bool hasConstantMultiplicity();

        bool isExcluded();

        virtual uchar getCount(const ushort) = 0;
        virtual uchar getInterclusterMultiplicity(const Utils::Sex) = 0;
        
        virtual bool isMulti() = 0;
        virtual bool isEmpty() = 0;

    protected:

        uchar has_max_count : 1, has_cluster_occ : 1, has_multicluster_occ : 1, has_incomplete_cluster_occ : 1, has_decoy_occ : 1, has_max_multiplicity : 1, has_multicluster_group_occ : 1, has_constant_multiplicity : 1;
};


class SampleKmerCounts : public KmerCounts {

    public:

        SampleKmerCounts(const ushort);
        ~SampleKmerCounts();

        void addCount(const ushort, const uchar);
        void addClusterMultiplicity(const uchar, const bool, const uint);
        void addInterclusterMultiplicity(const Utils::ChromosomeClass);

        uchar getCount(const ushort);
        uchar getInterclusterMultiplicity(const Utils::Sex);

        bool isMulti();
        bool isEmpty();

    private:

        uchar * counts;

        uchar max_cluster_multiplicity;
        uchar male_intercluster_multiplicity;
        uchar female_intercluster_multiplicity;

        uint variant_cluster_group_index;

        uchar addMultiplicity(const uchar, const uchar);
};


class MultiClusterKmerCounts : public KmerCounts {

    public:

        MultiClusterKmerCounts(const ushort, SampleKmerCounts *);
        ~MultiClusterKmerCounts();

        void reset(const ushort);

        uchar getMulticlusterMultiplicity(const ushort);
        void reduceMulticlusterMultiplicity(const ushort, const uchar);
        void addMulticlusterMultiplicity(const ushort, const uchar);

        uchar getCount(const ushort);
        uchar getInterclusterMultiplicity(const Utils::Sex);

        void setIndex(const ushort, const uint);
        uint getIndex(const ushort);

        bool isMulti();
        bool isEmpty();

    private:

        SampleKmerCounts * sample_kmer_counts;
        
        uchar * multiplicities;
        uint * indices;
};


class EmptyKmerCounts : public KmerCounts {

    public:

        EmptyKmerCounts();

        uchar getCount(const ushort);
        uchar getInterclusterMultiplicity(const Utils::Sex);

        bool isMulti();
        bool isEmpty();
};


#endif