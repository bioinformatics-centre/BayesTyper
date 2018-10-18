
/*
KmerCounts.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include "Utils.hpp"

class KmerCounts {

    public:

        KmerCounts();
        virtual ~KmerCounts() {};

        bool hasClusterOccurrence();
        bool hasMulticlusterOccurrence();
        bool hasMultigroupOccurrence();

        bool hasDecoyOccurrence();
        bool hasMaxMultiplicity();

        bool isParameter();
        void isParameter(const bool);

        bool isExcluded();

        void addInterclusterMultiplicity(const bool, const vector<Utils::Ploidy> &);
        uchar getInterclusterMultiplicity(const Utils::Gender);
        uchar getMaxInterclusterMultiplicity();

        void addClusterMultiplicity(const uchar, const bool);

        virtual void addSampleCount(const ushort, const uchar) = 0;
        virtual uchar getSampleCount(const ushort) = 0;

        virtual void resetSampleMultiplicity() = 0;
        virtual uchar getSampleMultiplicity(const ushort) = 0;
        virtual void reduceSampleMultiplicity(const ushort, const uchar) = 0;
        virtual void addSampleMultiplicity(const ushort, const uchar) = 0;
        
    protected:

        uchar has_cluster_occ : 1, has_multicluster_occ : 1, has_multigroup_occ : 1, has_decoy_occ : 1, has_max_multiplicity : 1, is_parameter : 1;

        uchar max_haploid_multiplicity;

        uchar female_intercluster_multiplicity;
        uchar male_intercluster_multiplicity;

        uchar updateMultiplicity(const uchar, const uchar);
};


template <uchar sample_bin>
class ObservedKmerCounts : public KmerCounts {

    public:

        ObservedKmerCounts();
        ~ObservedKmerCounts() {};

        void addSampleCount(const ushort, const uchar);
        uchar getSampleCount(const ushort);

        void resetSampleMultiplicity();
        uchar getSampleMultiplicity(const ushort);
        void reduceSampleMultiplicity(const ushort, const uchar);
        void addSampleMultiplicity(const ushort, const uchar);
        
    private:

        uchar counts[sample_bin] = {};
        uchar multiplicities[sample_bin] = {};
};


#endif