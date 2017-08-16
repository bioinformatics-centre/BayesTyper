
/*
Regions.cpp - This file is part of BayesTyper (v1.1)


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


#include <string>
#include <vector>
#include <unordered_map>

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"

#include "Regions.hpp"
#include "Utils.hpp"

Regions::Regions(const string & regions_string) {

    if (!(regions_string.empty())) {

        vector<string> regions_info_split;
        split(regions_info_split, regions_string, boost::is_any_of(","));

        for (auto &region_info_string: regions_info_split) {

            vector<string> region_info_split;
            split(region_info_split, region_info_string, boost::is_any_of(":"));

            if (region_info_split.size() == 1) {

                assert(regions.emplace(region_info_split.front(), vector<pair<uint, uint> >()).second);
            
            } else {

                assert(region_info_split.size() == 2);

                auto regions_emplace = regions.emplace(region_info_split.front(), std::vector<std::pair<uint, uint> >());
                assert(regions_emplace.second == regions_emplace.first->second.empty());

                vector<std::string> region_positions_split;
                split(region_positions_split, region_info_split.back(), boost::is_any_of("-"));

                assert(region_positions_split.size() == 2);

                regions_emplace.first->second.emplace_back(stoi(region_positions_split.front()), stoi(region_positions_split.back()));     
            }
        }
    }
}

bool Regions::overlaps(const string & chrom, const uint start_pos, const uint end_pos) const {

    if (regions.empty()) {

        return true;
    }

    auto regions_it = regions.find(chrom);

    if (regions_it != regions.end()) {

         if (!(regions_it->second.empty())) {

            for (auto &region: regions_it->second) {

                if ((end_pos >= region.first) and (start_pos <= region.second)) {

                    return true;
                }
            }

            return false;
         
        } else {

            return true;
        }
    }

    return false;
}
