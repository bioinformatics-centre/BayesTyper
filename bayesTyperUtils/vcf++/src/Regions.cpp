
/*
Regions.cpp - This file is part of BayesTyper (v0.9)


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
#include <assert.h>

#include "Regions.hpp"

Regions::Regions(const string & regions_string) {

    assert(!(regions_string.empty()));

    auto regions_info_split = Utils::splitString(regions_string, ':');

    for (auto &region_info_string: regions_info_split) {

        auto region_info = Utils::splitString(region_info_string, ',');

        if (region_info.size() == 1) {

            assert(regions.emplace(region_info.front(), vector<pair<uint, uint> >()).second);
        
        } else {

            assert(region_info.size() == 3);
            auto regions_emplace = regions.emplace(region_info.front(), std::vector<std::pair<uint, uint> >());

            assert(regions_emplace.second == regions_emplace.first->second.empty());
            regions_emplace.first->second.emplace_back(stoi(region_info.at(1)), stoi(region_info.at(2)));     
        }
    }
}

bool Regions::isIn(const Variant & cur_var) const {

    auto regions_it = regions.find(cur_var.chrom());

    if (regions_it != regions.end()) {

         if (!(regions_it->second.empty())) {

            bool is_in_region = false;

            for (auto &region: regions_it->second) {

                if ((cur_var.pos() >= region.first) and (cur_var.pos() <= region.second)) {

                    is_in_region = true;
                    break;
                }
            }

            return is_in_region;
         
         } else {

            return true;
         }
    
    } else {

        return false;
    } 
}

bool Regions::isNotIn(const Variant & cur_var) const {

    return !isIn(cur_var);
}


