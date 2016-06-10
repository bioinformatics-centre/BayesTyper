
/*
clean.cpp - This file is part of BayesTyper (v0.9)


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


#include <sstream>
#include <math.h>
#include <thread>
#include <map>

#include "vcf++/JoiningString.hpp"
#include "vcf++/Auxiliaries.hpp"

#include "clean.hpp"

namespace Clean {

	void clean(const string & in_vcf_filename, const string & output_prefix) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") clean ...\n" << endl;

        VcfFileReader vcf_reader(in_vcf_filename, false);

        vcf_reader.metaData().filterDescriptors().clear();
        vcf_reader.metaData().infoDescriptors().clear();
        vcf_reader.metaData().formatDescriptors().clear();
        vcf_reader.metaData().clearSamples();

        VcfFileWriter output_vcf(output_prefix + ".vcf", vcf_reader.metaData(), false);
        
        Variant * cur_var;

        uint num_variants = 0;
        uint num_alt_alleles = 0;

        while (vcf_reader.getNextVariant(&cur_var)) {

            num_variants++;
            num_alt_alleles += cur_var->numAlts();

            assert(!(cur_var->ref().isMissing()));

            vector<uint> missing_alt_indices;

            for (uint alt_idx = 0; alt_idx < cur_var->numAlts(); alt_idx++) {

                if (cur_var->alt(alt_idx).isMissing()) {

                    missing_alt_indices.push_back(alt_idx);
                }
            }

            cur_var->removeAlts(missing_alt_indices);
            cur_var->removeRedundantAlts();

            assert(cur_var->numAlls() > 1);            
            assert(cur_var->numAlts() > 0);
            
            Auxiliaries::rightTrimVariant(cur_var);
                
            cur_var->setIds({});
            cur_var->setQual({0, false});
            cur_var->setFilters({});

            output_vcf.write(cur_var);

            if ((num_variants % 100000) == 0) {

                std::cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variants" << endl;
            }

            delete cur_var;
        }

        cout << "\n[" << Utils::getLocalTime() << "] Cleaned " << num_variants << " variants and " << num_alt_alleles << " alternative alleles.\n" << endl;
    }
}
