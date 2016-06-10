
/*
remove.cpp - This file is part of BayesTyper (v0.9)


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

#include "vcf++/Regions.hpp"

#include "remove.hpp"

namespace Remove {

	void remove(const string & in_vcf_filename, const string & regions_string, const string & output_prefix) {

        cout << "[" << Utils::getLocalTime() << "] Running BayesTyperUtils (" << BTU_VERSION << ") remove ...\n" << endl;

        assert(in_vcf_filename != (output_prefix + ".vcf"));

        GenotypedVcfFileReader vcf_reader(in_vcf_filename, true);
        VcfFileWriter vcf_writer(output_prefix + ".vcf", vcf_reader.metaData(), true);

        Regions regions = Regions(regions_string);

        Variant * cur_var;
        uint num_variants = 0;
        uint num_removed_variants = 0;

        while (vcf_reader.getNextVariant(&cur_var)) {

            num_variants++;

            if (regions.isNotIn(*cur_var)) {

	            vcf_writer.write(cur_var);
            
            } else {

            	num_removed_variants++;
            }

            if (num_variants % 100000 == 0) {

                cout << "[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variant(s)" << endl;
            }

            delete cur_var;
        }

		cout << "\n[" << Utils::getLocalTime() << "] Parsed " << num_variants << " variant(s):\n" << endl;
		cout << "\t\t- Removed " << num_removed_variants << " variant(s)" << endl;
		cout << endl;
    }
}
