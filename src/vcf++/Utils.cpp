
/*
Utils.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#include "assert.h"
#include <sstream>
#include <iomanip>
#include <math.h>

#include "Utils.hpp"

namespace Utils {

  string boolToString(const bool in_bool) {

    if (in_bool) {

      return "1";

    } else {

      return "0";
    }
  }

  string floatToString(const float value, const uint precision) {

    stringstream float_ss;

    if (precision == 0) {

      float_ss << round(value);

    } else {

      float_ss << fixed << setprecision(precision) << value;
    }

    return float_ss.str();
  }

  vector<string> splitString(const string & str, const char delim) {

      stringstream ss(str);
      string item;
      vector<string> elems;
      elems.reserve(100);

      while (getline(ss, item, delim)) {

          elems.push_back(item);
      }

      return elems;
  }

  vector<string> splitStringEmptyIgnore(const string & str, const char delim) {

      stringstream ss(str);
      string item;
      vector<string> elems;
      elems.reserve(100);

      while (getline(ss, item, delim)) {

          if (!item.empty()) {

              elems.push_back(item);
          }
      }

      return elems;
  }

  vector<string> splitString(const string & str, const char delim, const uint max_elems) {

      assert(max_elems > 1);

      stringstream ss(str);
      string item;
      vector<string> elems;
      elems.reserve(max_elems);

      while (getline(ss, item, delim)) {

          elems.push_back(item);

          if (elems.size() == (max_elems - 1)) {

              getline(ss, item);
              elems.push_back(item);
              break;
          }
      }

      return elems;
  }
}
