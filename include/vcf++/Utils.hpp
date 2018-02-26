
/*
Utils.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef UTILS
#define UTILS

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>
#include <limits>
#include <cmath>

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Utils {

  using namespace std;

  static const double double_precision = numeric_limits<double>::epsilon();
  static const float float_precision = numeric_limits<float>::epsilon();

  static const char char_overflow = numeric_limits<char>::max();
  static const uchar uchar_overflow = numeric_limits<uchar>::max();
  static const ushort ushort_overflow = numeric_limits<ushort>::max();
  static const uint uint_overflow = numeric_limits<uint>::max();  
  static const ulong ulong_overflow = numeric_limits<ulong>::max();

  inline bool doubleCompare(const double a, const double b) {

      return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision * 100));
  }

  inline bool floatCompare(const float a, const float b) {

      return ((a == b) or (abs(a - b) < abs(min(a, b)) * float_precision * 100));
  }

  string boolToString(const bool);
  string floatToString(const float, const uint);

  vector<string> splitString(const string &, const char);
  vector<string> splitStringEmptyIgnore(const string &, const char);
  vector<string> splitString(const string &, const char, const uint);

  // Vector flush capability
  template < typename T >
  ostream& operator << (ostream& os, const vector<T>& v) {

      for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {

  		os << scientific << setprecision(6) << *it << "  ";
  	}

      return os;
  }

  // 2D-vector flush capability
  template < typename T >
  ostream& operator << (ostream& os, const vector<vector<T> >& v) {

      for (typename vector<vector<T> >::const_iterator it = v.begin(); it != v.end(); it++) {

          os << *it << "\n";
  	}

      return os;
  }


  inline string getLocalTime () {

      time_t now = time(NULL);
      struct tm * lt;
      lt = localtime (&now);

      stringstream output_string;

      if (lt->tm_mday < 10) {

          output_string << "0" << lt->tm_mday << "/";

      } else {

          output_string << lt->tm_mday << "/";

      }

      if ((1 + lt->tm_mon) < 10) {

          output_string << "0" << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";

      } else {

          output_string << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";

      }

      if (lt->tm_hour < 10) {

          output_string << "0" << lt->tm_hour << ":";

      } else {

          output_string << lt->tm_hour << ":";

      }

      if (lt->tm_min < 10) {

          output_string << "0" << lt->tm_min << ":";

      } else {

          output_string << lt->tm_min << ":";

      }

      if (lt->tm_sec < 10) {

          output_string << "0" << lt->tm_sec;

      } else {

          output_string << lt->tm_sec;

      }

      return output_string.str();

  }

};

#endif
