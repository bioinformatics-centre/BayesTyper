
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


#ifndef __bayesTyper__Utils_hpp
#define __bayesTyper__Utils_hpp

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <atomic>
#include <limits>
#include <bitset>
#include <math.h>
#include <list>
#include <assert.h>
#include <sys/resource.h>
#include <algorithm>

#include "Eigen/Dense"

using namespace std; 

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Utils {

    static const uint kmer_size = BT_KMER_SIZE;

    typedef Eigen::Matrix<uchar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> MatrixXuchar;
    typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixXbool;
    typedef Eigen::Matrix<bool,1,Eigen::Dynamic,Eigen::RowMajor> RowVectorXbool;
    typedef Eigen::Matrix<bool,Eigen::Dynamic,1> ColVectorXbool;

    static const double double_precision = numeric_limits<double>::epsilon();
    static const float float_precision = numeric_limits<float>::epsilon();

    static const uchar uchar_overflow = numeric_limits<uchar>::max();
    static const ushort ushort_overflow = numeric_limits<ushort>::max();
    static const uint uint_overflow = numeric_limits<uint>::max();  
    static const ulong ulong_overflow = numeric_limits<ulong>::max();

    static const uchar bit7_overflow = 127;

    static const ushort queue_size_thread_scaling = 2;

    enum class ChromClass : uchar {Autosomal = 0, X, Y, Decoy, CHROM_CLASS_SIZE};
    enum class Ploidy : uchar {Null = 0, Haploid, Diploid, PLOIDY_SIZE};
    enum class Gender : uchar {Male = 0, Female, GENDER_SIZE};


    inline bool doubleCompare(const double a, const double b) {

        assert(std::isfinite(a));
        assert(std::isfinite(b));

        return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision * 100));
    }

    inline bool floatCompare(const float a, const float b) {

        assert(std::isfinite(a));
        assert(std::isfinite(b));

        return ((a == b) or (abs(a - b) < abs(min(a, b)) * float_precision * 100));
    }

    inline bool floatLess(const float a, const float b) {

        assert(std::isfinite(a));
        assert(std::isfinite(b));

        return ((a < b) and !(floatCompare(a, b)));
    }

    inline double logAddition(const double log_summand1, const double log_summand2) {

        assert(std::isfinite(log_summand1));
        assert(std::isfinite(log_summand2));

        if (log_summand1 < log_summand2) {

            double log_sum = log_summand2 + log1p(exp(log_summand1 - log_summand2));
            assert(std::isfinite(log_sum));

            return log_sum;                
                        
        } else {

            double log_sum = log_summand1 + log1p(exp(log_summand2 - log_summand1)); 
            assert(std::isfinite(log_sum));

            return log_sum;
        }
    }
 
    inline string getMaxMemoryUsage() {

        stringstream output_string;

        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        output_string << "Maximum resident set size: " << usage.ru_maxrss/float(1e6) << " Gb";
        
        return output_string.str();
    }

    inline string getLocalTime () {
    
        time_t now = time(nullptr);
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
}


// pairs flush capability
template<typename T, typename T2>
std::ostream & operator << (std::ostream & os, const std::pair<T, T2> & std_pair) {
            
    os << "(" << to_string(std_pair.first) << ", " << to_string(std_pair.second) << ")";
    
    return os;
}

// Vector flush capability
template<typename T>
std::ostream & operator << (std::ostream & os, const std::vector<T> & std_vector) {
    
    for (typename std::vector<T>::const_iterator it = std_vector.cbegin(); it != std_vector.cend(); ++it) {
        
        if (it != std_vector.cbegin()) {

            os << "\t";
        }

		os << *it;
	}

    return os;
}

// 2D-vector flush capability
template<typename T>
std::ostream & operator << (std::ostream & os, const std::vector<std::vector<T> > & std_vector) {
    
    for (typename std::vector<std::vector<T> >::const_iterator it = std_vector.cbegin(); it != std_vector.cend(); ++it) {
		
		os << *it << "\n";
	}
	
    return os;
}

#endif