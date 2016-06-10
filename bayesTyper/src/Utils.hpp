
/*
Utils.hpp - This file is part of BayesTyper (v0.9)


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

#include "../Eigen/Dense"

#include "PerfectSet.hpp"
#include "HybridHash.hpp"

using namespace std; 

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Eigen {

    typedef Eigen::Matrix<uchar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> MatrixXuchar;
    typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> MatrixXbool;
    typedef Eigen::Matrix<uint,1,Eigen::Dynamic,Eigen::RowMajor> RowVectorXuint;
    typedef Eigen::Matrix<ushort,Eigen::Dynamic,1> ColVectorXushort;
    typedef Eigen::Matrix<uchar,1,Eigen::Dynamic,Eigen::RowMajor> RowVectorXuchar;
    typedef Eigen::Matrix<uchar,Eigen::Dynamic,1> ColVectorXuchar;
    typedef Eigen::Matrix<bool,1,Eigen::Dynamic,Eigen::RowMajor> RowVectorXbool;
    typedef Eigen::Matrix<bool,Eigen::Dynamic,1> ColVectorXbool;

}

namespace Utils {

    static const double double_underflow = numeric_limits<double>::min();
    static const double double_overflow = numeric_limits<double>::max();
    static const double double_precision = numeric_limits<double>::epsilon();
    static const double float_precision = numeric_limits<float>::epsilon();
    static const ulong ulong_overflow = numeric_limits<ulong>::max();
    static const uint uint_overflow = numeric_limits<uint>::max();
    static const ushort ushort_overflow = numeric_limits<ushort>::max();
    static const char char_overflow = numeric_limits<char>::max();
    static const uchar uchar_overflow = numeric_limits<uchar>::max();
   
    static const uchar bit7_overflow = 127;

    static const ushort queue_size_scaling = 3;
    static const uchar small_kmer_size = SMALLMERSIZE;
    
    typedef ThreadedPerfectSet<small_kmer_size * 2> SmallmerSet;

    enum class VariantType : uchar {SNP = 0, Insertion, Deletion, Complex, Mixture, Unsupported, VARIANT_TYPE_SIZE};
    static const vector<string> variant_type_strings = {"SNP", "Insertion", "Deletion", "Complex", "Mixture", "Unsupported"};

    enum class FilterStatus : uchar {PASS = 0, FILTER_STATUS_SIZE};
    static const vector<string> filter_status_strings = {"PASS"};

    enum class ChromosomeClass : uchar {Autosomal = 0, X, Y, Decoy, CHROMOSOME_CLASS_SIZE};
    static const vector<string> chromosome_class_strings = {"autosomal", "X", "Y", "decoy"};

    enum class Ploidy : uchar {Null = 0, Haploid, Diploid, PLOIDY_SIZE};
    enum class Sex : uchar {Male = 0, Female, SEX_SIZE};


    inline bool doubleCompare(const double a, const double b) {

        return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision * 100));
    }
 

    inline bool floatCompare(const float a, const float b) {

        return ((a == b) or (abs(a - b) < abs(min(a, b)) * float_precision * 100));
    }


    template <typename T>
    inline void hash_combine(size_t & seed, const T value) {

        hash<T> hasher;
        seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }


    inline double logAddition(const double log_summand1, const double log_summand2) {

        if (std::isfinite(log_summand1) and std::isfinite(log_summand2)) {

            if (log_summand1 < log_summand2) {

                double log_sum = log_summand2 + log(1 + exp(log_summand1 - log_summand2));
                assert(std::isfinite(log_sum));

                return log_sum;                
                            
            } else {

                double log_sum = log_summand1 + log(1 + exp(log_summand2 - log_summand1)); 
                assert(std::isfinite(log_sum));

                return log_sum;
            }

        } else if (std::isfinite(log_summand1)) {

            return log_summand1;
        
        } else if (std::isfinite(log_summand2)) {

            return log_summand2;

        } else {

            return -numeric_limits<double>::infinity();
        }
    }
 

    inline uint vectorBoolToUchar(const vector<bool> cur_vector) {

        uint integer_value = 0;

        for (uint i = 0; i < cur_vector.size(); i++) {

            if (cur_vector.at(i)) {
                
                integer_value += pow(2,(double) i);
            }
        }   

        return integer_value;
    }


    inline bool DNAtoBit(const char nt, bitset<2> * nt_bitkmer) {

        if ((nt == 'A') or (nt == 'a')) {

            nt_bitkmer->set(0, 0);
            nt_bitkmer->set(1, 0);

        } else if ((nt == 'C') or (nt == 'c')) {

            nt_bitkmer->set(0, 1);
            nt_bitkmer->set(1, 0);
        
        } else if ((nt == 'G') or (nt == 'g')) {

            nt_bitkmer->set(0, 0);
            nt_bitkmer->set(1, 1);
        
        } else if ((nt == 'T') or (nt == 't')) {

            nt_bitkmer->set(0, 1);
            nt_bitkmer->set(1, 1);
        
        } else {

            return false;
        }

        return true;
    }
    

    template <typename Type>
    inline vector<pair<typename Type::iterator, typename Type::iterator> > allocateToThreads(Type & container, const ushort num_threads) {

        vector<pair<typename Type::iterator, typename Type::iterator> > thread_object_allocations;

        int thread_object_base = container.size() / num_threads;
        int thread_object_modulo = container.size() % num_threads;

        if (container.size() <= num_threads) {

            thread_object_allocations.reserve(container.size());

            auto current_thread_start_object = container.begin();
            auto current_object = container.begin();
            uint object_idx = 0;
            
            while (current_object != container.end()) {

                current_object++;

                thread_object_allocations.emplace_back(current_thread_start_object, current_object);
                current_thread_start_object = current_object;

                object_idx++;
            }   

            assert(object_idx == container.size());
            assert(thread_object_allocations.size() == container.size());

        } else {

            thread_object_allocations.reserve(num_threads);

            auto current_object = container.begin();
            auto current_thread_start_object = container.begin();

            uint current_thread_end_object_idx = thread_object_base;
            
            if (thread_object_modulo > 0) {

                current_thread_end_object_idx++;
            }
            
            uint object_idx = 0;
            int current_thread = 0;

            while (current_object != container.end()) {

                if (object_idx == current_thread_end_object_idx) {

                    thread_object_allocations.emplace_back(current_thread_start_object, current_object);
                    current_thread_start_object = current_object;

                    current_thread_end_object_idx += thread_object_base;

                    if (current_thread < thread_object_modulo - 1) {

                        current_thread_end_object_idx++;
                    }

                    current_thread++;
                }

                current_object++;
                object_idx++;
            }
        
            thread_object_allocations.emplace_back(current_thread_start_object, container.end());
            assert(object_idx == container.size());
            assert(current_thread_end_object_idx == container.size());
            assert(thread_object_allocations.size() == num_threads);
        }
        
        return thread_object_allocations;
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
template < typename T, typename T2 >
std::ostream& operator << (std::ostream& os, typename std::pair<T, T2>& v) {
            
    os << "(" << to_string(v.first) << ", " << to_string(v.second) << ")";
    
    return os;
}


// Vector flush capability
template < typename T, typename T2 >
std::ostream& operator << (std::ostream& os, typename std::vector<std::pair<T, T2> >& v) {
    
    for (typename vector<std::pair<T, T2> >::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        
        os << "(" << to_string(ii->first) << ", " << to_string(ii->second) << ") ";
    }

    return os;
}

// Vector flush capability
template < typename T >
std::ostream& operator << (std::ostream& os, typename std::vector<T>& v) {
    
    for (typename vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        
		os << *ii << "  ";
	}

    return os;
}


// 2D-vector flush capability
template < typename T >
std::ostream& operator << (std::ostream& os, typename std::vector<std::vector<T> >& v) {
    
    for (typename vector<vector<T> >:: iterator ii = v.begin(); ii != v.end(); ++ii) {
		
		os << *ii << "\n";
	}
	
    return os;
}


// VariantType flush capability
inline std::ostream& operator << (std::ostream& os, const Utils::VariantType& vt) {
            
    os << Utils::variant_type_strings.at(static_cast<uchar>(vt));
    return os;
}


// ClusterClass flush capability
inline std::ostream& operator << (std::ostream& os, const Utils::FilterStatus& let) {
            
    os << Utils::filter_status_strings.at(static_cast<uchar>(let));
    return os;
}


// ChromosomeClass flush capability
inline std::ostream& operator << (std::ostream& os, const Utils::ChromosomeClass& cc) {
            
    os << Utils::chromosome_class_strings.at(static_cast<uchar>(cc));
    return os;
}

#endif