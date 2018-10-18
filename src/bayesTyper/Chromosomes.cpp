
/*
Chromosomes.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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

#include <iostream>
#include <fstream>

#include "boost/algorithm/string.hpp"

#include "Chromosomes.hpp"
#include "Utils.hpp"


Chromosomes::Chromosomes(const string & fasta_filename, const bool is_decoy) {

    total_length = 0;
    decoy_length = 0;

    parseFasta(fasta_filename, is_decoy);
}

void Chromosomes::addFasta(const string & fasta_filename, const bool is_decoy) {

    parseFasta(fasta_filename, is_decoy);
}

void Chromosomes::addSequence(const pair<string, string> & sequence, const bool is_decoy) {

    assert(!sequence.first.empty());

    _chromosomes.emplace_back(sequence); 
    
    if (!order.emplace(sequence.first, _chromosomes.size() - 1).second) {

        cerr << "\nERROR: Chromosome \"" << sequence.first << "\" appears multiple times in fasta file(s)\n" << endl;
        exit(1);
    }

    if (is_decoy) {

        assert(decoys.emplace(sequence.first).second);
        decoy_length += _chromosomes.back().second.size();
    }

    total_length += _chromosomes.back().second.size();
}

void Chromosomes::parseFasta(const string & fasta_filename, const bool is_decoy) {

    if (!fasta_filename.empty()) {

        ifstream fasta_infile(fasta_filename);

        if (!fasta_infile.is_open()) {

            cerr << "\nERROR: Unable to open file " << fasta_filename << "\n" << endl;
            exit(1);
        }

        string chrom_name = "";

        for (string fasta_line; getline(fasta_infile, fasta_line);) {

            if (fasta_line.substr(0,1) == ">") {

                vector<string> chrom_name_split;

                boost::split(chrom_name_split, fasta_line, boost::is_any_of("\t "));
                assert(chrom_name_split.front().at(0) == '>');

                chrom_name = chrom_name_split.front().substr(1);
                assert(!chrom_name.empty());

                addSequence(make_pair(chrom_name, ""), is_decoy);

            } else {

                assert(!chrom_name.empty());
                assert(!_chromosomes.empty());

                _chromosomes.back().second.append(fasta_line);
                total_length += fasta_line.size();

                 if (is_decoy) {

                    decoy_length += fasta_line.size();
                }                   
            }
        }

        fasta_infile.close();
    }
}

SequenceIterator Chromosomes::cbegin() const {

    return _chromosomes.cbegin();
}

SequenceIterator Chromosomes::cend() const {

    return _chromosomes.cend();
}

SequenceIterator Chromosomes::find(const string & name) const {

    auto order_it = order.find(name);

    if (order_it != order.end()) {

        return (_chromosomes.cbegin() + order_it->second);

    } else {

        return cend();
    }
}

bool Chromosomes::isDecoy(const string & name) const {

    assert(order.find(name) != order.end());
    return (decoys.count(name) > 0);
}

uint Chromosomes::getTotalCount() const {

    return _chromosomes.size();
}

uint Chromosomes::getDecoyCount() const {

    return decoys.size();
}

ulong Chromosomes::getTotalLength() const {

    return total_length;
}

ulong Chromosomes::getDecoyLength() const {

    return decoy_length;
}

void Chromosomes::convertToUpper() {

    for (auto chromosome: _chromosomes) {

        std::transform(chromosome.second.begin(), chromosome.second.end(), chromosome.second.begin(), ::toupper);
    }
}





