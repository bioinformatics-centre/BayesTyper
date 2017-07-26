
/*
VcfFile.hpp - This file is part of BayesTyper (v0.9)


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


#ifndef VCF_FILE
#define VCF_FILE

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>

#include "Variant.hpp"
#include "VcfMetaData.hpp"

typedef unsigned char uchar;
typedef unsigned int uint;

using namespace std;

class VcfFile {

    public:

      VcfFile(const bool);
      virtual ~VcfFile() {};

      VcfMetaData & metaData();

    protected:

        VcfMetaData meta_data;

        const bool is_sorted;

        string current_chromosome;
        long current_position;

        void checkChromosomeAndPositionOrder(const string &, const uint);
};

class VcfFileReaderBase : public VcfFile {

  public:

    VcfFileReaderBase(string, const bool);
    ~VcfFileReaderBase();

    virtual bool getNextVariant(Variant **) = 0;

  protected:

    ifstream vcf_file;

    uint num_cols;

    vector<string> cur_var_line;
    bool last_line_read;
};


class VcfFileReader : public VcfFileReaderBase {

  public:

    VcfFileReader(string, const bool);
    virtual bool getNextVariant(Variant **);

  private:

    void updateVariantLine();
};

class GenotypedVcfFileReader : public VcfFileReaderBase {

  public:

    GenotypedVcfFileReader(string, const bool);
    virtual bool getNextVariant(Variant **);

  private:

    void updateVariantLine();
};


class VcfFileWriter : public VcfFile {

  public:

    VcfFileWriter(string, const VcfMetaData &, const bool);
    ~VcfFileWriter();
    
    const VcfMetaData & metaData() const;

    void write(Variant *);

  protected:

      ofstream vcf_file;
};

#endif
