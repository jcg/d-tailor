/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef GENE_READER_H
#define GENE_READER_H
#include <iostream>
#include <fstream>
#include <string>
#include "seq.h"

using namespace std;

// abstract base class of all classes that can read annotation files
class GeneReader
{
public:
    virtual bool read_genes(Genome &) = 0;
    virtual ~GeneReader() {};
    virtual bool good() = 0;
};

// read a .coords file from TIGR CMR
class CoordsReader : public GeneReader
{
public:
    CoordsReader(const string &);
    virtual ~CoordsReader() {}
    bool read_genes(Genome &);
    bool good() { return _in.good(); }
private:
    ifstream _in;
};

// read a .ptt file from GenBank
class PTTReader : public GeneReader
{
public:
    PTTReader(const string &);
    virtual ~PTTReader() {}
    bool read_genes(Genome &);
    bool good() { return _in.good(); }
private:
    ifstream _in;
    string _id;
};

// return the correct reader class given a filename
GeneReader * gene_reader_factory(const string &);

#endif
