/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include "gene-reader.h"
#include "util.h"

// return the right kind of annotation reader given teh extn of the filename
// if none, return 0
GeneReader *
gene_reader_factory(const string & fn)
{
    string extn = fn.substr(fn.rfind('.')+1);
    if(extn == "ptt") return new PTTReader(fn);
    if(extn == "coords" || extn == "crd") return new CoordsReader(fn);
    return 0;
}


// make a new PTTReader.
PTTReader::PTTReader(const string & fn)
    : _in(fn.c_str())
{
    unsigned dotpos = fn.rfind('.');
    // remove the extn
    _id = (dotpos == 0)?fn:fn.substr(0, dotpos);

    // remove leading paths .../
    _id = _id.substr(_id.rfind('/')+1);
}


// read the genes and put them into the genome
bool
PTTReader::read_genes(Genome & g)
{
    string line;
    assert(good());

    // read the header
    getline(_in, line);
    getline(_in, line);
    getline(_in, line);

    while(getline(_in, line))
    {
        string loc, strand, pid, gene, syn, code, cog, desc, name;
        int len;

        istringstream iss(line);
        iss >> loc >> strand >> len >> pid >> gene >> syn >> code >> cog;
        getline(iss, desc);
        desc = trim_front(desc);

        if(gene != "-") name = gene;
        else if(syn != "-") name = syn;
        else if(pid != "-") name = pid;
        else name = "UNK";

        vector<string> locvec;
        split(loc, '.', locvec);
        if(locvec.size() != 3)
        {
            cerr << "Bad location format: " << loc << endl;
            exit(3);
        }

        unsigned long start, end;
        start = atol(locvec[0].c_str());
        end = atol(locvec[2].c_str());

        // make sure that the coords are given as left..right where left <= right
        if(end < start)
        {
            cerr << "WARNING: TransTerm does not handle genomes with genes that wrap around." 
                 << endl;
            continue;
        }

        if(strand == "-") 
        {
            swap(start, end);
        }
        else if(strand != "+")
        {
            cerr << "Unknown strand value: " << strand << endl;
            exit(3);
        }

        Seq * s = chrom_for_id(g, _id);
        
        if(s)
        {
            if(start > s->length || end > s->length || start < 0 || end < 0)
            {
                cerr << "Bad gene coordinates: " << start << " - " << end << endl;
                exit(3);
            }

            s->genes.push_back(
                new Region(name, s, s->dna + start - 1, s->dna + end - 1, desc));
        }
        else
        {
            cerr << "Can't find seq for id: " << _id << endl;
            exit(3);
        }
    }
    return true;
}


// make a new reader for .coords file
CoordsReader::CoordsReader(const string & fn)
    : _in(fn.c_str())
{
}


// Read the gene cords file a stream formated as:
//      gene_name   start   end chrom_id
// where if start > end the gene runs on the other strand.
// Start and end are 1-based
bool
CoordsReader::read_genes(Genome & g)
{
    assert(good());
    Seq * s;
    string line;

    while(getline(_in, line))
    {
        string name, chid;
        unsigned long startidx, endidx;

        istringstream iss(line);
        chid = "";
        iss >> name >> startidx >> endidx >> chid;

        // if there was no chid, we assume that its b/c we are missing
        // the gene names.
        if(chid == "")
        {
            istringstream iss(line);
            iss >> startidx >> endidx >> chid;
            name = "UNK";
        }

        s = chrom_for_id(g, chid);
        if(s)
        {
            if(startidx > s->length || endidx > s->length ||
               startidx <= 0 || endidx <= 0)
            {
                cerr << "Bad gene coordinates: " << startidx << " .. " 
                     << endidx << endl;
                exit(3);
            }

            s->genes.push_back(
                new Region(name, s, s->dna + startidx - 1, s->dna + endidx - 1));
        }
        else
        {
            cerr << "Unknown chromosome id: " << chid << endl;
            exit(3);
        }
    }
    return true;
}
