/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef DISTR_H
#define DISTR_H
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

// represents a "distribution" ---- really a binned vector so that float
// values go into the appropriate bin
struct Distribution
{
    // default constructor makes empty null distribution
    Distribution() : _min(0), _max(0), bs(0) {}

    // distributions have a fixed range [mmin, mmax] and a set # of bins
    // the range is broken into evenly sized bins
    Distribution(double mmin, double mmax, unsigned bins) :
        _min(mmin), _max(mmax), bs((_max-_min)/bins), d(bins) {}

    // the number of bins in the histogram
    unsigned size() const { return d.size(); }

    // the 'x' range of values represented by bin i
    double bin_lower(unsigned i) const { return _min + i*bs; }
    double bin_upper(unsigned i) const { return min(_max,_min + (i+1)*bs); }

    // the range of values represented by the histogram (as defined when
    // created --- does not mean that we've /seen/ a value that low)
    double low() const { return _min; }
    double high() const { return _max; }

    // the bin # (0-based) that the value v would go into
    unsigned binfor(double v) const; 

    // accessors: [i] gives the value of bin i (0-based) --- if the distr is
    // not const, can be used as an lvalue: dist[2] += 1.2; at(v) gives the
    // value of the bin that v would go into (also can be used as an lvalue)
    double operator[](unsigned i) const { return d[i]; }
    double & operator[](unsigned i) { return d[i]; }
    double & at(double v) { return d[binfor(v)]; }
    double at(double v) const { return d[binfor(v)]; }

    // interpolate a value for v
    double interp(double v) const;

protected:
    double _min, _max, bs;
    vector<double> d;
};

ostream & operator<<(ostream &, const Distribution &);

#endif
