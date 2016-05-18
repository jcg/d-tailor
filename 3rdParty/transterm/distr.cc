/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <cassert>
#include "distr.h"

// Distribution objects are used by the v1.0 confidence scheme

// output a distribution (surrounded by [])
ostream &
operator<<(ostream & out, const Distribution & dist)
{
    out << "[ ";
    for(unsigned i = 0; i< dist.size();i++)
    {
        out << dist[i] << " ";
    }
    out << "] ";
    return out;
}

// the bin # (0-based) that the value v would go into
unsigned 
Distribution::binfor(double v) const 
{
    if(v <= low()) return 0;
    if(v >= high()) return size()-1; 

    return min(unsigned((v-_min) / bs), size()-1);
}

// interpolate a value for v. Schematically:
//           *             ,
//           |      *      |
//    |      |      |      |
// (-----)(-----)(-----)(-----)
//    a      b      c      d
// if v is in [b,c] (the midpoints of two bins) then the line between the two
// * (the values for those two bins), is computed using the fact that the
// value of that line at x= (b+c)/2 is ([b]+[c])/2 (we use the bin value
// directly for values at the very begining and end)
double
Distribution::interp(double v) const
{
    unsigned b = binfor(v);

    unsigned otherb = 0;
    bool change = false;
    double x = 0.0;
    if(b > 0 && v <= (bin_lower(b) + bin_upper(b))/2.0)
    {
        otherb = b-1;
        x = bin_lower(b);
        change = true;
    }
    else if(b < size()-1)
    {
        otherb = b+1;
        x = bin_upper(b);
        change = true;
    }

    double E;
    if(change)
    {
        E = (d[otherb] + d[b])/2.0 + 
            (v - x) * (d[b] - d[otherb])/(bin_upper(b) - bin_upper(otherb));
    }
    else
    {
        E = d[b];
    }

    return E;
}
