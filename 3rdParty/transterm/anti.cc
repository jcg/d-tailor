/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

// The code to search for anti-terminators is currently in developement
// and does not work as well as one would hope

#include <iostream>
#include <cassert>
#include "conf.h"
#include "seq.h"
#include "transterm.h"

bool best_overlapping_hairpin(const Term &, Term &, int &);

class OutputAntiTerms : public EventResponder
{
public:
    OutputAntiTerms(ostream & out, Confidence & conf, int cutoff) :
        _out(out),
        _conf(conf),
        _cutoff(cutoff)
    {}

    virtual ~OutputAntiTerms() {}

    void 
    terminator(const Term * term)
    {
        int c = er_confidence(*this, _conf, *term);
        if(c >= _cutoff)
        {
            Term best_term;
            int overl;
            if(best_overlapping_hairpin(*term, best_term, overl) && 
               best_term.hp_energy < -2.0)
            {
                _out << c << " " << *term << " " 
                     <<  best_term.hp_energy << " " << overl;
                print_term_seq(_out, best_term);
                _out << endl;
            }
        }
    }
private:
    ostream & _out;
    const Confidence & _conf;
    int _cutoff;
};


/*
3215191 - 3215211 + R conf -5.5  -6.28967  -9.9 overlap
                        CCGTTTGTTCCCCGGctgttttttATCAAaaatcagttTTTTTTATTTCTA 
 TGACTTGATGACACaaaaaacagccgTTTGTTCCCcggctgttttttATC
*/

void
output_anti_terms(
    ostream & out, 
    const Genome & g, 
    Confidence & conf, 
    int cutoff)
{
    OutputAntiTerms oat(out, conf, cutoff);

    out << "Anti-terminators: " << endl;

    for(EVERY_CHROM_CONST(g, C))
    {
        out << "In " << (*C)->name << " :" << endl;
        scan_events(**C, oat, 100, 100); // XXX: think about what good cutoffs are
    }
}
