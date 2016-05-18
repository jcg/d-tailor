/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <fstream>
#include <cmath>
#include <cassert>
#include "transterm.h"
#include "seq.h"
#include "distr.h"
#include "conf.h"


// return the min and max values for the given energy kind among the terms
// high and low are /cumulative/ over successive calls
void
energy_range(
    Term::EnergyKind k, 
    const ConstTermVec & terms, 
    double & low, 
    double & high)
{
    for(EVERY_CTERM_CONST(terms, T))
    {
        low = min(low, (*T)->energy(k));
        high = max(high, (*T)->energy(k));
    }
}


// compute the energy distiribution of the given kind of energy
void
term_energy_dist(
    Term::EnergyKind k, 
    const ConstTermVec & terms, 
    Distribution & dist)
{
    for(EVERY_CTERM_CONST(terms, T))
    {
        dist.at((*T)->energy(k)) += 1.0;
    }

    for(unsigned i = 0; i < dist.size(); i++)
    {
        dist[i] /= terms.size();
    }
}


// return the signal to noise distribution 
Distribution
signal_to_noise(
    Term::EnergyKind k, 
    const ConstTermVec & signal, 
    const ConstTermVec & noise)
{
    const int NUM_BINS = 9;

    // get the range for the energy values
    Energy low = 1000, high = -1000;
    energy_range(k, signal, low, high);
    energy_range(k, noise, low, high);

    // compute the individual distributions
    Distribution sigdist(low, high, NUM_BINS), noidist(low, high, NUM_BINS);
    term_energy_dist(k, signal, sigdist);
    term_energy_dist(k, noise, noidist);

    // divide the noise by the signal
    Distribution dist(low, high, NUM_BINS); 
    for(int i = 0; i < NUM_BINS; i++)
    {
        dist[i] = (sigdist[i]>0)?(noidist[i]/sigdist[i]):0;
    }

    return dist;
}


// collect the terminators taht are in (and only in) the given regiontype. The
// list can be retrieved with terms()
class PureRegionTerms : public EventResponder
{
public:
    PureRegionTerms(RegionType r, bool c=true) : _region(r), _codir(c) {}
    virtual ~PureRegionTerms() {}

    void terminator(const Term * term)
    {
        bool save = false;
        switch(_region)
        {
            // terms in gene regions must not be in any kind of integenic
            // region and, if user wants codir only, the sense of the
            // terminator must match some enclosing gene
            case GENE:
                save = !in_t2t() && !in_h2t_fwd() && !in_h2t_rvs();
                if(save && _codir) 
                {
                    save = (fwd_gene_count() > 0 && term->dir() == FORWARD) ||
                           (rvs_gene_count() > 0 && term->dir() == REVERSE);
                }
                break;

            case HEAD2TAIL:
                save = gene_count() == 0 && (in_h2t_fwd() || in_h2t_rvs());
                if(save && _codir)
                {
                    save = (in_h2t_fwd() && term->dir() == FORWARD) ||
                           (in_h2t_rvs() && term->dir() == REVERSE);
                }
                break;

            case TAIL2TAIL:
                save = gene_count() == 0 && in_t2t();
                break;

            case HEAD2HEAD:
                save = gene_count() == 0 && in_h2h();
                break;
        }

        if(save) _tvec.push_back(term);
    }

    const ConstTermVec & terms() const { return _tvec; }

private:
    RegionType _region;
    ConstTermVec _tvec;
    bool _codir;
};


// count the # of AT characters between [start,end] (inclusive).  Also return
// the total # of characters
void
count_at(SeqPtr start, SeqPtr end, unsigned long & at, unsigned long & len)
{
    for(; start <= end; ++start)
    {
        if(*start == 'A' || *start == 'T') at++;
        len++;
    }
}


// count the AT content of the gene and intergenic regions
class ATContent : public EventResponder
{
public:
    ATContent() 
    {
        _gene_at = _gene_len = 0;
        _nongene_at = _nongene_len = 0;    
    }
    virtual ~ATContent() {}

    void start(const Seq & seq, Direction dir)
    {
        EventResponder::start(seq, dir);
        _nongene_start = seq.left();
        _gene_start = 0;
    }

    void enter_gene(const Event & e) 
    {
        EventResponder::enter_gene(e);
        // if we just left a nongene region
        if(gene_count() == 1)
        {
            if(_nongene_start && e.place > _nongene_start) 
            {
                count_at(_nongene_start, e.place-1, _nongene_at, _nongene_len);
//                cerr << "NONGENE: " << seqindex(*e.reg->seq, _nongene_start) << " " << seqindex(*e.reg->seq, e.place-1) << endl;
                _nongene_start = 0;
            }
            _gene_start = e.place;
        }
    }

    void leave_gene(const Event & e)
    {
        EventResponder::leave_gene(e);
//        cerr << "LEFT: " << e.reg->name << " at " << seqindex(*e.reg->seq, e.place) << endl;
        // if we just entered a non-gene region
        if(gene_count() == 0)
        {
            if(_gene_start) 
            {
                count_at(_gene_start, e.place, _gene_at, _gene_len);
//                cerr << "NONGENE: " << e.reg->name << " " << seqindex(*e.reg->seq, _gene_start) << " " << seqindex(*e.reg->seq, e.place-1) << endl;
                _gene_start = 0;
            }
            _nongene_start = e.place + 1;
        }
    }

    unsigned long gene_at() const { return _gene_at; }
    unsigned long gene_len() const { return _gene_len; }
    unsigned long nongene_at() const { return _nongene_at; }
    unsigned long nongene_len() const { return _nongene_len; }

private:
    SeqPtr _gene_start;
    SeqPtr _nongene_start;

    unsigned long _gene_at, _gene_len;
    unsigned long _nongene_at, _nongene_len;
};


// compute the total length of the intergenic regions of the given type
class RegionLength : public EventResponder
{
public:
    RegionLength(RegionType r) : _region(r), _len(0) {}
    virtual ~RegionLength() {}

    void start(const Seq & seq, Direction dir)
    {
        EventResponder::start(seq, dir);
        _ig_count = 0;
        _s = 0;
    }

    void enter_intergene(RegionType r, Direction d, const Event & e)
    {
        EventResponder::enter_intergene(r, d, e);

        if(r == _region)
        {
            _ig_count++;
            if(_ig_count == 1) _s = e.place + 1;
        }
    }

    void leave_intergene(RegionType r, Direction d, const Event & e)
    {
        EventResponder::leave_intergene(r, d, e);

        if(r == _region)
        {
            assert(_ig_count > 0);
            _ig_count--;
            if(_ig_count == 0)
            {
                assert(_s != 0);
                _len += e.place - _s; // correct without the +1
                _s = 0;
            }
        }
    }

    unsigned long length() const { return _len; }

private:
    RegionType _region;
    unsigned long _len;
    SeqPtr _s;
    int _ig_count;
};


// process the genome creating the statistic necessary for score() to assign a
// confidence value. 
void
ErmolaevaConfidence::prepare(const Genome & seqs)
{
    PureRegionTerms gene(GENE), h2t(HEAD2TAIL), t2t(TAIL2TAIL, false);
    RegionLength h2t_len(HEAD2TAIL), t2t_len(TAIL2TAIL);
    ATContent at;

    for(EVERY_CHROM_CONST(seqs, C))
    {
        // gene_er, h2t_er, and t2t_er are event responders that manintain
        // a running list of matching terms.
        scan_events(**C, gene, gene_start_cut, gene_end_cut);  //100
        scan_events(**C, h2t, gene_start_cut, gene_end_cut);   //50
        scan_events(**C, t2t, gene_start_cut, gene_end_cut);   //50

        // comptue the length of the T2T and H2T regions
        scan_events(**C, h2t_len, gene_start_cut, gene_end_cut); //50
        scan_events(**C, t2t_len, gene_start_cut, gene_end_cut); //50

        // compute the at content of the gene and integenic regions
        scan_events(**C, at, gene_start_cut, gene_end_cut); //50
    }

    // can't compute confidence if we have no gene_terms
    if(gene.terms().empty())
    {
        prepared = false;
        cout << "warning: no examples in genes; can't compute conf." << endl;
        return;
    }

    // compute K --- correction for AT content
    double at_in, at_not;
    at_in = double(at.gene_at()) / at.gene_len();
    at_not = double(at.nongene_at()) / at.nongene_len();

    K = (840*at_not*at_not - 1215.65*at_not + 448.9593) /
        (840*at_in*at_in   - 1215.65*at_in  + 448.9593);

    // output status messages
    cout << "Genes: " 
         << at_in << " %AT, " 
         << at.gene_len() << " nt, "
         << gene.terms().size() << " terms." << endl;

    cout << "Intergenic: "
         << at_not << " %AT, "
         << "H2T: " << h2t_len.length() << " nt, " 
         << h2t.terms().size() << " terms; "
         << "T2T: " << t2t_len.length() << " nt, " 
         << t2t.terms().size() << " terms. " 
         << endl;

    t2t_L = double(t2t_len.length()) / at.gene_len();
    h2t_L = double(h2t_len.length()) / at.gene_len();

    t2t_hp = signal_to_noise(Term::HAIRPIN, t2t.terms(), gene.terms());
    t2t_tail = signal_to_noise(Term::TAIL, t2t.terms(), gene.terms());

    h2t_hp = signal_to_noise(Term::HAIRPIN, h2t.terms(), gene.terms());
    h2t_tail = signal_to_noise(Term::TAIL, h2t.terms(), gene.terms());

    t2t_N = 2.0 * gene.terms().size() / t2t.terms().size();
    h2t_N = double(gene.terms().size()) / h2t.terms().size();

    prepared = true;
}


// compute the confidence score of the given terminator, which is in a region
// whose type is given by where (gene, head-to-tail or tail-to-tail). If the
// terminator has a 'partner' we modifiy the confidence score
int
ErmolaevaConfidence::score(const Term & term, RegionType where) const
{
    int c1 = score_one(term, where);

    // if this term has a partner (going in the other dir), and we're in a
    // TAIL2TAIL region (meaning that a bidirectional terminator makes sense)
    // we modify the confidence to take this into account
    if(term.partner && where == TAIL2TAIL) 
    {
        assert(term.partner->dir() != term.dir());
        int c2 = score_one(*term.partner, where);
        if(c1 >= 50 && c2 >= 50)
        {
            c1 = int((1.0 - (1.0 - c1/100.0)*(1.0 - c2/100.0))*100.0 + 0.5);
        }
    }
    return c1;
}


// compute the confidence score of the given terminator, which is in a region
// whose type is given by where (gene, head-to-tail or tail-to-tail)
int
ErmolaevaConfidence::score_one(const Term & term, RegionType where) const
{
    if(!prepared) return 0;

    const Distribution * taild, *hpd;
    double N, L;

    if(where == TAIL2TAIL)
    {
        taild = &t2t_tail;
        hpd = &t2t_hp;
        N = t2t_N;
        L = t2t_L;
    }
    else if(where == HEAD2TAIL)
    {
        taild = &h2t_tail;
        hpd = &h2t_hp;
        N = h2t_N;
        L = h2t_L;
    }
    else return 0;

    double qE, qT;
    qE = max(hpd->interp(term.hp_energy), 0.0);
    qT = max(taild->interp(term.tail_energy), 0.0);

    return unsigned(max(1-N*L*qE*qT*K, 0.0) * 100.0 + 0.5);
}


RandomPValueConfidence::RandomPValueConfidence(const string & fn)
    : RandomConfidence(fn)
{
    sum_exp_table();
}

int
RandomPValueConfidence::score(const Term & term, RegionType reg) const
{
    double at = _emp_at.find(reg)->second;
    return histvalue(get_table(at), term, at);
}

void
RandomPValueConfidence::sum_exp_table()
{
    for(map<int, Histogram2d>::iterator H = _exp_table.begin();
        H != _exp_table.end();
        ++H)
    {
        Histogram2d & hist = H->second;
        for(int i = 0; i < _nbins; i++)
        {
            for(int j = 0; j < _nbins; j++)
            {
                // In general, new_h[i][j] = new_h[i-1, j] + new_h[i,j-1] - new_h[i-1,j-1]
                // edge cases omit some of the terms
                if (i > 0) hist[i][j] += hist[i-1][j];
                if (j > 0) hist[i][j] += hist[i][j-1];
                if (i > 0 && j > 0) hist[i][j] -= hist[i-1][j-1];
            }
        }

//        double total_ex = double(hist[_nbins-1][_nbins-1]);
//CURRENT:        double max_log = log(1.0 / total_ex);

        double max_log = log(1.0 / _sample_size);
        double ss_log = log(double(_sample_size));

        // x is Rgc in the paper
        for(int i = 0; i < _nbins; i++)
        {
            for(int j = 0; j < _nbins; j++)
            {
//                cout << hist[i][j] << " " << total_ex << " " << max_log << " ";
                unsigned long x = max(hist[i][j], (unsigned long)1);
//CURRENT:                hist[i][j] = int(100.0 * log(x / total_ex) / max_log);
                hist[i][j] = int( 100.0 * (log((long double)x) - ss_log) / max_log );

//                cout << hist[i][j] << endl;
                assert(hist[i][j] >= 0 && hist[i][j] <= 100);

                //hist[i][j] = int(100.0 * (1.0 - (hist[i][j] / total_ex)));
            }
        }
    } 
}

//============================================================
// Random Confidence Scheme
//============================================================
void
debug_print_emp_table(RandomConfidence::Histogram2d & T)
{
    for(RandomConfidence::Histogram2d::iterator I = T.begin();
        I != T.end();
        ++I)
    {
        for(vector<unsigned long>::iterator J = I->begin();
            J != I->end();
            ++J)
        {
            cout << *J << " ";
        }
        cout << endl;
    }
}


// read the given exphist file
RandomConfidence::RandomConfidence(const string & fn)
    : _prepared(false)
{
    read_exp_table(fn);
}


// create the exp tables for each of the region types
void
RandomConfidence::prepare(const Genome & seqs)
{
    PureRegionTerms gene(GENE), h2t(HEAD2TAIL), t2t(TAIL2TAIL, false), 
                    h2h(HEAD2HEAD, false);
    RegionLength h2t_len(HEAD2TAIL), t2t_len(TAIL2TAIL), h2h_len(HEAD2HEAD);
    ATContent at;

    for(EVERY_CHROM_CONST(seqs, C))
    {
        // gene_er, h2t_er, and t2t_er are event responders that manintain
        // a running list of matching terms.
        scan_events(**C, gene, gene_start_cut, gene_end_cut);  // 100
        scan_events(**C, h2t, gene_start_cut, gene_end_cut);    // 50
        scan_events(**C, t2t, gene_start_cut, gene_end_cut);    // 50
        scan_events(**C, h2h, gene_start_cut, gene_end_cut);    // 50

        // comptue the length of the T2T and H2T regions
        scan_events(**C, h2t_len, gene_start_cut, gene_end_cut);  // 50
        scan_events(**C, t2t_len, gene_start_cut, gene_end_cut);  // 50
        scan_events(**C, h2h_len, gene_start_cut, gene_end_cut);  // 50

        // compute the at content of the gene and integenic regions
        scan_events(**C, at, gene_start_cut, gene_end_cut);       // 50
    }

    _emp_at[TAIL2TAIL] = double(at.nongene_at()) / at.nongene_len();
    _emp_at[HEAD2HEAD] = _emp_at[HEAD2TAIL] = _emp_at[TAIL2TAIL];
    _emp_at[GENE] = double(at.gene_at()) / at.gene_len();

    //_emp_at[GENE] = _emp_at[TAIL2TAIL]; // XXX: test using only intergenic AT content

    _emp_len[TAIL2TAIL] = 2*t2t_len.length(); 
    _emp_len[HEAD2TAIL] = h2t_len.length();
    _emp_len[HEAD2HEAD] = 2*h2h_len.length();
    _emp_len[GENE] = at.gene_len();

    // output status messages
    cout << "Genes: " 
         << _emp_at[GENE] << " %AT, " 
         << at.gene_len() << " nt, "
         << gene.terms().size() << " terms." << endl;

    cout << "Intergenic: "
         << _emp_at[TAIL2TAIL] << " %AT, "
         << "H2T: " << h2t_len.length() << " nt, " 
         << h2t.terms().size() << " terms; "
         << "T2T: " << t2t_len.length() << " nt, " 
         << t2t.terms().size() << " terms; " 
         << "H2H: " << h2h_len.length() << " nt, " 
         << h2h.terms().size() << " terms. " 
         << endl;
    
    fill_emp_table(GENE, gene.terms());
    fill_emp_table(HEAD2TAIL, h2t.terms());
    fill_emp_table(TAIL2TAIL, t2t.terms());
    fill_emp_table(HEAD2HEAD, h2h.terms());

    //debug_print_emp_table(_emp_table[TAIL2TAIL]);

    _prepared = true;
}


// return the score for the given terminator assuming its from the region reg
int
RandomConfidence::score(const Term & term, RegionType reg) const
{
    assert(_prepared);
    // NOTE: we must use x.find(reg)->second rather than x[reg] since
    // this function is const and operator[] is not

    unsigned long len = _emp_len.find(reg)->second;
    double at = _emp_at.find(reg)->second;

    double expv = double(histvalue(get_table(at), term, at));
    unsigned long empv = histvalue(_emp_table.find(reg)->second, term, at);


    if(empv == 0) return 0;
    //cout << expv << " " << _sample_size;

    expv *=  double(len) / _sample_size;
#if 0
    //cout << "C: " << reg << " " << len << " " << at << " " << expv << " " << empv << " " 
    cout << "C: " << reg << " " << len << " " << term.hp_energy << " " << term.tail_energy 
         << " " << expv << " " << empv << " " << " at=" << get_best_at(at) << " "
         <<  unsigned(100.0 * max(1.0 - expv / empv, 0.0)) << endl;
#endif

    return unsigned(100.0 * max(1.0 - expv / empv, 0.0));
}


// return teh bin index for a histogram between [low,hi] of n bins, for value
// v bins numbered 0 to n -1. if v is outside [low,hi], return the extream
// bins

int
RandomConfidence::hbin(double low, double hi, int n, double v) const
{
    if(v < low) return 0;
    if(v > hi) return n-1;    

    return unsigned((v-low) / ((hi-low)/n));
}


// return the value of the 2d histogram appropriate for term t
unsigned long 
RandomConfidence::histvalue(const Histogram2d & hist, const Term & t, double at) const
{
    int ati = get_best_at(at);
#if 0
    cout << ati << " " << _low_hp.find(ati)->second << " " << _high_hp.find(ati)->second << " "
         << _low_tail.find(ati)->second << " " << _high_tail.find(ati)->second << " "
         << t.hp_energy << " " << t.tail_energy << " ";
    int i = hbin(_low_hp.find(ati)->second, _high_hp.find(ati)->second,  _nbins, t.hp_energy);
    int j = hbin(_low_tail.find(ati)->second, _high_tail.find(ati)->second, _nbins, t.tail_energy);
    cout << i << " " << j << " | ";
#endif

    return hist[hbin(_low_hp.find(ati)->second, _high_hp.find(ati)->second,  _nbins, t.hp_energy)]
               [hbin(_low_tail.find(ati)->second, _high_tail.find(ati)->second, _nbins, t.tail_energy)];
}


// return the value of the 2d histogram appropriate for term t
unsigned long &
RandomConfidence::histvalue(Histogram2d & hist, const Term & t, double at)
{
    int ati = get_best_at(at);
    return hist[hbin(_low_hp[ati], _high_hp[ati], _nbins, t.hp_energy)]
               [hbin(_low_tail[ati], _high_tail[ati], _nbins, t.tail_energy)];
}


// given a list of terms that come from reg type reg, construct an emp
// 2d histogram
void
RandomConfidence::fill_emp_table(RegionType reg, const ConstTermVec & terms)
{
    _emp_table[reg].resize(_nbins);
    for(int i = 0; i < _nbins; i++) _emp_table[reg][i].resize(_nbins);

    for(ConstTermVec::const_iterator T = terms.begin();
        T != terms.end();
        ++T)
    {
        histvalue(_emp_table[reg], **T, _emp_at[reg])++;
    }
}

int
RandomConfidence::get_best_at(double at) const
{
    assert(0 < at && at < 1);

    double best_diff = 1000;
    int best_ati = -1;
    for(map<int, Histogram2d>::const_iterator H = _exp_table.begin();
        H != _exp_table.end();
        ++H)
    {
        if(fabs(int(H->first) - 100*at) < best_diff)
        {
            best_diff = fabs(int(H->first) - 100*at);
            best_ati = H->first;
        }
    }
    assert(best_ati > 0);
    return best_ati;
}


// find the histogram table that has an at content closest to at
RandomConfidence::Histogram2d &
RandomConfidence::get_table(double at)
{
    assert(0 < at && at < 1);

    //int ati = int(100*at);
    double best_diff = 1000;
    Histogram2d * hist = 0;
    for(map<int, Histogram2d>::iterator H = _exp_table.begin();
        H != _exp_table.end();
        ++H)
    {
        if(!hist || fabs(int(H->first) - 100*at) < best_diff)
        {
            hist = &H->second;
            best_diff = fabs(int(H->first) - 100*at);
        }
    }
    return *hist;
}


// find the histogram table that has an at content closest to at
const RandomConfidence::Histogram2d &
RandomConfidence::get_table(double at) const
{
    assert(0 < at && at < 1);

    //int ati = int(100*at);
    double best_diff = 1000;
    const Histogram2d * hist = 0;
    for(map<int, Histogram2d>::const_iterator H = _exp_table.begin();
        H != _exp_table.end();
        ++H)
    {
        if(!hist || fabs(int(H->first) - 100*at) < best_diff)
        {
            hist = &H->second;
            best_diff = fabs(int(H->first) - 100*at);
        }
    }
    return *hist;
}


// read the expterms.dat file given by fn the histogram matrix is written so
// that the tail scores go along the ROWS
// In the histograms, the organization is: H[hp_score][tail_score]
// in the file, rows are tail scores, columns are hp_scores
void
RandomConfidence::read_exp_table(const string & fn)
{
    ifstream in(fn.c_str());
    if(!in.good())
    {
        cerr << "Couldn't open data file: " << fn << endl;
        exit(3);
    }

    in >> _sample_size >> _nbins ;
       //>> _low_hp >> _high_hp >> _low_tail >> _high_tail;

    cout << _nbins << " bins; sample size is " << _sample_size << endl;
    //cout << "hp range = " << _low_hp << " to " << _high_hp 
    //     << ". tail range = " << _low_tail << " to " << _high_tail << endl;

    double at;
    while(in >> at)
    {
        if(at < 0.0) break;  // a negative AT value means stop

        int ati = int(at*100 + 0.001);

        in >> _low_hp[ati] >> _high_hp[ati] >> _low_tail[ati] >> _high_tail[ati];

        _exp_table[ati].resize(_nbins);
        for(int i = 0; i < _nbins; i++) _exp_table[ati][i].resize(_nbins);

        for(int i = 0; i < _nbins; i++)
        {
            // j is the column indx == hp scores
            for(int j = 0; j < _nbins; j++) 
            {
                in >> _exp_table[ati][j][i];
            }
        }
    }
}


// from inside an event responder, compute the best confidence score for the
// regions the event responder is currently in
int
er_confidence(
    const EventResponder & er, 
    const Confidence & conf, 
    const Term & term)
{
    int c_h2t, c_t2t, c_gene, c_h2h;
    c_h2t = c_t2t = c_gene = c_h2h = 0;

    if(er.in_t2t())
    {
        c_t2t = conf.score(term, TAIL2TAIL);
    }

#if 0
    if(er.in_h2t_fwd() && term.dir() == FORWARD ||
       er.in_h2t_rvs() && term.dir() == REVERSE)
#endif
    if(er.in_h2t_fwd() || er.in_h2t_rvs())
    {
        c_h2t = conf.score(term, HEAD2TAIL);
    }

    if(er.in_h2h())
    {
        c_h2h = conf.score(term, HEAD2HEAD);
    }

    if(er.gene_count() > 0) c_gene = conf.score(term, GENE);

    return max(max(max(c_h2t, c_t2t), c_gene), c_h2h);
}
