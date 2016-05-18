/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "map-output.h"
#include "conf.h"
#include "seq.h"
#include "util.h"
#include "ermolaeva-score.h"
#include "transterm.h"

class T2THits : public EventResponder
{
public:
    T2THits(int ms, Confidence & conf) : 
        _min_gene_span(ms),
        _conf(conf),
        _conf_histo(101), 
        _term_histo(101),
        _comp_histo(101) 
    { 
        init0(); 
        _t2t_count = 0;
        _good_term_count = 0;
    }
    virtual ~T2THits() {}

    void start(const Seq & seq, Direction dir)
    {
        init0();
        _dir = dir;
    }

    void init0()
    {
        _best_conf = 0;
        _good_region = false;
        _sense_gene_span = 0;
    }

    void enter_intergene(RegionType r, Direction d, const Event & e)
    {
        EventResponder::enter_intergene(r, d, e);
        bool gr = r == TAIL2TAIL && _sense_gene_span >= _min_gene_span;
        if(gr && !_good_region) _t2t_count++;
        _good_region = gr;
        _best_conf = 0;
        _region_hit = false;
    }

    void leave_intergene(RegionType r, Direction d, const Event & e)
    {
        EventResponder::leave_intergene(r, d, e);
        if(r == TAIL2TAIL && _good_region)
        {
            if(_region_hit) 
            {
                assert(_best_conf <= 100);
                _conf_histo[_best_conf]++;
            }

            _good_region = false;
        }
    }

    void leave_gene(const Event & e)
    {
        EventResponder::leave_gene(e);

        if((_dir == FORWARD && e.kind == Event::ForwardGeneEnd) ||
           (_dir == REVERSE && e.kind == Event::ReverseGeneEnd))
        {
            _sense_gene_span++;
        }
        else
        {
            _sense_gene_span = 0;
        }
    }

    void terminator(const Term * term)
    {
        EventResponder::terminator(term);
        if(_good_region && term->dir() != _dir)
        {
            int c = er_confidence(*this, _conf, *term);
            _best_conf = max(c, _best_conf);
            _good_term_count++;
            _region_hit = true;
            _term_histo[c]++;
        }
    }

    const vector<int> & term_histo() const { return _term_histo; }
    const vector<int> & histo() const { return _conf_histo; }
    const vector<int> & comp_histo() const { return _comp_histo; }
    virtual int t2tregion_count() const { return _t2t_count; }
    int good_terms() const { return _good_term_count; }

protected:
    Direction _dir;
    int _min_gene_span;
    Confidence & _conf;

    int _sense_gene_span, _best_conf;
    bool _good_region;
    bool _region_hit;

    int _t2t_count;
    int _good_term_count;
    vector<int> _conf_histo;
    vector<int> _term_histo;
    vector<int> _comp_histo;
};


// return true if (l1,r1) and (l2,r2) intersect
bool
interval_intersect(SeqPtr l1, SeqPtr r1, SeqPtr l2, SeqPtr r2)
{
    return (l1> l2 && l1< r2) || (r1> l2 && r1< r2) || (l1< l2 && r1> r2);
}

bool
by_left_side(const Region * a, const Region * b)
{
    return a->left() < b->left();
}


// from a small list of terminators, remove those that dominate something
// in the list
void
remove_dominating(vector<const Term *> & terms)
{
    vector<const Term*> out;
    bool keep_T;
    
    for(vector<const Term*>::iterator T = terms.begin();
        T != terms.end();
        ++T)
    {
        // check every term (!= T) to see if T dominates it
        keep_T = true;
        for(vector<const Term*>::iterator R = terms.begin();
            R != terms.end();
            ++R)
        {
            if (*R == *T) continue;

            // if R is dominated by T, we want to keep R (at this stage) and
            // not T --- the rationale is that since these all have the same
            // confidence, R probably includes less of the tail pairing
            if (dominates(**T, **R)) keep_T = false;
        }
        if(keep_T) out.push_back(*T);
    }

    // copy the new list into the old list
    terms = out;
}


/*
  When scanning in reverse, we're looking for terminators for FORWARD directed
genes; when scanning forward, we're looking for terminators for REVERSE
directed genes. We scan in the opposite direction so that we can count the # of
genes facing the gene of interest.
*/

class DeHoonT2THits : public T2THits
{
public:
    DeHoonT2THits(int ms, Confidence & conf, 
        int rs, ostream & out) : 
        T2THits(ms, conf), 
        _region_size(rs),
        _out(out),
        _t2t_region_count(0)
    {
        // init list
    }

    virtual ~DeHoonT2THits() {}

    int t2tregion_count() const { return _t2t_region_count; }

    void start(const Seq & seq, Direction dir)
    {
        T2THits::start(seq, dir);

        _terms.clear();
    }

    void terminator(const Term * term)
    {
        // intentailly bypass parent class
        EventResponder::terminator(term);

        if(term->dir() != _dir) _terms.push_back(term);
    }

    void event(const Event & e)
    {
        /* when scanning forward, e.g., looking for terminators for rvs genes,
           we don't look past starts of genes on the same strand: 
                ---->  *|*---> */
        EventResponder::event(e);
        if((_dir == FORWARD && e.kind == Event::ReverseGeneStart) ||
           (_dir == REVERSE && e.kind == Event::ForwardGeneStart))
        {
            _terms.clear();
        }
    }

    void leave_intergene(RegionType r, Direction d, const Event & e)
    {
        if(r == TAIL2TAIL && _good_region)
        {
            _t2t_region_count++;
            // copy the terms that are within the range
            vector<const Term*> regterms;
            for(int i = _terms.size()-1; i >= 0; i--)
            {
                if((_dir == FORWARD && _terms[i]->right() < e.place -  _region_size) ||
                   (_dir == REVERSE && _terms[i]->left() > e.place + _region_size)) break;
                regterms.push_back(_terms[i]);
            }

            _good_term_count += regterms.size();

            int best_conf = 0, comp_best = 0;
//            const Term * best = 0;

            // list of all the terminators that have the best confidence (handles ties)
            vector<const Term *> best_list;

            sort(regterms.begin(), regterms.end(), by_left_side);

            SeqPtr left = 0, right = 0;

            for(vector<const Term*>::const_iterator T = regterms.begin();
                T != regterms.end();
                ++T)
            {
                // track the best overall
                int c = er_confidence(*this, _conf, **T);
                if(c >= best_conf || best_list.empty()) 
                {
                    if(c > best_conf)
                    {
                        best_list.clear();
                    }
                    best_list.push_back(*T);
 //                   best = *T;
                    best_conf = c;
                }

                // track the histogram of confidences
                _term_histo[c]++;

                // track the best confidence in the current component
                // note: will be false when left == right == 0
                if(interval_intersect((*T)->left(), (*T)->right(), left, right))
                {
                    // we extend the current component
                    right = max((*T)->right(), right);
                    comp_best = max(comp_best, c);
                }
                else
                {
                    if(left != 0 && right != 0) _comp_histo[comp_best]++;

                    // starting a new component
                    left = (*T)->left();
                    right = (*T)->right();
                    comp_best = c;
                }
            }

            if(left != 0 && right != 0) _comp_histo[comp_best]++;
    
            if(!regterms.empty())
            {
                assert(best_conf <= 100);
                _conf_histo[best_conf]++;
            }

            // given the choice between A and B, if A domaininates B, we
            // output *B* (since they all have the same confidnece, by
            // definiation)
            remove_dominating(best_list);

            // output all the terminators with the best confidence (each on its own line)
            for(vector<const Term*>::const_iterator T = best_list.begin();
                T != best_list.end();
                ++T)
            {
                _out << setw(10) << e.reg->name 
                     << " " << setw(3) << regterms.size();
                //if(!regterms.empty())
                //{ 
                    SeqPtr tplace = (_dir==REVERSE)?(*T)->left():(*T)->right();
                    SeqPtr gene_end = (_dir==REVERSE)?e.reg->right():e.reg->left();
                    
                    _out << " " << setw(3) << best_conf 
                         << " " << setw(3) << (gene_end - tplace)*_dir 
                         << " " << **T;
                //}
                _out << endl;
            }

            if(best_list.empty())
            {
                _out << setw(10) << e.reg->name 
                     << " " << setw(3) << regterms.size() << " NONE" << endl;
            }

            _good_region = false;
            EventResponder::leave_intergene(r, d, e);
        }
    }

private:
    vector<const Term*> _terms;
    int _region_size;
    ostream & _out;
    int _t2t_region_count;
};




void
t2t_hitanal(
    ostream & out, 
    const Genome & g, 
    Confidence & conf, 
    int min_span,
    bool show_all)
{
    DeHoonT2THits hits(min_span, conf, 500 + gene_start_cut, out);

    for(EVERY_CHROM_CONST(g, C))
    {
        scan_events(**C, hits, gene_start_cut, gene_end_cut);
        reverse_scan_events(**C, hits, gene_start_cut, gene_end_cut);
    }

    int total, hit = 0, numterms = 0, numcomp = 0;
    total = hits.t2tregion_count();

    out << endl << "SUMMARY" << endl << endl;
    out << total << " putative operon ends contain " << hits.good_terms() 
        << " possible terminators." << endl;
    out << "Percent operon ends hit by terminators of confidence >= x:" << endl;

    out << "\t" << "x" << "\t" << "#>=x" << "\t" << "hits>=x" << "\t" 
        << "%hits>=x" << endl;

    for(int i = 100; i >= 0; i--)
    {
        hit += hits.histo()[i];
        numterms += hits.term_histo()[i];
        numcomp += hits.comp_histo()[i];
        if(i % 10 == 0 || show_all)
        {
            out << "\t" << i << "\t" << numterms << "\t" 
                << numcomp << "\t" << hit << "\t" 
                << int(float(hit)/total*100+0.5) 
                << ((i==100)?" %":"") << endl;
        }
    }
    out << endl;
}


class Tail2TailScores : public EventResponder
{
public:
    Tail2TailScores(Confidence & conf) :
        ttscores(102),
        allscores(102),
        _t2t(-1),
        _conf(conf)
    {}

    void start(const Seq & seq, Direction dir)
    {
        _t2t = -1;
    }

    void terminator(const Term * term)
    {
        int c = er_confidence(*this, _conf, *term);

        allscores[c+1]++;
        if(in_t2t())
        {
            _t2t = max(_t2t, int(c));
        }
    }

    void leave_intergene(RegionType r, Direction d, const Event & e)
    {
        EventResponder::leave_intergene(r, d, e);

        if(r == TAIL2TAIL)
        {
            ttscores[_t2t+1]++;
            _t2t = -1;
        }
    }

    vector<int> ttscores;
    vector<int> allscores;

private:
    int _t2t;
    Confidence & _conf;
};


// write  the data for a plot to out. For every confidnece, count the # of
// terms with confidence >= x, and the number of Tail-to-tail regions hit by
// those terms.
void
plot_tthits_vs_terms(ostream & out, Confidence & conf, Genome & g)
{
    Tail2TailScores tts(conf);
    for(EVERY_CHROM_CONST(g, C)) scan_events(**C, tts, gene_start_cut, gene_end_cut);

    unsigned long allsum = 0, ttsum = 0;
    for(int i=101; i>=0; i--)
    {
        allsum += tts.allscores[i];
        ttsum += tts.ttscores[i];
        out << allsum << " " << ttsum << " " << i-1 << endl;
    }
}


int 
count_starts_in_genes(const Seq & seq, Direction dir)
{
    int num_starts = 0;
    for(EVERY_REGION_CONST(seq.genes, G))
    {
        int in_window = 0;
        char letter = ((*G)->dir() == dir)?'T':'A';

        SeqPtr s;
        int i;
        for(i = 0, s = (*G)->left(); s<= (*G)->right() && i < UWINDOW_SIZE; s++, i++)
        {
            if(*s == letter) in_window++;
        }
        SeqPtr leaving = (*G)->left();

        for(; s <= (*G)->right(); s++, leaving++)
        {
            if(in_window >= UWINDOW_REQUIRE) num_starts++;
            
            if(*s == letter) in_window++;
            if(*leaving == letter) in_window--;
        }
    }
    return num_starts;
}

