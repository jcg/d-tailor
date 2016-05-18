/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include "map-output.h"
#include "conf.h"
#include "seq.h"
#include "util.h"
#include "transterm.h"


void
output_header(ostream & out)
{
    out << endl;
    out << "Each terminator entry starts in column 3 and is of the form:" << endl;
    out << left  << setw(11) << "  term #" << " "
        << right << setw(8) << "start" << " - " 
        << left  << setw(8) << "end" << " " 
        << "+/-" << " " << " region";

    out << right << setw(4) << "conf" << " "
        << setw(4) << "hp"   << " "
        << setw(8) << "tail" << " | notes";

    out << endl;
    out << "Followed by the sequence of the 5' tail, 5' stem, loop, 3' stem, and 3' tail." << endl;
    out << "Genes are interspersed, and start the first column." << endl;

    out << endl << endl;
}


// output a description of the gene
void
output_gene(ostream & out, const Region & gene)
{
    out << left << setw(12) << gene.name << " "
        << right << setw(8) << seqindex(*gene.seq, gene.start) << " - "
        << left << setw(8) << seqindex(*gene.seq, gene.end) << " "
        <<  dir_str(gene.dir()) 
        << " | " << gene.desc << endl;
}


// output the info about a terminator
void
output_terminator(
    ostream & out, 
    const Term & term, 
    int conf, 
    int count,
    bool print_seq,
    bool in_fwd_gene,
    bool in_rvs_gene,
    bool int2t,
    bool inh2tforward,
    bool inh2treverse,
    bool inh2h)
{
    //if(print_seq) out << endl;
    ostringstream oss;
    oss << "  TERM " << count;

    SeqPtr ll, rr;
    if(term.dir() == FORWARD)
    {
        ll = term.left();
        rr = term.right();
    }   
    else
    {
        ll = term.right();
        rr = term.left();
    }

    // output the type of region we are considered to be in
    string loc = "";
    if((in_fwd_gene && term.dir() == FORWARD) ||
       (in_rvs_gene && term.dir() == REVERSE)) loc += 'G'; // in a gene, coding strand

    if((in_fwd_gene && term.dir() == REVERSE) ||
       (in_rvs_gene && term.dir() == FORWARD)) loc += 'g'; // in a gene, noncoding

    if(int2t) loc += 'T';

    if(inh2tforward) loc += ((term.dir() == FORWARD)?'F':'f');
    if(inh2treverse) loc += ((term.dir() == REVERSE)?'R':'r');
    if(inh2h) loc += "H";

    if(loc == "") loc = "N"; 

    // the location & direction
    out << left  << setw(12) << oss.str() << " "
        << right << setw(8) << seqindex(*term.seq, ll) << " - " 
        << left  << setw(8) << seqindex(*term.seq, rr) << " " 
        << dir_str(term.dir()) << " ";
    out << setw(3) << loc << " ";

    // the scores
    out << right << setw(3) << conf << " "
        << setw(5) << term.hp_energy << " "
        << setw(8) << term.tail_energy << " | ";

    
    // the notes
    string comma = "";
#if 0
    if(term.partner) 
    {
        out << "bidir";
        comma = ", ";
    }
#endif

    if(term.gap != 0)
    {
        out << comma << "gap " << term.gaps.size();
        comma = ", ";
    }

    if(!term.opp_overlapping.empty())
    {
        out << comma << "opp_overlap";
        for(list<Term*>::const_iterator T = term.opp_overlapping.begin();
            T != term.opp_overlapping.end();
            ++T)
        {
            out << " " << seqindex(*(*T)->seq, (*T)->left()); 
        }
        comma = ", ";
    }

    if(!term.overlapping.empty())
    {
        out << comma << "overlap";
        for(list<Term*>::const_iterator T = term.overlapping.begin();
            T != term.overlapping.end();
            ++T)
        {
            out << " " << seqindex(*(*T)->seq, (*T)->left()); 
        }
        comma = ", ";
    }

#if 0
    // mark terminators that got 0 conf b/c their direction was inconsistent
    // with the region's direction
    if(!int2t && !(inh2tforward && inh2treverse) &&
       (inh2tforward && term.dir() == REVERSE ||
        inh2treverse && term.dir() == FORWARD)) out << comma << "wrong dir";
#endif

    out << endl;

    if(print_seq)
    {
        print_term_seq(out, term);
        out << endl;
    }
}


// used from output_map() in a call to scan_events(), which will send messages
// to this EventResponder when scanning a sequence from left to right.
class MapOutputer : public EventResponder
{
public:
    MapOutputer(ostream & out, Confidence & conf, 
        int co=90, bool print_seq=true, bool only_good_context=false) : 
        EventResponder(),
        _out(out),
        _conf(conf),
        _conf_cutoff(co), 
        _print_seq(print_seq),
        _only_good_context(only_good_context) {}

    virtual ~MapOutputer() {}

    // output the name of the seq and its length
    void start(const Seq & seq, Direction dir)
    {
        EventResponder::start(seq, dir);

        _out << "SEQUENCE " << seq.name << " " << seq.desc 
             << " (length " << seq.length << ")" << endl << endl;

        _last_was_term = false;
        _term_count = 0;
    }

    // when we see a terminator, compute its conf given its genomic context
    // and output it if it passes the cutoff
    void terminator(const Term * term)
    {
        EventResponder::terminator(term);

        int c = er_confidence(*this, _conf, *term);
        if(c >= _conf_cutoff && (!_only_good_context || good_context(term)))
        {
            output_terminator(_out, *term, c, ++_term_count, _print_seq, 
               fwd_gene_count()>0, rvs_gene_count()>0, in_t2t(), 
               in_h2t_fwd(), in_h2t_rvs(), in_h2h());
            _last_was_term = true;
        }
    }

    // for gene /ends/ (either -> or <-) we output the gene.
    void event(const Event & e)
    {
        if(e.kind == Event::ForwardGeneEnd ||
           e.kind == Event::ReverseGeneEnd)
        {
            //if(_last_was_term && _print_seq) _out << endl;
            _last_was_term = false;
            output_gene(_out, *e.reg);
        }
    }

    // output some newlines to end the printout
    void end()
    {
        _out << endl;
        // if(!_last_was_term || !_print_seq) _out << endl;
    }

    bool
    good_context(const Term * t)
    {
        return in_t2t() || 
            (t->dir() == FORWARD && in_h2t_fwd()) ||
            (t->dir() == REVERSE && in_h2t_rvs());
    }
    
private:
    ostream & _out;
    Confidence & _conf;
    int _conf_cutoff;
    bool _print_seq;
    bool _only_good_context;
    bool _last_was_term;
    int _term_count;
};


// output to out a map of each chrom seq in genome g
void
output_map(
    ostream & out, 
    const Genome & g, 
    Confidence & conf, 
    int conf_cutoff,
    bool print_seq,
    bool only_good_context)
{
    MapOutputer mo(out, conf, conf_cutoff, print_seq, only_good_context);
    output_header(out);
    for(EVERY_CHROM_CONST(g, C)) scan_events(**C, mo, gene_start_cut, gene_end_cut);
}

// |-----> ... |---> 
// |----->   <-----|

class BestAfterGene : public EventResponder
{
public:
    BestAfterGene(ostream & out, const Confidence & conf) : _out(out), _conf(conf) {}
    virtual ~BestAfterGene() {}
    void start(const Seq & seq, Direction dir)
    {
        EventResponder::start(seq, dir);
        _dir = dir;
        _best_terms.clear();
        _best_confs.clear();
        _best_dist.clear();
        _gene_names.clear();
        _name = "";
        _gene_end = 0;
        _in_intergene = false;
        //_best_term = 0;
        //_best_conf = 0;
    }

    void terminator(const Term * t)
    {
        EventResponder::terminator(t);

        SeqPtr term_pos = (_dir == FORWARD)?t->left():t->right();

        // if the terminator is facing the right way, if we have a gene that it
        // might be terminating, and its near enough to that gene
        // we continue until we leave the intergenic region or are 500bp from 
        // the end of the gene, whichever is FARTHER
        bool in_gene = ((_dir==FORWARD)?fwd_gene_count():rvs_gene_count())>0;
        if(t->dir() == _dir && _name != "" && 
           (abs(term_pos - _gene_end) <= 525) && !in_gene) //|| _in_intergene))
        {
            int c = er_confidence(*this, _conf, *t);
            if(c > _best_confs.back() || _best_terms.back() == 0)
            {
                _best_terms.back() = t;
                _best_confs.back() = c;
                _best_dist.back() = int(term_pos - _gene_end) * _dir;
            }
        }
    }

    void leave_gene(const Event & e)
    {
        EventResponder::leave_gene(e);
        if((_dir == FORWARD && e.kind == Event::ForwardGeneEnd) ||
           (_dir == REVERSE && e.kind == Event::ReverseGeneEnd))
        {
            _gene_end = e.place;

            // we append & + the position of the gene b/c gene names
            // are not unique
            ostringstream oss;
            oss << seqindex(*e.reg->seq, _gene_end);
            //oss << int(_gene_end);

            _name = e.reg->name + "&" + oss.str();

            _gene_names.push_back(_name);
            _best_terms.push_back(0);
            _best_confs.push_back(0);
            _best_dist.push_back(0);
            _in_intergene = true;
        }
    }

    void enter_gene(const Event & e)
    {
        EventResponder::enter_gene(e);
#if 0
        _name = "";
        _gene_end = 0;
#endif
        // replace the above two lines with the following commented lines
        // if you want to only stop when a co-directed gene start is encountered
        if((_dir == FORWARD && e.kind == Event::ForwardGeneStart) ||
           (_dir == REVERSE && e.kind == Event::ReverseGeneStart))
        {
            _name = "";
            _gene_end = 0;
        }
        else
        {
            _in_intergene = false;
        }
    }

    void end()
    {
        EventResponder::end();
        for(unsigned i = 0; i < _gene_names.size(); i++)
        {
            int best_conf = _best_confs[i];
            const Term * best_term = _best_terms[i];

            string name = _gene_names[i].substr(0,_gene_names[i].rfind('&'));

            _out << setw(12) << name << " "; // << dir_str(_dir) << " ";
            if(best_term)
            {
                _out << *best_term;
#if 0
                _out << setw(8) << seqindex(*best_term->seq, best_term->left()) << " - " 
                     << setw(8) << seqindex(*best_term->seq, best_term->right()) << " " 
                     << setw(3) << best_conf << " ";
                print_term_seq(_out, *best_term);  
#endif
                _out << " " << setw(3) << best_conf << " " << _best_dist[i];
            }
            else
            {
                _out << "NONE";
            }
            _out << endl;
        }
    }

private:
    Direction _dir;
    vector<string> _gene_names;
    vector<int> _best_confs;
    vector<int> _best_dist;
    vector<const Term*> _best_terms;

    SeqPtr _gene_end;
    string _name;
    ostream & _out;
    const Confidence & _conf;
    bool _in_intergene;
};

void
output_best_term(ostream & out, const Confidence & conf, const Genome & g)
{
    BestAfterGene bag(out, conf);    
    for(EVERY_CHROM_CONST(g, C))
    {
        scan_events(**C, bag, gene_start_cut, gene_end_cut);
        reverse_scan_events(**C, bag, gene_start_cut, gene_end_cut);
    }
}

