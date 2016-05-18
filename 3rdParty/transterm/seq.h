/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef SEQ_H
#define SEQ_H
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace std;

typedef double Energy;
enum Direction {REVERSE=-1, BIDIR=0, FORWARD=1};
enum RegionType { GENE, HEAD2TAIL, TAIL2TAIL, HEAD2HEAD };
typedef const char * SeqPtr;

class Seq;
class Confidence;

// Represents a region of a sequence [start,end] (both pointers are inclusive)
struct Region 
{
    string name;            // a description of the region
    SeqPtr start, end;      // pointers into seq at the start & end of the seq
    const Seq * seq;        // seq of which this is a region
    string desc;

    Region(const string & n, const Seq * seqq, SeqPtr s, SeqPtr e, 
        const string & d = "") : 
        name(n), start(s), end(e), seq(seqq), desc(d) {}

    Region() : name(""), start(0), end(0), seq(0) {}

    virtual ~Region() {}

    // Direction of the region. Regions always "run" from start to end
    // If end < start, the dreiction is "REVERSE"
    virtual Direction dir() const 
    { 
        return (end < start) ? REVERSE : FORWARD; 
    }

    // the 'left' of the region is always the 5' end of the sequence,
    // regardless of the direction of the region. The 'right', likewise, is
    // always the 3' end.
    SeqPtr left() const { return min(start,end); }
    SeqPtr right() const { return max(start,end); }

    int length() const { return abs(start-end) + 1; }
};


// represents a terminator region. the 'region' is defined as the hairpin
// sequence --- the tail 'region' is not included.
struct Term : public Region
{
    // terminator geometry
    int gap;
    int stem_len, loop_len;
    list<int> gaps;

    // terminator scores
    Energy hp_energy, tail_energy; 
    int conf;

    // a link to an 'equiv' term
    const Term * partner;
    list<Term*> opp_overlapping, overlapping;

    Term() : Region(), gap(0), stem_len(0), loop_len(0)
    {
        init0();
        lst = rst = 0;
    }

    Term(const Seq * s, Direction d, SeqPtr base, int sl, 
        int ll, int g) 
        : Region("", s, base - d*(geolength(sl, ll, g) - 1), base),
          gap(g), stem_len(sl), loop_len(ll), sense(d)
    {
        init0();
        rst = right() - stem_len + 1 - ((gap>0)?1:0);
        lst = right_stem_top() - loop_len - 1;
    }

    Term(const Seq * s, Direction d, SeqPtr lsb1, SeqPtr lst1, 
        SeqPtr rst1, SeqPtr rsb1, list<int> & glist, Energy hpe = 0) 
        : Region("", s, min(lsb1, rsb1), max(lsb1, rsb1))
    {
        init0();
        lst = min(lst1, rst1);
        rst = max(lst1, rst1);
        sense = d;
        hp_energy = hpe;

        loop_len = rst - lst - 1;
        stem_len = right() - rst + 1;

        // for backwards compatability, gap = true if there is a gap
        gap = (abs(lst - left() + 1) != stem_len) ? 1 : 0;

        gaps = glist;
    }

    virtual ~Term() 
    {
 //       if(gaps) 
 //       {
 //           delete gaps;
 //           gaps = 0;
 //       }
    }

    Direction dir() const { return sense; }
    Direction & dir() { return sense; }

    // accessors for the energy assigned to this terminator
    enum EnergyKind {TAIL, HAIRPIN};

    Energy energy(EnergyKind k) const 
    { 
        return (k==HAIRPIN)?hp_energy:tail_energy; 
    }

    // Pointers into the geometry of the terminator
    SeqPtr left_stem_base() const { return left(); }
    SeqPtr right_stem_base() const { return right(); }
    SeqPtr left_stem_top() const { return lst; }
    SeqPtr right_stem_top() const { return rst; }

private:
    // return the size of the hairpin with the given geometry
    int
    geolength(int stem_len, int loop_len, int gap)
    {
        return 2*stem_len + loop_len + ((gap!=0)?1:0);
    }

    // set the things that shd be 0 to 0
    void
    init0()
    {
        hp_energy = tail_energy = 0.0;
        conf = 0;
        partner = 0;
    }

    Direction sense;
    SeqPtr lst, rst;
};

typedef vector<const Term *> ConstTermVec;

// Represents a sequence
struct Seq 
{
    string name;            
    string desc;            
    unsigned long length;   // number of characters in the seq
    char * dna;             // pointer to the sequence data
    vector<Term*> terms;     // list of Term features
    vector<Region*> genes;   // list of gene features

    Seq() : name(""), desc(""), length(0), dna(0) {}

    ~Seq() { clear(); }

    void clear();
    
    SeqPtr left() const { return dna; }
    SeqPtr right() const { return dna + length - 1; }
};

typedef vector<Seq*> Genome;


// represents an event happening in the sequence.
struct Event
{
    const Region * reg;
    // for paired events, extent is the location of the other event
    SeqPtr place, extent; 
    enum Kind 
    { 
        Terminator,         
        ForwardGeneStart,  
        ForwardGeneEnd,  
        ReverseGeneStart,
        ReverseGeneEnd
    } kind;

    Event() : reg(0), place(0), extent(0) {}
    Event(const Region * r, SeqPtr p, Kind k, SeqPtr x=0) : 
        reg(r), place(p), extent(x), kind(k) {}
};

typedef vector<Event>::const_iterator event_iterator;

class EventResponder 
{
public:
    virtual ~EventResponder() {}

    virtual void start(const Seq & seq, Direction dir) 
    { 
        _fwd_gene = _rvs_gene = 0; 
    }
    virtual void end() {}
    virtual void event(const Event & e) {}
    virtual void terminator(const Term * term) {}

    // if these are called from the subclass, they'll manage gene_count()
    virtual void enter_gene(const Event & e);
    virtual void leave_gene(const Event & e); 

    virtual void 
    enter_intergene(RegionType r, Direction d, const Event & e)
    {
        set(r, d, true);
    }

    virtual void 
    leave_intergene(RegionType r, Direction d, const Event & e)
    {
        set(r, d, false);
    }

    friend int er_confidence(const EventResponder &, const Confidence &, const Term &);

protected:
    // you can't create a plain EventResponder --- must subclass
    EventResponder() 
    { 
        _fwd_gene = _rvs_gene = 0; 
        _t2t = _h2t_fwd = _h2t_rvs = _h2h = false; 
    } 

    bool in_t2t() const { return _t2t; }
    bool in_h2h() const { return _h2h; }
    bool in_h2t_fwd() const { return _h2t_fwd; }
    bool in_h2t_rvs() const { return _h2t_rvs; }

    int gene_count() const { return _fwd_gene + _rvs_gene; }
    int rvs_gene_count() const { return _rvs_gene; }
    int fwd_gene_count() const { return _fwd_gene; }

private:
    void set(RegionType r, Direction d, bool tf);

    int _fwd_gene, _rvs_gene;
    bool _t2t, _h2t_fwd, _h2t_rvs, _h2h;
};


#define EVERY_CHROM(ch, C) \
    Genome::iterator C = ch.begin(); C != ch.end(); ++C

#define EVERY_CHROM_CONST(ch, C) \
    Genome::const_iterator C = ch.begin(); C != ch.end(); ++C

#define EVERY_REGION(vec, R) \
    vector<Region*>::iterator R = vec.begin(); R != vec.end(); ++R

#define EVERY_REGION_CONST(vec, R) \
    vector<Region*>::const_iterator R = vec.begin(); R != vec.end(); ++R

#define EVERY_TERM(vec, T) \
    vector<Term*>::iterator T = vec.begin(); T != vec.end(); ++T

#define EVERY_TERM_CONST(vec, T) \
    vector<Term*>::const_iterator T = vec.begin(); T != vec.end(); ++T

#define EVERY_CTERM_CONST(vec, T) \
    vector<const Term*>::const_iterator T = vec.begin(); T != vec.end(); ++T

const char PADDING_CHAR = 'x';

int seqindex(const Seq &, SeqPtr);
string subseq(SeqPtr, SeqPtr);
void read_seqs(istream &, Genome &);
void pad_seqs(Genome &, int);
void pad_seq(Seq &, int);
bool region_isleftof(const Region *, const Region *);
bool hp_overlap(const Term &, const Term &);
bool dominates(const Term &, const Term &);
string dir_str(Direction);
void scan_events(const Seq &, EventResponder &, int, int);
void reverse_scan_events(const Seq &, EventResponder &, int, int);
void sort_genes(Genome &);
Seq * chrom_for_id(Genome &, const string &);
void print_term_seq(ostream &, const Term &);
ostream & operator<<(ostream &, const Term &);
#endif
