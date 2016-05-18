/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <cctype>
#include <vector>
#include <iomanip>
#include <algorithm>

#include "util.h"
#include "seq.h"

const size_t INIT_SEQ_ALLOC = 1000000;

// free the memory reserved for the sequence data. we clear the name
// and descriptions too
void 
Seq::clear()
{
    name = desc = "";
    if(dna) free(dna);
    dna = 0;
    length = 0;
    terms.clear();
    genes.clear();
}


// convert a direction into a string
string
dir_str(Direction dir)
{
    switch(dir)
    {
        case FORWARD: return "+";
        case REVERSE: return "-";
        case BIDIR: return "+/-";
    }
    return "";
}


// return the 1-based seqindex for the position cp
int
seqindex(const Seq & seq, SeqPtr cp)
{
    return cp - seq.left() + 1;
}


// extract a substring of a character array as a string
// inclusive: [cp1, cp2]
string
subseq(SeqPtr cp1, SeqPtr cp2)
{
    string str = "";
    for(; *cp1 && cp1 <= cp2; cp1++) str += *cp1;
    return str;
}


// used by read_seq to reallocate the space for the sequence
void
resize_buffer(Seq & seq, size_t bufsize)
{
    seq.dna = (char*)realloc(seq.dna, bufsize*sizeof(char));
    if(!seq.dna)
    {
        cerr << "Couldn't allocate enough memory for sequence." << endl;
        exit(3);
    }
}


// Read a single sequence from the input stream (fasta format)
// return a ptr to the new sequence object
Seq *
read_seq(istream & in)
{
    int i;

    Seq * seq = new Seq();

    // find first line starting with >
    while((i = in.get()) && i != EOF && i != '>') { }
    if(i == EOF) return false;

    // read the name
    while((i = in.get()) && i != EOF && !isspace(i))
    {
        seq->name += (char)i;
    }
    if(i == EOF) return false; 

    if(seq->name.length() == 0) 
    {
        cerr << "Sequence has no name on the '>' fasta line." << endl;
        exit(3);
    }

    // if i == '\n' then there is no description
    seq->desc = "";
    if(i != '\n')
    {
        // read the desc
        while((i = in.get()) && i != EOF && i != '\n')
        {
            seq->desc += (char)i;
        }
    }

    if(i == EOF) return false;

    // allocate initial buffer of 1Mb for seq
    size_t bufsize = INIT_SEQ_ALLOC;
    resize_buffer(*seq, bufsize);

    // read until EOF or '>'
    unsigned long cc = 0;
    while((i = in.get()) && i != EOF && i != '>')
    {
        // any whitespace is ignored
        if(isspace(i)) continue;
        if(i=='>')
        {
            cerr << "TransTermHP does not support multiple sequences in a single"
                 << " FASTA file. Spilt the sequences into multiple files." << endl;
            exit(3);
        }

        // if we've run out of space, double it
        if(cc > bufsize)
        {
            bufsize *= 2;
            resize_buffer(*seq, bufsize);
        }
        seq->dna[cc++] = (char)toupper(i);
    }

    // terminate sequence with '\0' && free unused memory
    seq->length = cc;
    seq->dna[cc++] = 0;
    resize_buffer(*seq, cc);

    // unget the stopping character
    if(i != EOF) in.putback(i);

    return seq;
}


// read all the sequences in the given file
void
read_seqs(istream & in, Genome & seqs)
{
    Seq * seq;
    while((seq = read_seq(in)))
    {
        seqs.push_back(seq);
    }
}


// add padding ' ' characters to each end of the sequence
void
pad_seq(Seq & dna, int padding)
{
    char * newseq = (char*)calloc(dna.length + 2*padding, sizeof(char));
    memset(newseq, PADDING_CHAR, padding);
    memcpy(newseq + padding, dna.dna, dna.length);
    memset(newseq + padding + dna.length, PADDING_CHAR, padding);
    free(dna.dna);
    dna.dna = newseq;
    dna.length += 2 * padding;
}


// add padding to each of the sequences in the genome
void
pad_seqs(Genome & chroms, int padding)
{
    for(Genome::iterator C = chroms.begin();
        C != chroms.end();
        ++C)
    {
        pad_seq(**C, padding);
    }
}


// for sorting: return true if region t1 starts to the left of t2
bool
region_isleftof(const Region * t1, const Region * t2)
{
    return t1->left() < t2->left();
}


// return true if the hairpin regions overlap
bool 
hp_overlap(const Term & t1, const Term & t2)
{
    SeqPtr b1 = t1.left_stem_base();
    SeqPtr e1 = t1.right_stem_base();
    SeqPtr b2 = t2.left_stem_base();
    SeqPtr e2 = t2.right_stem_base();
    return (b1>=b2 && b1 <= e2) || (e1>=b2 && e1<=e2) || (b1<=b2 && e1>=e2);
}

// return true if t1 dominates t2
bool
dominates(const Term & t1, const Term & t2)
{
    return t1.left() <= t2.left() && t1.right() >= t2.right();
}

//======================================================================
// TERMINATOR SEQUENCE OUTPUT
//======================================================================

string
concat_dir(const string & str1, char ch, Direction dir)
{
    if(dir == REVERSE)
    {
        return str1 + ch;
    }
    else
    {
        return ch + str1;
    }
}


void
popgaps(string & out, list<int> & gaps, int i, Direction dir)
{
    while(!gaps.empty() && abs(gaps.front()) == i)
    {
        out = concat_dir(out, '-', dir);
        gaps.pop_front();
    }
}

bool
abs_int_cmp(int a, int b)
{
    return abs(a) < abs(b);
}

string
addgaps(const string & str, const Term & term)
{
    if(term.gaps.size() == 0) return str;

#if 0
    cout << "ORIG: " << str << endl;
    for(list<int>::const_iterator G = term.gaps.begin();
        G != term.gaps.end();
        ++G)
    {
        cout << " " << *G;
    }
    cout << endl;
#endif

    // copy and sort the list of gaps
    list<int> gaps(term.gaps);
    gaps.sort(abs_int_cmp);

    Direction dir = term.dir();
    int starti, endi;
    if(dir == REVERSE)
    {
        starti = 0;
        endi = str.length()-1;
    }
    else
    {
        starti = str.length()-1;
        endi = 0;
    }

    // copy every character in the input str over to the output string
    // possibly inserting gaps as required by the gaps list
    int reali = 0;
    string out = "";
    for(int i = starti; ; i -= dir)
    {
        if(str[i] == ' ') 
        {
            out = concat_dir(out, ' ', dir);
            continue;
        }

        if(gaps.empty() || gaps.front() <= 0)
        {
            popgaps(out, gaps, reali, dir);
            out = concat_dir(out, str[i], dir);
        }
        else
        {
            out = concat_dir(out, str[i], dir);
            popgaps(out, gaps, reali, dir);
        }

        reali++;
        if(i == endi) break;
    }
    return out;
}


string
term_seq(const Term & term, bool withgaps)
{
    string lefttail, seq, righttail;

    lefttail = subseq(term.left_stem_base()-15, term.left_stem_base()-1);

    seq = subseq(term.left_stem_base(), term.left_stem_top()) + " " +
         subseq(term.left_stem_top()+1,term.right_stem_top()-1) + " " +
         subseq(term.right_stem_top(), term.right_stem_base());

    righttail =  subseq(term.right_stem_base()+1, term.right_stem_base()+15);

    if(withgaps) seq = addgaps(seq, term);
    
    return lefttail + " " + center(seq, 45) + " " + righttail;
}


// output the tails, stem and loop sequence of a terminator
void
print_term_seq(ostream & out, const Term & term)
{
    out << "  " << term_seq(term, true);
#if 0
    out << "  " 
      << subseq(term.left_stem_base()-15, term.left_stem_base()-1) << " "
      << center(subseq(term.left_stem_base(), term.left_stem_top()) + " " +
         subseq(term.left_stem_top()+1,term.right_stem_top()-1) + " " +
         subseq(term.right_stem_top(), term.right_stem_base()), 45) << " "
      << subseq(term.right_stem_base()+1, term.right_stem_base()+15);
#endif
}

// mostly for debugging, output a terminator to a stream
ostream &
operator<<(ostream & out, const Term & t)
{
    out << setw(7) << seqindex(*t.seq, t.left()) << " .. " 
        << setw(7) << seqindex(*t.seq, t.right()) 
        << " " << setw(1) << dir_str(t.dir()) << " " 
        << setw(5) << t.hp_energy << " " << setw(9) << t.tail_energy << " "; 
    print_term_seq(out, t);
    return out; // << endl;
}

//=============================================================================
// Event Functions
//=============================================================================

// true if event starts a region when scanning to the right
bool
is_gene_left_end(const Event & e)
{
    return e.kind == Event::ForwardGeneStart || 
           e.kind == Event::ReverseGeneEnd;
}


bool
is_gene_right_end(const Event & e)
{
    return e.kind == Event::ForwardGeneEnd ||
           e.kind == Event::ReverseGeneStart;
}


bool
is_rvs_gene_event(const Event & e)
{
    return e.kind == Event::ReverseGeneStart ||
           e.kind == Event::ReverseGeneEnd;
}


bool
is_fwd_gene_event(const Event & e)
{
    return e.kind == Event::ForwardGeneStart ||
           e.kind == Event::ForwardGeneEnd;
}


// track the number of fwd and rvs genes that we're in
void
EventResponder::enter_gene(const Event & e)
{
    if(is_fwd_gene_event(e)) 
    {
        _fwd_gene++;
    }
    else if(is_rvs_gene_event(e))
    {
        _rvs_gene++;
    }
    else assert(false);
}


// track the number of fwd and rvs genes that we're in
void
EventResponder::leave_gene(const Event & e)
{
    if(is_fwd_gene_event(e)) 
    {
        assert(_fwd_gene > 0);
        _fwd_gene--;
    }
    else if(is_rvs_gene_event(e))
    {
        assert(_rvs_gene > 0);
        _rvs_gene--;
    }
    else assert(false);
}


// (private) set the correct flags given the event type seen
void 
EventResponder::set(RegionType r, Direction d, bool tf)
{
    if(r == TAIL2TAIL)
    {
        _t2t = tf;
    } 
    else if(r == HEAD2TAIL && d == FORWARD)
    {
        _h2t_fwd = tf;
    }
    else if(r == HEAD2TAIL && d == REVERSE)
    {
        _h2t_rvs = tf;
    }
    else if(r == HEAD2HEAD)
    {
        _h2h = tf;
    }
}


// for sorting by event place
bool
by_event_location(const Event & e1, const Event & e2)
{
    return e1.place < e2.place;
}

// for sorting by event place
bool
reverse_by_event_location(const Event & e1, const Event & e2)
{
    return e1.place > e2.place;
}

// create an unsorted list of events
void
populate_events(
    const Seq & seq, 
    vector<Event> & events, 
    int start_cut_in,
    int end_cut_in,
    Direction dir)
{
    const int MIN_GENE_SIZE = 50; // never shrink genes smaller than this
    
    for(EVERY_REGION_CONST(seq.genes, G))
    {
        SeqPtr start, end;
        int gene_len;
        int start_cut, end_cut;

//        cut = max(0, (signed((*G)->length()) - signed(MIN_GENE_SIZE)) / 2);
//        cut = min(start_cut_in, cut);
        
        // first we make the full start cut (never making the gene less than MIN_GENE_SIZE)
        // then we make the end cut, never making the gene less than MIN_GENE_SIZE. In other
        // words, start_cut gets preference over end_cut
        gene_len = (*G)->length();
        start_cut = min(start_cut_in, max(0, gene_len - signed(MIN_GENE_SIZE)));
        gene_len -= start_cut;
        end_cut = min(end_cut_in, max(0, gene_len - signed(MIN_GENE_SIZE)));

//        cout << "CUT: " << (*G)->name << " " << (*G)->length() << " "
//             << cut << " " << start_cut << " " << end_cut << endl;

        if((*G)->dir() == FORWARD)
        {
            start = (*G)->start + start_cut;
            end = (*G)->end - end_cut;

            events.push_back(Event(*G, start, Event::ForwardGeneStart, end));
            events.push_back(Event(*G, end, Event::ForwardGeneEnd, start));
        }
        else
        {
            start = (*G)->start - start_cut;
            end = (*G)->end + end_cut;

            events.push_back(Event(*G, start, Event::ReverseGeneStart, end));
            events.push_back(Event(*G, end, Event::ReverseGeneEnd, start));
        }
    }

    // when scanning in reverse, we use the right() end of terminators
    for(EVERY_TERM_CONST(seq.terms, T))
    {
        events.push_back(Event(*T, 
            (dir==REVERSE) ? ((*T)->right()) : ((*T)->left()), 
            Event::Terminator));
    }

}



// return rightmost event of type k between E and end (exclusive) 
// that can be reached /without/ crossing a full gene.
event_iterator
rightmost_nocross(event_iterator E, event_iterator end, Event::Kind k)
{
    event_iterator rightmost = end;
    SeqPtr right = E->reg->seq->right();

    for(++E; E != end && E->place <= right; ++E)
    {
        if(is_gene_left_end(*E)) right = min(right, E->extent);
        if(E->kind == k) rightmost = E;
    }
    return rightmost;
}


// return leftmost event of type k between E and end (exclusive) that can be
// reached /without/ crossing a full gene. (events must be sorted from right to
// left)
event_iterator
leftmost_nocross(event_iterator E, event_iterator end, Event::Kind k)
{
    event_iterator leftmost = end;
    SeqPtr left = E->reg->seq->left();

    for(++E; E != end && E->place >= left; ++E)
    {
        if(is_gene_right_end(*E)) left = max(left, E->extent);
        if(E->kind == k) leftmost = E;
    }
    return leftmost;
}


// process the events for the sequence left to right 
void
scan_events(const Seq & seq, EventResponder & er, int start_cut, int end_cut)
{
    // create a event 'queue' that we'll process
    vector<Event> events;
    populate_events(seq, events, start_cut, end_cut, FORWARD);

    sort(events.begin(), events.end(), by_event_location);

    // these mark when the regions should end. If they equal events.end() then
    // we aren't in the region of the given type
    event_iterator h2t_fwd_end = events.end();
    event_iterator h2t_rvs_end = events.end();
    event_iterator t2t_end = events.end();
    event_iterator h2h_end = events.end();

    er.start(seq, FORWARD);

    for(event_iterator E = events.begin(); E != events.end(); ++E)
    {
        event_iterator e;

        er.event(*E); 

        if(E == h2t_fwd_end) 
        {
            er.leave_intergene(HEAD2TAIL, FORWARD, *h2t_fwd_end);
            h2t_fwd_end = events.end();
        }

        if(E == h2t_rvs_end) 
        {
            er.leave_intergene(HEAD2TAIL, REVERSE, *h2t_rvs_end);
            h2t_rvs_end = events.end();
        }

        if(E == t2t_end) 
        {
            er.leave_intergene(TAIL2TAIL, BIDIR, *t2t_end);
            t2t_end = events.end();
        }

        if(E == h2h_end)
        {
            er.leave_intergene(HEAD2HEAD, BIDIR, *h2h_end);
            h2h_end = events.end();
        }

        // enter and leave genes
        switch(E->kind)
        {
            case Event::ForwardGeneStart: er.enter_gene(*E); break;
            case Event::ReverseGeneEnd: er.enter_gene(*E); break;
            case Event::ForwardGeneEnd: er.leave_gene(*E); break;
            case Event::ReverseGeneStart: er.leave_gene(*E); break;
            default: break;
        }

        switch(E->kind)
        {
            case Event::Terminator: er.terminator(((Term*)E->reg)); break;

            case Event::ForwardGeneEnd:
                // Either a H2T forward or a t2t region can be begun with ->
                e = rightmost_nocross(E, events.end(), Event::ForwardGeneStart);
                if(h2t_fwd_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2TAIL, FORWARD, *E);
                }
                h2t_fwd_end = e;

                e = rightmost_nocross(E, events.end(), Event::ReverseGeneEnd);
                if(t2t_end == events.end() && e != events.end())
                {
                    er.enter_intergene(TAIL2TAIL, BIDIR, *E);
                }
                t2t_end = e;
                break;

            case Event::ReverseGeneStart:
                // only a H2T reverse or H2H can start at a --|
                e = rightmost_nocross(E, events.end(), Event::ReverseGeneEnd);
                if(h2t_rvs_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2TAIL, REVERSE, *E);
                }
                h2t_rvs_end = e;

                // head to heads start with --|     end with |--
                e = rightmost_nocross(E, events.end(), Event::ForwardGeneStart);
                if(h2h_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2HEAD, BIDIR, *E);
                }
                h2h_end = e;
                break;

            default: break;
        }
    }
    
    er.end();
}


// process the events for the sequence right to left
void
reverse_scan_events(const Seq & seq, EventResponder & er, int start_cut, int end_cut)
{
    // create a event 'queue' that we'll process
    vector<Event> events;
    populate_events(seq, events, start_cut, end_cut, REVERSE);

    sort(events.begin(), events.end(), reverse_by_event_location);

    // these mark when the regions should end. If they equal events.end() then
    // we aren't in the region of the given type
    event_iterator h2t_fwd_end = events.end();
    event_iterator h2t_rvs_end = events.end();
    event_iterator t2t_end = events.end();
    event_iterator h2h_end = events.end();

    er.start(seq, REVERSE);

    for(event_iterator E = events.begin(); E != events.end(); ++E)
    {
        event_iterator e;

        er.event(*E); 

        if(E == h2t_fwd_end) 
        {
            er.leave_intergene(HEAD2TAIL, FORWARD, *h2t_fwd_end);
            h2t_fwd_end = events.end();
        }

        if(E == h2t_rvs_end) 
        {
            er.leave_intergene(HEAD2TAIL, REVERSE, *h2t_rvs_end);
            h2t_rvs_end = events.end();
        }

        if(E == t2t_end) 
        {
            er.leave_intergene(TAIL2TAIL, BIDIR, *t2t_end);
            t2t_end = events.end();
        }

        if(E == h2h_end)
        {
            er.leave_intergene(HEAD2HEAD, BIDIR, *h2h_end);
            h2h_end = events.end();
        }

        switch(E->kind)
        {
            case Event::ForwardGeneEnd: er.enter_gene(*E); break;
            case Event::ReverseGeneStart: er.enter_gene(*E); break;
            case Event::ForwardGeneStart: er.leave_gene(*E); break;
            case Event::ReverseGeneEnd: er.leave_gene(*E); break;
            default: break;
        }

        switch(E->kind)
        {
            case Event::Terminator: er.terminator(((Term*)E->reg)); break;

            case Event::ForwardGeneStart:
                // Either a H2T forward or a t2t region can be begun with ->
                e = leftmost_nocross(E, events.end(), Event::ForwardGeneEnd);
                if(h2t_fwd_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2TAIL, FORWARD, *E);
                }
                h2t_fwd_end = e;

                e = leftmost_nocross(E, events.end(), Event::ReverseGeneStart);
                if(h2h_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2HEAD, BIDIR, *E);
                }
                h2h_end = e;
                break;

            case Event::ReverseGeneEnd:
                // only a H2T reverse can start at a --|
                e = leftmost_nocross(E, events.end(), Event::ReverseGeneStart);
                if(h2t_rvs_end == events.end() && e != events.end())
                {
                    er.enter_intergene(HEAD2TAIL, REVERSE, *E);
                }
                h2t_rvs_end = e;

                e = leftmost_nocross(E, events.end(), Event::ForwardGeneEnd);
                if(t2t_end == events.end() && e != events.end())
                {
                    er.enter_intergene(TAIL2TAIL, BIDIR, *E);
                }
                t2t_end = e;
                break;

            default: break;
        }
    }
    
    er.end();
}


// given an id, try to find the rigth chromosome (seq) in the genome
// b/c of poor consistency (aka no consistency) in the naming schemes, 
// we try a bunch of heuristics.
Seq *
chrom_for_id(Genome & g, const string & id)
{
    vector<string> vec;

    // try exact match first
    for(EVERY_CHROM(g, C))
    {
        if((*C)->name == id) return *C;
    }

    // we also try to find a chrom with junk|cmr:ID|junk
    for(EVERY_CHROM(g, C))
    {
        split((*C)->name, '|', vec);
        for(unsigned i = 0; i < vec.size(); ++i)
        {
            if(vec[i].substr(0,4) == "cmr:" &&
               vec[i].substr(4) == id) return *C;
        }
    }

    // next we try to find gi|ID or gb|ID or gb|ID.junk
    for(EVERY_CHROM(g, C))
    {
        split((*C)->name, '|', vec);
        for(unsigned i = 0; i < vec.size() - 1; ++i)
        {
            if((vec[i] == "gi" && vec[i+1] == id) ||
               ((vec[i] == "gb" || vec[i] == "ref") && 
                (vec[i+1] == id || vec[i+1].substr(0,vec[i+1].rfind('.')) == id)))
            {
                return *C;
            }
        }
    }

    // finally, we check the first part of x|y|z to see if x = junk:ID
    for(EVERY_CHROM(g, C))
    {
        split((*C)->name, '|', vec);
        if(!vec.empty() && vec[0].substr(vec[0].rfind(':')+1) == id)
        {
            return *C;
        }
    }

    return 0;
}

bool
same_coords(const Region * r1, const Region * r2)
{
    return r1->start == r2->start && r1->end == r2->end;
}


// sort all the genes in the genome by their left end point
void
sort_genes(Genome & g)
{
    for(EVERY_CHROM(g, C))
    {
        sort((*C)->genes.begin(), (*C)->genes.end(), region_isleftof);

        // remove any duplicate genes
        vector<Region*>::iterator e = unique((*C)->genes.begin(), 
            (*C)->genes.end(), same_coords);
        (*C)->genes.erase(e, (*C)->genes.end());
    } 
}

