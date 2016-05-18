/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <iostream>
#include <vector>
#include <sstream>
#include <list>
#include <queue>
#include <algorithm>
#include "transterm.h"
#include "seq.h"
#include "ermolaeva-score.h"
#include "util.h"

int UWINDOW_SIZE = 6;
int UWINDOW_REQUIRE = 3;
int MIN_STEM    = 4;
int MIN_LOOP    = 3;
Energy ENERGY_CUTOFF = -2; 
Energy TAIL_CUTOFF   = -2.5;

// MAX_STEM is only used for the version 1 search scheme --- 
// it should not be changed or used in new code.
int MAX_STEM = 22;

int MAX_LOOP = 13;
int MAX_HP = 2*MAX_STEM + MAX_LOOP + 2;

void 
set_max_len_loop(int len, int loop)
{
    MAX_LOOP = loop;
    MAX_HP = len;
    assert(MAX_HP < REALLY_MAX_HP);
    assert(MAX_LOOP < MAX_HP);
}


// return true if the hp represented by t is a candidate hp
bool
is_candidate_hp(HPScoreFcn score, const Term & t)
{
    return score(*t.left_stem_base(), *t.right_stem_base()) < AU &&
           score(*t.left_stem_top(),  *t.right_stem_top()) < 0 &&
           (score(*(t.left_stem_top() + 1), *(t.right_stem_top() - 1)) > 0 ||
                t.loop_len == MIN_LOOP || t.loop_len == MIN_LOOP + 1);
}


// Check if there is at least three consecutive nucleotides of appropriate type
// that letter among 7 characters after cp
inline
bool
check_tail(char letter, SeqPtr cp)
{
    return (cp[0] == letter && cp[1] == letter && cp[2] == letter) ||
           (cp[1] == letter && cp[2] == letter && cp[3] == letter) ||
           (cp[2] == letter && cp[3] == letter && cp[4] == letter) ||
           (cp[3] == letter && cp[4] == letter && cp[5] == letter) ||
           (cp[4] == letter && cp[5] == letter && cp[6] == letter);
}


// pair terminators that have the same coordinates going in the opposite
// directions. (the 'partner' member points to the other member of the pair)
void
pair_bidirect(vector<Term*> & in)
{
    Term * prev = 0;
    for(EVERY_TERM(in, T))
    {
        if(prev && 
           (*T)->left() == prev->left() && 
           (*T)->right() == prev->right() &&
           (*T)->dir() != prev->dir())
        {
            assert(!(*T)->partner && !prev->partner);
            (*T)->partner = prev;
            prev->partner = *T;
            prev = 0;
        }
        else 
        {
            prev = *T;
        }
    }
}

struct by_ascending_rightend
{
    bool
    operator()(const Term * t1, const Term * t2)
    {
        return t1->right() < t2->right();
    }
};

void
insert_by_rightend(list<Term*> & pq, Term * term)
{
    list<Term*> x;
    x.push_back(term);
    pq.merge(x, by_ascending_rightend());
}


// precond: terminators are sorted by left endpoint
void
find_same_overlapping(vector<Term*> & in)
{
    typedef list<Term*> TermPQ;
    TermPQ forward_queue, reverse_queue;

    for(EVERY_TERM(in, T))
    {
        TermPQ * que = ((*T)->dir() == FORWARD) ? &forward_queue : &reverse_queue;

        while(!que->empty() && que->front()->right() <= (*T)->left())
        {
            que->pop_front();
        }

        copy(que->begin(), que->end(), back_inserter((*T)->overlapping));
        for(TermPQ::iterator R = que->begin();
            R != que->end();
            ++R)
        {
            (*R)->overlapping.push_back(*T);
        }

        insert_by_rightend(*que, *T);
    }
}


// precond: terminators are sorted by left endpoint
void
find_opp_overlapping(vector<Term*> & in)
{
    typedef list<Term*> TermPQ;
    TermPQ forward_queue, reverse_queue;
    TermPQ *mydirQ, *oppdirQ;

    // for every terminator in the list sorted by leftend point
    for(EVERY_TERM(in, T)) 
    {
        // use direction to decide which queues are which
        if((*T)->dir() == FORWARD)
        {
            mydirQ = &forward_queue;
            oppdirQ = &reverse_queue;
        }
        else
        {
            mydirQ = &reverse_queue;
            oppdirQ = &forward_queue;
        }

        // save me in my direction's queue and the global queue
        insert_by_rightend(*mydirQ, *T);
        
        // remove all the guys that we've passed
        while(!oppdirQ->empty() && oppdirQ->front()->right() <= (*T)->left()) 
        {
            oppdirQ->pop_front();
        }

        // everyone still in the queue starts before me and ends after my
        // left end, and thus overlaps me
        copy(oppdirQ->begin(), oppdirQ->end(), back_inserter((*T)->opp_overlapping));

        // everyone that overlaps me, I overlap
        for(TermPQ::iterator R = oppdirQ->begin();
            R != oppdirQ->end();
            ++R)
        {
            (*R)->opp_overlapping.push_back(*T);
        }
    }

    // this is just a big assertion statement --- it checks that what we
    // tagged as overlapping above really does
    for(EVERY_TERM_CONST(in, T))
    {
        for(list<Term*>::const_iterator R = (*T)->opp_overlapping.begin();
            R != (*T)->opp_overlapping.end();
            ++R)
        {
            // all the opposite overlapping genes must be on the other strand
            assert((*T)->dir() != (*R)->dir());

            // the terminators must overlap
            assert(hp_overlap(**T, **R));
        }
    }
}




// handles the dynamic programmign table in a memory efficient way because we
// are looking for hairpins of bounded size (<MAX_HP), we need only keep a
// MAX_HP x MAX_HP table. as we move down the sequence we rotate the table,
// reusing the entries that we can.
class HPDPTable
{
public:
    HPDPTable(SeqPtr, Direction);

    void rotate(); 
    void update();

    // return the score between chars at i and j
    Energy score(int i, int j) const
    {
        return sc(*ptr_for(i), *ptr_for(j));
    }

    // return a pointer into the seq for the given index
    SeqPtr ptr_for(int i) const { return (cp - dir*i); }

    // the tracebacj 'arrows' in the DP algorithm
    enum Reason { LOOP, MISMATCH, MATCH, I_GAP, J_GAP, BAD };

    // access the tables
    float & stbl(int i, int j) { return s[idx(i)][idx(j)]; }
    Reason & rtbl(int i, int j) { return r[idx(i)][idx(j)]; }
    float stbl(int i, int j) const { return s[idx(i)][idx(j)]; }
    Reason rtbl(int i, int j) const { return r[idx(i)][idx(j)]; }

private:
    // translate indices to internal format
    int idx(int i) const { return (i+shift) % MAX_HP; }

    // access the table given an internal index (From idx())
    // same as stbl() and rtbl() except faster b/c we assume you
    // are giving us an internal index
    float & stbl_idx(int i, int j) { return s[i][j]; }
    Reason & rtbl_idx(int i, int j) { return r[i][j]; }
    float stbl_idx(int i, int j) const { return s[i][j]; }
    Reason rtbl_idx(int i, int j) const { return r[i][j]; }


    int rots_since_update;
    int shift;
    float s[REALLY_MAX_HP][REALLY_MAX_HP];  // changed to make big table
    Reason r[REALLY_MAX_HP][REALLY_MAX_HP];

    SeqPtr cp;
    Direction dir;
    HPScoreFcn sc; 
};


// make a new HPDPTable for the sequence that c points into
// we'll be scanning in direction d
HPDPTable::HPDPTable(SeqPtr c, Direction d) : 
    rots_since_update(MAX_HP), 
    shift(MAX_HP),
    cp(c),
    dir(d) 
{
    sc = (dir==FORWARD)?&forward_pair:&reverse_pair;
    for(int i = 0; i < MAX_HP; i++)
    {
        for(int j = 0; j < MAX_HP; j++)
        {
            s[i][j] = 1000;
            r[i][j] = BAD;
        }
    }
}


// Each time cp is incremented, we must rotate the table
void
HPDPTable::rotate()
{
    rots_since_update = min(rots_since_update+1, MAX_HP);
    shift = (shift==0)?MAX_HP-1:(shift-1);
    cp += dir;
}


void
HPDPTable::update()
{
    for(int i = rots_since_update - 1; i >= 0; i--)
    {
        int ii = idx(i);
        int ip1 = idx(i+1);

        if(i+MIN_LOOP-1 < MAX_HP)
        {
            int iml = idx(i+MIN_LOOP-1);
            stbl_idx(ii, iml) = loop_penalty(MIN_LOOP);
            rtbl_idx(ii, iml) = LOOP;
        }
        if(i+MIN_LOOP < MAX_HP)
        {
            int iml = idx(i+MIN_LOOP);
            stbl_idx(ii, iml) = loop_penalty(MIN_LOOP+1);
            rtbl_idx(ii, iml) = LOOP;
        }

        for(int j = i + MIN_LOOP+1; j < MAX_HP; j++)
        {
            int jj = idx(j);
            int jm1 = idx(j-1);

            Reason & rr = rtbl_idx(ii,jj);
            float & ss = stbl_idx(ii,jj);
            float y;

            rr = LOOP;
            ss = loop_penalty(j-i+1);

            Energy sij = score(i,j);
            bool mm_open = sij >= MM && rtbl_idx(ip1, jm1)!=MISMATCH;

            if((y = sij + stbl_idx(ip1,jm1) + (mm_open?MM_OPEN:0.0)) <= ss)
            {
                rr = (sij<MM)?MATCH:MISMATCH;
                ss = y;
            }

            if((y = GAP + stbl_idx(ip1, jj)) <= ss)
            {
                rr = I_GAP;
                ss = y;
            }

            if((y = GAP + stbl_idx(ii,jm1)) <= ss)
            {
                rr = J_GAP;
                ss = y;
            }
        }
    }
    rots_since_update = 0;
}

// find the lowest scoring hp anchored at the current position (cp)
// its other end will be returned (in index-space) in best_j and its
// energy will be returned. Use make_best_term() to follow the traceback
// arrows to get the actual terminator.
Energy
find_best_hp(const HPDPTable & dp, int & best_j)
{
    float hpe = 1000;
    best_j = 0;
    for(int j = MAX_HP-1; j >= 2*MIN_STEM + MIN_LOOP - 1; j--)
    {
        // first clause requires at least a weak pairing for first base pair in
        // the stem
        if(//dp.score(0,j) < MM  &&  // MM was AU in transterm 1.0
           (j == MAX_HP-1 || dp.stbl(0,j) < hpe))
        {
            hpe = dp.stbl(0, j);
            best_j = j;
        }
    }
    return hpe;
}


// follow the traceback arrows to get the best terminator with endpoints [cp,
// best_j].
Term
make_best_term(
    const Seq & seq, 
    Direction dir, 
    const HPDPTable & dp, 
    int best_j,
    Energy hpe)
{
    int i,j;
    i = 0;
    j = best_j;

    list<int> gaps;

    while(dp.rtbl(i, j) != HPDPTable::LOOP)
    {
        switch(dp.rtbl(i, j))
        {
            case HPDPTable::LOOP: break;
            case HPDPTable::I_GAP: i++; gaps.push_back(j); break;
            case HPDPTable::J_GAP: j--; gaps.push_back(-i); break;
            case HPDPTable::MATCH: // fall through
            case HPDPTable::MISMATCH: i++; j--; break;
            default: 
                cerr << "No value for: " << i << " " << j << endl;
                assert(false);
        }
    }

/*
    char tail_proximal = *ptr_for(i);
    char tail_distal = *ptr_for(j);

    if (dp.rtbl(0, best_j) == HPDPTable::I_GAP || 
        dp.rtbl(0, best_j) == HPDPTable::J_GAP || 
        (dir == FORWARD && tail_proximal != '
*/
    return Term(&seq, dir, dp.ptr_for(best_j), dp.ptr_for(j+1), 
                           dp.ptr_for(i-1), dp.ptr_for(0), gaps, hpe);
}


// look at last added term, if this one dominates, keep the one with the
// better tail score. Because of intervening overlaping (but non-dominating)
// terminators, dominating terms may appear in output
void
add_greedy_nodominating(const Term & t, vector<Term*> & terms)
{
    // if last term is inside t, keep the one with the best tail score.
    if(!terms.empty() && 
       ((t.dir() == FORWARD && terms.back()->left() >= t.left()) ||
       (t.dir() ==  REVERSE && terms.back()->right() <= t.right())))
    {
        if(t.tail_energy < terms.back()->tail_energy)
        {
            *terms.back() = t;
        }
    }
    else
    {
        terms.push_back(new Term(t)); 
    }
}


// keep all overlapping terminators and let the confidence sort things out
void
add_all_terminators(const Term & t, vector<Term*> & terms)
{
    // make sure that we extend the tail as far as possible
//    if((dir == FORWARD && first_pair(t) != "AT") ||
//       (dir == REVERSE && first_pair(t) != "TA"))
//    {
        terms.push_back(new Term(t));
//    }
}




// search all previous terms, if this one dominates any, keep the one with the
// best tail score. If nooverlaps == true, then we remove all /overlapping/
// (not just dominating) --- this may not output the /best/ set of
// non-overlapping hp, however
void
add_nodominating(
    const Term & t, 
    vector<Term*> & terms, 
    bool nooverlaps = false
    )
{
    // among the terminators in a domination chain involving t, find the one
    // with the best tail score (may be t)
    const Term * best = &t;
    for(vector<Term*>::reverse_iterator T = terms.rbegin();
        T != terms.rend();
        ++T)
    {
        if(!hp_overlap(t, **T)) break;

        if(nooverlaps || dominates(t, **T))
        {
            best = (t.tail_energy <= (*T)->tail_energy) ? &t : *T;
        }
    }

    // if the best one is t, then we have to remove the ones that it
    // dominates. Otherwise, the best one is already in the list and t is
    // redundant.
    if(best == &t)
    {
        // (note: reverse_iterator was segfaulting for no apparent reason
        // when we erase() --- thus, we count a forward iterator backward)
        for(vector<Term*>::iterator T = terms.end();
            T != terms.begin();)
        {
            T--;
            if(!hp_overlap(t, **T)) break;

            // invalidates T and those T we've already checked
            if(nooverlaps || dominates(t, **T)) terms.erase(T); 
        }

        terms.push_back(new Term(t));
    }
}


// return true if there's a run of >= 5 A within 3 of the left side and >= 5 T within 3 of the
// right side and the terminators stem is >= 10 bases
bool
has_bad_tails(const Term & t)
{
//    01234567
//    ...AAAAA
//    ..AAAAA.
//    .AAAAA..
//    AAAAA...

    int stem_len = min(t.left_stem_top() - t.left_stem_base(),
                       t.right_stem_base() - t.right_stem_top());

    const int REQUIRE_AT = 5;
    const int LEADING = stem_len / 2 - REQUIRE_AT;


    // if both stems are at least 10 bases
    if (stem_len > 12)
    {
        // check to see if there is a stretch of As on the left
        int A = 0;
        for(SeqPtr cp = t.left_stem_base(); 
            cp < t.left_stem_base() + LEADING + REQUIRE_AT &&
            cp <= t.left_stem_top(); 
            ++cp)
        {
            A = (*cp == 'A') ? A+1 : 0;
            if(A >= REQUIRE_AT) break;
        }

        if (A < REQUIRE_AT) return false;

        // check to see if there is a stretch of Ts on the right
        int T = 0;
        for(SeqPtr cp = t.right_stem_base();
            cp > t.right_stem_base() - LEADING - REQUIRE_AT &&
            cp >= t.right_stem_top();
            --cp)
        {
            T = (*cp == 'T') ? T+1 : 0;
            if(T >= REQUIRE_AT) break;
        }

        if (T < REQUIRE_AT) return false;
        return true;
    }
    
    return false;
}


// we add t (with scores (H,T)) if there are no terminators that completely dominate
// t, where a terminator A completely dominates B if A hairpin is a supersequence of
// Bs and H(B) > H(A) and T(B) > T(A). 
void
add_non_completely_dominating(
    const Term & t, 
    vector<Term*> & terms)
{
    list<Term*> remove_these;

    // we assume we're going to keep t
    bool add_t = true;

    //if(abs(seqindex(*t.seq, t.right()) - 253468) < 1000) cerr << "XX: " << t << " ";

    // walk backwards untill we no longer overlap t
    for(vector<Term*>::reverse_iterator T = terms.rbegin();
        T != terms.rend();
        ++T)
    {
        if(!hp_overlap(t, **T)) break;

        // if t dominates one of the prevous ones:
        if(dominates(t, **T))
        {
            // if there is already something there that t dominates and t looks like
            // it might be a bad extension, don't add it
            if(has_bad_tails(t)) return;
            //if(abs(seqindex(*t.seq, t.right()) - 253468) < 1000) cerr << "good tails ";

            // if new terminator t dominates this previous terminator and uniformly
            // has better hp and tail scores: then there is no reason to keep the old one
            if(t.tail_energy <= (*T)->tail_energy && t.hp_energy <= (*T)->hp_energy)
            {
                remove_these.push_back(*T);
            }

            // convesely, if the old one is better uniformly better, then there's no
            // need to keep this new one
            if(remove_these.empty() && 
               t.tail_energy > (*T)->tail_energy && t.hp_energy > (*T)->hp_energy)
            {
                add_t = false;
            }
        }
    }

    //if(abs(seqindex(*t.seq, t.right()) - 253468) < 1000) cerr << add_t << endl;

    // below depends on the ordering of the terms in the sequences. In
    // remove_these they are ordered in REVERSE order so that the first
    // terminator is the lastest one in terms.

    // if we have to remove some things that t completely domianates
    if(!remove_these.empty())
    {
        // (note: reverse_iterator was segfaulting for no apparent reason
        // when we erase() --- thus, we count a forward iterator backward)
        for(vector<Term*>::iterator T = terms.end();
            T != terms.begin();)
        {
            T--;
            if(!hp_overlap(t, **T)) break;

            if(*T == remove_these.front())
            {
                // invalidates T and those T we've already checked
                terms.erase(T);
                remove_these.pop_front();
            }
        }
    }

    if(add_t) terms.push_back(new Term(t));
}


// return the number of letter in between cp and last (inclusive)
int
count_letters(char letter, SeqPtr cp, SeqPtr last)
{
    SeqPtr end = max(cp, last);
    int count = 0;
    for(SeqPtr start = min(cp, last); start != end; start++)
    {
        if(*start == letter) count++;
    }
    if(*end == letter) count++;
    return count;
}


// find the best hairpin A overlapping the given terminator on the left for
// FORWARD terminators, or the right for REVERSE terminators
//                ----------.....---------  Term
// --------------------- Antiterm
//
// return true if anything found
bool
best_overlapping_hairpin(
    const Term & term, 
    Term & best_term, 
    int & best_overlap)
{
    Direction dir = term.dir();
    SeqPtr stem_end, stem_start;
    if(dir == FORWARD)
    {
        stem_end = term.left_stem_top()+1;
        stem_start = term.left_stem_base();
    }
    else
    {
        stem_end = term.right_stem_top()-1;
        stem_start = term.right_stem_base();
    }

    SeqPtr cp = stem_start - dir * MAX_HP + 2;

    int best_j = 0, j;
    Energy best_hpe = 1000, hpe;

    HPDPTable dp(cp, dir);

    for(; cp != stem_end; cp += dir)
    {
        // cp is a valid starting point for an antiterm
        if((dir == FORWARD && cp >= stem_start) ||
           (dir == REVERSE && cp <= stem_start))
        {
            dp.update();
            if((hpe = find_best_hp(dp, j)) < best_hpe)
            {
                best_j = j;
                best_hpe = hpe;
                best_overlap = abs(stem_start - cp) + 1;
                best_term = make_best_term(*term.seq, dir, dp, best_j, hpe);
            }
        }
        dp.rotate();
    }

    return best_j > 0;
}


// search for terminators along the 'dir' strand using an efficient dynamic
// programming algorithm. puts the terms into seq.terms. 
void
find_terms_dp(Seq & seq, Direction dir)
{
    Energy hpe;
    int best_j;
    SeqPtr cp, end_cp;
    int tailoffset = 0;
    char letter = 'T';
    
    if(dir == FORWARD)
    {
        cp = seq.left() + MAX_HP;
        end_cp = seq.right() - 15;
        tailoffset = 1;
        letter = 'T';
    }
    else
    {
        cp = seq.right() - MAX_HP;
        end_cp = seq.left() + 15;
        tailoffset = -7;
        letter = 'A';
    }

    int in_window = count_letters(letter, cp+1*dir, cp+UWINDOW_SIZE*dir);

    HPDPTable dp(cp, dir);

    // for every character that has a plausible tail
    for (; cp != end_cp; cp += dir)
    {
        if(in_window >= UWINDOW_REQUIRE && *cp != letter) //check_tail(letter, cp + tailoffset)) 
        {
            dp.update();

            if((hpe = find_best_hp(dp, best_j)) < ENERGY_CUTOFF)
            {
                Term t = make_best_term(seq, dir, dp, best_j, hpe);
                //t.hp_energy /= t.stem_len;  // uncomment to use per-base energies
                // filter "hairpins" with too small a stem
                if(t.stem_len >= MIN_STEM) 
                {
                    t.tail_energy = tail_score(t, dir); 
                    if(t.tail_energy < TAIL_CUTOFF)
                    {
                        //add_greedy_nodominating(t, seq.terms);
                        //add_all_terminators(t, seq.terms);
                        add_non_completely_dominating(t, seq.terms);
                    }
                }
            }
        }

        dp.rotate();

        // update the sliding window of size 6
        if(*(cp+1*dir) == letter) in_window--;
        if(*(cp+(UWINDOW_SIZE+1)*dir) == letter) in_window++;
    }
}


// find the terminators using a efficient dynamic programmign algorithm
void
find_terms_dp(Seq & seq)
{
    find_terms_dp(seq, FORWARD);
    find_terms_dp(seq, REVERSE);
    sort(seq.terms.begin(), seq.terms.end(), region_isleftof);
    pair_bidirect(seq.terms);
    find_opp_overlapping(seq.terms);
    find_same_overlapping(seq.terms);
}

//============================================================
// 2ndscore.cc code: used to output hp energies for every posn
//============================================================

// return a vector scores[] of length = len(seq) where scores[i] == the best
// hairpin score for a hairpin anchored at position i. (the positions are
// 1-based and entry 0 is not used.)
void
every_hairpin_energy(Seq & seq, Direction dir, vector<Term> & scores)
{
    // compute the start and end of the ranges
    SeqPtr cp, end_cp;
    if(dir == FORWARD)
    {
        cp = seq.left() + MAX_HP;
        end_cp = seq.right() - MAX_HP;
    }
    else
    {
        cp = seq.right() - MAX_HP;
        end_cp = seq.left() + MAX_HP;
    }

    // this variable stores the complete dynamic programming table
    HPDPTable dp(cp, dir);

    // we assume all energies are < MAX_ENERGY (= 10000) [XXX: this assignment
    // is needed only b/c we don't yet compute the scores for the very start
    // and ends of the sequences]
    scores.resize(seq.length+1);

    for(; cp != end_cp; cp += dir)
    {
        int best_j;  // not used, but must be passed to find_best_hp()

        dp.update();

        Energy hpe = find_best_hp(dp, best_j);
        Term t = make_best_term(seq, dir, dp, best_j, hpe);
        t.tail_energy = 0.0;
        scores[seqindex(seq, cp)] = t;
        
        dp.rotate();
    }
}


//============================================================
// Version 1.0 Search Scheme
//============================================================

// find the terminators going in a particular direction. We split the directions
// so that we can bail out of the loop early to save time
void
find_terms_ermolaeva(Seq & seq, Direction dir)
{
    // Ermolaeva, et al algorithm: brute force. for every position, stem
    // length, gap position, and loop length, check both strands to see if we
    // get a terminator with good score
    HPScoreFcn score;
    int tailoffset;
    char letter;
    int start_buf, end_buf;

    if(dir == FORWARD)
    {
        score = &forward_pair;
        tailoffset = 1;
        letter = 'T';
        start_buf = 2*MAX_STEM + MAX_LOOP + 15 + ((no_gaps)?0:1);
        end_buf = 15;
    }
    else
    {
        score = &reverse_pair;
        tailoffset = -7;
        letter = 'A';
        start_buf = 15;
        end_buf = 2*MAX_STEM + MAX_LOOP + 15 + ((no_gaps)?0:1);
    }

    Term last_term;
    bool have_last_term = false;
    for (SeqPtr cp = seq.left() + start_buf; cp < seq.right()-end_buf; cp++)
    {
        if(!check_tail(letter, cp + tailoffset)) continue;
        

        for(int stem_len = MAX_STEM; stem_len >= MIN_STEM; stem_len--)
        {
            // gap is a location of the gap in the stem (maximum one gap per
            // hairpin is allowed). If gap = 0 than there is no gaps in athe
            // stem.  If gap < 0 than the gap is in the right side of the stem
            // (so a nucleotide that don't have pair is in the left side. If
            // gap > 0 than the gap is located in the left side.  The absolute
            // value of the gap is a number of nucleotides located in the stem
            // before the gap (counting from the end opposite to loop)

            int min_gap = 0, max_gap = 0;
            if(!no_gaps)
            {
                min_gap = -stem_len+2;
                max_gap = stem_len-2;
            }

            for(int loop_len = MAX_LOOP; loop_len >= MIN_LOOP; loop_len--)
            {
                for(int gap = min_gap; gap <= max_gap; gap++)
                {
                    Term term(&seq, dir, cp, stem_len, loop_len, gap);

                    if(term.left_stem_base() < seq.left()) continue;
                    
                    if(!is_candidate_hp(score, term)) continue;

                    term.hp_energy = hairpin_score(score, term);

                    if(term.hp_energy < ENERGY_CUTOFF)
                    {
                        term.tail_energy = tail_score(term, dir);

                        // these rules for printing are a little strange ---
                        // but they duplicate what Ermolaeva et al. use for
                        // their paper 
                        if(have_last_term && !hp_overlap(term, last_term))
                        {
                            if(last_term.tail_energy < TAIL_CUTOFF) 
                            {
                                seq.terms.push_back(new Term(last_term));
                            }
                            last_term = term;
                            have_last_term = true;
                        }
                        else if(!have_last_term || 
                                term.hp_energy < last_term.hp_energy)
                        {
                            last_term = term;
                            have_last_term = true;
                        }
                    }
                }
            }
        }
    }

    if(last_term.tail_energy < TAIL_CUTOFF) 
    {
        seq.terms.push_back(new Term(last_term));
    }
}


// find the terminators. will place the term objects into the sequence's
// term vector. The list of terms will be softed by their right() endpoint
// and their partners will be identified.
void
find_terms_ermolaeva(Seq & seq)
{
    find_terms_ermolaeva(seq, REVERSE);
    find_terms_ermolaeva(seq, FORWARD);
    sort(seq.terms.begin(), seq.terms.end(), region_isleftof);
    pair_bidirect(seq.terms);
}
