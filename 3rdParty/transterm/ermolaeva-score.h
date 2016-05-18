/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef ERMOLAEVA_SCORE_H
#define ERMOLAEVA_SCORE_H
#include "seq.h"

typedef Energy (*HPScoreFcn)(char, char);

Energy tail_score(const Term &, Direction);
Energy hairpin_score(HPScoreFcn, const Term & );
Energy reverse_pair(char, char);
Energy forward_pair(char, char);
Energy loop_penalty(int);
void set_loop_pen(const string &);



//==== simple scores ===
/*
const Energy GC            = -3;
const Energy AU            = -2;
const Energy GU            = -1;
const Energy MM            = 4;
*/

//=== SVM model scores ===

/*
const Energy GC            = -1.99;
const Energy AU            = -0.99;
const Energy GU            = 1;
const Energy MM            = 4;
const float LOOP_PENALTY[20]={1000,1000,1000,
    -0.75,0.25,1.24,
    1.25, 2.50, 3.75, 5.00, 6.25, 7.50, 8.75, 
    10.00, 11.25, 12.50, 13.75, 15.00, 16.25, 17.5};
*/

//=== SVM model scores ===
/*
const Energy GC            = -1;
const Energy AU            = -0.67;
const Energy GU            = 0.46;
const Energy MM            = 1.83;
const float LOOP_PENALTY[20]={1000,1000,1000,
    0.18, 0.0, 0.21,
    0.33,  0.66,  0.99,  1.32,  1.65,  1.98,  2.31,  
    2.64,  2.97,  3.30,  3.63,  3.96,  4.29 };
*/
//=== Ermolaeva et al Scores ===

extern Energy ENERGY_CUTOFF;
extern Energy TAIL_CUTOFF;

extern Energy LOOP_PENALTY[20];

// stem energies
extern Energy GC;
extern Energy AU;
extern Energy GU;
extern Energy MM;
extern Energy GAP;
extern Energy MM_OPEN;

// stem size limits
//const unsigned MAX_STEM    = 22;
//const unsigned MAX_LOOP    = 13;
//const int MAX_STEM    = 22;
//const int MAX_STEM    = 200;
//const int MAX_LOOP    = 13;
//const int MAX_HP      = 2*MAX_STEM + MAX_LOOP + 2;
extern int MIN_STEM;
extern int MIN_LOOP;
extern int MAX_STEM;
extern int MAX_LOOP;
extern int MAX_HP;
const int REALLY_MAX_HP = 1000;  // size of the allocated table

// required # of U's in the window next to the tail
extern int UWINDOW_SIZE;
extern int UWINDOW_REQUIRE;

// the energy of a hairpin pair on the forward strand
inline
Energy
forward_pair(char ch1, char ch2)
{
    if((ch1 == 'G' && ch2 == 'C') || (ch1 == 'C' && ch2 == 'G')) return GC;
    if((ch1 == 'T' && ch2 == 'A') || (ch1 == 'A' && ch2 == 'T')) return AU;
    if((ch1 == 'T' && ch2 == 'G') || (ch1 == 'G' && ch2 == 'T')) return GU;
    if(ch1 == PADDING_CHAR || ch2 == PADDING_CHAR) return 1000.0;
    return MM;
}


// the energy of a hairpin pair on the reverse strand
inline
Energy
reverse_pair(char ch1, char ch2)
{
    if((ch1 == 'G' && ch2 == 'C') || (ch1 == 'C' && ch2 == 'G')) return GC;
    if((ch1 == 'T' && ch2 == 'A') || (ch1 == 'A' && ch2 == 'T')) return AU;
    if((ch1 == 'A' && ch2 == 'C') || (ch1 == 'C' && ch2 == 'A')) return GU;
    if(ch1 == PADDING_CHAR || ch2 == PADDING_CHAR) return 1000.0;
    return MM;
}

// calculate the loop penalty
inline
Energy
loop_penalty(int len)
{
    if(len > MAX_LOOP) return 1000.0;
//    if(len >= sizeof(LOOP_PENALTY)/sizeof(float)) return 1000.0;
    return LOOP_PENALTY[len];
}


#endif
