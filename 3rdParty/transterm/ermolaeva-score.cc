/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <cassert>
#include "transterm.h"
#include "ermolaeva-score.h"
#include "util.h"

Energy MM_OPEN = 0.0;  // this must be 0.0 (do not change)


// the params from the v1.0 paper (they have not yet been updated)


Energy GC  = -2.3;
Energy AU  = -0.9;
Energy GU  = 1.3;
Energy MM  = 3.5;
Energy GAP = 6.0;

Energy LOOP_PENALTY[20]={1000,1000,1000,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};



/*
// mfold-derived parameters
Energy GC = -2.1;
Energy AU = -1.3;
Energy GU = -0.8;
//Energy MM = 0.8;
Energy MM = 3.5;
Energy GAP = 3.9;
//Energy GAP = 6;

Energy LOOP_PENALTY[20] = {1000,1000,1000, 4.10, 4.90, 4.40, 4.70, 5.00, 5.10, 5.20, 5.30, 5.40,
                           5.50, 5.60, 5.70, 5.80, 5.80, 5.90, 5.90, 6.00};
*/

// format of str: f1,f2,f3,f4,f5
// f1 is the cost of a loop of length MIN_LOOP
// theere are fewer terms than needed to get up to MAX_LOOP, then
// the last term is repeated so "0,2" would give cost 0 to any loop
// of length MIN_LOOP, and cost 2 to any larger loop.
// extra terms are ignored
void
set_loop_pen(const string & str)
{
    for(int i = 0; i < MIN_LOOP;i++) LOOP_PENALTY[i] = 1000;

    vector<string> pen;
    split(str, ',', pen);
    for(int i = MIN_LOOP; i <= MAX_LOOP; i++)
    {
        int j = i - MIN_LOOP;
        LOOP_PENALTY[i] = atof(((unsigned(j) < pen.size())?pen[j]:pen.back()).c_str());
    }
}

// score the hairpin in the terminator t
Energy
hairpin_score(HPScoreFcn score, const Term & t)
{
    Energy e = loop_penalty(t.loop_len) + ((t.gap != 0) ? GAP : 0.0);

    for(int j = 0; j < t.stem_len; j++)
    {
        // adjust for the single gap if we pass it
        int left, right;
        left = (t.gap < 0 && j >= abs(t.gap)) ? 1 : 0;
        right = (t.gap > 0 && j >= abs(t.gap)) ? 1 : 0;

        e += score(*(t.left_stem_base() + j + left), 
                   *(t.right_stem_base() - j - right));
    }
    return e;
}


// calculate the tail score (at the given end of the Term)
Energy
tail_score(const Term & ter, Direction dir)
{
    int inc = 1;
    const char * start = ter.right_stem_base() + 1;
    char letter = 'T';

    if(dir == REVERSE)
    {
        inc = -1;
        start = ter.left_stem_base() - 1;
        letter = 'A';
    }
        
    Energy sum = 0.0, prev = 1.0;

    for(int i = 0; i < 15; i++)
    {
        prev *= ((*start == letter) ? 0.9 : 0.6);
        start += inc;
        sum += prev;
    }
    return -sum;
}
