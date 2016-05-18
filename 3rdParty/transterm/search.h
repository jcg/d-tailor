#ifndef SEARCH_H

#include <vector>
#include "seq.h"

using namespace std;

void every_hairpin_energy(Seq &, Direction, vector<Term> &);
int set_max_len_loop(int, int);

extern bool no_gaps;

#endif
