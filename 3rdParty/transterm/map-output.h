/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef MAP_OUTPUT_H
#define MAP_OUTPUT_H
#include <iostream>
#include "seq.h"
#include "conf.h"

void output_map(ostream &, const Genome &, Confidence &, int=90, bool=true, bool=true);
void output_best_term(ostream &, const Confidence &, const Genome &);

#endif
