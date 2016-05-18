/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef ANALYSIS_H
#define ANALYSIS_H

void plot_tthits_vs_terms(ostream &, Confidence &, Genome &);
void t2t_hitanal(ostream &, const Genome &, Confidence &, int, bool );

unsigned count_starts_in_genes(const Seq &, Direction dir);

#endif
