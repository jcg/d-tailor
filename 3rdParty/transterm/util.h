/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace std;

void print_status(ostream &, unsigned long, unsigned long);
void split(const string &, char, vector<string> &);
string center(const string &, int);
string trim_front(const string &);

#endif
