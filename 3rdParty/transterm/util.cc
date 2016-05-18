/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <sstream>
#include "util.h"


// print a status percent
void
print_status(ostream & out, unsigned long cur, unsigned long max)
{
    ostringstream oss;
    oss << (int)(100 * ((float)cur) / max) << "%"; 
    out << oss.str();
    for(unsigned i = 0; i < oss.str().length(); i++) out << "\b";
}

// split a string into fields separated by a single character
void
split(const string & s, char sep, vector<string> & out)
{
    unsigned i = 0;
    out.clear();
    out.push_back("");
    for(unsigned j = 0; j < s.length(); j++)
    {
        if(s[j] == sep)
        {
            i++;
            out.push_back("");
        }
        else
        {
            out[i] += s[j];
        }
    }
}  

// remove whitespace at the front of s
string
trim_front(const string & s)
{
    unsigned i;
    for(i = 0; i<s.length(); i++)
    {
        if(!isspace(s[i])) break;
    }
    return s.substr(i);
}

 
// center a string 
string
center(const string & s, int fieldsize)
{
    string news = "";
    int ex = fieldsize - s.length();
    if(ex <= 0) return s;
    int pad = ex/2;
    for(int i = 0; i < pad; i++) news += ' ';
    news += s;
    for(int i = 0; i < pad; i++) news += ' ';
    if(ex % 2 != 0) news += ' ';
    return news;
}
