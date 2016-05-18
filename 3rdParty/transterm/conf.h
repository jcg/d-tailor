/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#ifndef CONF_H
#define CONF_H
#include <map>
#include "distr.h"
#include "seq.h"

// Abstract base class for a confidence algorithm. Users must call prepare()
// first and then can call score() to get  confidence score for a terminator
class Confidence
{
public:

    virtual void prepare(const Genome &) = 0;
    virtual int score(const Term &, RegionType) const = 0;
    virtual ~Confidence() {}
};

// compute the confidence as described in the Ermolaeva et al 2000 paper.
// Note: the terms and genes must be sorted by their leftmost point before
// the call to prepare.
class ErmolaevaConfidence : public Confidence
{
public:
    ErmolaevaConfidence() : prepared(false) {}
    virtual ~ErmolaevaConfidence() {}
    void prepare(const Genome &);
    int score(const Term &, RegionType) const;

protected:
    double K;
    double t2t_L, h2t_L;
    double t2t_N, h2t_N;

    Distribution h2t_hp, h2t_tail;
    Distribution t2t_hp, t2t_tail;

    bool prepared;

    int score_one(const Term &, RegionType) const;
};


// compute the confidence using the background distribution (usually random)
// read in from a given file.
class RandomConfidence : public Confidence
{
public:
    RandomConfidence(const string &);
    virtual ~RandomConfidence() {}
    void prepare(const Genome &);
    int score(const Term &, RegionType) const;

    typedef vector<vector<unsigned long> > Histogram2d;

protected:
    Histogram2d & get_table(double);
    const Histogram2d & get_table(double) const;
    void read_exp_table(const string &);
    void fill_emp_table(RegionType, const ConstTermVec &);
    unsigned long histvalue(const Histogram2d & hist, const Term &, double) const;
    unsigned long & histvalue(Histogram2d & hist, const Term &, double);
    int hbin(double, double, int, double) const;
    int get_best_at(double) const;

    unsigned long _sample_size;
    //double _low_hp, _high_hp;
    //double _low_tail, _high_tail;
    map<int, double> _low_hp, _high_hp, _low_tail, _high_tail;
    int _nbins;
    map<int, Histogram2d> _exp_table;
    map<RegionType, Histogram2d> _emp_table;
    map<RegionType, unsigned long> _emp_len;
    map<RegionType, double> _emp_at;

    bool _prepared;
};

class RandomPValueConfidence : public RandomConfidence
{
public:
    RandomPValueConfidence(const string &);
    virtual ~RandomPValueConfidence() {}
    virtual int score(const Term &, RegionType) const ;
    
protected:
    void sum_exp_table();
};

int er_confidence(const EventResponder &, const Confidence &, const Term &);
Distribution signal_to_noise(Term::EnergyKind, const ConstTermVec &, const ConstTermVec &);

#endif
