/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. */

#include <vector>
#include <iostream>
#include "seq.h"
#include "conf.h"
#include "util.h"


//==========================================================================
// Version 1.0's confidence & output scheme (deprecated)
// the code below is /not/ efficient /or/ elegant. 
//==========================================================================

class ErmolaevaConfVer1 : public ErmolaevaConfidence
{
public:
    ErmolaevaConfVer1() : ErmolaevaConfidence() {}
    virtual ~ErmolaevaConfVer1() {}
    void prepare(const Genome &);
};

// count the # of at and gc in the string. Counts are added to at and gc
void
count_atgc(const string & s, unsigned long * at, unsigned long * gc)
{
    for(unsigned i = 0; i < s.length(); i++)
    {
        if(s[i] == 'A' || s[i] == 'T') (*at)++;
        else (*gc)++;
    }
}


// the fraction of ingene bases that are A or T
// does not count the first or last 100 bases of the gene
double
at_percent_ingenes(const Genome & seqs, int buf)
{
    string gene;
    unsigned long at = 0, gc = 0;

    for(EVERY_CHROM_CONST(seqs, C))
    {
        for(EVERY_REGION_CONST((*C)->genes, G))
        {
            gene = subseq((*G)->left()+buf, (*G)->right()-buf);
            count_atgc(gene, &at, &gc);
        }
    }
    return ((float)at)/(at+gc);
}


// assumes that the genes are sorted
double
at_percent_notgenes(const Genome & seqs, int buf)
{
    string gene;
    unsigned long at = 0, gc = 0;

    for(EVERY_CHROM_CONST(seqs, C)) 
    {
        SeqPtr max_seen_pos = (*C)->left();

        for(EVERY_REGION_CONST((*C)->genes, G)) 
        {
            SeqPtr nextgene = (*G)->left();
            
            string intergene = subseq(
                max(max_seen_pos-buf, (*C)->left()), 
                nextgene+buf);
            count_atgc(intergene, &at, &gc);

            max_seen_pos = (*G)->right(); 
        }
    }
    return ((float)at)/(at + gc);
}


// return a list of tail-to-tail regions
void
tail_to_tail_regions(const vector<Region*> & genes, vector<Region*> & reg)
{
    Region * prev = 0;

    reg.clear();

    for(EVERY_REGION_CONST(genes, G))
    {
        // Look for conseqative genes of the form ---> <---.
        // If they overlap: 
        //          ---->
        //            <-----
        // then the 'tail to tail' region is the overlap region.
        if(G != genes.begin() && prev->dir() != (*G)->dir() && 
           abs(prev->end - (*G)->end) < abs(prev->start - (*G)->start))
        {
            reg.push_back(
                new Region(prev->name + "><" + (*G)->name,
                       (*G)->seq,
                       min(prev->end, (*G)->end),
                       max(prev->end, (*G)->end)));
        }

        prev = *G;
    }
}


// return a list of head-to-tail regions
void
head_to_tail_regions(const vector<Region*> & genes, vector<Region*> & reg)
{
    Region * prev = 0;
    reg.clear();

    for(EVERY_REGION_CONST(genes, G))
    {
        // look for two genes in the same dir: -->  -->  or <-- <--
        if(G != genes.begin() && prev->dir() == (*G)->dir())
        {
            string name;
            SeqPtr begin, end;

            // If they overlap, then the region is the overlap region,
            // going in the same direction as the genes direction
            if(prev->dir() == FORWARD)
            {
                name = prev->name + "->" + (*G)->name;
                begin = min(prev->end, (*G)->start);
                end = max(prev->end, (*G)->start);
                if(begin==end) end++;
            }
            else
            {
                name = prev->name + "<-" + (*G)->name;
                begin = max((*G)->end, prev->start);
                end = min((*G)->end, prev->start);
                if(begin == end) begin++;
            }
            reg.push_back(new Region(name, (*G)->seq, begin, end));
        }

        prev = *G;
    }
}


// output new terms that are in the given regions. The terms must be sorted
// the terms will be copies of the old terms
unsigned long
copy_terms_in_regions(
    const ConstTermVec & terms, 
    const vector<Region*> & reg, 
    ConstTermVec & out, 
    int cut, 
    bool require_codirect = false)
{
    unsigned long len = 0;
    unsigned i = 0;
    for(EVERY_REGION_CONST(reg, R))
    {
        if(abs((*R)->start - (*R)->end) > cut*2)
        {
            SeqPtr start, end, reg_right, reg_left;

            // cut the sequences by the given amount
            int dir = ((*R)->dir() == FORWARD) ? 1 : -1;
            start = (*R)->start + cut * dir;
            end = (*R)->end - cut * dir;
            len += abs(start - end);
            reg_right = max(start, end);
            reg_left = min(start, end);

            for(; i < terms.size() && 
                  terms[i]->right_stem_base() <= reg_right; i++)
            {
                if(terms[i]->left_stem_base() >= reg_left && 
                   (!require_codirect || terms[i]->dir() == (*R)->dir()))
                {
                    Term * t = new Term(*terms[i]);
                    t->name = (*R)->name;
                    out.push_back(t);
                }
            }
        }
    }
    return len;
}

const int GENE_CUT = 100;
const int H2T_CUT = -50;
const int T2T_CUT = -50;

// following two functions used /ONLY/ to wedge old code into the new scheme to
// support the version 1.0 way of doing this. Version 1.0 is only supported for
// comparison, debuggin, and completness -- do not use these two functions
void to_const_term_vec(vector<Term*> vec, ConstTermVec & out)
{
    out.clear();
    copy(vec.begin(), vec.end(), back_inserter(out));
}

void force_remove_const(const ConstTermVec & in, vector<Term*> vec)
{
    for(EVERY_CTERM_CONST(in, T))
    {
        vec.push_back(const_cast<Term*>(*T));
    }
}

// compute statistics on the genome object to prepare for assessing the
// confidence of a terminator with score(). This must be called before score
// and both genes, and terminators must be sorted by their leftmost point
void
ErmolaevaConfVer1::prepare(const Genome & seqs)
{
    ConstTermVec gene_terms, h2t_terms, t2t_terms;
    unsigned long gene_len = 0, h2t_len = 0, t2t_len = 0;

    // get the terms in tail-to-tail and tail-to-head and gene regions
    for(EVERY_CHROM_CONST(seqs, C))
    {
        // hack to conver vector<term*> to vector<const Term *>
        ConstTermVec terms;
        to_const_term_vec((*C)->terms, terms);

        // "true" means only co-directional terms
        gene_len += copy_terms_in_regions(terms, (*C)->genes, 
            gene_terms, GENE_CUT, true); 

        vector<Region*> reg;

        head_to_tail_regions((*C)->genes, reg);
        h2t_len += copy_terms_in_regions(terms, reg, h2t_terms, H2T_CUT, true);
        tail_to_tail_regions((*C)->genes, reg);
        t2t_len += copy_terms_in_regions(terms, reg, t2t_terms, T2T_CUT);
    }

    // can't compute confidence if we have no gene_terms
    if(gene_terms.empty())
    {
        prepared = false;
        cout << "warning: no examples in genes; can't compute conf." << endl;
        return;
    }

    // compute K --- correction for AT content
    double at_in, at_not;
    at_in = at_percent_ingenes(seqs, 100);
    at_not = at_percent_notgenes(seqs, 100);

    K = (840*at_not*at_not - 1215.65*at_not + 448.9593) /
        (840*at_in*at_in   - 1215.65*at_in  + 448.9593);

    cout << "Genes: " 
         << at_in << " %AT, " 
         << gene_len << " nt, "
         << gene_terms.size() << " terms." << endl;

    cout << "Intergenic: "
         << at_not << " %AT, "
         << "H2T: " << h2t_len << " nt, " << h2t_terms.size() << " terms; "
         << "T2T: " << t2t_len << " nt, " << t2t_terms.size() << " terms. " 
         << endl;

    
    t2t_L = double(t2t_len) / gene_len;
    h2t_L = double(h2t_len) / gene_len;

    t2t_hp = signal_to_noise(Term::HAIRPIN, t2t_terms, gene_terms);
    t2t_tail = signal_to_noise(Term::TAIL, t2t_terms, gene_terms);

    h2t_hp = signal_to_noise(Term::HAIRPIN, h2t_terms, gene_terms);
    h2t_tail = signal_to_noise(Term::TAIL, h2t_terms, gene_terms);

    t2t_N = 2.0 * gene_terms.size() / t2t_terms.size();
    h2t_N = double(gene_terms.size()) / h2t_terms.size();

    prepared = true;
}


// patch up the scores for terms in out to account for possible bidirected
// terminators. (this is a hack to duplicate version 1.0's scheme)
void
pair_bidirect(
    ConstTermVec & in,
    ConstTermVec & out)
{
    const Term * prev = 0;

    for(EVERY_CTERM_CONST(in, T))
    {
        if(prev && 
           (*T)->left() == prev->left() && 
           (*T)->right() == prev->right())
        {
            Term * t = new Term(*prev);
            t->conf = int((1.0 - (1.0 - 
                (*T)->conf/100.0)*(1.0 - t->conf/100.0))*100.0 + 0.5);
            t->dir() = BIDIR;
            out.push_back(t);
            prev = 0;
        }
        else if(prev) 
        {
            Term * t = new Term(*prev);
            out.push_back(t);
            prev = 0;
        }
        else
        {
            prev = *T;
        }
    }
}



// called by confidence_ermolaeva() to copy and add the confidence
// this is used only to duplicate the version 1.0 scheme
void
add_confidence(
    ConstTermVec & out,
    const ConstTermVec & in,
    RegionType where,
    Confidence & conf)
{
    for(EVERY_CTERM_CONST(in, T))
    {
        Term * t = new Term(**T);
        t->conf = conf.score(**T, where);
        out.push_back(t);
    }
}


// for backwards compatibility, this function will compute something similar to
// TransTerm version 1.0
void
confidence_ermolaeva(Genome & seqs, ConstTermVec & out)
{
    ErmolaevaConfidence conf;
    conf.prepare(seqs);

    ConstTermVec gene_terms, h2t_terms, t2t_terms;
    vector<Region*> reg;
    for(EVERY_CHROM_CONST(seqs, C))
    {
        ConstTermVec terms;
        to_const_term_vec((*C)->terms, terms);

        head_to_tail_regions((*C)->genes, reg);
        copy_terms_in_regions(terms, reg, h2t_terms, H2T_CUT, true);
        
        tail_to_tail_regions((*C)->genes, reg);
        copy_terms_in_regions(terms, reg, t2t_terms, T2T_CUT);
    }

    ConstTermVec tmp;
    add_confidence(tmp, t2t_terms, TAIL2TAIL, conf);
    pair_bidirect(tmp, out);
    add_confidence(out, h2t_terms, HEAD2TAIL, conf);
}

// print the terminators
void
print_terms(ostream & out, const vector<Term*> & vec)
{
    for(EVERY_TERM_CONST(vec, T))
    {
        Term & ter = **T;

        out << setw(15) << ((ter.name=="")?"n/a":ter.name) << " "
            << setw(3) << ter.conf << " "
            << setw(5) << ter.hp_energy << " " 
            << setw(8) << ter.tail_energy << " " 
            << setw(7) << seqindex(*ter.seq, ter.right_stem_base()) << " "
            << setw(2) << ter.dir() << " ";
 
        out << subseq(ter.left_stem_base()-15, ter.left_stem_base()-1) << " "
            << center(subseq(ter.left_stem_base(), ter.left_stem_top()) + " " +
               subseq(ter.left_stem_top()+1,ter.right_stem_top()-1) + " " +
               subseq(ter.right_stem_top(), ter.right_stem_base()), 45) << " "
            << subseq(ter.right_stem_base()+1, ter.right_stem_base()+15) << " ";

        out << ter.gap << " " << ter.seq->name << endl;
    }
}
