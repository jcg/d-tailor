/* This file is part of TransTerm v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. 
 *
 * Author: Carl Kingsford, (c) 2005-2006 
 */

// 2ndscore will read in a fasta sequence and assign two secondary structure
// scores to every position: the score of the best hairpin to the left the
// score of the best hairpin on the complement-strand, anchored at the right

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <getopt.h>
#include <cassert>

#include "seq.h"
#include "util.h"
#include "ermolaeva-score.h"
#include "search.h"

bool no_gaps = false; // needed for legacy code
bool print_seq = true; // print the hairpin sequence
bool print_pos_strand = true;
bool print_neg_strand = true;

void
usage()
{
    cerr << "Usage: 2ndscore in.fasta" << endl;
    exit(3);
}

// print out the values for the essential options
void
print_options(ostream & out)
{
    out << endl 
        << "--gc=" << GC << " "
        << "--au=" << AU << " "
        << "--gu=" << GU << " "
        << "--mm=" << MM << " "
        << "--gap=" << GAP << endl;
    out << "--max-len=" << MAX_HP << " "
        << "--max-loop=" << MAX_LOOP << " "
        << "--min-loop=" << MIN_LOOP << endl;

    out << "--loop-penalty=";
    for(int i = MIN_LOOP; i <= MAX_LOOP; i++)
    {
        if(i != MIN_LOOP) out << ",";
        out << LOOP_PENALTY[i];
    }

    out << endl << endl;
}

int
process_options(int argc, char * argv[])
{
    const char * OPTIONS = "hS";

    enum {GC_OPT, AU_OPT, GU_OPT, MM_OPT, GAP_OPT, MINSTEM_OPT, 
          MINLOOP_OPT, NOPOS_OPT, NONEG_OPT, MAXLEN_OPT, MAXLOOP_OPT,
          LOOPPEN_OPT };

    static struct option long_options[] = {
        {"gc", 1, 0, GC_OPT},
        {"au", 1, 0, AU_OPT},
        {"gu", 1, 0, GU_OPT},
        {"mm", 1, 0, MM_OPT},
        {"gap", 1, 0, GAP_OPT},
        {"loop-penalty", 1, 0, LOOPPEN_OPT}, 

        {"min-loop", 1, 0, MINLOOP_OPT},
        {"max-len", 1, 0, MAXLEN_OPT},
        {"max-loop", 1, 0, MAXLOOP_OPT},

        {"no-fwd", 0, 0, NOPOS_OPT},
        {"no-rvs", 0, 0, NONEG_OPT},

        {"help", 0, 0, 'h'},
        {0,0,0,0}
    };

    int len = MAX_HP;
    int loop = MAX_LOOP;

//    opterr = 0; // don't print msg --- we'll do it
    int a;
    while((a=getopt_long(argc, argv, OPTIONS, long_options, 0)) != -1)
    {
        switch(a)
        {
            case 'h': usage(); break;

            // energy function options
            case GC_OPT: GC = atof(optarg); break;
            case AU_OPT: AU = atof(optarg); break;
            case GU_OPT: GU = atof(optarg); break;
            case MM_OPT: MM = atof(optarg); break;
            case GAP_OPT: GAP = atof(optarg); break;
            case LOOPPEN_OPT: set_loop_pen(optarg); break;

            // filtering options
            case MINSTEM_OPT: MIN_STEM = atoi(optarg); break;
            case MINLOOP_OPT: MIN_LOOP = atoi(optarg); break;
            case MAXLEN_OPT: len = atoi(optarg); break;
            case MAXLOOP_OPT: loop = atoi(optarg); break;

            case 'S': print_seq = false; break;
            case NOPOS_OPT: print_pos_strand = false; break;
            case NONEG_OPT: print_neg_strand = false; break;

            default:
                cerr << "Error: unknown option. " << endl;
                exit(3);
        }
    }

    if(len > REALLY_MAX_HP)
    {
        cerr << "Error: must search for hairpins with total length smaller than " 
             << REALLY_MAX_HP << endl;
        cerr << "(recompile after changing REALLY_MAX_HP to increase)" << endl;
        exit(3);
    }

    if (loop >= len)
    {
        cerr << "Error: max-loop must be less than max-len" << endl;
        exit(3);
    }

    // set the global constants MAX_STEM, MAX_LOOP, MAX_HP
    set_max_len_loop(len, loop);

    if(MAX_STEM < MIN_STEM || MIN_STEM < 1)
    {
        cerr << "Error: min-stem must be <= max-stem and > 0" << endl;
        exit(3);
    }

    if(MAX_LOOP < MIN_LOOP)
    {
        cerr << "Error: max-loop must be >= min-loop" << endl;
        exit(3);
    }

    if(optind >= argc) usage();

    return optind;
}


// run every_hairpin_energy() above going both forward and backward
void
every_hairpin_energy(
    Seq & seq, 
    vector<Term> & fwd_strand, 
    vector<Term> & rvs_strand)
{
    every_hairpin_energy(seq, FORWARD, fwd_strand);
    every_hairpin_energy(seq, REVERSE, rvs_strand);
}


// output the list of scores to the stream given by out
void
print_hp_energies(
    ostream & out, 
    const Seq & seq, 
    vector<Term> & scores, 
    Direction dir,
    int padding)
{
    // format is header line starting with > followed by a description and the
    // last word on the header line will be FORWARD or REVERSE indicating to
    // which strand the scores apply
    out << ">" << seq.name 
        << " " << seq.desc 
        << (dir==FORWARD?" FORWARD":" REVERSE") << endl;

    // XXX: remember that scores.size() = 1 mroe than the actual size b/c the
    // zero entry is not used
    for(unsigned i = MAX_HP + 1; i < scores.size() - MAX_HP; i++)
    {
        // i is one based:
        assert(seq.dna[i-1] != PADDING_CHAR);

        //Energy hpe = (scores[i].stem_len == 0)?10.0:scores[i].hp_energy;
        int s = seqindex(seq, scores[i].start) - padding;
        int e = seqindex(seq, scores[i].end) - padding;

        if(s > 0 && e < (int)scores.size() - padding)
        {
            // print the energy (or None if no stem was found)
            if(scores[i].stem_len == 0)
            {
                out << setw(4) << "None" << " ";
            }
            else
            {
                out << setw(4) << setprecision(4) << scores[i].hp_energy << " ";
            }

            // print the coordinates
            out << setw(7) << s << " .. " << setw(7) << e << " "; 

            // print the hairpin
            if(print_seq) print_term_seq(out, scores[i]);
            out << endl;
        }
    }
}


int
main(int argc, char * argv[])
{
    cerr << "2ndscore (" << __DATE__ << ")" << endl;
    int first_file_index = process_options(argc, argv);

    print_options(cerr);

    // holds all the sequences read
    Genome dna;

    // read the all given fasta file into the dna variable
    for(int i = first_file_index; i < argc; i++)
    {
        string seq_filename = argv[first_file_index];
        ifstream seq_file(seq_filename.c_str());
        if(!seq_file)
        {
            cerr << "Error: couldn't read file: " << seq_filename << endl;
            exit(3);
        }
        read_seqs(seq_file, dna);
    }
    
    // for every sequence, output the scores for both the positive and
    // negative strands
    int i = 0;
    for(EVERY_CHROM(dna, C))
    {
        vector<Term> terms;
        cerr << ++i << ". Seq: " << (*C)->name << " (length " << (*C)->length << ")";
        pad_seq(**C, MAX_HP);

        if (print_pos_strand)
        {
            every_hairpin_energy(**C, FORWARD, terms);
            print_hp_energies(cout, **C, terms, FORWARD, MAX_HP);
        }

        if(print_neg_strand)
        {
            every_hairpin_energy(**C, REVERSE, terms);
            print_hp_energies(cout, **C, terms, REVERSE, MAX_HP);
        }
        cerr << endl;
    }
}
