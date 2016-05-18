/* This file is part of TransTermHP v2.0 BETA and is covered by the GNU GPL
 * License version 2.0. See file LICENSE.txt for more details. 
 *
 * Author: Carl Kingsford, (c) 2005-2006 
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <getopt.h>
#include <time.h>

#include "transterm.h"
#include "seq.h"
#include "util.h"
#include "map-output.h"
#include "gene-reader.h"
#include "analysis.h"
#include "ermolaeva-score.h"
#include "search.h"

void find_terms_dp(Seq & );
void output_anti_terms(ostream &, const Genome &, Confidence &, int);

// OPTIONS
bool no_gaps = false;
int conf_cutoff = 76;
bool print_seq = true;
bool show_all_t2t_roc = false;
bool show_gaps = false;
bool only_good_context = true;
string t2tperf_file = "";
string tthitfile = "";
string bag_file = "";
string antifile = "";
Confidence * conf;
int gene_start_cut = 0;
int gene_end_cut = 25;


// output the usage information & exit
void
usage(void)
{
    cerr << "usage: transterm [options] *.fasta *.coords" << endl;
    cerr << "See the USAGE.txt file for available options" << endl;
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
//    out << "--max-stem=" << MAX_STEM << " "
    out << "--max-len=" << MAX_HP << " "
        << "--min-stem=" << MIN_STEM << " "
        << "--max-loop=" << MAX_LOOP << " "
        << "--min-loop=" << MIN_LOOP << endl;
    out << "--uwin-length=" << UWINDOW_SIZE << " "
        << "--uwin-require=" << UWINDOW_REQUIRE << " "
        << "--max-hp-score=" << ENERGY_CUTOFF << " "
        << "--max-tail-score=" << TAIL_CUTOFF << endl; 

    out << "--loop-penalty=";
    for(int i = MIN_LOOP; i <= MAX_LOOP; i++)
    {
        if(i != MIN_LOOP) out << ",";
        out << LOOP_PENALTY[i];
    }

    out << endl << "--start-cut=" << gene_start_cut << " --end-cut=" << gene_end_cut;
    out << endl << endl;
}


// parse the command line
int
process_options(int argc, char * argv[])
{
    const char * OPTIONS = "hc:Sr:p:";

    enum {GC_OPT, AU_OPT, GU_OPT, MM_OPT, GAP_OPT, MINSTEM_OPT, 
          MINLOOP_OPT, LOOPPEN_OPT, UWINLEN_OPT, UWINREQ_OPT, 
          OVERLAP_OPT, MAXHP_OPT, MAXTAIL_OPT, V1CONF_OPT,
          RANDCONF_OPT, BAGOUTPUT_OPT, T2TPERF_OPT, ANTITERMS_OPT,
          SHOW_ALL_T2T_ROC_OPT, SHOW_GAPS_OPT, STARTCUT_OPT, ENDCUT_OPT, 
          MAXLEN_OPT, MAXLOOP_OPT, ALLCONTEXT_OPT, OLDRANDCONF_OPT};

    static struct option long_options[] = {
        {"gc", 1, 0, GC_OPT},
        {"au", 1, 0, AU_OPT},
        {"gu", 1, 0, GU_OPT},
        {"mm", 1, 0, MM_OPT},
        {"gap", 1, 0, GAP_OPT},
        {"min-stem", 1, 0, MINSTEM_OPT},
        {"min-loop", 1, 0, MINLOOP_OPT},
        {"max-loop", 1, 0, MAXLOOP_OPT},
        {"loop-penalty", 1, 0, LOOPPEN_OPT}, 
        {"uwin-size", 1, 0, UWINLEN_OPT},
        {"uwin-require", 1, 0, UWINREQ_OPT},
        {"overlap", 1, 0, OVERLAP_OPT}, // NYI
        {"max-hp-score", 1, 0, MAXHP_OPT},
        {"max-tail-score", 1, 0, MAXTAIL_OPT},
        {"max-len", 1, 0, MAXLEN_OPT},

        {"start-cut", 1, 0, STARTCUT_OPT},
        {"end-cut", 1, 0, ENDCUT_OPT},
        
        {"v1-conf", 0, 0, V1CONF_OPT},
        {"old-rand-conf", 1, 0, OLDRANDCONF_OPT},
        {"rand-conf", 1, 0, 'r'},
        {"pval-conf", 1, 0, 'p'},

        {"t2t-perf", 1, 0, T2TPERF_OPT},
        {"full-t2t-roc", 0, 0, SHOW_ALL_T2T_ROC_OPT},
        {"show-gaps", 0, 0, SHOW_GAPS_OPT},
        {"antiterms", 1, 0, ANTITERMS_OPT},
        {"all-context", 0, 0, ALLCONTEXT_OPT},

        {"help", 0, 0, 'h'},
        {"min-conf", 1, 0, 'c'},
        {"bag-output", 1, 0, BAGOUTPUT_OPT},
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
            case UWINLEN_OPT: UWINDOW_SIZE = atoi(optarg); break;
            case UWINREQ_OPT: UWINDOW_REQUIRE = atoi(optarg); break;
            case MAXHP_OPT: ENERGY_CUTOFF = atof(optarg); break;
            case MAXTAIL_OPT: TAIL_CUTOFF = atof(optarg); break;
            case MAXLEN_OPT: len = atoi(optarg); break;
            case MAXLOOP_OPT: loop = atoi(optarg); break;

            case OVERLAP_OPT: 
                cerr << "Error: not yet implemented" << endl;
                exit(3);

            case STARTCUT_OPT: gene_start_cut = atoi(optarg); break;
            case ENDCUT_OPT: gene_end_cut = atoi(optarg); break;

            // confidence options
            case V1CONF_OPT: 
                conf = new ErmolaevaConfidence(); 
                cerr << "WARNING: USING VERSION 1.0 CONFIDENCE." << endl
                     << "         (use -p expterm.dat to use v2.0 confidence)" << endl
                     << "Version 1 confidence scheme exists only for debugging and"
                     << "historical reasons. Use the updated scoring system." << endl;
                break;

            case OLDRANDCONF_OPT:
                conf = new RandomConfidence(optarg); 
                cout << "--rand-conf=" << string(optarg) << endl;
                break;

            case 'p': case 'r':
                conf = new RandomPValueConfidence(optarg);
                cout << "--pval-conf=" << string(optarg) << endl;
                break;

            // output options
            case BAGOUTPUT_OPT: bag_file = optarg; break;
            case T2TPERF_OPT: t2tperf_file = optarg; break;
            case SHOW_ALL_T2T_ROC_OPT: show_all_t2t_roc = true; break;
            case SHOW_GAPS_OPT: show_gaps = true; break;
            case ALLCONTEXT_OPT: only_good_context = false; break;
            case ANTITERMS_OPT: antifile = optarg; break;
            case 'c': conf_cutoff = atoi(optarg); break;
            case 'S': print_seq = false; break;

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

    // by default, use v1 confidence (since we don't know where to find
    // the data file necessary for rand-conf
    if(!conf) 
    {
        cerr << "You must specify a background distribution file with '-p expterm.dat'." << endl;
        exit(3);
    }

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

    if(optind+1 >= argc) usage();

    return optind;
}


// the main program
int
main(int argc, char * argv[])
{
    cout << "TransTermHP v2.08 (built on " << __DATE__ << ")" << endl;
    time_t start_time = time(NULL);

    // read commandline options
    int first_fileindex = process_options(argc, argv);
    print_options(cout);

    Genome chroms;

    // the rest of the args are filenames
    // NOTE: the seq file must come before the annotation file
    for(int i = first_fileindex; i < argc; i++)
    {
        // if the reader factor can find a reader, then we assume this is an
        // annotation file
        GeneReader * reader = gene_reader_factory(argv[i]);
        if(reader)
        {
            if(reader->good())
            {
                reader->read_genes(chroms);
                delete reader;
            }
            else
            {
                cerr << "Couldn't read: " << argv[i] << endl;
                exit(3);
            }
        }
        else
        {
            ifstream seq_file(argv[i]);
            if(!seq_file)
            {
                cerr << "Error: couldn't read file: " << argv[i]<< endl;
                exit(3);
            }
            read_seqs(seq_file, chroms);
        }
    }

    sort_genes(chroms);

    // for every sequence, find the terms
    for(EVERY_CHROM(chroms, C))
    {
        cerr << "Seq: " << (*C)->name << " (length " << (*C)->length << ", " 
             << (*C)->genes.size() << " genes) ";
        if((*C)->length <= unsigned(MAX_HP + 15)) 
        {
            cerr << endl 
                 << "Error: input sequences must have length > " << MAX_HP + 15
                 << endl;
            exit(3);
        }
        find_terms_dp(**C);
        cerr << endl;
    }
    cerr << endl;

    // analize the found terminators to make the confidence function
    conf->prepare(chroms);

    // output summary of the % of T2T regions we hit
    if(t2tperf_file != "")
    {
        ofstream out(t2tperf_file.c_str());
        if(!out)
        {
            cerr << "Couldn't open output file: " << t2tperf_file << endl;
            exit(3);
        }
        t2t_hitanal(out, chroms, *conf, 2, show_all_t2t_roc);
    }

    // output the main map of terminators & genes
    output_map(cout, chroms, *conf, conf_cutoff, print_seq, only_good_context); 

#if 0
    if(tthitfile != "")
    {
        ofstream out(tthitfile.c_str());
        plot_tthits_vs_terms(out, *conf, chroms);
    }
#endif

    if(antifile != "")
    {
        ofstream out(antifile.c_str());
        output_anti_terms(out, chroms, *conf, 90);
    }

    if(bag_file != "")
    {
        ofstream out(bag_file.c_str());
        output_best_term(out, *conf, chroms);
    }

    // print out the elapsed time
    time_t seconds = time(NULL) - start_time;
    cerr << "Wall clock time = " << seconds << " seconds." << endl;

    return 0;
}

