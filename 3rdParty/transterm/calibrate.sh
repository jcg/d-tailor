#!/bin/sh

#=============================
# Global Constants
#=============================

#HP_RANGE="-30 -2"      # range of possible hp scores
#TAIL_RANGE="-7.15 -2"    # range of possible tail scoress
NUM_BINS=25           # number of bins in the histogram
SAMPLESIZE=20000000   # length of random dna sequence to generate

#HP_RANGE="-25 -2"      # range of possible hp scores
#TAIL_RANGE="-7 -2"    # range of possible tail scoress
#NUM_BINS=10           # number of bins in the histogram
#SAMPLESIZE=200000   # length of random dna sequence to generate

# values for %AT that we generate data for:
AT="0.26 0.28 0.30 0.32 0.34 0.36 0.38"
AT="$AT 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58"
AT="$AT 0.60 0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78" 
AT="$AT 0.80 0.82 0.84"

#=============================
# Process commandline options
#=============================

if [ ${#} -lt 1 ] ; then
    echo "usage: `basename $0` outputfile.dat [transterm options]" > /dev/stderr
    exit 1
fi

output=$1
shift

#===========================
# Find helper scripts
#===========================

findprog() {
    if [ -x "./$1" ] ; then
        echo "./$1"
    elif [ -x "`which $1`" ] ; then
        echo "$1"
    else
        echo "`basename $0`: $1 must be in the current directory or on your PATH." > /dev/stderr
        exit 2
    fi
}

RAND=`findprog random_fasta.py`
TRANSTERM=`findprog transterm`
MAKE_EXPTERM=`findprog make_expterm.py`

#==============================
# Ensure python is avialable
#==============================

if [ ! -x "`which python`" ] ; then
    echo "`basename $0`: python must be installed to run calibration scripts" > /dev/stderr
    exit 3
fi

#=========================================
# Generate and calibrate on random data
#=========================================

rm -f random_terms.dat

echo "NOTE: warnings about 'using version 1.0 confidence' are expected and OK." > /dev/stderr
echo "TransTerm Options = " $*

$TRANSTERM $* --v1-conf -S -c 0 --all-context /dev/null /dev/null > /dev/stderr

for at in $AT ; do
    echo "`basename $0`: Running TransTerm on random sequence with %AT = $at" > /dev/stderr

    #================================
    # Generate random DNA data
    #================================
    $RAND $at $SAMPLESIZE "random$at.fasta" "random$at.coords"

    #==================================================
    # Run Transterm on that random sequence 
    # output all terminators regardless of confidence
    # in a striped down format
    #==================================================
    $TRANSTERM $* --v1-conf -S -c 0 "random$at.fasta" "random$at.coords" \
        | awk '/TERM/ { print at, $9, $10 }' at=$at >> random_terms.dat

    #=================================================
    # Clean everything up
    #=================================================
    rm random$at.coords random$at.fasta
done

#===================================
# Make the actual expterms.dat file
#===================================
echo "Creating distribution file in $output" > /dev/stderr
$MAKE_EXPTERM random_terms.dat $SAMPLESIZE $NUM_BINS > $output

echo "-1" >> $output
$TRANSTERM $* -S -c 0 /dev/null /dev/null | grep -- "--" >> $output

echo "Done. You can run transterm with this background distribution using '-r $output'"

