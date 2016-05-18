#!/usr/bin/env python

import sys

def make_zero_matrix(n):
    A = []
    for i in range(0, n):
        A.append([0] * n)
    return A

def dxdy(at, num_bins):
    dx = (max_hp[at] - min_hp[at]) / float(num_bins)
    dy = (max_tail[at] - min_tail[at]) / float(num_bins)
    return dx, dy


def set_ranges_from_file(termfile, num_bins):
    global min_hp, max_hp, min_tail, max_tail

    min_hp = {}
    max_hp = {}
    min_tail = {}
    max_tail = {}

    infile = open(termfile)
    for line in infile:
        s = line[:-1].split()
        at, hp, tail = s[0], float(s[1]), float(s[2])
        if at not in min_hp or hp < min_hp[at]:
            min_hp[at] = hp
        if at not in max_hp or hp > max_hp[at]:
            max_hp[at] = hp
        if at not in min_tail or tail < min_tail[at]:
            min_tail[at] = tail
        if at not in max_tail or tail > max_tail[at]:
            max_tail[at] = tail

    # add in a buffer in the range = to the size of 1 bin on the low end
    # and 0.1 on the high end (the later is for numerical problems)
    for at in min_hp:
        # NOTE: b/c we change the ranges here, this dx dy is only an estimate...
        dx, dy = dxdy(at, num_bins)
        min_hp[at] -= 1.5*dx
        min_tail[at] -= 1.5*dy
        max_hp[at] += 0.1
        max_tail[at] += 0.1

    infile.close()


def hist2d_from_file(termfile, num_bins): # hp_range, tail_range):
    # compute the bin sizes
    #dx = (hp_range[1] - hp_range[0]) / float(num_bins)
    #dy = (tail_range[1] - tail_range[0]) / float(num_bins)

    D = {}

    infile = open(termfile)
    for line in infile:
        s = line[:-1].split()
        at, hp, tail = s[0], float(s[1]), float(s[2])

    #    if warn_if_out_of_range(hp, hp_range, "hairpin") or \
    #       warn_if_out_of_range(tail, tail_range, "tail"):
    #        continue

        if at not in D:
            D[at] = make_zero_matrix(num_bins)

        dx, dy = dxdy(at, num_bins)

        i = int((hp - min_hp[at])/dx)
        j = int((tail - min_tail[at])/dy)
        
        if i == 0 or j == 0:
            print >> sys.stderr, at, i, j, hp, tail, min_hp[at], min_tail[at], dx, dy
        
        if not (0 <= i < num_bins and 0 <= j < num_bins):
            print >> sys.stderr, "WARNING: out of range values:", i, j, at, hp, tail
            print >> sys.stderr, "Ranges=", min_hp[at], max_hp[at], min_tail[at], max_tail[at]
            continue

        D[at][i][j] += 1

    infile.close()
    return D


warned = {}
def warn_if_out_of_range(value, rng, title):
    if value <= rng[0] or value >= rng[1]:
        if title not in warned:
            print >> sys.stderr, "@" * 60
            print >> sys.stderr, """WARNING: random %s generated with energy lower than supplied range.
            Such examples are ignored. I suggest you re-run calibrate.sh after changing
            the lowerbound in the range variable therein.""" % (title)
            print >> sys.stderr, "Range = ", rng, "Seen = ", value
            print >> sys.stderr, "@" * 60
            warned[title] = True
        return True
    return False


def print_matrix(A):
    """Write out the matrix (transposed)"""
    for j in range(len(A)):
        for i in range(len(A)):
            print A[i][j],
        print


def main():
    # read the input
    infile = sys.argv[1]
    seqlen = int(sys.argv[2])
    num_bins = int(sys.argv[3])
    #hp_range = (float(sys.argv[4]), float(sys.argv[5]))
    #tail_range = (float(sys.argv[6]), float(sys.argv[7]))

    # print the header
    print seqlen, num_bins # hp_range[0], hp_range[1], tail_range[0], tail_range[1]

    set_ranges_from_file(infile, num_bins)
    D = hist2d_from_file(infile, num_bins) #, hp_range, tail_range)

    # for every at value, compute and print the matrix
    for at in sorted(D):
        print at, min_hp[at], max_hp[at], min_tail[at], max_tail[at]
        print_matrix(D[at])

if __name__ == '__main__': main()


def read_random_terms(infile, at):
    """Returns a 3-tuple: (D, HPRange, TailRange), where D is a dict maping %at values
    to a list of (hp, tail) tuples, and the ranges are pairs giving the (min, max) seen"""

    L = []
    inf=open(infile)
    first = True
    for line in inf:
        line = line[:-1]
        at, hp, tail = line.split()
        at, hp, tail = float(at), float(hp), float(tail)
        if at not in D: D[at] = []
        D[at].append((hp, tail))

        # track the global ranges of the hp and tail scores
        if first:
            minhp = maxhp = hp
            mintail = maxtail = tail
            first = False
        else:
            minhp, maxhp = min(minhp, hp), max(maxhp, hp)
            mintail, maxtail = min(mintail, tail), max(maxtail, tail)
        
    inf.close()
    return (D, (minhp, maxhp), (mintail, maxtail))


def hist2d(terms, num_bins, hp_range, tail_range):
    """Build a matrix that is the 2d histogram"""

    # compute the bin sizes
    dx = (hp_range[1] - hp_range[0]) / float(num_bins)
    dy = (tail_range[1] - tail_range[0]) / float(num_bins)

    # make a num_bins by num_bins zero matrix
    A = []
    for i in range(0,num_bins):
        A.append([0] * num_bins)

    # fill in the matrix
    total = 0
    for i in range(0, num_bins):
        hp_slice = [(hp, tail) for (hp, tail) in terms
                               if hp >= i * dx + hp_range[0] and hp < (i+1)*dx + hp_range[0]]

        for j in range(0, num_bins):
            A[i][j] = len([1 for (hp, tail) in hp_slice
                             if tail >= j * dy + tail_range[0] and tail < (j+1)*dy + tail_range[0]])
            total += A[i][j]

    return A

