#! /usr/bin/perl -w

use strict;
use warnings;

unless (@ARGV) {
    print "Usage: $0 file.plot file.ct\n";
    exit;
}

my $plot = shift;
my $ct = shift or die 'Must specify .ct file';

my @prob;

open IN, "<$plot" or die $!;
scalar <IN>;
while (<IN>) {
    chomp;
    my @fields = split;
    my ($i, $j, $p) = splice @fields, $#fields - 2;
    $prob[$i][$j] = $prob[$j][$i] = $p;
}
close IN or die $!;

my @ss = (1) x (1 + @prob);
for (my $i = 1; $i < @prob; ++$i) {
    for (my $j = 1; $j < @prob; ++$j) {
	$ss[$i] -= $prob[$i][$j] if defined $prob[$i][$j];
    }
}

open IN, "<$ct" or die $!;
scalar <IN>;
while (<IN>) {
    chomp;
    my ($base, undef, undef, undef, $pair) = split;
    if ($pair) {
	$prob[$base][$pair] ||= 0;
	print "$base\t$prob[$base][$pair]\n";
    } else {
	$ss[$base] ||= 1;
	print "$base\t$ss[$base]\n";
    }
}
close IN or die $!;
