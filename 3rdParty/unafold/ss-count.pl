#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure 'gnu_getopt', 'no_auto_abbrev', 'no_ignore_case';

use constant R => 0.0019872;

sub version ($) {
    print "$_[0] (UNAFold) 3.6\n";
    print "By Nicholas R. Markham and Michael Zuker\n";
    print "Copyright (C) 2006\n";
    print "Rensselaer Polytechnic Institute\n";
    print "Troy, NY 12810-3590 USA\n";
    exit;
}

sub usage () {
    print <<EOF;
Usage: ss-count.pl [options] file.ct [file.ext]

Options:
-n, --normalized
-w, --weighted
-t, --temperature=<number>

EOF
    print 'Report bugs to markhn@rpi.edu', "\n";
    exit 0;
}

my ($weighted, $normalized, $temp) = (0, 0, 37);
GetOptions 'version|V' => sub { version('ss-count.pl.pl') }, 'help|h' => sub { usage() }, 'normalized|n' => \$normalized, 'weighted|w' => \$weighted, 'temperature|t=f' => \$temp;

unless (@ARGV) {
    print STDERR "Error: file not specified\nRun 'ss-count.pl -h' for help\n";
    exit 1;
}
my $ctFile = shift;
my $extFile = shift;

my ($bases, @ss, @seq, @prob);
my $structure = 0;
my $sum = 0.0;

if ($extFile) {
    open IN, "<$extFile" or die $!;
    my $header = <IN>;
    if ($header =~ /sequence/) {
	my $len = 0;
	while (<IN>) {
	    my ($s, $i, $p) = split;
	    if ($s == 1 and $i > $len) {
		$len = $i;
	    } elsif ($s == 2) {
		$i += $len;
	    }
	    $prob[$i] = $p;
	}
    } else {
	while (<IN>) {
	    my ($i, $p) = split;
	    $prob[$i] = $p;
	}
    }
    close IN or die $!;
}

my $weight = 1;
open IN, "<$ctFile" or die $!;
while (defined(my $header = <IN>)) {
    $bases = (split /\s+/, $header)[0];
    if ($weighted) {
	my ($dG) = $header =~ /dG = (.+?)\s/;
	$weight = exp(-$dG / R / (273.15 + $temp));
	$sum += $weight;
    }
    for (my $i = 0; $i < $bases; ++$i) {
	my $line = <IN>;
	my ($num, $chr, undef, undef, $pair) = split /\s+/, $line;
	$seq[$num] = $chr;
	$ss[$num] += $weight if $pair == 0;
    }
    ++$structure;
}
close IN or die $!;
$sum = $structure unless $weighted;

print "$structure\n" unless $normalized or $weighted;
for (my $i = 1; $i <= $bases; ++$i) {
    $ss[$i] or $ss[$i] = 0;
    if ($extFile) {
	my $sd = sqrt($prob[$i] * (1.0 - $prob[$i]));
	if ($sd == 0) {
	    printf "%d\t%.3f\t%.3f\t%.3f\t%s\n", $i, $prob[$i], 0, $ss[$i] / $sum, $ss[$i] / $sum == $prob[$i] ? '+0.000' : 'inf';
	} else {
	    printf "%d\t%.3f\t%.3f\t%.3f\t%+.3f\n", $i, $prob[$i], $sd, $ss[$i] / $sum, ($ss[$i] / $sum - $prob[$i]) / $sd;
	}
    } else {
	if ($weighted or $normalized) {
	    $ss[$i] = sprintf '%.3f', $ss[$i] / $sum;
	}
	print "$i\t$ss[$i]\t$seq[$i]\n";
    }
}
