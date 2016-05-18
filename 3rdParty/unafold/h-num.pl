#! /usr/bin/perl -w

use strict;
use warnings;

sub byHnum ($$) { $_[1][4] <=> $_[0][4] }

unless (@ARGV) {
    print "Usage: $0 <file prefix>\n";
    exit;
}

my ($prefix) = @ARGV;
my (@p_num, @h_num);

open IN, "<$prefix.plot" or die "Couldn't read $prefix.plot: $!";

scalar <IN>;

while (<IN>) {
    my ($level, $length, $istart, $jstart) = split;
    for (my $k = 0; $k < $length; ++$k) {
	++$p_num[$istart + $k];
	++$p_num[$jstart - $k];
    }
    push @h_num, [$level, $length, $istart, $jstart];
}

close IN or die $!;

foreach (@h_num) {
    my ($level, $length, $istart, $jstart) = @$_;
    my $h_num = 0;
    for (my $k = 0; $k < $length; ++$k) {
	$h_num += $p_num[$istart + $k] + $p_num[$jstart - $k];
    }
    push @$_, $h_num / $length;
}

open OUT, ">$prefix.h-num" or die "Couldn't write $prefix.h-num: $!";
print OUT "level\tlength\tistart\tjstart\th-num\n";
foreach my $line (sort byHnum @h_num) {
    printf OUT "%d\t%d\t%d\t%d\t%g\n", @$line;
}
close OUT or die $!;
