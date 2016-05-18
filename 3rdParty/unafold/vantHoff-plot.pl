#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure 'gnu_getopt', 'no_auto_abbrev', 'no_ignore_case';

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
Usage: vantHoff-plot.pl [options] file1 file2

Options:
-n, --NA=(RNA | DNA) (defaults to RNA)
-C, --Ct=<total strand concentration>
-x, --exclude=(A|B|AA|BB)
    --fraction=<fraction of ensemble enthalpy> (defaults to 0.1)
    --nofraction
    --Tmelt=<melting temp. for single strands> (defaults to 50)
-o, --output=<file>

EOF
    print 'Report bugs to markhn@rpi.edu', "\n";
    exit;
}

my ($NA, $Ct, $fraction, $tMelt) = ('RNA', 2e-5, 0.1, 50);
my (@exclude, $prefix);

GetOptions 'version|V' => sub { version('vantHoff-plot.pl') }, 'help|h' => sub { usage() }, 'NA|n=s' => \$NA, 'fraction=f' => \$fraction, 'nofraction' => sub { $fraction = undef; }, 'Tmelt=f' => \$tMelt, 'Ct|C=f' => \$Ct, 'exclude|x=s' => \@exclude, 'output|o=s' => \$prefix or die $!;
my ($exclA, $exclB, $exclAA, $exclBB);
foreach my $x (@exclude) {
    if ($x eq 'A') {
	++$exclA;
    } elsif ($x eq 'B') {
	++$exclB;
    } elsif ($x eq 'AA') {
	++$exclAA;
    } elsif ($x eq 'BB') {
	++$exclBB;
    }
}

unless (@ARGV) {
    print STDERR "Error: data not specified\nRun 'vantHoff-plot.pl -h' for help\n";
    exit 1;
}

my ($file1, $file2) = @ARGV;
my $same = not defined $file2;
my ($prefix1, $prefix2);
($prefix1 = $file1) =~ s/\.seq$//;
($prefix2 = $file2) =~ s/\.seq$// unless $same;
defined $prefix or $prefix = $same ? "$prefix1-$prefix1.vh" : "$prefix1-$prefix2.vh";

my @concArgs = ('--output', $prefix);
push @concArgs, '--exclude', 'A' if $exclA;
push @concArgs, '--exclude', 'B' if $exclB and not $same;
push @concArgs, '--exclude', 'AA' if $exclAA;
push @concArgs, '--exclude', 'BB' if $exclBB and not $same;

$exclA or -e "$prefix1.dG" or die "$prefix1.dG not found";
$exclAA or -e "$prefix1-$prefix1.dG" or die "$prefix1-$prefix1.dG not found";
unless ($same) {
    $exclB or -e "$prefix2.dG" or die "$prefix2.dG not found";
    $exclBB or -e "$prefix2-$prefix2.dG" or die "$prefix2-$prefix2.dG not found";
    -e "$prefix1-$prefix2.dG" or die "$prefix1-$prefix2.dG not found";
}

if (defined $fraction) {
    chomp(my $dHa = `sbs -n$NA $file1`);
    $dHa *= $fraction / 2.0;
    my $dSa = $dHa / ($tMelt + 273.15) * 1000.0;
    if ($same) {
	push @concArgs, '--enthalpy', $dHa, '--entropy', $dSa;
    } else {
	chomp(my $dHb = `sbs -n$NA $file2`);
	$dHb *= $fraction / 2.0;
	my $dSb = $dHb / ($tMelt + 273.15) * 1000.0;
	push @concArgs, '--enthalpy', "$dHa,$dHb", '--entropy', "$dSa,$dSb";
    }
}

open OUT, ">$prefix" or die $!;
print OUT "#Ct\t1/TmConc\t1/TmExt1\t1/TmExt2\n";
for (my $c = log($Ct) - log(10); $c <= log($Ct) + log(10); $c += log(10) / 25) {
    my $conc = exp($c);
    if ($same) {
	system('concentration-same', '--A0', $conc, @concArgs, $prefix1) == 0 or systemError('concentration-same');
	system('concentrations-same.pl', $prefix) == 0 or systemError('concentrations-same.pl');
	system('ensemble-ext-same', '--NA', $NA, '--points', 5, @concArgs, $prefix1) == 0 or systemError('ensemble-ext-same');
    } else {
	$conc /= 2;
	system('concentration', '--A0', $conc, '--B0', $conc, @concArgs, $prefix1, $prefix2) == 0 or systemError('concentration');
	system('concentrations.pl', $prefix) == 0 or systemError('concentrations.pl');
	system('ensemble-ext', '--NA', $NA, '--points', 5, @concArgs, $prefix1, $prefix2) == 0 or systemError('ensemble-dg');
    }
    my ($TmConc, $TmExt1, $TmExt2);
    open IN, "<$prefix.TmConc" or die $!;
    chomp($_ = <IN>);
    ($TmConc) = split;
    close IN or die $!;
    open IN, "<$prefix.ens.TmExt1" or die $!;
    while (<IN>) {
	chomp;
	($TmExt1) = split;
    }
    close IN or die $!;
    open IN, "<$prefix.ens.TmExt2" or die $!;
    chomp($_ = <IN>);
    ($TmExt2) = split;
    close IN or die $!;
    printf OUT "%g\t%g\t%g\t%g\n", $same ? $conc : 2 * $conc, 1 / (273.15 + $TmConc), 1 / (273.15 + $TmExt1), 1 / (273.15 + $TmExt2)
}
close OUT or die $!;

my $LINEWIDTH = 2;
open PLOT, ">$prefix.gp" or die $!;
print PLOT "set terminal postscript enhanced color\n";
print PLOT "set output '$prefix.ps'\n";
print PLOT "set data style lines\n";
print PLOT "set key right\n";
print PLOT "set xlabel 'C_T (M)'\n";
print PLOT "set ylabel '1 / T_m'\n";
print PLOT "set logscale x\n";
print PLOT "plot '$prefix' u 1:2 t 'T_m(Conc)' lw $LINEWIDTH, '$prefix' u 1:3 t 'T_m(Ext1)' lw $LINEWIDTH, '$prefix' u 1:4 t 'T_m(Ext2)' lw $LINEWIDTH\n";

