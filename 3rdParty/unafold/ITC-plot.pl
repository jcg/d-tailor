#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure 'gnu_getopt', 'no_auto_abbrev', 'no_ignore_case';

sub max ($$) {
    $_[0] > $_[1] ? $_[0] : $_[1];
}

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
Usage: ITC-plot.pl [options] file1 file2

Options:
-n, --NA=(RNA | DNA) (defaults to RNA)
-t, --temperature=<number> (defaults to 37)
-A, --A0=<total A>
-B, --B0=<total B>
-x, --exclude=(A|B|AA|BB)
    --fraction=<fraction of ensemble enthalpy> (defaults to 0.1)
    --nofraction
    --Tmelt=<melting temp. for single strands> (defaults to 50)
-o, --output=<file>
-r, --reuse
    --hyperbolic

EOF
print 'Report bugs to markhn@rpi.edu', "\n";
    exit;
}

my @t;
my ($NA, $fraction, $tMelt, $reuse, $hyper) = ('RNA', 0.1, 50, 0, 0);
my ($A0, $B0, @exclude, $prefix);

GetOptions 'version|V' => sub { version('ITC-plot.pl') }, 'help|h' => sub { usage() }, 'temperature|t=s' => \@t, 'NA|n=s' => \$NA, 'fraction=f' => \$fraction, 'nofraction' => sub { $fraction = undef; }, 'Tmelt=f' => \$tMelt, 'A0|A=f' => \$A0, 'B0|B=f' => \$B0, 'exclude|x=s' => \@exclude, 'reuse|r' => \$reuse, 'output|o=s' => \$prefix, 'hyperbolic' => \$hyper or die $!;
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
    print STDERR "Error: data not specified\nRun 'ITC-plot.pl -h' for help\n";
    exit 1;
}
$ARGV[0] and $ARGV[1] or die "Must specify files";
my ($file1, $file2) = @ARGV;
my ($prefix1, $prefix2);
($prefix1 = $file1) =~ s/\.seq$//;
($prefix2 = $file2) =~ s/\.seq$//;
defined $prefix or $prefix = "$prefix1-$prefix2.itc";

@t = split ',', join ',', @t;

my ($max, $change);
if (defined $A0 and defined $B0) {
    die 'Can\'t specify [A0] and [B0]';
} elsif (defined $A0) {
    $change = \$B0;
    $max = 2 * $A0;
} elsif (defined $B0) {
    $change = \$A0;
    $max = 2 * $B0;
} else {
    die 'Mus\'t specify [A0] or [B0]';
}

my @concArgs = ('--output', $prefix);
push @concArgs, '--exclude', 'A' if $exclA;
push @concArgs, '--exclude', 'B' if $exclB;
push @concArgs, '--exclude', 'AA' if $exclAA;
push @concArgs, '--exclude', 'BB' if $exclBB;

$exclA or -e "$prefix1.dG" or die "$prefix1.dG not found";
$exclB or -e "$prefix2.dG" or die "$prefix2.dG not found";
$exclAA or -e "$prefix1-$prefix1.dG" or die "$prefix1-$prefix1.dG not found";
$exclBB or -e "$prefix2-$prefix2.dG" or die "$prefix2-$prefix2.dG not found";
-e "$prefix1-$prefix2.dG" or die "$prefix1-$prefix2.dG not found";

my %index;
if ($reuse) {
    my ($i, $old) = (-1, -300);
    open IN, "<$prefix" or die $!;
    scalar <IN>;
    while (<IN>) {
	my ($t) = split or next;
	unless ($t == $old) {
	    $index{$t} = ++$i;
	    $old = $t;
	}
    }
    close IN or die $!;
} else {
    if (defined $fraction) {
	chomp(my $dHa = `sbs -n$NA $file1`);
	chomp(my $dHb = `sbs -n$NA $file2`);
	$dHa *= $fraction / 2.0;
	$dHb *= $fraction / 2.0;	
	my $dSa = $dHa / ($tMelt + 273.15) * 1000.0;
	my $dSb = $dHb / ($tMelt + 273.15) * 1000.0;
	push @concArgs, '--enthalpy', "$dHa,$dHb", '--entropy', "$dSa,$dSb";
    }

    my %data;
    for (my $c = $max / 50; $c <= 1.01 * $max; $c += $max / 50) {
	$$change = $c;
	system('concentration', '--A0', $A0, '--B0', $B0, @concArgs, $prefix1, $prefix2) == 0 or systemError('concentration');
	system('ensemble-dg', @concArgs, $prefix1, $prefix2) == 0 or systemError('ensemble-dg');
	system('dG2dH', '--points', 5, "$prefix.ens.dG") == 0 or systemError('dG2dH');
	open IN, "<$prefix.ens.H" or die $!;
	while (<IN>) {
	    next if /^#/;
	    chomp;
	    my ($T, $H) = split;
	    push @{$data{$T}}, $H * max($A0, $B0);
	}
	close IN or die $!;
    }

    open OUT, ">$prefix" or die $!;
    print OUT "#T\tC\tH\tdH/dc\n";
    my $i = -1;
    foreach my $T (sort {$a <=> $b} keys %data) {
	$index{$T} = ++$i;
	my $old = shift @{$data{$T}};
	for (my $c = $max / 25; @{$data{$T}}; $c += $max / 50) {
	    my $H = shift @{$data{$T}};
	    printf OUT "%g\t%g\t%g\t%g\n", $T, $c, $H, ($H - $old) / $max * 50;
	    $old = $H;
	}
	print OUT "\n";
    }
    close OUT or die $!;
}

my @temps;
foreach my $temp (@t) {
    push @temps, $temp if defined $index{$temp};
}
@temps or push @temps, (keys %index)[0];

my $LINEWIDTH = 2;
open PLOT, ">$prefix.gp" or die $!;
print PLOT "set terminal postscript enhanced color\n";
print PLOT "set output '$prefix.ps'\n";
print PLOT "set data style lines\n";
print PLOT "set key bottom right\n";
printf PLOT "set xrange [0:%g]\n", $max;
if ($change == \$A0) {
    print PLOT "set xlabel '[A_0] (M)'\n";
    print PLOT "set ylabel '{/Symbol DD}H / {/Symbol D}[A_0] (kcal/mol)'\n";
} else {
    print PLOT "set xlabel '[B_0] (M)'\n";
    print PLOT "set ylabel '{/Symbol DD}H / {/Symbol D}[B_0] (kcal/mol)'\n";
}
my $command = "'$prefix' every :::$index{$temps[0]}::$index{$temps[0]} u 2:" . ($hyper ? 3 : 4 ) . " t '$temps[0]\{/Symbol \\260\}' lw $LINEWIDTH";
for (my $i = 1; $i < @temps; ++$i) {
    $command .= ", '$prefix' every :::$index{$temps[$i]}::$index{$temps[$i]} u 2:" . ($hyper ? 3 : 4 ) . " t '$temps[$i]\{/Symbol \\260\}' lw $LINEWIDTH";
}
print PLOT "plot $command\n";
close PLOT or die $!;

