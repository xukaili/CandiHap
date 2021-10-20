#!/usr/bin/perl
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
sub usage {
    die(
        qq!
CandiHap: An R Platform for haplotype analysis on variation study

Usage:    perl  $0  -f genome.gff  -m ann.hmp  -p Phenotype.txt   -l LDkb  -c Chr:position
  e.g.    perl  $0  -f test.gff  -m haplotypes.hmp  -p Phenotype.txt  -l 50kb  -c 9:54583294

Command:  -m    input hmp file name (Must)
          -p    input phenotype file name (Must)
          -f    input gff file name (Must)
          -l    LD length by kb
          -c    Chr and position
          -h    this (help) message

Author:   Xukai Li, xukai_li\@sxau.edu.cn
Version:  V1.2.0
Update:   2021/04/26
Notes:    bioRxiv, 2020.02.27.967539.
          doi: https://doi.org/10.1101/2020.02.27.967539
\n!
    )
}

my %opts;
getopts('m:p:f:l:c:h', \%opts);
&usage unless ( exists $opts{m} && exists $opts{p} && exists $opts{f}  && exists $opts{l} && exists $opts{c});
&usage if $opt{h};
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "\nInput files are:\n";
foreach my $key ( keys %opts){
    print $key, "--> ",$opts{$key},"\n";
}

basename($opts{p}) =~ /(\S+)\.txt/;
$Phenotype = $1;
$opts{l} =~ /(\d+)/;
$kb = $1 * 1000;
print "\n\nThe LD is $opts{l} --> $kb bp\n\n";

open GFF ,"$opts{f}" or die "$!";
while (<GFF>) {
    if ($_ =~ /\tgene\t/) {
        $_ =~ s/\n*|\r*//g;
        @F = split(/\t/,$_);
        $F[8] =~ /ID=(\w+);?(.*)/;
        $lxk{$1} = "$F[0]\t$F[3]\t$F[4]\t$1\t$2";
    }
}
close GFF;
@F='';

print "The number of $i genes in $opts{l} LD are going to haplotype analysis ......\n\n\n";
@b = split(/:/,$opts{c});
$start = $b[1] - $kb;
$end   = $b[1] + $kb;
open RES,">>LD_$opts{l}_gene-$b[0]_$b[1].txt" or die $!;
for (keys %lxk) {
    @F = split(/\t/,$lxk{$_});
    if (($F[0] eq $b[0]) and ($F[1] > $start) and ($F[2] < $end)) {
        print  RES "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\n";
        system("perl  CandiHap.pl  -m $opts{m}  -p $opts{p}  -f $opts{f}  -g $F[3]");
        $i ++;
    }
}

print "Haplotype analysis of $i genes done.\n\n\n";
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
exit(1);