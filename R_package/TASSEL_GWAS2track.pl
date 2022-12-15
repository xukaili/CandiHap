#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  Phenotype.gwas.txt
e.g.
     cd  Your_dir
     perl  GWAS2track.pl  Phenotype.gwas.txt
     
     cd /Users/lxk/GitHub-lxk/test/
     perl  GWAS2track.pl  Phenotype.gwas.txt
USAGE
print $usage if(@ARGV==0);
exit if(@ARGV==0);

open OUT,">track-$ARGV[0]" or die $!;
print OUT "seqnames\tstart\twidth\tstrand\tscores\tr2\n";
open GWAS ,"$ARGV[0]" or die "$!";
while (<GWAS>) {
    next if $_ =~ /^Trait/;
    $_ =~ s/\n*|\r*//g;
    @F = split(/\t/,$_);
    $log = -log($F[6])/log(10);
    print OUT "$F[2]\t$F[3]\t1\t*\t$log\t$F[14]\n";
}
close GWAS;
#seqnames	start	width	strand	scores	r2
#9	54506145	1	*	11.8079357806654	0.191464850728128

