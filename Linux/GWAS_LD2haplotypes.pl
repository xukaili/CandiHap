#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  ./genome.gff  ./ann.hmp  ./Phenotype.txt   LDkb  Chr:position
e.g. perl  GWAS_LD2haplotypes.pl  ./test.gff  ./haplotypes.hmp  ./Phenotype.txt  50kb  9:54583294
perl  GWAS_LD2haplotypes.pl  ./Si.Ann_RNA.gff   ./haplotypes.hmp4   ./Millet-b-BLUP.mGWAS.txt   50kb  4:37876478

USAGE
print $usage if(@ARGV==0);
exit if(@ARGV==0);
$ARGV[2] =~ /.*\/(\S+)\.txt/;
$Phenotype = $1;
$ARGV[3] =~ /(\d+)/;
$kb = $1 * 1000;
print "The LD is $ARGV[3] --> $kb bp\n\n\n";

open GFF ,"$ARGV[0]" or die "$!";
while (<GFF>) {
    if ($_ =~ /\tgene\t/) {
        $_ =~ s/\n*|\r*//g;
        @F = split(/\t/,$_);
        $F[8] =~ /ID=(\w+);?(.*)/;
        $lxk{$1} = "$F[0]\t$F[3]\t$F[4]\t$1\t$2";
        $i ++;
    }
}
close GFF;
@F='';

print "The number of $i genes in $ARGV[3] LD are going to haplotype analysis ......\n\n\n";
open RES,">>LD_$ARGV[3]_gene_$ARGV[4].txt" or die $!;
@b = split(/:/,$ARGV[4]);
$start = $b[1] - $kb;
$end   = $b[1] + $kb;
for (keys %lxk) {
    @F = split(/\t/,$lxk{$_});
    if (($F[0] eq $b[0]) and ($F[1] > $start) and ($F[2] < $end)) {
        print  RES "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\n";
        system("perl  CandiHap.pl  -m $ARGV[1]   -f $ARGV[0] -p $ARGV[2] -g  $F[3]");
    }
}
print "Haplotype analysis of $i genes done.\n\n\n";

