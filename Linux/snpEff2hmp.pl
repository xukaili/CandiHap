#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  snpEff_Your.vcf
e.g.
     perl  snpEff2hmp.pl  snpEff_test.vcf
USAGE

%dic=('AC' => "M",'CA' => 'M',"AG" => "R","GA" => "R","AT" => "W","TA" => "W","GC" => "S","CG" => "S","CT" => "Y","TC" => "Y","GT" => "K","TG" => "K","--" => "-","++" => "+","+-" => "H","-+" => "H","AA" => "A","TT" => "T","GG" => "G","CC" => "C","NA" => "N");

print $usage if(@ARGV==0);
exit if(@ARGV==0);

open OUT,">haplotypes_$ARGV[0].hmp" or die $!;

print "Reading snpEff file ......\n";
open SNPEFF,$ARGV[0] or die $!;
#open(SNPEFF,"gzip -dc  $ARGV[0]|") or die ("can not open $infilename\n");
while(<SNPEFF>){
    $_ =~ s/\r*\n*//g;
    if ($_ =~ /CHROM/) {
        @names = split(/\t/,$_);
        print OUT  "$names[0]\t$names[1]\t$names[3]\t$names[4]\t$names[7]\tFrequency\t" . join ("\t",@names[9..$#names]) . "\n";
    }
    next if $_ =~ /^#/;
    @a = split(/\t/,$_);
    $ref =$a[3];
    $alt = $a[4];
    $name = "$a[0]_$a[1]_$a[3]_$a[4]";
    $N{$name} = 0; $N0{$name} = 0; $N1{$name} =0;
    $ann = "$a[7]";
    $h{$chrpos} = $ann;
    for ($i = 9; $i <= $#a; $i++) {
        if ($a[$i] =~ /(\.)\/(\.)/) {
            $a[$i] = 'N';
            $N{$name} ++  if ($1 eq '.');
            $N{$name} ++  if ($2 eq '.');
        }
        elsif ($a[$i] =~ /(.)(.)(.)/ ) {
            $Fir = $1;
            $Sen = $3;
            $SeqRef = length($a[3]);
            $SeqAlt = length($a[4]);
            $N0{$name} ++ if ($Fir eq 0);
            $N0{$name} ++ if ($Sen eq 0);
            $N1{$name} ++ if ($Fir >  0);
            $N1{$name} ++ if ($Sen >  0);
            @alt = split(/,/,$a[4]);
            unshift( @alt, $a[3]);
            if  ($SeqRef == $SeqAlt) {
                $dd = $alt[$Fir].$alt[$Sen];
                $a[$i]  = $dic{$dd};
            }
            
            if  ($SeqRef != $SeqAlt) {
            $a[3] = '+' if  ($SeqRef > $SeqAlt);
            $a[3] = '-' if  ($SeqRef < $SeqAlt);
            $a[4] = '-' if  ($SeqRef > $SeqAlt);
            $a[4] = '+' if  ($SeqRef < $SeqAlt);
            @alt = split(/,/,$a[4]);
            unshift( @alt, $a[3]);
            $dd2 = $alt[$Fir].$alt[$Sen];
            $a[$i]  = $dic{$dd};
            }
        }
    }
   # $maf = $N1{$name} / ($N0{$name} + $N1{$name} + $N{$name});
    $Frequency = "Ref: $N0{$name}, Alt: $N1{$name}, N: $N{$name}, Maf: $maf";
    $l = "$a[0]\t$a[1]\t$ref\t$alt\t$a[7]\t$Frequency\t" . join ("\t",@a[9..$#a]);
    print OUT "$l\n";
}
close SNPEFF;
print "All done.\n"; 



