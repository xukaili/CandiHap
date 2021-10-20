#!/usr/bin/perl
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
use threads;
sub usage {
    die(
        qq!
CandiHap: An R Platform for haplotype analysis on variation study

Usage:    perl  $0  -m Your.hmp  -f Genome.gff  -p Phenotype.txt  -g Your_gene_ID
  e.g.    perl  CandiHap.pl  -m haplotypes.hmp  -f test.gff  -p Phenotype.txt  -g Si9g49990 -n 1
  e.g.    perl  CandiHap.pl  -m haplotypes.hmp1  -f Si.Ann_RNA.gff -p Hull-a-BLUP.mGWAS.txt  -g Si1g06520 -n 1

Command:  -m    input hmp file name (Must)
          -p    input phenotype file name (Must)
          -f    input gff file name (Must)
          -g    Your gene ID  (Must)
          -k    keek all tmp files
          -h    this (help) message
          -u    gene upstream. default is  2000 bp.
          -d    gene downstream. default is  500 bp.
          -n    input pop file name and plot haploNet figure. default is 0. require R package "pegas".

Author:   Xukai Li, xukai_li\@sxau.edu.cn
Version:  V1.2.0
Update:   2021/09/06
Notes:    CandiHap: a haplotype analysis toolkit for natural variation study.
          bioRxiv, 2020.02.27.967539.
          doi: https://doi.org/10.1101/2020.02.27.967539
\n!
    )
}

my %opts;
getopts('m:p:f:g:n:kh', \%opts);
&usage unless ( exists $opts{m} && exists $opts{p} && exists $opts{f}  && exists $opts{g});
&usage if $opt{h};
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input files are:\n";
foreach my $key ( keys %opts){
    print $key, "--> ",$opts{$key},"\n";
}
$opts{k}=0 unless defined($opts{k});
$opts{n}=0 unless defined($opts{n});
$opts{u}=2000 unless defined($opts{u});
$opts{d}=500 unless defined($opts{d});

basename($opts{p}) =~ /(\S+)\.txt/;
$Phenotype = $1;
$mychr = $mystart = $myend =0;
print "\n\nReading $opts{f} GFF file for gene $opts{g} ......\n";
system("grep  '$opts{g}' $opts{f} > tmp-gene_gff-$opts{g}-$Phenotype.txt");
open GFF ,"tmp-gene_gff-$opts{g}-$Phenotype.txt" or die "$!";
while (<GFF>) {
    $_ =~ s/\n*|\r*//g;
    if($_ =~ /\texon\t/){
        @F = split(/\t/,$_);
        $F[8]=~/Parent=(.*?);/;
        $cds_pos{$1}{"$F[3] $F[4] $F[7]"}=1;
        $contig{$1}=$F[0];
        $stream{$1}=$F[6];
    }
    if ($_ =~ /\tgene\t/){
        $_ =~ /;Note=(\d+) /;
        $gene_n = $1;
        @gff = split(/\t/,$_);
        $mychr = $gff[0];
        if ($gff[6] eq '+') {
            $mystart = $gff[3] - $opts{u};
            $myend   = $gff[4] + $opts{d};
        }
        elsif ($gff[6] eq '-') {
            $mystart = $gff[3] - $opts{d};
            $myend   = $gff[4] + $opts{u};
        }
    }
}
close GFF;
$gene_n = `grep 'mRNA' tmp-gene_gff-$opts{g}-$Phenotype.txt |wc -l`;
$gene_n =~ s/\s+//g;
print "Read $opts{f} file done.\n\n\n";
print "Query:  $opts{g} $gene_n --> $mychr  $mystart  $myend\n\n";

print "Reading $opts{p} file for phenotype ......\n";
open PHENOTYPE ,"$opts{p}" or die "$!";
while(<PHENOTYPE>){
    $_ =~ s/\n*|\r*//g;
    @F = split(/\t/,$_);
    $name{$F[0]} = $F[1] if $_ ne 'NaN';
}
close PHENOTYPE;
print "Read $opts{p} file done.\n\n\n";

print "Reading $opts{m} file for gene $opts{g} ......\n";
#system("grep -E '$opts{g}|CHROM' $opts{m}$mychr > tmp-hmp-$opts{g}-$Phenotype.txt");
system("grep -E '$opts{g}|CHROM' $opts{m} > tmp-hmp-$opts{g}-$Phenotype.txt");

open HMP, "tmp-hmp-$opts{g}-$Phenotype.txt" or die $!;
while (<HMP>) {
    $_ =~ s/\r*\n*//g;
    if ($_ =~ /CHROM/) {
        @names = split(/\t/,$_);
        $n  = "$names[1]\t$names[2]\t$names[3]\t$names[4]\t$names[5]\t" . join ("\t",@names[6..$#names]);
    }
    next if $_ !~ /^$mychr\t/;
    #next if ($_ =~ /:: synonymous/);
    @a = split(/\t/,$_);
    next if ($a[1] < $mystart or $a[1] > $myend);
    $BP{$a[1]} = $a[2] . '->'. $a[3];
    $INFO{$a[1]} = $a[4];
    #$a[7] =~ /EFF=(.*)/;
    for ($i = 6; $i <= $#a; $i++){
        if (length($a[2]) ==1  and length($a[3]) ==1 ) {
            $a[$i] =~ /(.)(.)(.)/;
            $seq{$names[$i]} .=  $1 if ($opts{n});
        }
        if (($name{$names[$i]} > 0 or $name{$names[$i]} eq '0') and ($a[$i] !~ /N/)) {
            $a[$i] =~ /(.*)\/(.*)/;
            next if ($1 ne $2);
            next if ($a[$i] !~ /[ATCG+-]/);
            $TT = "$a[1]\t$a[$i]";
            $Ttest{$TT} .= "$TT\t$name{$names[$i]}\t$names[$i]\n";
        }
    }
    $l = "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t" . join ("\t",@a[6..$#a]);
    push @lxk , $l;
}
push @lxk , $n;
print "Read $opts{m} file done.\n\n\n";

open OUT,">tmp-wilcox-$opts{g}-$Phenotype.txt" or die $!;
foreach  $k (sort {$a<=>$b} keys %Ttest ){
    print OUT "$Ttest{$k}";
}

open OUTR,">wilcox-$opts{g}-$Phenotype.R" or die $!;
print OUTR <<EOF;
data <- read.table("tmp-wilcox-$opts{g}-$Phenotype.txt", header=F,sep="\\t")
for (i in unique(data\$V1)) {
    test =''
    pvalue =1
    merge=''
    mydata=data[which(data\$V1 == i), ]
    if (length(unique(mydata\$V2)) ==1)
    merge = cbind(i,pvalue)
    if (length(unique(mydata\$V2)) ==1)
    write.table(merge, file = "tmp-wilcox_result-$opts{g}-$Phenotype.txt", sep = "\\t", quote = FALSE, col.names = FALSE, row.names = FALSE, append = T)
    if (length(unique(mydata\$V2)) ==1)
    next
    test <- wilcox.test(V3 ~ V2, data = mydata)
    pvalue <- test\$p.value
    #fdr=p.adjust(pvalue, method="fdr", n=nrow(mydata))
    merge = cbind(i,pvalue)
    write.table(merge, file = "tmp-wilcox_result-$opts{g}-$Phenotype.txt", sep = "\\t", quote = FALSE, col.names = FALSE, row.names = FALSE, append = T)
}
EOF
system("Rscript wilcox-$opts{g}-$Phenotype.R");

@a ='';
open TEST, "tmp-wilcox_result-$opts{g}-$Phenotype.txt" or die $!;
while (<TEST>) {
    $_ =~ s/\r*\n*//g;
    @a = split(/\t/,$_);
    $wilcox{$a[0]} = $a[1];
    push (@b , $a[1]) if ($a[1] > 0);
}
@b = grep(/[\d]/, @b);
@b=sort {$a<=>$b}@b;
$wilcox_max = $b[$#b];
$wilcox_min = $b[0];
$wilcox_min = 0.01  if ($wilcox_min > 0.01);
@b='';

#$threshold = 10 ** (log($wilcox_min)/log(10) * 0.6);
$threshold = 1;
print "The minimum wilcox test is: $wilcox_min\nThe wilcox test threshold is: $threshold\n\n";

open OUT,">tmp-$opts{g}-$Phenotype.txt" or die $!;
for (@lxk) {
    my @a = split /\t/, $_;
    for (my $i=0;$i<@a;$i++){
        $out{$i}{$a[0]} = $a[$i];
    }
}
foreach my $k (sort {$a<=>$b} keys %out ){
    my $out;
    foreach $d (sort {$a<=>$b} keys %{$out{$k}}){
        $out .= "$out{$k}{$d}\t"  if ($wilcox{$d} < $threshold);
    }
    $out =~ s/\t$//g;
    print OUT "$out\n";
}
close IN;
close OUT;


open OUTR,">LDheatmap-$opts{g}-$Phenotype.R" or die $!;
print OUTR <<EOF;
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (! require("LDheatmap")) BiocManager::install("LDheatmap")
if (! require("genetics")) install.packages("genetics")
data <- read.table("tmp-$opts{g}-$Phenotype.txt", header=F,sep="\t",row.names=1)
NAs <- data == "N/N"
is.na(data)[NAs] <- TRUE
data[NAs] <- NA
Dist = as.numeric(data[1,])
mydata = data[-1,]
mydata = mydata[-1,]
mydata = mydata[-1,]
mydata = mydata[-1,]
snp = mydata[-1,]
num<-ncol(snp) 
 for(i in 1:num){
   snp[,i]<-as.genotype(snp[,i]) 
 }
pdf(file="Plot-LDheatmap-$opts{g}.pdf", width=6, height=4)
rgb.palette <- colorRampPalette(rev(c("#1F77B4", "#AEC7E8", "#FF7F0E")), space = "rgb")
LDheatmap(snp, genetic.distances = Dist, flip = TRUE, color=rgb.palette(40))
dev.off()
EOF
system("Rscript LDheatmap-$opts{g}-$Phenotype.R");
print "Plotted Gene LDheatmap to file Plot-LDheatmap-$opts{g}.pd\n\n";

@a ='';
open TEST, "tmp-wilcox_result-$opts{g}-$Phenotype.txt" or die $!;
while (<TEST>) {
    $_ =~ s/\r*\n*//g;
    @a = split(/\t/,$_);
    $wilcox{$a[0]} = $a[1];
    push (@b , $a[1]) if $a[1] > 0;
}
@b = grep(/[\d]/, @b);
@b=sort {$a<=>$b}@b;
$wilcox_max = $b[$#b];
$wilcox_min = $b[0];
@b='';

@a ='';
open FILE, "tmp-$opts{g}-$Phenotype.txt" or die $!;
while (<FILE>) {
    if ($_ =~ /^POS/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $POS = "SNP positions\t". join ("\t",@a[1..$#a]) . "\tNumber of varieties\tVarieties ID\tAverage\tStdev\tSignificant Difference";
        $test = "Wilcox test\t";
       # $test2 = "Normalization\t";
        #$score = "Score\t";
        for (1..$#a) {
            $test .= "$wilcox{$a[$_]}\t";
            #$normalization = ($wilcox{$a[$_]} - $wilcox_min) / ($wilcox_max - $wilcox_min);
            #$test2 .= "$normalization\t";
            #$S = (1- ($wilcox{$a[$_]} / $wilcox_max))*100;
            #$S = (1- ($wilcox{$a[$_]} - $wilcox_min) / ($wilcox_max - $wilcox_min))*100;
            #$score .= "$S\t";
        }
        next;
    }

    if ($_ =~ /^REF/  ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $ref = "References allele\t". join ("\t",@a[1..$#a]);
        next;
    }
   if ($_ =~ /^ALT/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $alt = "Alternative allele\t". join ("\t",@a[1..$#a]);
        next;
    }
    if ($_ =~ /^INFO/ ) {
        $_ =~ s/\r*\n*//g;
        $stop = $1 if $_ =~ /(stop\w+)/;
        @a = split(/\t/,$_);
        $EFF = "SNP annotation\t". join ("\t",@a[1..$#a]);
        next;
    }
    if ($_ =~ /^Allele/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $All = "Allele Frequency\t". join ("\t",@a[1..$#a]);
        next;
    }
    $_ =~ s/\r*\n*//g;
    @a = split(/\t/,$_);
    $n = join ("\t",@a[1..$#a]);
    $lxk{$n} .= "$a[0] : $name{$a[0]}, ";
    $Line{$a[0]} = "$n!!\t$a[0]\t$name{$a[0]}";
    $ave{$n} .= "$name{$a[0]} " if $name{$a[0]} =~ /\d/;
    $lxk2{$n} ++;
}



open BOX,">tmp-Boxplot-$Phenotype-$opts{g}.txt" or die $!;
open TMP,">tmp-Hap-$opts{g}-$Phenotype.txt" or die $!;
print TMP "$ref\n$alt\n$All\n$EFF\n$test\n$POS\n";
#print TMP "$ref\n$alt\n$All\n$EFF\n$test\n$test2\n$score\n$POS\n";
$i =1 ;
foreach my $k (sort {$lxk2{$b} <=> $lxk2{$a} } keys %lxk2 ){
    $hap = Hap_ . $i;
    $lxk{$k} =~ s/, $//g;
    @ave2 = split(/ /,$ave{$k});
    if (@ave2 > 1) {
        $ave{$k} = &average(\@ave2);
        $std{$k} = &stdev(\@ave2);
    }
    if (@ave2 > 2) {
        for (@ave2) {
            print BOX "$hap\t$_\n";
        }
    }
    if (@ave2 < 2) {
        for (@ave2) {
            $lxk{$k} =~ /(.*) :/;
            $seqX{$1} = $1;
        }
    }
    $Line2{$k} = $hap;
    print TMP "$hap\t$k\t$lxk2{$k}\t$lxk{$k}\t$ave{$k}\t$std{$k}\t$diff{$hap}\n";
    $i ++;
}
if ($opts{n}) {
    open SEQ, ">tmp-seq-$opts{g}.txt" or die $!;
    foreach my $k (sort {$a<=>$b} keys %seq ){
        next if $seq{$k} =~ /N/;
        next if exists $seqX{$k};
        print SEQ ">$k\n$seq{$k}\n" ;
    }
}
open LINE,">Line_hap-$opts{g}-$Phenotype.txt" or die $!;
$POS =~ s/Average\tStdev\tSignificant Difference//g;
$POS =~ s/ /_/g;
print LINE "$POS\n";
foreach my $k (sort keys %Line ) {
    @a = split(/!!/,$Line{$k});
    $Line{$k} =~ s/!!//g;
    print LINE "$Line2{$a[0]}\t$Line{$k}\n";
}

##$gene_n
for ($i =1; $i <= $gene_n; $i++) {
  $new_gff = $opts{g} . '.' . $i;
  #system("grep  'gene'     tmp-gene_gff-$opts{g}-$Phenotype.txt >  tmp-gene_gff-$new_gff-$Phenotype.txt");
  system("grep  '$new_gff' tmp-gene_gff-$opts{g}-$Phenotype.txt > tmp-gene_gff-$new_gff-$Phenotype.txt");
}

foreach $k (keys %BP ){
    if ($wilcox{$k} >= $threshold) {
        delete($BP{$k});
    }
}
@hash_size = keys %BP;
$hash_size = @hash_size;


sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

print "Plotting Gene Structure of $opts{g} ......\nAnd plotting Histogram of $Phenotype ......\n";
open OUTR,">Plot_gene-$opts{g}-$Phenotype.R" or die $!;
print OUTR <<EOF;
x = read.table("$opts{p}")
y = x\$V2
Pheno = y[2:length(y)]
Phenoty =  as.numeric(as.character(Pheno))
Phenotype = Phenoty[!is.nan(Phenoty)]
P=shapiro.test(Phenotype)\$p.value
h = hist(Phenotype)
m=max(h\$counts) * 1.2
pdf(file="Plot_hist-$Phenotype.pdf", width=10, height=5)
par(mfrow=c(1,2))
h = hist(Phenotype, col='gray' ,xlab='Phenotype' , main="Histogram of $Phenotype", labels = TRUE, ylim = c(0,m))
rug(jitter(Phenotype))
xfit<-seq(min(Phenotype), max(Phenotype), length=100)
yfit<-dnorm(xfit, mean=mean(Phenotype), sd=sd(Phenotype))
yfit<-yfit*diff(h\$mids[1:2]) * length(Phenotype)
lines(xfit, yfit, col='darkblue', lwd=1)
lines(density(Phenotype)\$x, density(Phenotype)\$y * max(h\$counts)/max(h\$density), col='red',lwd=1)
text(min(Phenotype),m, expression(paste(italic("P =  "))))
text(min(Phenotype),m, labels=round(P,3), pos=4)

d=max(h\$density) * 1.2
hist(Phenotype, freq=FALSE, col='gray', xlab='Phenotype', main="Histogram of $Phenotype", ylim = c(0,d))
rug(jitter(Phenotype))
lines(density(Phenotype),col='red',lwd=1)
curve(dnorm(x, mean=mean(Phenotype), sd=sd(Phenotype)), add=TRUE, col="darkblue", lwd=1)
text(min(Phenotype),d, expression(paste(italic("P =  "))))
text(min(Phenotype),d, labels=P, pos=4)
dev.off()

mutation_plot<-function(start, stop, text="", drop=-0.15, col="red") {
   rect(start, -0.03, stop, drop, col=col, border=col)
   text(stop+180, drop-0.15, text, cex=0.5, col="black", pos=4, offset=-1)
}
genemodel_plot<-function(model, xaxis=TRUE, drop=0) {
  par(mar=c(1,1,3,1), cex=1)
  colnames(model) =c('chr','db','feature','start','end','tmp1','orientation','tmp3','gene')
  orientation = model\$orientation[1]
  start <- min(c(model\$start,model\$end))
  end <- max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model\$start=model\$start-tmp_min+1
  model\$end=model\$end-tmp_min+1
  tmp_max=max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model<-cbind(as.character(model[,3]), as.numeric(model[,4]), as.numeric(model[,5]) )
  model<-as.data.frame(model)
  colnames(model)<-c("feature", "start", "bpstop")
  model\$start<-as.numeric(as.character(model\$start));model\$bpstop<-as.numeric(as.character(model\$bpstop))
  length<-tmp_max-tmp_min
  if (orientation=="-") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2010, end+2500), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="mRNA") {
        arrows(x0=model\$newstop[i]+2000,y0=0.1 - drop,x1=model\$newstart[i]-500,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
        #arrows(x0=model\$newstop[i]+.2*length,y0=.1,x1=model\$newstart[i]-.2*length,y1=.1, lwd=1, col="dodgerblue4", length = 0.1)
      }
      if (type=="three_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
        #x<-c(model\$newstop[i], model\$newstop[i], model\$newstart[i]+.02*length, start, model\$newstart[i]+.02*length)
        #y<-c(0 - drop,.2 - drop,.2 - drop,.1 - drop,0 - drop)
        #polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (orientation=="+") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2500, end+2010), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="mRNA") {
        arrows(x0=model\$newstart[i]-2000,y0=0.1 - drop,x1=model\$newstop[i]+500,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
        #segments(x0=model\$newstart[i],y0=.1,x1=model\$newstop[i],y1=.1, lwd=1, col="dodgerblue4")
      }
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="three_prime_UTR") {
         rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
         #x<-c(model\$newstart[i], model\$newstart[i], model\$newstop[i]-.02*length, end, model\$newstop[i]-.02*length)
         #y<-c(0 - drop,.2 - drop,.2 - drop,.1 - drop,0 - drop)
         #polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (xaxis==T)   {
    Axis(side=3, labels=T, cex.axis=0.7)
  }
}
pdf(file="Plot_gene-$opts{g}-$Phenotype.pdf") 
EOF
for ($i = 1; $i <= $gene_n; $i++) {
  $new_gff = $opts{g} . '.' . $i;
  $new_gff = $opts{g} if (-z "tmp-gene_gff-$new_gff-$Phenotype.txt");
  #$j = ($i)*0.5;
  print OUTR "$new_gff <- read.table('tmp-gene_gff-$new_gff-$Phenotype.txt',stringsAsFactors = F, header = F,comment.char = \"#\",sep = '\\t')\n";
  print OUTR "genemodel_plot(model=$new_gff, xaxis=T, drop=0)\n";
  print OUTR "text($mystart + 1000, 0.5, \"$new_gff\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
  $w =1;$old=0;
  foreach $k (sort {$a<=>$b} keys %BP ){
    $w++ if ($k < $old + 0.1* ($myend - $mystart));
    $w=1 if ($k > $old + 0.15* ($myend - $mystart));
    $d= $w * (-0.2);
    #print OUTR "mutation_plot($k, $k, \"$BP{$k}\", $d, col=\"black\")\n";
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"#E69F00\")\n" if $INFO{$k} =~ /upstream/;
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"#0072B2\")\n" if $INFO{$k} =~ /UTR5/;
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"red\")\n" if $INFO{$k} =~ /exonic/;
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"#009E73\")\n" if $INFO{$k} =~ /intronic/;
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"#56B4E9\")\n" if $INFO{$k} =~ /UTR3/;
    print OUTR "mutation_plot($k, $k, \">$k\", $d, col=\"#CC79A7\")\n" if $INFO{$k} =~ /downstream/;
    $old = $k;
  }
}
print OUTR "dev.off()\n";
system("Rscript Plot_gene-$opts{g}-$Phenotype.R");
print "Plotted Gene Structure to file Plot_gene-$opts{g}-$Phenotype.pdf\n\n";
print "Plotted Histogram to file Plot_boxplot-$Phenotype-$opts{g}.pdf\n\n\n\n";

system("touch tmp-Diff-$opts{g}.txt");
print "Plotting Boxplot of $opts{g} - $Phenotype ......\n";
open OUTR,">Plot_hist-$Phenotype.R" or die $!;
print OUTR <<EOF;
if (! require("ggplot2")) install.packages("ggplot2")
df<-read.delim("tmp-Boxplot-$Phenotype-$opts{g}.txt",sep="\t",header=F)
names(df)<-c("hap","Phenotype")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (! require("agricolae")) BiocManager::install("agricolae")
tmp <- aov(df\$Phenotype ~ df\$hap)
res <- LSD.test(tmp, 'df\$hap', p.adj = 'bonferroni',alpha =0.001)
diff <- res\$groups
write.table(diff, file = "tmp-Diff-$opts{g}.txt", sep = "\t", quote = FALSE, col.names = FALSE)
diff[,'hap']<-factor(rownames(diff))
df\$hap<-factor(df\$hap,levels=unique(df\$hap))
if (nrow(diff)-3 > 3) {
    mywidth = (nrow(diff)-3)
} else {mywidth = 4 }
if (! require("ggbeeswarm")) install.packages("ggbeeswarm")
pdf(file="Plot_boxplot-$Phenotype-$opts{g}.pdf", width=mywidth, height=3)
ggplot(data = df,aes(x = hap, y = Phenotype, group= hap)) + geom_boxplot(fill=rainbow(length(levels(factor(df\$hap)))), outlier.colour = NA, alpha=.1, color=rainbow(length(levels(factor(df\$hap)))))  + geom_beeswarm(size = 0.05, cex = 0.6, alpha=0.5) + labs(title="Boxplot of $opts{g} Gene's Haplotypes", x="Haplotypes", y="Value of Phenotype") + geom_text(data=diff,aes(y=-0.1,label=groups)) + theme(legend.position="none")
ggplot(data = df,aes(x = hap, y = Phenotype, group= hap)) + geom_boxplot(fill=rainbow(length(levels(factor(df\$hap)))), outlier.colour = NA, alpha=.1, color=rainbow(length(levels(factor(df\$hap)))))  +  geom_point(size=0.1,position = position_jitter(width=0.15, height = 0), alpha=.5,shape=20) +labs(title="Boxplot of $opts{g} Gene's Haplotypes", x="Haplotypes", y="Value of Phenotype") + geom_text(data=diff,aes(y=-0.1,label=groups))
ggplot(data = df,aes(x = hap, y = Phenotype, group= hap)) + geom_boxplot(fill=rainbow(length(levels(factor(df\$hap)))), outlier.colour = NA) +  geom_point(size=0.1,position = position_jitter(width=0.15, height = 0), alpha=.5,shape=20) +labs(title="Boxplot of $opts{g} Gene's Haplotypes", x="Haplotypes", y="Value of Phenotype") + geom_text(data=diff,aes(y=-0.1,label=groups))
#boxplot(Phenotype~hap,data=df, col=rainbow(length(levels(factor(df\$hap)))), xlab = "Haplotypes", ylab = "Value of Phenotype", main = "Boxplot of $opts{g} Gene's Haplotypes")
dev.off()
EOF
system("Rscript Plot_hist-$Phenotype.R");
print "Plotted Boxplot to file Plot_boxplot-$Phenotype-$opts{g}.pdf\n\n\n\n";

if ($opts{n}) {
print "Plotting Haplotype Network of $opts{g} ......\n";
open OUTR,">Plot_Network-$opts{g}.R" or die $!;
print OUTR <<EOF;
if (! require("sf")) install.packages('sf', repos='https://cran.rstudio.com/')
if (! require("pegas")) install.packages('pegas', repos='https://cran.rstudio.com/')
fa <- read.FASTA("tmp-seq-$opts{g}.txt")
haps <- haplotype(fa)
(network <- haploNet(haps))
pdf(file="Plot-haploNet-$opts{g}.pdf")
plot(network, size = log2(attr(network, "freq")),cex = 0.8,show.mutation=2)
dev.off()
cat("\n\nThe sequence of each haplotype:")
as.data.frame(diffHaplo(haps,1:nrow(haps)))
EOF
system("Rscript Plot_Network-$opts{g}.R");
print "Plotted Haplotype Network to file Plot-haploNet-$opts{g}.pdf\n\n\n";
}

open DIFF,"tmp-Diff-$opts{g}.txt" or die $!;
while (<DIFF>) {
    $_ =~ s/\n*|\r*//g;
    @D = split(/\t/,$_);
    $diff{$D[0]}=$D[2];
}
if ($stop) {
    open HAP,">Hap_result-$opts{g}-$Phenotype-$stop.txt" or die $!;
    print "Print Haplotypes results to file: Hap_result-$opts{g}-$Phenotype-$stop.txt\n\n";
} else {
    open HAP,">Hap_result-$opts{g}-$Phenotype.txt" or die $!;
    print "Print Haplotypes results to file: Hap_result-$opts{g}-$Phenotype.txt\n\n";
}
open HAP1,"tmp-Hap-$opts{g}-$Phenotype.txt" or die $!;
while (<HAP1>) {
    $_ =~ s/\s+$//g;
    ($head,$info) = split(/\t/,$_,2);
    print HAP "$head\t$info\t$diff{$head}\n";
}
if  ($opts{k} == 0) {
    system("rm -rf tmp*");
    system("rm -rf *.R");
}
system("rm -rf Rplots.pdf");
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
exit(1);