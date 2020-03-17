# CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_mac.gif" width="100" height="100">

# For Linux system (command lines)
First of all, please install the R software environment (https://www.r-project.org), and three packages.</br>

## Download:
Software: https://github.com/xukaili/CandiHap/raw/master/Linux/CandiHap-1.0.1.zip</br>
      Or: https://pan.baidu.com/s/1qGYtzhH7YAwOiro9oqbKSg            access_code: e4mv</br></br>
Test Data: https://github.com/xukaili/CandiHap/raw/master/test_data.zip</br>
       Or: https://pan.baidu.com/s/1X4Tu1ha6d1caC518CBSHVA           access_code: i2sr</br></br>

## License
Academic users may download and use the application free of charge according to the accompanying license. Commercial users must obtain a commercial license from Xukai Li. If you have used the program to obtain results, please cite the following paper:</br>

> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang (王兴春)*. CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## To Install __`R`__ for Linux system
```
   1. Open an internet browser and go to link: https://www.r-project.org
   2. Click the "download R" link in the middle of the page under "Getting Started."
   3. Select a CRAN location (a mirror site) and click the corresponding link.
   4. Click on the "Download R for Linux" link at the top of the page.
   5. Click on Download R-3.5.0 (or a newer version).
   6. Install R. Leave all default settings in the installation options.
   7. Open R and install three packages by command： 
      install.packages(c("ggplot2", "agricolae" , "pegas")
```
</br>

## Getting started
There are mainly three steps included in the CandiHap analytical through command lines, and the test data files can freely download at __`test_data.zip`__.</br>
Put __`vcf2hmp.pl`__  test.gff, test.vcf, and genome.fa files in a same dir, then run:</br>
```
     # 1. To annotate the vcf by ANNOVAR: 
     gffread  test.gff   -T -o test.gtf
     gtfToGenePred -genePredExt test.gtf  si_refGene.txt
     retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
     table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test --buildver  si --protocol refGene --operation g -remove

     # 2. To convert the txt result of annovar to hapmap format, 0.1 means the minor allele frequency (MAF):
     perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt  0.1
```
</br>

Put __`CandiHap.pl`__ and Phenotype.txt, Your.hmp, genome.gff files in a same dir, then run:</br>
```
     # 3. To run CandiHaplotypes
     perl  CandiHap.pl  ./Your.hmp  ./Phenotype.txt  ./genome.gff  Your_gene_ID
e.g. perl  CandiHap.pl  ./haplotypes.hmp   ./Phenotype.txt  ./test.gff  Si9g49990
```
</br>

By the way, if you want do All gene in LD region of a position, then run:</br>
```
     perl  GWAS_LD2haplotypes.pl   ./genome.gff  ./ann.hmp  ./Phenotype.txt  50kb  Chr:position
e.g. perl  GWAS_LD2haplotypes.pl   ./test.gff  ./haplotypes.hmp   ./Phenotype.txt   50kb  9:54583294
```
</br>

## Contact information
In the future, CandiHap will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
