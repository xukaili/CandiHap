## CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_mac.gif" width="100" height="100">

# For Linux system (command lines)
First of all, please install the [**R software environment**](https://www.r-project.org), and three packages.</br>

## Download:
**Software and Test Data**:      [BioCode: https://ngdc.cncb.ac.cn/biocode/tools/7080/releases/V1.3.0](https://ngdc.cncb.ac.cn/biocode/tools/7080/releases/V1.3.0)

## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Jianhua Gao (高建华), Xingchun Wang (王兴春), kai Guo* (郭凯). CandiHap: a haplotype analysis toolkit for natural variation study. Molecular Breeding, 2023, 43:21. doi: https://doi.org/10.1007/s11032-023-01366-4</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## To Install __`R`__ for Linux and packages
      1. Open an internet browser and go to link: https://www.r-project.org</br>
      2. Click the '__`download R`__' link in the middle of the page under '__`Getting Started`__'.</br>
      3. Select a CRAN location (a __`mirror site`__) and click the corresponding link.</br>
      4. Click on the '__`Download R for Linux`__' link at the top of the page.</br>
      5. Click on Download '__`R-3.5.0`__' (or a newer version).</br>
      6. Install R and leave all default settings in the installation options.</br>
      7. Open R and install three packages by command: </br>
          __`install.packages(c("ggplot2", "agricolae", "pegas"))`__</br>
</br>

## Getting started
There are mainly three steps included in the CandiHap analytical through command lines, and the test data files can freely download at __`test_data.zip`__.</br>
Put __`vcf2hmp.pl`__  test.gff, test.vcf, and genome.fa files in a same dir, then run:</br>
```sh
     # 1. To annotate the vcf by ANNOVAR (Version: 2019-10-24 00:05:27 -0400): 
     gffread  test.gff   -T -o test.gtf
     gtfToGenePred -genePredExt test.gtf  si_refGene.txt
     retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
     table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test  --buildver  si  --protocol refGene --operation g -remove

     # 2. To convert the txt result of annovar to hapmap format:
     perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt
```
</br>

Put __`CandiHap.pl`__ and Phenotype.txt, Your.hmp, genome.gff files in a same dir, then run:</br>
```sh
     # 3. To run CandiHaplotypes
     perl  CandiHap.pl  -m Your.hmp  -f Genome.gff  -p Phenotype.txt  -g Your_gene_ID
e.g. perl  CandiHap.pl  -m haplotypes.hmp  -f test.gff  -p Phenotype.txt  -g Si9g49990
     perl  CandiHap.pl  -m haplotypes.hmp  -f test.gff -p Phenotype.txt  -g Si9g49990 -s 0.5 -u 2000 -d 500 -l 1 -n Structure.txt
```
</br>

If you want do analysis __`All gene in LD region of a position`__, please run:</br>
```sh
     perl  GWAS_LD2haplotypes.pl  -f genome.gff  -m ann.hmp  -p Phenotype.txt   -l LDkb  -c Chr:position
e.g. perl  GWAS_LD2haplotypes.pl  -f test.gff  -m haplotypes.hmp  -p Phenotype.txt  -l 50kb  -c 9:54583294
```
</br>

## Contact information
In the future, **CandiHap** will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@sxau.edu.cn </br>
