## CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_win.gif" width="100" height="100">

# For Windows
The installation package **integrates all the necessary modules** for **running independently**, meaning no more software installation required.</br>

## Download:
**Software**:      [BioCode: CandiHap_win64_V1.3.0.zip](https://ngdc.cncb.ac.cn/biocode/tools/7080/releases/V1.3.0)


## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Jianhua Gao (高建华), Xingchun Wang (王兴春), kai Guo* (郭凯). CandiHap: a haplotype analysis toolkit for natural variation study. Molecular Breeding, 2023, 43:21. doi: https://doi.org/10.1007/s11032-023-01366-4</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## Getting started
To annotate the vcf by ANNOVAR (Version: 2019-10-24 00:05:27 -0400):</br>
```sh
     gffread  test.gff   -T -o test.gtf
     gtfToGenePred -genePredExt test.gtf  si_refGene.txt
     retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
     table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test  --buildver  si  --protocol refGene --operation g -remove
```
</br>

## Contact information
In the future, **CandiHap** will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
