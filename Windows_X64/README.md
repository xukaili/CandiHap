## CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_win.gif" width="100" height="100">

# For Windows
The installation package **integrates all the necessary modules** for **running independently**, meaning no more software installation required.</br>

## Download:
**Software**:      [BioCode: CandiHap_win64_V1.3.0.zip](https://ngdc.cncb.ac.cn/biocode/tools/7080/releases/V1.3.0)

**Test Data**:     [GitHub](https://github.com/xukaili/CandiHap/raw/master/test_data.zip)                [Google Drive](https://drive.google.com/file/d/1L2FTr1ktxU5Jgkuk4QXIJIMHHSzri9l4/view?usp=sharing)                [Baidu Pan](https://pan.baidu.com/s/1X4Tu1ha6d1caC518CBSHVA):   access_code: **i2sr**</br></br>

## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang* (王兴春). CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## Getting started
To annotate the vcf by ANNOVAR:</br>
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
