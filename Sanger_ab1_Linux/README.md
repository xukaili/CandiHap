# CandiHap: a haplotype analysis toolkit for natural variation study.
# Haplotype analysis in Sanger __`.ab1`__ files (Linux system)
First of all, please install GATK (GenomeAnalysisTK.jar), Picard (picard.jar), bwa, samtools, bcftools, bgzip, java, perl and R (with sangerseqR)</br></br>

## Download:
**Software**:   [GitHub](https://github.com/xukaili/CandiHap/raw/master/Sanger_ab1_Linux/sanger_CandiHap-1.0.1.zip)             [Google Drive](https://drive.google.com/file/d/1QgVMbSYx27_j_OYzu5OWWRshSehsoZkJ/view?usp=sharing)             [Baidu Pan](https://pan.baidu.com/s/1ShMMMNogsJNdx3GOyMIj_w) access_code: **zuqy**</br></br>

**Test Data**: [GitHub](https://github.com/xukaili/CandiHap/raw/master/Sanger_ab1_Linux/sanger_teat_data.zip)             [Google Drive](https://drive.google.com/file/d/1ZiFGuFG01b4r_zsnbhqpaIZAZ-yffUCC/view?usp=sharing)             [Baidu Pan](https://pan.baidu.com/s/1Y-Ohg-Q8AiXLDJDoQxkMzA) access_code: **cp6r**</br></br>

## License
Academic users may download and use the application free of charge according to the accompanying license. Commercial users must obtain a commercial license from Xukai Li. If you have used the program to obtain results, please cite the following paper:</br>

> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang* (王兴春). CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## Dependencies
install __`sangerseqR`__ packages in R:</br>
```
if (! requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (! require("sangerseqR")) BiocManager::install("sangerseqR")
```

## Getting started
Put __`sanger_CandiHap.sh`__, __`Gene_VCF2haplotypes.pl`__, __`ab1-fastq.pl`__ and and all __`.ab1`__ files in a same dir, then run:</br>
```
     sh  sanger_CandiHap.sh  Gene_ref.fa
e.g. sh  sanger_CandiHap.sh  PHYC.txt
```

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/Sanger_Figure.png">

**Fig. 1 | Overview of the sanger_CandiHap process. a,** General scheme of the process from Sanger ab1 files. **b,** PeakTrace of ab1 images of three main genotypes. **c,** The statistics of haplotypes. </br>

## Contact information
In the future, CandiHap will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
