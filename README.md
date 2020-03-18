# CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_mac.gif" width="100" height="100"> <img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_win.gif" width="100" height="100">

__`CandiHap`__ is a user-friendly local software, that can fast preselect candidate causal SNPs from __`Sanger or next-generation sequencing data`__, and report results in table and exquisite vector-graphs within a minute. Investigators can use CandiHap to specify a gene or linkage sites based on GWAS and explore favourable haplotypes of candidate genes for target traits. CandiHap can be run on computers with __`Windows`__, __`Mac OS X`__, or __`Linux`__ platforms in graphical user interface or command lines, and applied to __`any species`__ of plant, animal and microbial. The analysis of CandiHap can do as the followings: </br></br>
    **1).** Convert the VCF file to the hapmap format for CandiHap (__`vcf2hmp`__);</br>
    **2).** Haplotype analysis for a gene (__`CandiHap`__);</br>
    **3).** Haplotype analysis for all genes in the LD regions of a significant SNP one by one (__`GWAS_LD2haplotypes`__);</br>
    **4).** Haplotype analysis for Sanger sequencing data of population variation (__`sanger_CandiHap.sh`__).</br>

## Download:
**Download All Files**:      [GitHub](https://github.com/xukaili/CandiHap/archive/master.zip)                [Google Drive](https://drive.google.com/drive/folders/1rkm__0jlPrWvs8Cy_T46d3sKBf6UIVhQ?usp=sharing)                [Baidu Pan](https://pan.baidu.com/s/1migiq6x5v5CgYlBntY8aYw):   access_code: **4d4q**</br></br>

## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang* (王兴春). CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## Dependencies
__`perl 5`__, __`R ≥ 3.2`__ (with ggplot2, agricolae, pegas and sangerseqR), and __`electron`__. </br></br>

## Figures
![CandiHap](Figures/CandiHap.png)
**Fig. 1 | Overview of the CandiHap process. a,** A GWAS result. **b,** General scheme of the process. **c,** The histogram of phenotype. **d,** The statistics of haplotypes and significant differences haplotypes are highlighted by color boxes. **e,** Gene structure and SNPs of a critical gene. **f,** Boxplot of a critical gene’s haplotypes. </br>

![HaploNet](Figures/HaploNet.png)
**Fig. 2 | Haplotype network analysis for *Si9g49990*. a,** The difference of haplotypes. **b,** Haplotype network. Note: only the SNPs and haplotypes found in ≥2 accessions were used to construct the haplotype network. The value of circle size had converted into log2. </br>

![Rice-2018_Nat_Commun_9_735](Figures/Rice-2018_Nat_Commun_9_735.png)
**Fig. 3 | Haplotype analysis of the *ARE1* gene in rice compared with the results by [Wang *et al*. 2018, *Nat*. *Commun*. 9, 735](https://www.nature.com/articles/s41467-017-02781-w/figures/5). a,** Gene structure and SNPs of *ARE1*. **b,** Major haplotypes of SNPs in the *ARE1* coding region of 2747 rice varieties. **c,** The haplotype results of *ARE1* coding region of 3023 rice varieties using CandiHap (SNPs data were downloaded from [RFGB](http://www.rmbreeding.cn)). Major SNP haplotypes and casual variations in the encoded amino acid residues are shown. The five more SNPs is due to the fact that there are 276 more rice varieties used in our study (highlighted by blue boxes), and two errors highlighted by red boxes. </br></br>

## Contact information
In the future, **CandiHap** will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
