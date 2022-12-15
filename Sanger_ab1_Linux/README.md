## CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_mac.gif" width="100" height="100">

# Haplotype analysis for Sanger __`.ab1`__ files on Linux system
First of all, please install __`GATK`__ (GenomeAnalysisTK.jar), __`Picard`__ (picard.jar), __`bwa`__, __`samtools`__, __`bcftools`__, __`bgzip`__, __`java`__ and __`R`__ (with sangerseqR).</br></br>

## Download:
**Software**:       [BioCode](https://ngdc.cncb.ac.cn/biocode/tools/BT007080/releases/V1.3.0)                [Google Drive](https://drive.google.com/file/d/1QgVMbSYx27_j_OYzu5OWWRshSehsoZkJ/view?usp=sharing)

**Test Data**:     [BioCode](https://ngdc.cncb.ac.cn/biocode/tools/BT007080/releases/V1.3.0)                [Google Drive](https://drive.google.com/file/d/1ZiFGuFG01b4r_zsnbhqpaIZAZ-yffUCC/view?usp=sharing)
</br></br>

## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang* (王兴春). CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## To Install __`R`__ and __`sangerseqR`__ package
      1. Open an internet browser and go to link: https://www.r-project.org</br>
      2. Click the '__`download R`__' link in the middle of the page under '__`Getting Started`__'.</br>
      3. Select a CRAN location (a __`mirror site`__) and click the corresponding link.</br>
      4. Click on the '__`Download R for Linux`__' link at the top of the page.</br>
      5. Click on Download '__`R-3.5.0`__' (or a newer version).</br>
      6. Install R and leave all default settings in the installation options.</br>
      7. Open R and install the package by command: </br>
          __`if (! requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")`__</br>
          __`if (! require("sangerseqR")) BiocManager::install("sangerseqR")`__</br>
</br>

## Getting started
Put __`sanger_CandiHap.sh`__, __`Gene_VCF2haplotypes.pl`__, __`ab1-fastq.pl`__ and and all __`.ab1`__ files in a same dir, then run:</br>
```sh
     sh  sanger_CandiHap.sh  Gene_ref.fa
e.g. sh  sanger_CandiHap.sh  PHYC.txt
```

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/Sanger_Figure.png">

**Fig. 1 | Overview of the sanger_CandiHap process.** __`a,`__ General scheme of the process from Sanger **.ab1** files. __`b,`__ PeakTrace of ab1 images of three main genotypes. __`c,`__ The statistics of haplotypes. </br></br>

## Contact information
In the future, **CandiHap** will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
