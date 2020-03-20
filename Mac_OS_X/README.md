## CandiHap: a haplotype analysis toolkit for natural variation study.

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/logo_mac.gif" width="100" height="100">

# For Mac OS X
First of all, please install the [**R software environment**](https://www.r-project.org), and three packages.</br>

## Download:
**Software**:      [GitHub](https://github.com/xukaili/CandiHap/raw/master/Mac_OS_X/CandiHap-1.0.1.dmg)                [Google Drive](https://drive.google.com/file/d/1BluXrsNxBgQksamGH-prb-BOWv9wcotr/view?usp=sharing)                [Baidu Pan](https://pan.baidu.com/s/1Cd4luIxElLHRkiay9isFmQ):   access_code: **xx7n**</br>

**Test Data**:     [GitHub](https://github.com/xukaili/CandiHap/raw/master/test_data.zip)                [Google Drive](https://drive.google.com/file/d/1L2FTr1ktxU5Jgkuk4QXIJIMHHSzri9l4/view?usp=sharing)                [Baidu Pan](https://pan.baidu.com/s/1X4Tu1ha6d1caC518CBSHVA):   access_code: **i2sr**</br></br>

## License
__`Academic users`__ may download and use the application free of charge according to the accompanying license.</br>
__`Commercial users`__ must obtain a commercial license from Xukai Li.</br>
**If you have used the program to obtain results, please cite the following paper:**</br>
> Xukai Li☯* (李旭凯), Zhiyong Shi☯ (石志勇), Qianru Qie (郄倩茹), Jianhua Gao (高建华), Yiwei Jiang (姜亦巍), Yuanhuai Han (韩渊怀) & Xingchun Wang* (王兴春). CandiHap: a haplotype analysis toolkit for natural variation study. bioRxiv 2020.02.27.967539. doi: https://doi.org/10.1101/2020.02.27.967539</br>
> （☯ Equal contributors; * Correspondence）</br>
</br>

## To Install __`R`__ for Mac OS X and packages
      1. Open an internet browser and go to link: https://www.r-project.org</br>
      2. Click the '__`download R`__' link in the middle of the page under '__`Getting Started`__'.</br>
      3. Select a CRAN location (a __`mirror site`__) and click the corresponding link.</br>
      4. Click on the '__`Download R for (Mac) OS X`__' link at the top of the page.</br>
      5. Click on Download '__`R-3.5.0.pkg`__' (or a newer version).</br>
      6. Install R and leave all default settings in the installation options.</br>
      7. Open R and install three packages by command: </br>
          __`install.packages(c("ggplot2", "agricolae", "pegas"))`__</br>
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

## To Install __`CandiHap.app`__ for Mac OS X
If you attempt to open CandiHap.app and macOS stops you from doing so, that doesn't necessarily mean there is something wrong with the app. But it will indicate that the app is from an '__`unidentified developer`__'.</br>
You can open the app and override the block. Here's how:</br>
      1. Open '__`System Preferences`__'.</br>
      2. Go to '__`Security & Privacy`__' and select the '__`General`__' tab.</br>
      3. Click on the button '__`Open Anyway`__'.</br>
      4. You’ll be asked one more time, and clicking '__`Open`__'.</br>

<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/Mac_ReadMe_1.png"  width="458" height="180">
<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/Mac_ReadMe_2.png"  width="652" height="460">
<img src="https://github.com/xukaili/CandiHap/blob/master/Figures/Mac_ReadMe_3.png"  width="404" height="200">
</br>

## Contact information
In the future, **CandiHap** will be regularly updated, and extended to fulfill more functions with more user-friendly options.</br>
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com </br>
