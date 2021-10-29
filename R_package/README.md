```sh
# install the CandiHap package:
if (! require("devtools ")) install.packages("devtools", dependencies = T)
install_github("guokai8/CandiHap")
if (! require("eoffice")) install.packages("eoffice")

# download the test data from: https://bigd.big.ac.cn/biocode/tools/7080/releases/v2.0.15 

# run CandiHap
library(CandiHap)
if (! require("eoffice")) install.packages("eoffice")

gff <- importGFF("test.gff",format="gff3")
gr<-preGranges(gff,gene="Si9g49990",cds = T)

#hmp <- read_vcf('test.vcf')           # or if you use vcf file
hmp <- read_hmp("haplotypes.hmp")
ovl <- findover(gr,hmp)

pheno <- read_pheno("Phenotype.txt",sep="\t")
hap <- snp2hap(pheno,ovl, hapname="Hap")

## want to extract results
Hap_results <- results(hap, gene="Si9g49990")
write.table(Hap_results, file = "Hap_results.txt", sep = "\t", quote = FALSE, col.names = T, row.names = F)

## Plot Gene LDheatmap
if (! require("LDheatmap")) install.packages("LDheatmap")
if (! require("genetics")) install.packages("genetics")
LDsnp <- Hap_results[, !names(Hap_results) %in% c("sample","Info")]
LDsnp = LDsnp[,-ncol(LDsnp)]
NAs <- LDsnp == "N/N"
is.na(LDsnp)[NAs] <- TRUE
LDsnp[NAs] <- NA
LDsnp = LDsnp[-1,]
Dist = as.numeric(LDsnp[1,])
LDsnp = LDsnp[-c(1:6),]
num<-ncol(LDsnp) 
for(i in 1:num){
    LDsnp[,i]<-as.genotype(LDsnp[,i]) 
}
rgb.palette <- colorRampPalette(rev(c("#1F77B4", "#AEC7E8", "#FF7F0E")), space = "rgb")
LDheatmap(LDsnp, genetic.distances = Dist, flip = TRUE, color=rgb.palette(40))
topptx(filename ="LDheatmap.pptx")

# gene figures
snplot(hap,gene="Si9g49990",side=F)
snplot(hap,gene="Si9g49990",side=F,random = F)
topptx(filename ="gene_figure.pptx")
snplot(hap,gene="Si9g49990",side=T,random = F,hapname="haplotype3")
snplot(hap,gene="Si9g49990",side=F,random = F,hapname="haplotype15",mutateOnly=TRUE)

## boxplot
snboxplot(hap,gene="Si9g49990",feature = "Height")
topptx(filename ="snboxplot.pptx")

hapnet(hap,gene="Si9g49990",feature = "Height")
topptx(filename ="hapnet.pptx")


# plot gene track with snp
library(CandiHap)
dat <- read_data("track-Phenotype.gwas.txt",sep="\t")
# notice that the dat should have the chromosome name in the first column, position in the second column and the values in the following columns 

## id is the gene name you want to display, in gff3 file should be 'Parent'
#show snp locate in gene only
snptrack(gff,dat=dat,id="Parent",geneOnly=F, color='r2', low='green',high='red', point.size = 0.5, chr = 9,  region = c(54520000, 54620000))
topptx(filename ="snptrack.pptx")

## show some genes
snptrack(gff,dat=dat,id="Parent",gene='Si9g49990',color='r2',upstream=2000)
topptx(filename ="gene_snptrack.pptx")

     
```
