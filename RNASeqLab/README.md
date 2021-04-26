---
title: 'RNASeq Lab: Parts I-IV'
author: "Sid"
date: "4/25/2021"
output:
  word_document: default
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=12, fig.height=8) 
```

# Part I: Alignment 

**Goals:** 

  1) Answer all the Questions
  2) Submit shell scripts that you used
  3) Run star alignment (hopper:/opt/STAR-master/bin/Linux_x86_64_static/)
  
### 1) Answer all the Questions

**Question 1:** How do you save the alignment summary? 
**Answer:** Pipe into a new file in the command line 

### 2) Submit shell scripts that you used 

In order to align the RNA seq reads, I used hisat2 and the following **align.sh** shell script. 

```{r eval=FALSE}
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188044_chrX_1.fastq.gz -2 samples/ERR188044_chrX_2.fastq.gz -S hisat_out/ERR188044_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188104_chrX_1.fastq.gz -2 samples/ERR188104_chrX_2.fastq.gz -S hisat_out/ERR188104_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188234_chrX_1.fastq.gz -2 samples/ERR188234_chrX_2.fastq.gz -S hisat_out/ERR188234_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188245_chrX_1.fastq.gz -2 samples/ERR188245_chrX_2.fastq.gz -S hisat_out/ERR188245_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188257_chrX_1.fastq.gz -2 samples/ERR188257_chrX_2.fastq.gz -S hisat_out/ERR188257_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188273_chrX_1.fastq.gz -2 samples/ERR188273_chrX_2.fastq.gz -S hisat_out/ERR188273_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188337_chrX_1.fastq.gz -2 samples/ERR188337_chrX_2.fastq.gz -S hisat_out/ERR188337_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188383_chrX_1.fastq.gz -2 samples/ERR188383_chrX_2.fastq.gz -S hisat_out/ERR188383_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188401_chrX_1.fastq.gz -2 samples/ERR188401_chrX_2.fastq.gz -S hisat_out/ERR188401_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188428_chrX_1.fastq.gz -2 samples/ERR188428_chrX_2.fastq.gz -S hisat_out/ERR188428_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR188454_chrX_1.fastq.gz -2 samples/ERR188454_chrX_2.fastq.gz -S hisat_out/ERR188454_chrX.sam
hisat2 -p 8 --dta -x indexes/chrX_tran -1 samples/ERR204916_chrX_1.fastq.gz -2 samples/ERR204916_chrX_2.fastq.gz -S hisat_out/ERR204916_chrX.sam
```

Unfortunately, hisat2 output each of the aligned files in SAM format, while future processing required BAM format. Therefore, I used the provided **sam_to_bam.sh** shell script to run the conversion. 

```{r eval=FALSE}
samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam
samtools sort -@ 8 -o ERR188104_chrX.bam ERR188104_chrX.sam
samtools sort -@ 8 -o ERR188234_chrX.bam ERR188234_chrX.sam
samtools sort -@ 8 -o ERR188245_chrX.bam ERR188245_chrX.sam
samtools sort -@ 8 -o ERR188257_chrX.bam ERR188257_chrX.sam
samtools sort -@ 8 -o ERR188273_chrX.bam ERR188273_chrX.sam
samtools sort -@ 8 -o ERR188337_chrX.bam ERR188337_chrX.sam
samtools sort -@ 8 -o ERR188383_chrX.bam ERR188383_chrX.sam
samtools sort -@ 8 -o ERR188401_chrX.bam ERR188401_chrX.sam
samtools sort -@ 8 -o ERR188428_chrX.bam ERR188428_chrX.sam
samtools sort -@ 8 -o ERR188454_chrX.bam ERR188454_chrX.sam
samtools sort -@ 8 -o ERR204916_chrX.bam ERR204916_chrX.sam
```

# Part II: Expression Calculation 

**Goals:** 

  1) Answer for the Questions
  2) Submit shell scripts that you used
  3) Run the steps using STAR alignment results

### 1) Answer for the Questions

No questions were provided for these slides. 

### 2) Submit shell scripts that you used

First I ran **stringtie** to estimate the expression levels from the BAM files. Dr. Ahn provided a script, **stringtie.sh** which I have provided below. 

```{r eval=FALSE}
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188044_chrX.gtf -l ERR188044 hisat_out/ERR188044_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188104_chrX.gtf -l ERR188104 hisat_out/ERR188104_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188234_chrX.gtf -l ERR188234 hisat_out/ERR188234_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188245_chrX.gtf -l ERR188245 hisat_out/ERR188245_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188257_chrX.gtf -l ERR188257 hisat_out/ERR188257_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188273_chrX.gtf -l ERR188273 hisat_out/ERR188273_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188337_chrX.gtf -l ERR188337 hisat_out/ERR188337_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188383_chrX.gtf -l ERR188383 hisat_out/ERR188383_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188401_chrX.gtf -l ERR188401 hisat_out/ERR188401_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188428_chrX.gtf -l ERR188428 hisat_out/ERR188428_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR188454_chrX.gtf -l ERR188454 hisat_out/ERR188454_chrX.bam
stringtie -p 8 -G genes/chrX.gtf -o stringtie_out/ERR204916_chrX.gtf -l ERR204916 hisat_out/ERR204916_chrX.bam
```

Finally, I merged each individual stringtie run into one master file using **stringtie_merge.sh**, which required just one linux command execution. 

```{r eval=FALSE}
stringtie --merge -p 8 -G ../genes/chrX.gtf -o stringtie_merged.gtf ../mergelist.txt
```

At this point, we wanted to determine the number of reads per gene, which we approached using two different techniques. The first technique we used was **HTSeq**. No script was provided to run this program so I had to create my own, **htseq.sh**, which I have listed below. 

```{r eval=FALSE}
htseq-count -f bam hisat_out/ERR188044_chrX.bam genes/chrX.gtf > htseq_out/ERR188044_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188104_chrX.bam genes/chrX.gtf > htseq_out/ERR188104_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188234_chrX.bam genes/chrX.gtf > htseq_out/ERR188234_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188245_chrX.bam genes/chrX.gtf > htseq_out/ERR188245_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188257_chrX.bam genes/chrX.gtf > htseq_out/ERR188257_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188273_chrX.bam genes/chrX.gtf > htseq_out/ERR188273_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188337_chrX.bam genes/chrX.gtf > htseq_out/ERR188337_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188383_chrX.bam genes/chrX.gtf > htseq_out/ERR188383_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188401_chrX.bam genes/chrX.gtf > htseq_out/ERR188401_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188428_chrX.bam genes/chrX.gtf > htseq_out/ERR188428_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR188454_chrX.bam genes/chrX.gtf > htseq_out/ERR188454_chrX_htseq_count.txt
htseq-count -f bam hisat_out/ERR204916_chrX.bam genes/chrX.gtf > htseq_out/ERR204916_chrX_htseq_count.txt
```

Now that I had gene count data for each of the reads, I applied the second approach, **featureCounts**, to solve this problem. I ran the command given below, which produced a summary file. 

```{r eval=FALSE}
featureCounts -p -T 10 -a genes/chrX.gtf -o featurecounts_out/all_featureCountx.txt hisat_out/ERR*.bam
```

```{r eval=FALSE}
Status	hisat_out/ERR188044_chrX.bam	hisat_out/ERR188104_chrX.bam	hisat_out/ERR188234_chrX.bam	hisat_out/ERR188245_chrX.bam	hisat_out/ERR188257_chrX.bam	hisat_out/ERR188273_chrX.bam	hisat_out/ERR188337_chrX.bam	hisat_out/ERR188383_chrX.bam	hisat_out/ERR188401_chrX.bam	hisat_out/ERR188428_chrX.bam	hisat_out/ERR188454_chrX.bam	hisat_out/ERR204916_chrX.bam
Assigned	651625	646894	820514	429846	477388	307547	595829	477717	651752	429821	501095	545909
Unassigned_Unmapped	5174	4237	5910	2765	3602	2173	4493	2501	4150	2721	3026	3969
Unassigned_Read_Type	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_Singleton	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_MappingQuality	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_Chimera	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_FragmentLength	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_Duplicate	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_MultiMapping	48140	45750	52632	26925	35745	15246	52393	28739	47310	25055	43431	37664
Unassigned_Secondary	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_NonSplit	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_NoFeatures	195821	201832	228150	116791	127190	74297	188470	135194	172421	105864	162209	158032
Unassigned_Overlapping_Length	0	0	0	0	0	0	0	0	0	0	0	0
Unassigned_Ambiguity	450052	421469	493625	296407	339642	189903	470959	337776	470798	295233	360768	375652
```


# Part III: Expression Analysis with R 

**Goals:**

  1) Answer for the Questions (from slides)
  2) Complete the "Data visualization" section (Procedure 17-21) from <https://www.nature.com/articles/nprot.2016.095#an1> 
  3) Submit your R script and plots

### 1) Answer for the Questions (from slides)

First, we will load the necessary R libraries required to complete our first goal. Each of these libraries must be pre-installed locally in order for successful compilation. We will then set the working directory to the location in which the RNAseq files are stored in. This is the same director in which all output files will be placed. 

```{r, message=FALSE}
# Load necessary libraries (must be pre-installed)
library(ballgown)
library(RSkittleBrewer) 
library(genefilter)
library(dplyr)
library(devtools)
# set working directory appropriately 
setwd("~/Documents/2020-21/BCB 5250/Labs/RNAseq") 
```

We can now begin analyzing the data after we load it using Ball gown. A sample csv (geuvadis_phenodata) containing phenotype data was provided by the authors of this paper, which is the same one we will use for our analysis. A Large ballgown data structure will be created by inputting the data directory containing the stringtie output files, a pattern contained within the stringtie files, and the phenotype information loaded in the prior line. The final preprocessing step is to filter out low abundance reads. The authors defined "low abundance" as any genes with less than 1% variance so I chose to maintain this threshold for consistency. We can now begin running data analysis on the Large ballgown data. 

```{r, message=FALSE} 
# load phenotype data (ID, sex, population)
pheno_data = read.csv("geuvadis_phenodata.csv")
# great Large ballgown data structure by running ballgown  
bg_chrX = ballgown(dataDir = "ballgown_out", samplePattern = "ERR", pData=pheno_data)
# filter out reads with low abundance (<1%)
bg_chrX_filt = subset(bg_chrX, "rowVars(texpr(bg_chrX))> 1", genomesubset=TRUE)
```

We have the choice between ID'ing statistically significant divergent *genes* or *transcripts*. The authors chose to run both using *sex* as the covariate. For now we will keep this consistent with the authors but later we will circle back to this and use *population* as a covariate. 

```{r} 
# ID transcripts using sex as independent variable 
results_transcripts = stattest(bg_chrX_filt,feature="transcript", covariate="sex", adjustvars = c("population"), getFC=TRUE, meas="FPKM")
# ID genes located on transcripts from above 
results_genes = stattest(bg_chrX_filt, feature="gene", covariate="sex", adjustvars = c("population"), getFC=TRUE, meas="FPKM")
# add gene names and IDs to transcripts (from line 38)
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)
# sort output transcripts and genes by p values 
results_transcripts = arrange(results_transcripts,pval) 
results_genes = arrange(results_genes,pval)
# save sorted output as csv (in working directory)
write.csv(results_transcripts, "chrX_transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv", row.names=FALSE)
```

Since we have conducted a large number of statistical testing, we will have to consider the high amount of false positive errors. Some false positive correction methods, such as the Bonferroni method, are too conservative - they reduce the false positive rate but also remove true positives. Therefore, we will attempt to manage the false positive error by controlling the false-discovery rate (FDR) through adjusted p values, which we call "q values." The q values are found by using characteristics of the p value distribution to adjust for multiple testing. Hence, the q value has higher power for finding true positives than conservative corrections like the Bonferroni method. 

Below we will use a q value threshold of 0.05 to identify statistically significant transcripts and genes. q = 0.05 means that 5% of *significant* tests will result in false positives, in comparison, p = 0.05 means that 5* of *all* tests will result in false positives. The true positives will be displayed as an output below each command. 

``` {r}
# ID (statistically significant) transcripts by q value < 0.05
subset(results_transcripts,results_transcripts$qval<0.05)
# ID (statistically significant) genes by q value < 0.05
subset(results_genes,results_genes$qval<0.05)
```

### 2) Complete the "Data visualization" section (Procedure 17-21) from <https://www.nature.com/articles/nprot.2016.095#an1> 

At this point, we have completed our first goal of answering the questions and following the guide provided in the Lecture 16 slides. So we will transition into the second goal of following the Data Visualization steps (17-21) used by the authors of the original paper (<https://www.nature.com/articles/nprot.2016.095#an1>). However, many of the commands the authors provided in this section do not run correctly; I suspect this is likely due to changes in R syntax over the past 5 years. Accordingly, I have provided the authors' original commands (commented out) and the modified commands which do result in the proper visualization. 

**Figure 3** of the paper shows the distribution of FPKM values for 12 samples, stratified by sex. I managed to get the boxplots to appear but could not manage to fill them by sex as the authors did. 

```{r}
# The provided command by the authors does not run because they used single quotes instead of double
# tropical= c(′darkorange′, ′dodgerblue′, ′hotpink′, ′limegreen′, ′yellow′)
# When replaced with double quotes the palette can be correctly created. 
tropical= c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
# The authors used the non-filtered Ballgown data and should have used single quotes about FPKM
# fpkm = texpr(bg_chrX,meas=″FPKM″)
# I believe the correct command should instead be: 
fpkm = texpr(bg_chrX_filt, meas='FPKM')
# Calculate log2 values of FPKM
fpkm = log2(fpkm+1)
# The boxplot command works but does not fill them by sex 
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')
```

I received one error "NAs introduced by coercion" which means some of the data was not formatted correctly for plotting. However the program did manage to overlook this to plot the majority of the data so this was not a critical point. 

**Figure 4** plots the FPKMs distributions in males and females for the 12th transcript. We begin by printing the transcript name and the gene containing it. When I tried to run the plot function, I received the following error: "unexpected input in "plot(fpkm[12,] ∼". Since I am not as familiar with R, I decided to run just the **> plot(fpkm[12,])** command to see if I could at least get the FPKMs for the 12th transcript to appear. This worked successfully. I had the same error message appear for the next command, points, so I again removed everything except **> points(fpkm[12,])** which ran correctly. In the future, I hope to be able to solve this problem once I am more familiar with the R syntax. 

```{r}
# print transcript name
ballgown::transcriptNames(bg_chrX)[12]
# print gene name 
ballgown::geneNames(bg_chrX)[12]
# Provided code for boxplots and points does not function 
# plot(fpkm[12,] ∼ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(bg_chrX)[12],' : ', ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab=″Sex″, ylab='log2(FPKM+1)')
# points(fpkm[12,] ∼ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))
# So instead I used a simpler verison
plot(fpkm[12,])
points(fpkm[12,])
```

**Figure 5** plots all transcript structure and expression levels which come from the same gene. The given code plots all transcripts within gene containing the first transcript, NR_001564. There are some errors in the first command: the authors' command plots data from the non-filtered Ballgown data structure and use single instead of double quotes. The correct code is provided below the commented-out authors' code. The second command requires similar modifications too, the first arg must be encased in double quotes while the groupvar arg must be in single quotes. 

```{r}
# authors' code =
# plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234′), sample=c('ERR188234'))
plotTranscripts(ballgown::geneIDs(bg_chrX_filt)[1729], bg_chrX, main=c("Gene XIST in sample ERR188234"), sample=c("ERR188234"))
# plotMeans(′MSTRG.56', bg_chrX_filt,groupvar=″sex″,legend=FALSE)
plotMeans("MSTRG.56", bg_chrX_filt, groupvar = 'sex', legend=FALSE) 
```

My figure has only one minor difference to the figure presented in the paper. The genomic position numbers are horizontal rather than diagonal which makes the boxes on my figure appear much smaller relative to the paper's **Figure 5**. However, this is just a minor difference. 

I also plotted a new figure which the authors' did not provide in this section. This figure shows the difference in isoform expression stratified by sex. 

# Part IV: Expression Analysis 2 

**Goals:**

  1) Run DESeq2 for two comparisons (sex, population) separately from HISAT2+featureCounts results.
  2) Run DESeq2 for two comparisons (sex, population) separately from STAR+featureCounts results.
  3) Generate meaning plots (e.g., volcano plot or MA-plot) and compare the above results.

Dr. Ahn provided a script, **DESeq2.r** to run this section of the lab. This script has been pre-commented to describe each command's function. 

### 1) Run DESeq2 for two comparisons (sex, population) separately from HISAT2+featureCounts results.
and 
### 2) Run DESeq2 for two comparisons (sex, population) separately from STAR+featureCounts results.

From our Lecture 17 notes, and in-class discussion, I attempted to modify the count_matrix.txt file which contained the RNA expression data to allow for sex- or population-based DE analysis. However, I was unsuccessful in both attempts because I was unsure how to modify the count matrix file to reflect sex or population information. 

One idea I had was to manually change the "ctrl" or "treat" designators to "male" or "female" at random to create 2 sham groups I could run for the purpose of analysis. I have provided the images of the sham figures (MA plot, heatmap) at the end of the file for viewing. 

### 3) Generate meaning plots (e.g., volcano plot or MA-plot) and compare the above results.

Returning to the given data, I ran the deseq2.r script Dr. Ahn provided to generate an MA plot and heatmap. 
```{r, message=FALSE}
# working directory
getwd() 
# read in count matrix
countData <- read.csv("count_matrix.txt", header=T, row.names=1, sep="\t") 
dim(countData)
head(countData) 
# basic QC
par(mar = rep(2, 4))
barplot(colSums(countData)*1e-6,
        names=colnames(countData),
        ylab="Library size (millions)")
```

Above we can see some **elementary QC** run on the RNA expression data. Each of the samples from control and treatment groups have roughly the same number of reads, or are at least on the same order of magnitude (15-20e6). Therefore, since we are assured our data is of comparable scale, we may proceed. 

```{r, message=FALSE}
# load library
library(DESeq2)
# create experiment labels (two conditions)
colData <- DataFrame(condition=factor(c("ctrl","ctrl","ctrl", "treat", "treat", "treat")))
# create DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
# run DEseq
dds <- DESeq(dds)
# visualize differentially expressed genes
par(mar = rep(2, 4))
plotMA(dds)
```
Above we can see the log-2 fold change (Y axis) against the mean number of normalized counts (X axis). The statistically significant DE genes are presented in blue while the non-significant are in gray. The normalization occurs by applying a negative binomial distribution to the raw number of mean counts, which flattens the central axis of the data. 

This data alone has no functional interpretation due to the complexity of gene networks. Therefore, to gain a biological interpretation, we must consider their known functions, or annotate it ourselves, through GO or KEGG pathways. 

```{r, message=FALSE}
# get differentially expressed genes
res <- results(dds)
# order by BH adjusted p-value
resOrdered <- res[order(res$padj),]
# top of ordered matrix
head(resOrdered)
# how many differentially expressed genes ? FDR=10%, |fold-change|>2 (up and down)
# get differentially expressed gene matrix
sig <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj<0.10 &
                    abs(resOrdered$log2FoldChange)>=1,]
# top of the differentially expressed genes
head(sig)
# how to create a heat map
# select genes
selected <- rownames(sig);selected
# load libraries for the heat map
library("RColorBrewer")
library("gplots")
# colors of the heat map
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
# heatmap
heatmap.2( log2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,]),
           col = hmcol, scale="row",
           Rowv = TRUE, Colv = FALSE,
           dendrogram="row",
           trace="none",
           margin=c(4,6), cexRow=0.5, cexCol=1, keysize=1 )
