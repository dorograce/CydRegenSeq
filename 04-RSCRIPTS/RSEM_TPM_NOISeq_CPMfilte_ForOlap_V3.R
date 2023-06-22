## Dorothy Mitchell ##
## Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
#Copyright (C) 2021 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin17.0 (64-bit)

## TPMS = Transcripts per million 

## Using gene counts (TPM) from RSEM to run differential gene expression (DGE) with NOISeqbio 

## loading required packages ##



## NOISeq ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("NOISeq")
library(NOISeq)


# Loading gene quantification results # 
# Uncut Replicate 1 #
unct_i_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/uncut-i-rep1_S21_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
unct_i_Rp1_B2_Trm


# Uncut Replicate 2 # 
unct_i_Rp2_B2_Trm <-read.delim("../03-DATA_FILES/uncut-i-rep2_S22_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
unct_i_Rp2_B2_Trm


# Uncut Replicate 3 # 
unct_i_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/uncut-i-rep3_S23_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
unct_i_Rp3_B2_Trm

# Uncut Replicate 3 ii # 
unct_ii_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/uncut-ii-rep3_S24_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
unct_ii_Rp3_B2_Trm

# 0Hab Replicate 1 #
hab0_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/0hab-rep1_S25_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab0_Rp1_B2_Trm

# 0Hab Replicate 2 # 
hab0_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/0hab-rep2_S26_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab0_Rp2_B2_Trm

# 0Hab Replicate 3 # 
hab0_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/0hab-rep3_S27_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab0_Rp3_B2_Trm

# 1Hab Replicate 1 # 
hab1_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/1hab-rep3_S30_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab1_Rp1_B2_Trm

# 1Hab Replicate 2 # 
hab1_Rp2_B2_Trm <- read.delim("../03-DATA_FILES/1hab-rep2_S29_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab1_Rp2_B2_Trm

# 1Hab Replicate 3 # 
hab1_Rp3_B2_Trm <- read.delim("../03-DATA_FILES/1hab-rep3_S30_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab1_Rp3_B2_Trm

# 3Hab Replicate 1 #
hab3_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/3hab-rep1_S31_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab3_Rp1_B2_Trm

# 3Hab Replicate 2 #
hab3_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/3hab-rep2_S32_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab3_Rp2_B2_Trm

# 3Hab Replicate 3 # 
hab3_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/3hab-rep3_S33_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab3_Rp3_B2_Trm

# 6Hab Replicate 1 #
hab6_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/6hab-rep1_S34_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab6_Rp1_B2_Trm

# 6Hab Replicate 2 # 
hab6_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/6hab-rep2_S35_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab6_Rp2_B2_Trm

# 6Hab Replicate 3 # 
hab6_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/6hab-rep3_S36_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab6_Rp3_B2_Trm

# 12Hab Replicate 1 # 
hab12_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/12hab-rep1_S37_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab12_Rp1_B2_Trm

# 12Hab Replicate 2 # 
hab12_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/12hab-rep2_S38_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab12_Rp2_B2_Trm

# 12Hab Replicate 3 #
hab12_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/12hab-rep3_S39_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab12_Rp3_B2_Trm

# 24Hab Replicate 1 # 
hab24_i_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/24hab-i-rep1_S40_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab24_i_Rp1_B2_Trm

# 24Hab Replicate 2 # 
hab24_i_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/24hab-i-rep2_S41_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab24_i_Rp2_B2_Trm

# 24Hab Replicate 3 #
hab24_i_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/24hab-i-rep3_S42_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab24_i_Rp3_B2_Trm


# 24Hab Replicate 3 ii #
hab24_ii_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/24hab-ii-rep3_S43_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab24_ii_Rp3_B2_Trm

# 48Hab Replicate 1 # 
hab48_i_Rp1_B2_Trm<-read.delim("../03-DATA_FILES/48hab-i-rep1_S44_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab48_i_Rp1_B2_Trm

# 48Hab Replicate 2 # 
hab48_i_Rp2_B2_Trm<-read.delim("../03-DATA_FILES/48hab-i-rep2_S45_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab48_i_Rp2_B2_Trm

# 48Hab Replicate 3 # 
hab48_i_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/48hab-i-rep3_S46_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab48_i_Rp3_B2_Trm

# 48Hab Replicate 3 ii # 
hab48_ii_Rp3_B2_Trm<-read.delim("../03-DATA_FILES/48hab-ii-rep3_S47_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
hab48_ii_Rp3_B2_Trm





## preparing the data for Noiseqbio 

RSEMcounts_B2_Trm_TPM <- data.frame(row.names = hab48_ii_Rp3_B2_Trm$gene_id,
                                    uncut_1_i=unct_i_Rp1_B2_Trm$TPM,
                                    uncut_2_i=unct_i_Rp2_B2_Trm$TPM,
                                    uncut_3_i=unct_i_Rp3_B2_Trm$TPM,
                                    hab00_1=hab0_Rp1_B2_Trm$TPM,
                                    hab00_2=hab0_Rp2_B2_Trm$TPM,
                                    hab00_3=hab0_Rp3_B2_Trm$TPM,
                                    hab01_1=hab1_Rp1_B2_Trm$TPM,
                                    hab01_2=hab1_Rp2_B2_Trm$TPM,
                                    hab01_3=hab1_Rp3_B2_Trm$TPM,
                                    hab03_1=hab3_Rp1_B2_Trm$TPM,
                                    hab03_2=hab3_Rp2_B2_Trm$TPM,
                                    hab03_3=hab3_Rp3_B2_Trm$TPM,
                                    hab06_1=hab6_Rp1_B2_Trm$TPM,
                                    hab06_2=hab6_Rp2_B2_Trm$TPM,
                                    hab06_3=hab6_Rp3_B2_Trm$TPM,
                                    hab12_1=hab12_Rp1_B2_Trm$TPM,
                                    hab12_2=hab12_Rp2_B2_Trm$TPM,
                                    hab12_3=hab12_Rp3_B2_Trm$TPM,
                                    hab24_1_i=hab24_i_Rp1_B2_Trm$TPM,
                                    hab24_2_i=hab24_i_Rp2_B2_Trm$TPM,
                                    hab24_3_i=hab24_i_Rp3_B2_Trm$TPM,
                                    hab48_1_i=hab48_i_Rp1_B2_Trm$TPM,
                                    hab48_2_i=hab48_i_Rp2_B2_Trm$TPM,
                                    hab48_3_i=hab48_i_Rp3_B2_Trm$TPM)

length(RSEMcounts_B2_Trm_TPM$uncut_1_i) #16548



### Differential Expression ###
### with Pre-Normalized Data ###
########################


### NOISeq ## - non parametric ## trying with filtering option CPM
## Noiseqbio accepts pre-normalized data -- TPM

length(RSEMcounts_B2_Trm_TPM) # 24 samples(column headers)

## Comparing Uncut with 10min
Noiseqbio_U_00_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(uncut_1_i,uncut_2_i,uncut_3_i,hab00_1,hab00_2,hab00_3))
mygroups_nois_U_00 = rep(c("2_uncut", "1_hab00"), each = 3)
mynoiseqbio_U_00 = NOISeq::readData(Noiseqbio_U_00_counts, factors = data.frame("group" = mygroups_nois_U_00))
mynoiseqbio_U_00 = noiseqbio(mynoiseqbio_U_00, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_U_00_down = degenes(mynoiseqbio_U_00, q = 0.95,M="down") #34
mynoiseqbio_U_00_up = degenes(mynoiseqbio_U_00, q = 0.95,M="up") #36
mynoiseqbio_DE_U_00 = degenes(mynoiseqbio_U_00, q = 0.95) #V.2.38.0=70 


## Comparing 10min with 01
Noiseqbio_00_01_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab00_1,hab00_2,hab00_3,hab01_1,hab01_2,hab01_3))
mygroups_nois_00_01 = rep(c("2_hab00", "1_hab01"), each = 3)
mynoiseqbio_00_01 = NOISeq::readData(Noiseqbio_00_01_counts, factors = data.frame("group" = mygroups_nois_00_01))
mynoiseqbio_00_01 = noiseqbio(mynoiseqbio_00_01, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_00_01_down = degenes(mynoiseqbio_00_01, q = 0.95,M="down") #  26
mynoiseqbio_00_01_up = degenes(mynoiseqbio_00_01, q = 0.95,M="up") # 34
mynoiseqbio_DE_00_01 = degenes(mynoiseqbio_00_01, q = 0.95) #V.2.38.0=60 


## Comparing 01 with 03
Noiseqbio_01_03_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab01_1,hab01_2,hab01_3,hab03_1,hab03_2,hab03_3))
mygroups_nois_01_03 = rep(c("2_hab01", "1_hab03"), each = 3)
mynoiseqbio_01_03 = NOISeq::readData(Noiseqbio_01_03_counts, factors = data.frame("group"= mygroups_nois_01_03))
mynoiseqbio_01_03 = noiseqbio(mynoiseqbio_01_03, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_01_03_down = degenes(mynoiseqbio_01_03, q = 0.95,M="down") #1257
mynoiseqbio_01_03_up = degenes(mynoiseqbio_01_03, q = 0.95,M="up") #1561
mynoiseqbio_DE_01_03 = degenes(mynoiseqbio_01_03, q = 0.95) #V.2.38.0=2818

## Comparing 03 with 06
Noiseqbio_03_06_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab03_1,hab03_2,hab03_3,hab06_1,hab06_2,hab06_3))
mygroups_nois_03_06 = rep(c("2_hab03", "1_hab06"), each = 3)
mynoiseqbio_03_06 = NOISeq::readData(Noiseqbio_03_06_counts, factors = data.frame("group"= mygroups_nois_03_06))
mynoiseqbio_03_06 = noiseqbio(mynoiseqbio_03_06, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_03_06_down = degenes(mynoiseqbio_03_06, q = 0.95,M="down") #2223
mynoiseqbio_03_06_up = degenes(mynoiseqbio_03_06, q = 0.95,M="up") #3699
mynoiseqbio_DE_03_06 = degenes(mynoiseqbio_03_06, q = 0.95)#V.2.38.0=5922


## Comparing 06 with 12
Noiseqbio_06_12_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab06_1,hab06_2,hab06_3,hab12_1,hab12_2,hab12_3))
mygroups_nois_06_12 = rep(c("2_hab06", "1_hab12"), each = 3)
mynoiseqbio_06_12 = NOISeq::readData(Noiseqbio_06_12_counts, factors = data.frame("group"= mygroups_nois_06_12))
mynoiseqbio_06_12 = noiseqbio(mynoiseqbio_06_12, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_06_12_down = degenes(mynoiseqbio_06_12, q = 0.95,M="down") #69
mynoiseqbio_06_12_up = degenes(mynoiseqbio_06_12, q = 0.95,M="up") #99
mynoiseqbio_DE_06_12 = degenes(mynoiseqbio_06_12, q = 0.95)#V.2.38.0=168

## Comparing 12 with 24
Noiseqbio_12_24_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab12_1,hab12_2,hab12_3,hab24_1_i,hab24_2_i,hab24_3_i))
mygroups_nois_12_24 = rep(c("2_hab12", "1_hab24"), each = 3)
mynoiseqbio_12_24 = NOISeq::readData(Noiseqbio_12_24_counts, factors = data.frame("group"= mygroups_nois_12_24))
mynoiseqbio_12_24 = noiseqbio(mynoiseqbio_12_24, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_12_24_down = degenes(mynoiseqbio_12_24, q = 0.95,M="down") # 67
mynoiseqbio_12_24_up = degenes(mynoiseqbio_12_24, q = 0.95,M="up") #45
mynoiseqbio_DE_12_24 = degenes(mynoiseqbio_12_24, q = 0.95) #V.2.38.0=112

## Comparing 24 hab with 48hab
Noiseqbio_24_48_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab24_1_i,hab24_2_i,hab24_3_i,hab48_1_i,hab48_2_i,hab48_3_i))
mygroups_nois_24_48 = rep(c("2_hab24", "1_hab48"), each = 3)
mynoiseqbio_24_48 = NOISeq::readData(Noiseqbio_24_48_counts, factors = data.frame("group"= mygroups_nois_24_48))
mynoiseqbio_24_48 = noiseqbio(mynoiseqbio_24_48, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_24_48_down = degenes(mynoiseqbio_24_48, q = 0.95,M="down") #119
mynoiseqbio_24_48_up = degenes(mynoiseqbio_24_48, q = 0.95,M="up") #127
mynoiseqbio_DE_24_48 = degenes(mynoiseqbio_24_48, q = 0.95) #V.2.38.0=246



#Make a list of all of the result data frames 

Noiseqbio_TPM_res_list <- list(mynoiseqbio_DE_U_00,
                               mynoiseqbio_DE_00_01,
                               mynoiseqbio_DE_01_03,
                               mynoiseqbio_DE_06_12,
                               mynoiseqbio_DE_12_24,
                               mynoiseqbio_DE_24_48)


Noiseqbio_TPM_res_list_matrix <- lapply(Noiseqbio_TPM_res_list, function(x) as.matrix(x))

Noiseqbio_TPM_res_list_abs_filt_matrix_bind <- do.call(rbind,Noiseqbio_TPM_res_list_matrix)

Noiseqbio_TPM_res_list_abs_filt_matrix_bind <- as.matrix(rownames(Noiseqbio_TPM_res_list_abs_filt_matrix_bind))


length(Noiseqbio_TPM_res_list_abs_filt_matrix_bind)#3474

# Want to collect all of the DE genes, but take away any repeats so a gene does not appear more than once in a list 

Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups <-  as.matrix(Noiseqbio_TPM_res_list_abs_filt_matrix_bind[!duplicated(Noiseqbio_TPM_res_list_abs_filt_matrix_bind)])
length(Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups)#3092


write.csv(Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups,row.names=F,file = "NOISeqbio_sig_genes_TPM_nodups_v2.38.0_2_V3.csv")

# For figure in main text
#For BLAST

HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)


## Uncut - 10 min
#DOWN 
length(mynoiseqbio_U_00_down$`1_hab00_mean`) #34

#BLAST 
mynoiseqbio_U_00_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_U_00_down,by='row.names')
#order log2FC 
mynoiseqbio_U_00_down_blastp_merge_order<- mynoiseqbio_U_00_down_blastp_merge[order(mynoiseqbio_U_00_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_U_00_down_blastp_merge_order$Row.names) #8

#UP
length(mynoiseqbio_U_00_up$`1_hab00_mean`) #36
#BLAST 
mynoiseqbio_U_00_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_U_00_up,by='row.names')

#order log2FC
mynoiseqbio_U_00_up_blastp_merge_order<- mynoiseqbio_U_00_up_blastp_merge[order(mynoiseqbio_U_00_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_U_00_up_blastp_merge_order$Row.names) #13

#10 min - 01
#DOWN
length(mynoiseqbio_00_01_down$`1_hab01_mean`) #26
#BLAST 
mynoiseqbio_00_01_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_00_01_down,by='row.names')
#order log2FC
mynoiseqbio_00_01_down_blastp_merge_order<- mynoiseqbio_00_01_down_blastp_merge[order(mynoiseqbio_00_01_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_00_01_down_blastp_merge_order$Row.names)#5

#UP
length(mynoiseqbio_00_01_up$`1_hab01_mean`) #34
#BLAST 
mynoiseqbio_00_01_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_00_01_up,by='row.names')

#order log2FC 
mynoiseqbio_00_01_up_blastp_merge_order <- mynoiseqbio_00_01_up_blastp_merge[order(mynoiseqbio_00_01_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_00_01_up_blastp_merge_order$Row.names) #12

## 01-03 
#DOWN 
length(mynoiseqbio_01_03_down$`1_hab03_mean`) #1257

#BLAST 
mynoiseqbio_01_03_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_01_03_down,by='row.names')
#order log2FC
mynoiseqbio_01_03_down_blastp_merge_order<- mynoiseqbio_01_03_down_blastp_merge[order(mynoiseqbio_01_03_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_01_03_down_blastp_merge_order$Row.names) #769


#UP
length(mynoiseqbio_01_03_up$`1_hab03_mean`)#1561
#BLAST 
mynoiseqbio_01_03_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_01_03_up,by='row.names')

#order log2FC 
mynoiseqbio_01_03_up_blastp_merge_order <- mynoiseqbio_01_03_up_blastp_merge[order(mynoiseqbio_01_03_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_01_03_up_blastp_merge_order$Row.names) #1081

## 03-06
#DOWN 
length(mynoiseqbio_03_06_down$`1_hab06_mean`)#2223

#BLAST 
mynoiseqbio_03_06_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_03_06_down,by='row.names')
#order log2FC
mynoiseqbio_03_06_down_blastp_merge_order<- mynoiseqbio_03_06_down_blastp_merge[order(mynoiseqbio_03_06_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_03_06_down_blastp_merge_order$Row.names) #1458

#UP
length(mynoiseqbio_03_06_up$`1_hab06_mean`)#3699
#BLAST 
mynoiseqbio_03_06_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_03_06_up,by='row.names')

#order log2FC 
mynoiseqbio_03_06_up_blastp_merge_order <- mynoiseqbio_03_06_up_blastp_merge[order(mynoiseqbio_03_06_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_03_06_up_blastp_merge_order$Row.names)#2226

## 06-12
#DOWN 
length(mynoiseqbio_06_12_down$`1_hab12_mean`) #69

#BLAST 
mynoiseqbio_06_12_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_06_12_down,by='row.names')
#order log2FC
mynoiseqbio_06_12_down_blastp_merge_order<- mynoiseqbio_06_12_down_blastp_merge[order(mynoiseqbio_06_12_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_06_12_down_blastp_merge_order$Row.names) #32

#UP
length(mynoiseqbio_06_12_up$`1_hab12_mean`)#99
#BLAST 
mynoiseqbio_06_12_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_06_12_up,by='row.names')

#order log2FC
mynoiseqbio_06_12_up_blastp_merge_order <- mynoiseqbio_06_12_up_blastp_merge[order(mynoiseqbio_06_12_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_06_12_up_blastp_merge_order$Row.names) #41
## 12-24
#DOWN 
length(mynoiseqbio_12_24_down$`1_hab24_mean`)#67

#BLAST 
mynoiseqbio_12_24_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_12_24_down,by='row.names')
#order log2FC 
mynoiseqbio_12_24_down_blastp_merge_order<- mynoiseqbio_12_24_down_blastp_merge[order(mynoiseqbio_12_24_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_12_24_down_blastp_merge_order$Row.names)#32

#UP
length(mynoiseqbio_12_24_up$`1_hab24_mean`)#45
#BLAST 
mynoiseqbio_12_24_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_12_24_up,by='row.names')

#order log2FC 
mynoiseqbio_12_24_up_blastp_merge_order <- mynoiseqbio_12_24_up_blastp_merge[order(mynoiseqbio_12_24_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_12_24_up_blastp_merge_order$Row.names) #14

## 24-48
#DOWN 
mynoiseqbio_24_48_down
length(mynoiseqbio_24_48_down$`1_hab48_mean`) #119

#BLAST 
mynoiseqbio_24_48_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_24_48_down,by='row.names')

#order log2FC
mynoiseqbio_24_48_down_blastp_merge_order<- mynoiseqbio_24_48_down_blastp_merge[order(mynoiseqbio_24_48_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_24_48_down_blastp_merge_order$Row.names)#56

#UP
length(mynoiseqbio_24_48_up$`1_hab48_mean`)#127
#BLAST 
mynoiseqbio_24_48_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_24_48_up,by='row.names')

#order log2FC
mynoiseqbio_24_48_up_blastp_merge_order <- mynoiseqbio_24_48_up_blastp_merge[order(mynoiseqbio_24_48_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_24_48_up_blastp_merge_order$Row.names)#61

#Merge all of the sig genes in each interval with their corresponding p-values and LogFC into a list 
NOIS_DE_log2FC <- list(sig_U_00_up = mynoiseqbio_U_00_up_blastp_merge_order ,
                       sig_U_00_down = mynoiseqbio_U_00_down_blastp_merge_order,
                             sig_00_01_down= mynoiseqbio_00_01_down_blastp_merge_order ,
                             sig_00_01_up=mynoiseqbio_00_01_up_blastp_merge_order ,
                             sig_01_03_down =mynoiseqbio_01_03_down_blastp_merge_order ,
                             sig_01_03_up =mynoiseqbio_01_03_up_blastp_merge_order ,
                             sig_03_06_down = mynoiseqbio_03_06_down_blastp_merge_order,
                             sig_03_06_up =mynoiseqbio_03_06_up_blastp_merge_order ,
                             sig_06_12_down = mynoiseqbio_06_12_down_blastp_merge_order,
                             sig_06_12_up = mynoiseqbio_06_12_up_blastp_merge_order,
                      sig_12_24_down = mynoiseqbio_12_24_down_blastp_merge_order,
                      sig_12_24_up = mynoiseqbio_12_24_up_blastp_merge_order,
                      sig_24_48_down = mynoiseqbio_24_48_down_blastp_merge_order,
                      sig_24_48_up = mynoiseqbio_24_48_up_blastp_merge_order)


write.csv(NOIS_DE_log2FC, file = "NOIS_DE_log2FC_V3.csv")

