## Dorothy Mitchell ##
## Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

#R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
#Copyright (C) 2023 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin20 (64-bit)

## TPMS = Transcripts per million 

## Using gene counts (TPM) from RSEM to run differential gene expression (DGE) with NOISeqbio 

## loading required packages ##



## NOISeq ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("NOISeq")
library(NOISeq)
sessionInfo()
v2.44.0

# Loading gene quantification results # 
# Uncut Replicate 1 #
unct_i_Rp1_B2_Trm<-read.delim("03-RSEM/uncut-i-rep1_S21.RSEM.genes.results")
unct_i_Rp1_B2_Trm


# Uncut Replicate 2 # 
unct_i_Rp2_B2_Trm <-read.delim("03-RSEM/uncut-i-rep2_S22.RSEM.genes.results")
unct_i_Rp2_B2_Trm


# Uncut Replicate 3 # 
unct_i_Rp3_B2_Trm<-read.delim("03-RSEM/uncut-i-rep3_S23.RSEM.genes.results")
unct_i_Rp3_B2_Trm


# 0Hab Replicate 1 #
hab0_Rp1_B2_Trm<-read.delim("03-RSEM/0hab-rep1_S25.RSEM.genes.results")
hab0_Rp1_B2_Trm

# 0Hab Replicate 2 # 
hab0_Rp2_B2_Trm<-read.delim("03-RSEM/0hab-rep2_S26.RSEM.genes.results")
hab0_Rp2_B2_Trm

# 0Hab Replicate 3 # 
hab0_Rp3_B2_Trm<-read.delim("03-RSEM/0hab-rep3_S27.RSEM.genes.results")
hab0_Rp3_B2_Trm

# 1Hab Replicate 1 # 
hab1_Rp1_B2_Trm<-read.delim("03-RSEM/1hab-rep1_S28.RSEM.genes.results")
hab1_Rp1_B2_Trm

# 1Hab Replicate 2 # 
hab1_Rp2_B2_Trm <- read.delim("03-RSEM/1hab-rep2_S29.RSEM.genes.results")
hab1_Rp2_B2_Trm

# 1Hab Replicate 3 # 
hab1_Rp3_B2_Trm <- read.delim("03-RSEM/1hab-rep3_S30.RSEM.genes.results")
hab1_Rp3_B2_Trm

# 3Hab Replicate 1 #
hab3_Rp1_B2_Trm<-read.delim("03-RSEM/3hab-rep1_S31.RSEM.genes.results")
hab3_Rp1_B2_Trm

# 3Hab Replicate 2 #
hab3_Rp2_B2_Trm<-read.delim("03-RSEM/3hab-rep2_S32.RSEM.genes.results")
hab3_Rp2_B2_Trm

# 3Hab Replicate 3 # 
hab3_Rp3_B2_Trm<-read.delim("03-RSEM/3hab-rep3_S33.RSEM.genes.results")
hab3_Rp3_B2_Trm

# 6Hab Replicate 1 #
hab6_Rp1_B2_Trm<-read.delim("03-RSEM/6hab-rep1_S34.RSEM.genes.results")
hab6_Rp1_B2_Trm

# 6Hab Replicate 2 # 
hab6_Rp2_B2_Trm<-read.delim("03-RSEM/6hab-rep2_S35.RSEM.genes.results")
hab6_Rp2_B2_Trm

# 6Hab Replicate 3 # 
hab6_Rp3_B2_Trm<-read.delim("03-RSEM/6hab-rep3_S36.RSEM.genes.results")
hab6_Rp3_B2_Trm

# 12Hab Replicate 1 # 
hab12_Rp1_B2_Trm<-read.delim("03-RSEM/12hab-rep1_S37.RSEM.genes.results")
hab12_Rp1_B2_Trm

# 12Hab Replicate 2 # 
hab12_Rp2_B2_Trm<-read.delim("03-RSEM/12hab-rep2_S38.RSEM.genes.results")
hab12_Rp2_B2_Trm

# 12Hab Replicate 3 #
hab12_Rp3_B2_Trm<-read.delim("03-RSEM/12hab-rep3_S39.RSEM.genes.results")
hab12_Rp3_B2_Trm

# 24Hab Replicate 1 # 
hab24_i_Rp1_B2_Trm<-read.delim("03-RSEM/24hab-i-rep1_S40.RSEM.genes.results")
hab24_i_Rp1_B2_Trm

# 24Hab Replicate 2 # 
hab24_i_Rp2_B2_Trm<-read.delim("03-RSEM/24hab-i-rep2_S41.RSEM.genes.results")
hab24_i_Rp2_B2_Trm

# 24Hab Replicate 3 #
hab24_i_Rp3_B2_Trm<-read.delim("03-RSEM/24hab-i-rep3_S42.RSEM.genes.results")
hab24_i_Rp3_B2_Trm

# 48Hab Replicate 1 # 
hab48_i_Rp1_B2_Trm<-read.delim("03-RSEM/48hab-i-rep1_S44.RSEM.genes.results")
hab48_i_Rp1_B2_Trm

# 48Hab Replicate 2 # 
hab48_i_Rp2_B2_Trm<-read.delim("03-RSEM/48hab-i-rep2_S45.RSEM.genes.results")
hab48_i_Rp2_B2_Trm

# 48Hab Replicate 3 # 
hab48_i_Rp3_B2_Trm<-read.delim("03-RSEM/48hab-i-rep3_S46.RSEM.genes.results")
hab48_i_Rp3_B2_Trm


## preparing the data for Noiseqbio 

RSEMcounts_B2_Trm_TPM <- data.frame(row.names = hab0_Rp1_B2_Trm$gene_id,
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
mynoiseqbio_U_00_down = degenes(mynoiseqbio_U_00, q = 0.95,M="down") #57
mynoiseqbio_U_00_up = degenes(mynoiseqbio_U_00, q = 0.95,M="up") #91
mynoiseqbio_DE_U_00 = degenes(mynoiseqbio_U_00, q = 0.95) #148 v2.44.0 


## Comparing 10min with 01
Noiseqbio_00_01_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab00_1,hab00_2,hab00_3,hab01_1,hab01_2,hab01_3))
mygroups_nois_00_01 = rep(c("2_hab00", "1_hab01"), each = 3)
mynoiseqbio_00_01 = NOISeq::readData(Noiseqbio_00_01_counts, factors = data.frame("group" = mygroups_nois_00_01))
mynoiseqbio_00_01 = noiseqbio(mynoiseqbio_00_01, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_00_01_down = degenes(mynoiseqbio_00_01, q = 0.95,M="down") #16
mynoiseqbio_00_01_up = degenes(mynoiseqbio_00_01, q = 0.95,M="up") #50
mynoiseqbio_DE_00_01 = degenes(mynoiseqbio_00_01, q = 0.95) #66 v2.44.0


## Comparing 01 with 03
Noiseqbio_01_03_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab01_1,hab01_2,hab01_3,hab03_1,hab03_2,hab03_3))
mygroups_nois_01_03 = rep(c("2_hab01", "1_hab03"), each = 3)
mynoiseqbio_01_03 = NOISeq::readData(Noiseqbio_01_03_counts, factors = data.frame("group"= mygroups_nois_01_03))
mynoiseqbio_01_03 = noiseqbio(mynoiseqbio_01_03, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_01_03_down = degenes(mynoiseqbio_01_03, q = 0.95,M="down") #340
mynoiseqbio_01_03_up = degenes(mynoiseqbio_01_03, q = 0.95,M="up") #94
mynoiseqbio_DE_01_03 = degenes(mynoiseqbio_01_03, q = 0.95) #434 v2.44.0

## Comparing 03 with 06
Noiseqbio_03_06_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab03_1,hab03_2,hab03_3,hab06_1,hab06_2,hab06_3))
mygroups_nois_03_06 = rep(c("2_hab03", "1_hab06"), each = 3)
mynoiseqbio_03_06 = NOISeq::readData(Noiseqbio_03_06_counts, factors = data.frame("group"= mygroups_nois_03_06))
mynoiseqbio_03_06 = noiseqbio(mynoiseqbio_03_06, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_03_06_down = degenes(mynoiseqbio_03_06, q = 0.95,M="down") #2250
mynoiseqbio_03_06_up = degenes(mynoiseqbio_03_06, q = 0.95,M="up") #3843
mynoiseqbio_DE_03_06 = degenes(mynoiseqbio_03_06, q = 0.95)#6093 v2.44.0


## Comparing 06 with 12
Noiseqbio_06_12_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab06_1,hab06_2,hab06_3,hab12_1,hab12_2,hab12_3))
mygroups_nois_06_12 = rep(c("2_hab06", "1_hab12"), each = 3)
mynoiseqbio_06_12 = NOISeq::readData(Noiseqbio_06_12_counts, factors = data.frame("group"= mygroups_nois_06_12))
mynoiseqbio_06_12 = noiseqbio(mynoiseqbio_06_12, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_06_12_down = degenes(mynoiseqbio_06_12, q = 0.95,M="down") #79
mynoiseqbio_06_12_up = degenes(mynoiseqbio_06_12, q = 0.95,M="up") #58
mynoiseqbio_DE_06_12 = degenes(mynoiseqbio_06_12, q = 0.95)#137 v2.44.0

## Comparing 12 with 24
Noiseqbio_12_24_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab12_1,hab12_2,hab12_3,hab24_1_i,hab24_2_i,hab24_3_i))
mygroups_nois_12_24 = rep(c("2_hab12", "1_hab24"), each = 3)
mynoiseqbio_12_24 = NOISeq::readData(Noiseqbio_12_24_counts, factors = data.frame("group"= mygroups_nois_12_24))
mynoiseqbio_12_24 = noiseqbio(mynoiseqbio_12_24, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_12_24_down = degenes(mynoiseqbio_12_24, q = 0.95,M="down") #78
mynoiseqbio_12_24_up = degenes(mynoiseqbio_12_24, q = 0.95,M="up") #45
mynoiseqbio_DE_12_24 = degenes(mynoiseqbio_12_24, q = 0.95) #123 v2.44.0

## Comparing 24 hab with 48hab
Noiseqbio_24_48_counts = subset(RSEMcounts_B2_Trm_TPM, select = c(hab24_1_i,hab24_2_i,hab24_3_i,hab48_1_i,hab48_2_i,hab48_3_i))
mygroups_nois_24_48 = rep(c("2_hab24", "1_hab48"), each = 3)
mynoiseqbio_24_48 = NOISeq::readData(Noiseqbio_24_48_counts, factors = data.frame("group"= mygroups_nois_24_48))
mynoiseqbio_24_48 = noiseqbio(mynoiseqbio_24_48, norm = "n", k = NULL, factor = "group", r = 50,filter=1)
mynoiseqbio_24_48_down = degenes(mynoiseqbio_24_48, q = 0.95,M="down") #62
mynoiseqbio_24_48_up = degenes(mynoiseqbio_24_48, q = 0.95,M="up") #108
mynoiseqbio_DE_24_48 = degenes(mynoiseqbio_24_48, q = 0.95) #170 v2.44.0



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


length(Noiseqbio_TPM_res_list_abs_filt_matrix_bind)#1078

# Want to collect all of the DE genes, but take away any repeats so a gene does not appear more than once in a list 

Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups <-  as.matrix(Noiseqbio_TPM_res_list_abs_filt_matrix_bind[!duplicated(Noiseqbio_TPM_res_list_abs_filt_matrix_bind)])
length(Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups)#834


write.csv(Noiseqbio_TPM_res_list_abs_filt_matrix_bind_nodups,row.names=F,file = "NOISeq_sig_v4.csv")

# For NOISeq summary figure in supplement
#For BLAST
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

# For consensus gene set annotation 
con_genenames <- read.csv(file="consensus_DE_v4.csv",1) # ML gene names of the genes identified as DE in all three methods
con_genenames_df <- data.frame(row.names= con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111...,V1 = con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111... )

## Uncut - 10 min DOWN
length(mynoiseqbio_U_00_down$`1_hab00_mean`) #57

#BLAST DOWN
mynoiseqbio_U_00_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_U_00_down,by='row.names')

#order log2FC 
mynoiseqbio_U_00_down_blastp_merge_order<- mynoiseqbio_U_00_down_blastp_merge[order(mynoiseqbio_U_00_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_U_00_down_blastp_merge_order$Row.names) #18

#Consensus DOWN order
mynoiseqbio_U_00_down_blastp_con_merge <- merge(mynoiseqbio_U_00_down_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_U_00_down_blastp_con_merge_order<- mynoiseqbio_U_00_down_blastp_con_merge[order(mynoiseqbio_U_00_down_blastp_con_merge$log2FC, decreasing = F),]

## Uncut - 10 min UP
length(mynoiseqbio_U_00_up$`1_hab00_mean`) #91

#BLAST UP
mynoiseqbio_U_00_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_U_00_up,by='row.names')

#order log2FC
mynoiseqbio_U_00_up_blastp_merge_order<- mynoiseqbio_U_00_up_blastp_merge[order(mynoiseqbio_U_00_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_U_00_up_blastp_merge_order$Row.names) #36

#Consensus UP order 
mynoiseqbio_U_00_up_blastp_con_merge <- merge(mynoiseqbio_U_00_up_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_U_00_up_blastp_con_merge_order <- mynoiseqbio_U_00_up_blastp_con_merge[order(mynoiseqbio_U_00_up_blastp_con_merge$log2FC, decreasing = T),]

##10 min - 01 DOWN
# DOWN
length(mynoiseqbio_00_01_down$`1_hab01_mean`) #16
#BLAST DOWN 
mynoiseqbio_00_01_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_00_01_down,by='row.names')
#order log2FC
mynoiseqbio_00_01_down_blastp_merge_order<- mynoiseqbio_00_01_down_blastp_merge[order(mynoiseqbio_00_01_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_00_01_down_blastp_merge_order$Row.names)#5

#Consensus DOWN
mynoiseqbio_00_01_down_blastp_con_merge <- merge(mynoiseqbio_00_01_down_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_00_01_down_blastp_con_merge_order <- mynoiseqbio_00_01_down_blastp_con_merge[order(mynoiseqbio_00_01_down_blastp_con_merge$log2FC, decreasing = F),]


##10 min - 01 UP
length(mynoiseqbio_00_01_up$`1_hab01_mean`) #50

#BLAST UP
mynoiseqbio_00_01_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_00_01_up,by='row.names')

#order log2FC 
mynoiseqbio_00_01_up_blastp_merge_order <- mynoiseqbio_00_01_up_blastp_merge[order(mynoiseqbio_00_01_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_00_01_up_blastp_merge_order$Row.names) #27

#Consensus UP
mynoiseqbio_00_01_up_blastp_con_merge<- merge(mynoiseqbio_00_01_up_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_00_01_up_blastp_con_merge_order <- mynoiseqbio_00_01_up_blastp_con_merge[order(mynoiseqbio_00_01_up_blastp_con_merge$log2FC, decreasing = T),]


## 01-03 DOWN
#DOWN 
length(mynoiseqbio_01_03_down$`1_hab03_mean`) #340

#BLAST DOWN
mynoiseqbio_01_03_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_01_03_down,by='row.names')
length(mynoiseqbio_01_03_down_blastp_merge$V1)#227

#order log2FC
mynoiseqbio_01_03_down_blastp_merge_order<- mynoiseqbio_01_03_down_blastp_merge[order(mynoiseqbio_01_03_down_blastp_merge$log2FC, decreasing = F),]


#Consensus DOWN
mynoiseqbio_01_03_down_blastp_con_merge<- merge(mynoiseqbio_01_03_down_blastp_merge,con_genenames_df,by='V1')
length(mynoiseqbio_01_03_down_blastp_con_merge$V1) #37
mynoiseqbio_01_03_down_blastp_con_merge_order<- mynoiseqbio_01_03_down_blastp_con_merge[order(mynoiseqbio_01_03_down_blastp_con_merge$log2FC, decreasing = F),]
length(mynoiseqbio_01_03_down_blastp_con_merge_order$Row.names) #37


## 01-03 UP
length(mynoiseqbio_01_03_up$`1_hab03_mean`)#94

#BLAST 
mynoiseqbio_01_03_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_01_03_up,by='row.names')

#order log2FC 
mynoiseqbio_01_03_up_blastp_merge_order <- mynoiseqbio_01_03_up_blastp_merge[order(mynoiseqbio_01_03_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_01_03_up_blastp_merge_order$Row.names) #61

#Consensus UP 
mynoiseqbio_01_03_up_blastp_con_merge<- merge(mynoiseqbio_01_03_up_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_01_03_up_blastp_con_merge_order <- mynoiseqbio_01_03_up_blastp_con_merge[order(mynoiseqbio_01_03_up_blastp_con_merge$log2FC, decreasing = T),]


## 03-06 DOWN
#DOWN 
length(mynoiseqbio_03_06_down$`1_hab06_mean`)#2250

#BLAST 
mynoiseqbio_03_06_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_03_06_down,by='row.names')
#order log2FC
mynoiseqbio_03_06_down_blastp_merge_order<- mynoiseqbio_03_06_down_blastp_merge[order(mynoiseqbio_03_06_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_03_06_down_blastp_merge_order$Row.names) #1483

#Consensus DOWN 
mynoiseqbio_03_06_down_blastp_con_merge<- merge(mynoiseqbio_03_06_down_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_03_06_down_blastp_con_merge_order<- mynoiseqbio_03_06_down_blastp_con_merge[order(mynoiseqbio_03_06_down_blastp_con_merge$log2FC, decreasing = F),]


## 03-06 UP
length(mynoiseqbio_03_06_up$`1_hab06_mean`)#3843
#BLAST 
mynoiseqbio_03_06_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_03_06_up,by='row.names')

#order log2FC 
mynoiseqbio_03_06_up_blastp_merge_order <- mynoiseqbio_03_06_up_blastp_merge[order(mynoiseqbio_03_06_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_03_06_up_blastp_merge_order$Row.names)#2305

#Consensus UP
mynoiseqbio_03_06_up_blastp_con_merge <- merge(mynoiseqbio_03_06_up_blastp_merge,con_genenames_df,by='V1')
mynoiseqbio_03_06_up_blastp_con_merge_order <- mynoiseqbio_03_06_up_blastp_con_merge[order(mynoiseqbio_03_06_up_blastp_con_merge$log2FC, decreasing = T),]



## 06-12
#DOWN 
length(mynoiseqbio_06_12_down$`1_hab12_mean`) #79

#BLAST 
mynoiseqbio_06_12_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_06_12_down,by='row.names')
#order log2FC
mynoiseqbio_06_12_down_blastp_merge_order<- mynoiseqbio_06_12_down_blastp_merge[order(mynoiseqbio_06_12_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_06_12_down_blastp_merge_order$Row.names) #39

#Consensus Down 
mynoiseqbio_06_12_down_blastp_con_merge_order<- merge(mynoiseqbio_06_12_down_blastp_merge_order,con_genenames_df,by='V1')


#UP
length(mynoiseqbio_06_12_up$`1_hab12_mean`)#58
#BLAST 
mynoiseqbio_06_12_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_06_12_up,by='row.names')

#order log2FC
mynoiseqbio_06_12_up_blastp_merge_order <- mynoiseqbio_06_12_up_blastp_merge[order(mynoiseqbio_06_12_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_06_12_up_blastp_merge_order$Row.names) #23

#Consensus up 
mynoiseqbio_06_12_up_blastp_con_merge<- merge(mynoiseqbio_06_12_up_blastp_merge_order,con_genenames_df,by='V1')
mynoiseqbio_06_12_up_blastp_con_merge_order <- mynoiseqbio_06_12_up_blastp_con_merge[order(mynoiseqbio_06_12_up_blastp_con_merge$log2FC, decreasing = T),]


## 12-24
#DOWN 
length(mynoiseqbio_12_24_down$`1_hab24_mean`)#78

#BLAST 
mynoiseqbio_12_24_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_12_24_down,by='row.names')
#order log2FC 
mynoiseqbio_12_24_down_blastp_merge_order<- mynoiseqbio_12_24_down_blastp_merge[order(mynoiseqbio_12_24_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_12_24_down_blastp_merge_order$Row.names)#38

#Consensus down 
mynoiseqbio_12_24_down_blastp_con_merge<- merge(mynoiseqbio_12_24_down_blastp_merge_order,con_genenames_df,by='V1')
mynoiseqbio_12_24_down_blastp_con_merge_order <- mynoiseqbio_12_24_down_blastp_con_merge[order(mynoiseqbio_12_24_down_blastp_con_merge$log2FC, decreasing = F),]



#UP
length(mynoiseqbio_12_24_up$`1_hab24_mean`)#45
#BLAST 
mynoiseqbio_12_24_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_12_24_up,by='row.names')

#order log2FC 
mynoiseqbio_12_24_up_blastp_merge_order <- mynoiseqbio_12_24_up_blastp_merge[order(mynoiseqbio_12_24_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_12_24_up_blastp_merge_order$Row.names) #15

#Consensus UP
mynoiseqbio_12_24_up_blastp_con_merge<- merge(mynoiseqbio_12_24_up_blastp_merge_order,con_genenames_df,by='V1')
mynoiseqbio_12_24_up_blastp_con_merge_order <- mynoiseqbio_12_24_up_blastp_con_merge[order(mynoiseqbio_12_24_up_blastp_con_merge$log2FC, decreasing = T),]


## 24-48
#DOWN 
mynoiseqbio_24_48_down
length(mynoiseqbio_24_48_down$`1_hab48_mean`) #62

#BLAST 
mynoiseqbio_24_48_down_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_24_48_down,by='row.names')

#order log2FC
mynoiseqbio_24_48_down_blastp_merge_order<- mynoiseqbio_24_48_down_blastp_merge[order(mynoiseqbio_24_48_down_blastp_merge$log2FC, decreasing = F),]
length(mynoiseqbio_24_48_down_blastp_merge_order$Row.names)#26

#Consensus down 
mynoiseqbio_24_48_down_blastp_con_merge<- merge(mynoiseqbio_24_48_down_blastp_merge_order,con_genenames_df,by='V1')
mynoiseqbio_24_48_down_blastp_con_merge_order <- mynoiseqbio_24_48_down_blastp_con_merge[order(mynoiseqbio_24_48_down_blastp_con_merge$log2FC, decreasing = F),]



#UP
length(mynoiseqbio_24_48_up$`1_hab48_mean`)#108
#BLAST 
mynoiseqbio_24_48_up_blastp_merge <- merge(HumRef2020ML_df, mynoiseqbio_24_48_up,by='row.names')

#order log2FC
mynoiseqbio_24_48_up_blastp_merge_order <- mynoiseqbio_24_48_up_blastp_merge[order(mynoiseqbio_24_48_up_blastp_merge$log2FC, decreasing = T),]
length(mynoiseqbio_24_48_up_blastp_merge_order$Row.names)#52

#Consensus up 
mynoiseqbio_24_48_up_blastp_con_merge<- merge(mynoiseqbio_24_48_up_blastp_merge_order,con_genenames_df,by='V1')
mynoiseqbio_24_48_up_blastp_con_merge_order <- mynoiseqbio_24_48_up_blastp_con_merge[order(mynoiseqbio_24_48_up_blastp_con_merge$log2FC, decreasing = T),]


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

options(max.print=999999)
capture.output(NOIS_DE_log2FC, file = "NOISeq_sig_log2FC_v4.csv")


#Merge all of the sig genes also in the consensus in each interval with their corresponding p-values and LogFC into a list 
NOIS_DE_log2FC_consensus <- list(sig_U_00_up = mynoiseqbio_U_00_up_blastp_con_merge_order ,
                       sig_U_00_down = mynoiseqbio_U_00_down_blastp_con_merge_order,
                       sig_00_01_down= mynoiseqbio_00_01_down_blastp_con_merge_order ,
                       sig_00_01_up=mynoiseqbio_00_01_up_blastp_con_merge_order ,
                       sig_01_03_down =mynoiseqbio_01_03_down_blastp_con_merge_order ,
                       sig_01_03_up =mynoiseqbio_01_03_up_blastp_con_merge_order ,
                       sig_03_06_down = mynoiseqbio_03_06_down_blastp_con_merge_order,
                       sig_03_06_up =mynoiseqbio_03_06_up_blastp_con_merge_order ,
                       sig_06_12_down = mynoiseqbio_06_12_down_blastp_con_merge_order,
                       sig_06_12_up = mynoiseqbio_06_12_up_blastp_con_merge_order,
                       sig_12_24_down = mynoiseqbio_12_24_down_blastp_con_merge_order,
                       sig_12_24_up = mynoiseqbio_12_24_up_blastp_con_merge_order,
                       sig_24_48_down = mynoiseqbio_24_48_down_blastp_con_merge_order,
                       sig_24_48_up = mynoiseqbio_24_48_up_blastp_con_merge_order)

options(max.print=999999)
capture.output(NOIS_DE_log2FC_consensus, file = "NOISeq_sig_log2FC_consensus_v4.csv")




#TopGO on all of the NOISeq results for each time interval - 
## Mnemiopsis GO terms and TopGO code adopted from Melissa DeBiasse 
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")=
BiocManager::install("rrvgo")
BiocManager::install("org.Hs.eg.db")
library("topGO")
library(rrvgo)
#rrvgo v1.12.2
library("org.Hs.eg.db")
library(ggplot2)

#v.2.52.0
geneID2GO <- readMappings(file = "ML2.2.aa_goterms_2.txt") #Mnemiopsis genes and their GO annotations generated from interproscan 
length(geneID2GO)#8040
geneUniverse <- names(geneID2GO) #sets gene names from annotation file as the gene universe

#A function to make y-axis values integer only
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

## GO Uncut-10mins UP
U_10m_UP_char <- as.character(row.names(mynoiseqbio_U_00_up))
length(U_10m_UP_char) #91
geneList_U_10m_UP_char <- factor(as.integer(geneUniverse %in% U_10m_UP_char)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_U_10m_UP_char) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_U_10m_UP_char)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_U_10m_UP_char))#DEG with GO annotation 


## Molecular Function 
GOdata_MF_U_10m_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_U_10m_UP_char, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_U_10m_UP <- runTest(GOdata_MF_U_10m_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_U_10m_UP <- GenTable(GOdata_MF_U_10m_UP, classicFisher =  resultFisher_MF_U_10m_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar=1000)
goEnrich_MF_U_10m_UP$classicFisher <- as.numeric(goEnrich_MF_U_10m_UP$classicFisher)
goEnrich_MF_U_10m_UP <-goEnrich_MF_U_10m_UP[order(goEnrich_MF_U_10m_UP$classicFisher, decreasing = T),]
goEnrich_MF_U_10m_UP <- goEnrich_MF_U_10m_UP[goEnrich_MF_U_10m_UP$classicFisher <=0.05,]
goEnrich_MF_U_10m_UP$Term <- factor(goEnrich_MF_U_10m_UP$Term,levels = goEnrich_MF_U_10m_UP$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_U_10m_UP <- calculateSimMatrix(goEnrich_MF_U_10m_UP$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_U_10m_UP <- setNames(-log10(goEnrich_MF_U_10m_UP$classicFisher), goEnrich_MF_U_10m_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_U_10m_UP <- reduceSimMatrix(simMatrix_MF_U_10m_UP ,scores_MF_U_10m_UP,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_U_10m_UP <- unique(reducedTerms_MF_U_10m_UP$parentTerm)
reducedMatrix_MF_U_10m_UP <- subset.data.frame(goEnrich_MF_U_10m_UP,(goEnrich_MF_U_10m_UP$Term %in% parentTerms_MF_U_10m_UP))                                    
reducedMatrix_MF_U_10m_UP <- data.frame(reducedMatrix_MF_U_10m_UP,group= "up")                             


## GO Uncut-10mins DOWN
U_10m_D_char <- as.character(row.names(mynoiseqbio_U_00_down))
length(U_10m_D_char) #57
geneList_U_10m_D_char <- factor(as.integer(geneUniverse %in% U_10m_D_char)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_U_10m_D_char) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_U_10m_D_char)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_U_10m_D_char))#DEG with GO annotation 

## Molecular Function 
GOdata_MF_U_10m_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_U_10m_D_char, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_U_10m_D <- runTest(GOdata_MF_U_10m_D, algorithm="classic", statistic="fisher")
goEnrich_MF_U_10m_D <- GenTable(GOdata_MF_U_10m_D, classicFisher =  resultFisher_MF_U_10m_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_U_10m_D$classicFisher <- as.numeric(goEnrich_MF_U_10m_D$classicFisher)
goEnrich_MF_U_10m_D <-goEnrich_MF_U_10m_D[order(goEnrich_MF_U_10m_D$classicFisher, decreasing = T),]
goEnrich_MF_U_10m_D <- goEnrich_MF_U_10m_D[goEnrich_MF_U_10m_D$classicFisher <=0.05,]
goEnrich_MF_U_10m_D$Term <- factor(goEnrich_MF_U_10m_D$Term,levels = goEnrich_MF_U_10m_D$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_U_10m_D <- calculateSimMatrix(goEnrich_MF_U_10m_D$GO.ID,
                                            orgdb="org.Hs.eg.db",
                                            ont="MF",
                                            method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_U_10m_D <- setNames(-log10(goEnrich_MF_U_10m_D$classicFisher), goEnrich_MF_U_10m_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_U_10m_D <- reduceSimMatrix(simMatrix_MF_U_10m_D ,scores_MF_U_10m_D,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_U_10m_D <- unique(reducedTerms_MF_U_10m_D$parentTerm)
reducedMatrix_MF_U_10m_D <- subset.data.frame(goEnrich_MF_U_10m_D,(goEnrich_MF_U_10m_D$Term %in% parentTerms_MF_U_10m_D))                                    
reducedMatrix_MF_U_10m_D <- data.frame(reducedMatrix_MF_U_10m_D,group= "down") 


#Combine Up and Down 
MF_U_10m_Cmb <- rbind(reducedMatrix_MF_U_10m_D,reducedMatrix_MF_U_10m_UP)
MF_U_10m_Cmb$group <- factor(MF_U_10m_Cmb$group, levels = c("up","down"))

ggplot(data = MF_U_10m_Cmb, aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity",size= 0.5)+
  labs( x = "# DEG", y= "Gene Ontology",title = "Uncut-10m") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space="free_y")


## GO 10mins-1hr UP
NOIS_10m_01h_char_UP <- as.character(row.names(mynoiseqbio_00_01_up))
length(NOIS_10m_01h_char_UP) #50
geneList_NOIS_10m_01h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_10m_01h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_10m_01h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_10m_01h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_10m_01h_char_UP))#DEG with GO annotation 


## Molecular Function 
GOdata_MF_10m_01h_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_10m_01h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_10m_01h_UP <- runTest(GOdata_MF_10m_01h_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_10m_01h_UP <- GenTable(GOdata_MF_10m_01h_UP, classicFisher =  resultFisher_MF_10m_01h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar= 1000)
goEnrich_MF_10m_01h_UP$classicFisher <- as.numeric(goEnrich_MF_10m_01h_UP$classicFisher)
goEnrich_MF_10m_01h_UP <-goEnrich_MF_10m_01h_UP[order(goEnrich_MF_10m_01h_UP$classicFisher, decreasing = T),]
goEnrich_MF_10m_01h_UP <- goEnrich_MF_10m_01h_UP[goEnrich_MF_10m_01h_UP$classicFisher <=0.05,]
goEnrich_MF_10m_01h_UP$Term <- factor(goEnrich_MF_10m_01h_UP$Term,levels = goEnrich_MF_10m_01h_UP$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_10m_01h_UP <- calculateSimMatrix(goEnrich_MF_10m_01h_UP$GO.ID,
                                           orgdb="org.Hs.eg.db",
                                           ont="MF",
                                           method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_10m_01h_UP <- setNames(-log10(goEnrich_MF_10m_01h_UP$classicFisher), goEnrich_MF_10m_01h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_10m_01h_UP <- reduceSimMatrix(simMatrix_MF_10m_01h_UP ,scores_MF_10m_01h_UP,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_10m_01h_UP <- unique(reducedTerms_MF_10m_01h_UP$parentTerm)
reducedMatrix_MF_10m_01h_UP <- subset.data.frame(goEnrich_MF_10m_01h_UP,(goEnrich_MF_10m_01h_UP$Term %in% parentTerms_MF_10m_01h_UP))                                    
reducedMatrix_MF_10m_01h_UP <- data.frame(reducedMatrix_MF_10m_01h_UP,group= "up") 
reducedMatrix_MF_10m_01h_UP$Term  <- factor(reducedMatrix_MF_10m_01h_UP$Term, levels =reducedMatrix_MF_10m_01h_UP$Term)




## GO 10mins-1hr Down
NOIS_10m_01h_char_D <- as.character(row.names(mynoiseqbio_00_01_down))
length(NOIS_10m_01h_char_D) #16
geneList_NOIS_10m_01h_char_D <- factor(as.integer(geneUniverse %in% NOIS_10m_01h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_10m_01h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_10m_01h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_10m_01h_char_D))#DEG with GO annotation 

## Molecular Function 
GOdata_MF_10m_01h_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_10m_01h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_10m_01h_D <- runTest(GOdata_MF_10m_01h_D, algorithm="classic", statistic="fisher")
goEnrich_MF_10m_01h_D <- GenTable(GOdata_MF_10m_01h_D, classicFisher =  resultFisher_MF_10m_01h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar = 1000)
goEnrich_MF_10m_01h_D$classicFisher <- as.numeric(goEnrich_MF_10m_01h_D$classicFisher)
goEnrich_MF_10m_01h_D <-goEnrich_MF_10m_01h_D[order(goEnrich_MF_10m_01h_D$classicFisher, decreasing = T),]
goEnrich_MF_10m_01h_D <- goEnrich_MF_10m_01h_D[goEnrich_MF_10m_01h_D$classicFisher <=0.05,]
goEnrich_MF_10m_01h_D$Term <- factor(goEnrich_MF_10m_01h_D$Term,levels = goEnrich_MF_10m_01h_D$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_10m_01h_D <- calculateSimMatrix(goEnrich_MF_10m_01h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_10m_01h_D <- setNames(-log10(goEnrich_MF_10m_01h_D$classicFisher), goEnrich_MF_10m_01h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_10m_01h_D <- reduceSimMatrix(simMatrix_MF_10m_01h_D ,scores_MF_10m_01h_D,
                                              threshold=0.7,
                                              orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_10m_01h_D <- unique(reducedTerms_MF_10m_01h_D$parentTerm)
reducedMatrix_MF_10m_01h_D <- subset.data.frame(goEnrich_MF_10m_01h_D,(goEnrich_MF_10m_01h_D$Term %in% parentTerms_MF_10m_01h_D))                                    
reducedMatrix_MF_10m_01h_D  <- data.frame(reducedMatrix_MF_10m_01h_D,group= "down") 


#Combine Up and Down 
MF_10m_01h_Cmb <- rbind(reducedMatrix_MF_10m_01h_UP,reducedMatrix_MF_10m_01h_D)
MF_10m_01h_Cmb$group <- factor(MF_10m_01h_Cmb$group, levels = c("up","down"))

ggplot(data =  MF_10m_01h_Cmb, aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.5)+
  labs( x = "# DEG", y= "Gene Ontology", title = "10m-1h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")

## GO 1hr-3hr UP
NOIS_01h_03h_char_UP <- as.character(row.names(mynoiseqbio_01_03_up))
length(NOIS_01h_03h_char_UP) #94
geneList_NOIS_01h_03h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_01h_03h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_01h_03h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_01h_03h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_01h_03h_char_UP))#DEG with GO annotation 

## Molecular Function 
GOdata_MF_01h_03h_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_01h_03h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_01h_03h_UP <- runTest(GOdata_MF_01h_03h_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_01h_03h_UP <- GenTable(GOdata_MF_01h_03h_UP, classicFisher =  resultFisher_MF_01h_03h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar=100)
goEnrich_MF_01h_03h_UP$classicFisher <- as.numeric(goEnrich_MF_01h_03h_UP$classicFisher)
goEnrich_MF_01h_03h_UP <-goEnrich_MF_01h_03h_UP[order(goEnrich_MF_01h_03h_UP$classicFisher, decreasing = T),]
goEnrich_MF_01h_03h_UP <- goEnrich_MF_01h_03h_UP[goEnrich_MF_01h_03h_UP$classicFisher <=0.05,]
goEnrich_MF_01h_03h_UP$Term <- factor(goEnrich_MF_01h_03h_UP$Term,levels = goEnrich_MF_01h_03h_UP$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_01h_03h_UP <- calculateSimMatrix(goEnrich_MF_01h_03h_UP$GO.ID,
                                            orgdb="org.Hs.eg.db",
                                            ont="MF",
                                            method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_01h_03h_UP <- setNames(-log10(goEnrich_MF_01h_03h_UP$classicFisher), goEnrich_MF_01h_03h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_01h_03h_UP <- reduceSimMatrix(simMatrix_MF_01h_03h_UP ,scores_MF_01h_03h_UP,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_01h_03h_UP <- unique(reducedTerms_MF_01h_03h_UP$parentTerm)
reducedMatrix_MF_01h_03h_UP <- subset.data.frame(goEnrich_MF_01h_03h_UP,(goEnrich_MF_01h_03h_UP$Term %in% parentTerms_MF_01h_03h_UP))                                    
reducedMatrix_MF_01h_03h_UP <- data.frame(reducedMatrix_MF_01h_03h_UP,group= "up")                             




## GO 1hr-3hr Down
NOIS_01h_03h_char_D <- as.character(row.names(mynoiseqbio_01_03_down))
length(NOIS_01h_03h_char_D) #340
geneList_NOIS_01h_03h_char_D <- factor(as.integer(geneUniverse %in% NOIS_01h_03h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_01h_03h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_01h_03h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_01h_03h_char_D))#DEG with GO annotation 


## Molecular Function 
GOdata_MF_01h_03h_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_01h_03h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_01h_03h_D <- runTest(GOdata_MF_01h_03h_D, algorithm="classic", statistic="fisher")
goEnrich_MF_01h_03h_D <- GenTable(GOdata_MF_01h_03h_D, classicFisher =  resultFisher_MF_01h_03h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_01h_03h_D$classicFisher <- as.numeric(goEnrich_MF_01h_03h_D$classicFisher)
goEnrich_MF_01h_03h_D <-goEnrich_MF_01h_03h_D[order(goEnrich_MF_01h_03h_D$classicFisher, decreasing = T),]
goEnrich_MF_01h_03h_D <- goEnrich_MF_01h_03h_D[goEnrich_MF_01h_03h_D$classicFisher <=0.05,]
goEnrich_MF_01h_03h_D$Term <- factor(goEnrich_MF_01h_03h_D$Term,levels = goEnrich_MF_01h_03h_D$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_01h_03h_D <- calculateSimMatrix(goEnrich_MF_01h_03h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_01h_03h_D <- setNames(-log10(goEnrich_MF_01h_03h_D$classicFisher), goEnrich_MF_01h_03h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_01h_03h_D <- reduceSimMatrix(simMatrix_MF_01h_03h_D ,scores_MF_01h_03h_D,
                                              threshold=0.7,
                                              orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_01h_03h_D <- unique(reducedTerms_MF_01h_03h_D$parentTerm)
reducedMatrix_MF_01h_03h_D <- subset.data.frame(goEnrich_MF_01h_03h_D,(goEnrich_MF_01h_03h_D$Term %in% parentTerms_MF_01h_03h_D))                                    
reducedMatrix_MF_01h_03h_D <- data.frame(reducedMatrix_MF_01h_03h_D,group= "down")                             

#Combine Up and Down 
MF_01h_03h_Cmb <- rbind(reducedMatrix_MF_01h_03h_UP,reducedMatrix_MF_01h_03h_D)
MF_01h_03h_Cmb$group <- factor(MF_01h_03h_Cmb$group, levels = c("up","down"))

ggplot(data = MF_01h_03h_Cmb , aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.9)+
  labs( x = "# DEG", y= "Gene Ontology",title = "1h-3h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")


## GO 3hr-6hr UP
NOIS_03h_06h_char_UP <- as.character(row.names(mynoiseqbio_03_06_up))
length(NOIS_03h_06h_char_UP) #3843
geneList_NOIS_03h_06h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_03h_06h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_03h_06h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_03h_06h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_03h_06h_char_UP))#DEG with GO annotation 

## Molecular Function 
GOdata_MF_03h_06h_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_03h_06h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_03h_06h_UP <- runTest(GOdata_MF_03h_06h_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_03h_06h_UP <- GenTable(GOdata_MF_03h_06h_UP, classicFisher =  resultFisher_MF_03h_06h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=100)
goEnrich_MF_03h_06h_UP$classicFisher <- as.numeric(goEnrich_MF_03h_06h_UP$classicFisher)
goEnrich_MF_03h_06h_UP <-goEnrich_MF_03h_06h_UP[order(goEnrich_MF_03h_06h_UP$classicFisher, decreasing = T),]
goEnrich_MF_03h_06h_UP <- goEnrich_MF_03h_06h_UP[goEnrich_MF_03h_06h_UP$classicFisher <=0.05,]
goEnrich_MF_03h_06h_UP$Term <- factor(goEnrich_MF_03h_06h_UP$Term,levels = goEnrich_MF_03h_06h_UP$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_03h_06h_UP <- calculateSimMatrix(goEnrich_MF_03h_06h_UP$GO.ID,
                                             orgdb="org.Hs.eg.db",
                                             ont="MF",
                                             method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_03h_06h_UP <- setNames(-log10(goEnrich_MF_03h_06h_UP$classicFisher), goEnrich_MF_03h_06h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_03h_06h_UP <- reduceSimMatrix(simMatrix_MF_03h_06h_UP ,scores_MF_03h_06h_UP,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_03h_06h_UP <- unique(reducedTerms_MF_03h_06h_UP$parentTerm)
reducedMatrix_MF_03h_06h_UP <- subset.data.frame(goEnrich_MF_03h_06h_UP,(goEnrich_MF_03h_06h_UP$Term %in% parentTerms_MF_03h_06h_UP))                                    
reducedMatrix_MF_03h_06h_UP <- data.frame(reducedMatrix_MF_03h_06h_UP,group= "up")                             


## GO 3hr-6hr DOWN
NOIS_03h_06h_char_D <- as.character(row.names(mynoiseqbio_03_06_down))
length(NOIS_03h_06h_char_D) #2250
geneList_NOIS_03h_06h_char_D <- factor(as.integer(geneUniverse %in% NOIS_03h_06h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_03h_06h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_03h_06h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_03h_06h_char_D))#DEG with GO annotation 

## Molecular Function 
GOdata_MF_03h_06h_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_03h_06h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_03h_06h_D <- runTest(GOdata_MF_03h_06h_D, algorithm="classic", statistic="fisher")
goEnrich_MF_03h_06h_D <- GenTable(GOdata_MF_03h_06h_D, classicFisher =  resultFisher_MF_03h_06h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar = 1000)
goEnrich_MF_03h_06h_D$classicFisher <- as.numeric(goEnrich_MF_03h_06h_D$classicFisher)
goEnrich_MF_03h_06h_D <-goEnrich_MF_03h_06h_D[order(goEnrich_MF_03h_06h_D$classicFisher, decreasing = T),]
goEnrich_MF_03h_06h_D <- goEnrich_MF_03h_06h_D[goEnrich_MF_03h_06h_D$classicFisher <=0.05,]
goEnrich_MF_03h_06h_D$Term <- factor(goEnrich_MF_03h_06h_D$Term,levels = goEnrich_MF_03h_06h_D$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_03h_06h_D <- calculateSimMatrix(goEnrich_MF_03h_06h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_03h_06h_D<- setNames(-log10(goEnrich_MF_03h_06h_D$classicFisher), goEnrich_MF_03h_06h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_03h_06h_D <- reduceSimMatrix(simMatrix_MF_03h_06h_D ,scores_MF_03h_06h_D,
                                              threshold=0.7,
                                              orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_03h_06h_D <- unique(reducedTerms_MF_03h_06h_D$parentTerm)
reducedMatrix_MF_03h_06h_D <- subset.data.frame(goEnrich_MF_03h_06h_D,(goEnrich_MF_03h_06h_D$Term %in% parentTerms_MF_03h_06h_D))                                    
reducedMatrix_MF_03h_06h_D <- data.frame(reducedMatrix_MF_03h_06h_D,group= "down")                             

#Combine Up and Down 
MF_03h_06h_Cmb <- rbind(reducedMatrix_MF_03h_06h_UP,reducedMatrix_MF_03h_06h_D)
MF_03h_06h_Cmb$group <- factor(MF_03h_06h_Cmb$group, levels = c("up","down"))

ggplot(data = MF_03h_06h_Cmb  , aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.9)+
  labs( x = "# DEG", y= "Gene Ontology",title= "3h-6h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")


## GO 6hr-12hr UP
NOIS_06h_12h_char_UP <- as.character(row.names(mynoiseqbio_06_12_up))
length(NOIS_06h_12h_char_UP) #58
geneList_NOIS_06h_12h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_06h_12h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_06h_12h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_06h_12h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_06h_12h_char_UP))#DEG with GO annotation

## Molecular Function 
GOdata_06h_12h_MF_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_06h_12h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_06h_12h_UP <- runTest(GOdata_06h_12h_MF_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_06h_12h_UP <- GenTable(GOdata_06h_12h_MF_UP, classicFisher =  resultFisher_MF_06h_12h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_06h_12h_UP$classicFisher <- as.numeric(goEnrich_MF_06h_12h_UP$classicFisher)
goEnrich_MF_06h_12h_UP <-goEnrich_MF_06h_12h_UP[order(goEnrich_MF_06h_12h_UP$classicFisher, decreasing = T),]
goEnrich_MF_06h_12h_UP <- goEnrich_MF_06h_12h_UP[goEnrich_MF_06h_12h_UP$classicFisher <=0.05,]
goEnrich_MF_06h_12h_UP$Term <- factor(goEnrich_MF_06h_12h_UP$Term,levels = goEnrich_MF_06h_12h_UP$Term)

# using REVIGO, calculate similarity matrix
simMatrix_MF_06h_12h_UP <- calculateSimMatrix(goEnrich_MF_06h_12h_UP$GO.ID,
                                             orgdb="org.Hs.eg.db",
                                             ont="MF",
                                             method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_06h_12h_UP<- setNames(-log10(goEnrich_MF_06h_12h_UP$classicFisher), goEnrich_MF_06h_12h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_06h_12h_UP <- reduceSimMatrix(simMatrix_MF_06h_12h_UP ,scores_MF_06h_12h_UP,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_06h_12h_UP<- unique(reducedTerms_MF_06h_12h_UP$parentTerm)
reducedMatrix_MF_06h_12h_UP <- subset.data.frame(goEnrich_MF_06h_12h_UP,(goEnrich_MF_06h_12h_UP$Term %in% parentTerms_MF_06h_12h_UP))                                    
reducedMatrix_MF_06h_12h_UP <- data.frame(reducedMatrix_MF_06h_12h_UP,group= "up")                             


## GO 6-12hr DOWN
NOIS_06h_12h_char_D <- as.character(row.names(mynoiseqbio_06_12_down))
length(NOIS_06h_12h_char_D) #79
geneList_NOIS_06h_12h_char_D <- factor(as.integer(geneUniverse %in% NOIS_06h_12h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_06h_12h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_06h_12h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_06h_12h_char_D))#DEG with GO annotation

## Molecular Function 
GOdata_06h_12h_MF_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_06h_12h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_06h_12h_D <- runTest(GOdata_06h_12h_MF_D, algorithm="classic", statistic="fisher")
goEnrich_MF_06h_12h_D <- GenTable(GOdata_06h_12h_MF_D, classicFisher =  resultFisher_MF_06h_12h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10, numChar = 1000)
goEnrich_MF_06h_12h_D$classicFisher <- as.numeric(goEnrich_MF_06h_12h_D$classicFisher)
goEnrich_MF_06h_12h_D <-goEnrich_MF_06h_12h_D[order(goEnrich_MF_06h_12h_D$classicFisher, decreasing = T),]
goEnrich_MF_06h_12h_D <- goEnrich_MF_06h_12h_D[goEnrich_MF_06h_12h_D$classicFisher <=0.05,]
goEnrich_MF_06h_12h_D$Term <- factor(goEnrich_MF_06h_12h_D$Term,levels = goEnrich_MF_06h_12h_D$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_06h_12h_D <- calculateSimMatrix(goEnrich_MF_06h_12h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_06h_12h_D<- setNames(-log10(goEnrich_MF_06h_12h_D$classicFisher), goEnrich_MF_06h_12h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_06h_12h_D <- reduceSimMatrix(simMatrix_MF_06h_12h_D ,scores_MF_06h_12h_D,
                                              threshold=0.7,
                                              orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_06h_12h_D<- unique(reducedTerms_MF_06h_12h_D$parentTerm)
reducedMatrix_MF_06h_12h_D <- subset.data.frame(goEnrich_MF_06h_12h_D,(goEnrich_MF_06h_12h_D$Term %in% parentTerms_MF_06h_12h_D))                                    
reducedMatrix_MF_06h_12h_D <- data.frame(reducedMatrix_MF_06h_12h_D,group= "down")                             


#Combine Up and Down 
MF_06h_12h_Cmb <- rbind(reducedMatrix_MF_06h_12h_UP,reducedMatrix_MF_06h_12h_D)
MF_06h_12h_Cmb$group <- factor(MF_06h_12h_Cmb$group, levels = c("up","down"))

ggplot(data = MF_06h_12h_Cmb , aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.9)+
  labs( x = "# DEG", y= "Gene Ontology",title= "6h-12h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")



## GO 12hr-24hr UP
NOIS_12h_24h_char_UP <- as.character(row.names(mynoiseqbio_12_24_up))
length(NOIS_12h_24h_char_UP) #45
geneList_NOIS_12h_24h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_12h_24h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_12h_24h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_12h_24h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_12h_24h_char_UP))#DEG with GO annotation

## Molecular Function 
GOdata_12h_24h_MF_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_12h_24h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_12h_24h_UP <- runTest(GOdata_12h_24h_MF_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_12h_24h_UP <- GenTable(GOdata_12h_24h_MF_UP, classicFisher =  resultFisher_MF_12h_24h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_12h_24h_UP$classicFisher <- as.numeric(goEnrich_MF_12h_24h_UP$classicFisher)
goEnrich_MF_12h_24h_UP <-goEnrich_MF_12h_24h_UP[order(goEnrich_MF_12h_24h_UP$classicFisher, decreasing = T),]
goEnrich_MF_12h_24h_UP <- goEnrich_MF_12h_24h_UP[goEnrich_MF_12h_24h_UP$classicFisher <=0.05,]
goEnrich_MF_12h_24h_UP$Term <- factor(goEnrich_MF_12h_24h_UP$Term,levels = goEnrich_MF_12h_24h_UP$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_12h_24h_UP <- calculateSimMatrix(goEnrich_MF_12h_24h_UP$GO.ID,
                                             orgdb="org.Hs.eg.db",
                                             ont="MF",
                                             method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_12h_24h_UP<- setNames(-log10(goEnrich_MF_12h_24h_UP$classicFisher), goEnrich_MF_12h_24h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_12h_24h_UP<- reduceSimMatrix(simMatrix_MF_12h_24h_UP ,scores_MF_12h_24h_UP,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_12h_24h_UP<- unique(reducedTerms_MF_12h_24h_UP$parentTerm)
reducedMatrix_MF_12h_24h_UP <- subset.data.frame(goEnrich_MF_12h_24h_UP,(goEnrich_MF_12h_24h_UP$Term %in% parentTerms_MF_12h_24h_UP))                                    
reducedMatrix_MF_12h_24h_UP <- data.frame(reducedMatrix_MF_12h_24h_UP,group= "up")                             



## GO 12hr-24hr DOWN
NOIS_12h_24h_char_D <- as.character(row.names(mynoiseqbio_12_24_down))
length(NOIS_12h_24h_char_D) #78
geneList_NOIS_12h_24h_char_D <- factor(as.integer(geneUniverse %in% NOIS_12h_24h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_12h_24h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_12h_24h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_12h_24h_char_D))#DEG with GO annotation

## Molecular Function 
GOdata_12h_24h_MF_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_12h_24h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_12h_24h_D <- runTest(GOdata_12h_24h_MF_D, algorithm="classic", statistic="fisher")
goEnrich_MF_12h_24h_D <- GenTable(GOdata_12h_24h_MF_D, classicFisher =  resultFisher_MF_12h_24h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_12h_24h_D$classicFisher <- as.numeric(goEnrich_MF_12h_24h_D$classicFisher)
goEnrich_MF_12h_24h_D <-goEnrich_MF_12h_24h_D[order(goEnrich_MF_12h_24h_D$classicFisher, decreasing = T),]
goEnrich_MF_12h_24h_D <- goEnrich_MF_12h_24h_D[goEnrich_MF_12h_24h_D$classicFisher <=0.05,]
goEnrich_MF_12h_24h_D$Term <- factor(goEnrich_MF_12h_24h_D$Term,levels = goEnrich_MF_12h_24h_D$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_12h_24h_D <- calculateSimMatrix(goEnrich_MF_12h_24h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_12h_24h_D <- setNames(-log10(goEnrich_MF_12h_24h_D$classicFisher), goEnrich_MF_12h_24h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_12h_24h_D <- reduceSimMatrix(simMatrix_MF_12h_24h_D ,scores_MF_12h_24h_D,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_12h_24h_D<- unique(reducedTerms_MF_12h_24h_D$parentTerm)
reducedMatrix_MF_12h_24h_D <- subset.data.frame(goEnrich_MF_12h_24h_D,(goEnrich_MF_12h_24h_D$Term %in% parentTerms_MF_12h_24h_D))                                    
reducedMatrix_MF_12h_24h_D <- data.frame(reducedMatrix_MF_12h_24h_D,group= "down")                             

#Combine Up and Down 
MF_12h_24h_Cmb <- rbind(reducedMatrix_MF_12h_24h_D,reducedMatrix_MF_12h_24h_UP)
MF_12h_24h_Cmb$group <- factor(MF_12h_24h_Cmb$group, levels = c("up","down"))

ggplot(data = MF_12h_24h_Cmb , aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.5)+
  labs( x = "# DEG", y= "Gene Ontology", title = "12h-24h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")



## GO 24hr-48hr UP
NOIS_24h_48h_char_UP <- as.character(row.names(mynoiseqbio_24_48_up))
length(NOIS_24h_48h_char_UP) #108
geneList_NOIS_24h_48h_char_UP <- factor(as.integer(geneUniverse %in% NOIS_24h_48h_char_UP)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_24h_48h_char_UP) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_24h_48h_char_UP)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_24h_48h_char_UP))#DEG with GO annotation

## Molecular Function 
GOdata_24h_48h_MF_UP <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_24h_48h_char_UP, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_24h_48h_UP <- runTest(GOdata_24h_48h_MF_UP, algorithm="classic", statistic="fisher")
goEnrich_MF_24h_48h_UP <- GenTable(GOdata_24h_48h_MF_UP, classicFisher =  resultFisher_MF_24h_48h_UP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_24h_48h_UP$classicFisher <- as.numeric(goEnrich_MF_24h_48h_UP$classicFisher)
goEnrich_MF_24h_48h_UP <-goEnrich_MF_24h_48h_UP[order(goEnrich_MF_24h_48h_UP$classicFisher, decreasing = T),]
goEnrich_MF_24h_48h_UP <- goEnrich_MF_24h_48h_UP[goEnrich_MF_24h_48h_UP$classicFisher <=0.05,]
goEnrich_MF_24h_48h_UP$Term <- factor(goEnrich_MF_24h_48h_UP$Term,levels = goEnrich_MF_24h_48h_UP$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_24h_48h_UP <- calculateSimMatrix(goEnrich_MF_24h_48h_UP$GO.ID,
                                             orgdb="org.Hs.eg.db",
                                             ont="MF",
                                             method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_24h_48h_UP <- setNames(-log10(goEnrich_MF_24h_48h_UP$classicFisher), goEnrich_MF_24h_48h_UP$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_24h_48h_UP<- reduceSimMatrix(simMatrix_MF_24h_48h_UP ,scores_MF_24h_48h_UP,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_24h_48h_UP<- unique(reducedTerms_MF_24h_48h_UP$parentTerm)
reducedMatrix_MF_24h_48h_UP <- subset.data.frame(goEnrich_MF_24h_48h_UP,(goEnrich_MF_24h_48h_UP$Term %in% parentTerms_MF_24h_48h_UP))                                    
reducedMatrix_MF_24h_48h_UP <- data.frame(reducedMatrix_MF_24h_48h_UP,group= "up")                             


## GO 24hr-48hr DOWN
NOIS_24h_48h_char_D <- as.character(row.names(mynoiseqbio_24_48_down))
length(NOIS_24h_48h_char_D) #62
geneList_NOIS_24h_48h_char_D <- factor(as.integer(geneUniverse %in% NOIS_24h_48h_char_D)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_NOIS_24h_48h_char_D) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)#8040
length(Filter(function(x) x == 0, geneList_NOIS_24h_48h_char_D)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_NOIS_24h_48h_char_D))#DEG with GO annotation

## Molecular Function 
GOdata_24h_48h_MF_D <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_NOIS_24h_48h_char_D, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_24h_48h_D <- runTest(GOdata_24h_48h_MF_D, algorithm="classic", statistic="fisher")
goEnrich_MF_24h_48h_D <- GenTable(GOdata_24h_48h_MF_D, classicFisher =  resultFisher_MF_24h_48h_D, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10,numChar=1000)
goEnrich_MF_24h_48h_D$classicFisher <- as.numeric(goEnrich_MF_24h_48h_D$classicFisher)
goEnrich_MF_24h_48h_D <-goEnrich_MF_24h_48h_D[order(goEnrich_MF_24h_48h_D$classicFisher, decreasing = T),]
goEnrich_MF_24h_48h_D <- goEnrich_MF_24h_48h_D[goEnrich_MF_24h_48h_D$classicFisher <=0.05,]
goEnrich_MF_24h_48h_D$Term <- factor(goEnrich_MF_24h_48h_D$Term,levels = goEnrich_MF_24h_48h_D$Term)


# using REVIGO, calculate similarity matrix
simMatrix_MF_24h_48h_D <- calculateSimMatrix(goEnrich_MF_24h_48h_D$GO.ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="MF",
                                              method="Rel")

#Set scores from classic Fisher p-value, transforming with -log10
scores_MF_24h_48h_D <- setNames(-log10(goEnrich_MF_24h_48h_D$classicFisher), goEnrich_MF_24h_48h_D$GO.ID)

#reduce terms using Sim matrix
reducedTerms_MF_24h_48h_D<- reduceSimMatrix(simMatrix_MF_24h_48h_D ,scores_MF_24h_48h_D,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")

#Vector with unique parent terms to subset reduced matrix
parentTerms_MF_24h_48h_D<- unique(reducedTerms_MF_24h_48h_D$parentTerm)
reducedMatrix_MF_24h_48h_D <- subset.data.frame(goEnrich_MF_24h_48h_D,(goEnrich_MF_24h_48h_D$Term %in% parentTerms_MF_24h_48h_D))                                    
reducedMatrix_MF_24h_48h_D <- data.frame(reducedMatrix_MF_24h_48h_D,group= "down")                             

#Combine Up and Down 
MF_24h_48h_Cmb <- rbind(reducedMatrix_MF_24h_48h_D,reducedMatrix_MF_24h_48h_UP)
MF_24h_48h_Cmb$group <- factor(MF_24h_48h_Cmb$group, levels = c("up","down"))

ggplot(data = MF_24h_48h_Cmb , aes(x = Significant , y = Term, fill= log10(classicFisher))) + 
  scale_fill_gradient(low = "black", high = "darkgrey") +
  geom_bar(stat = "identity", width = 0.9)+
  labs( x = "# DEG", y= "Gene Ontology",title = "24h-48h") +
  scale_x_continuous(breaks = integer_breaks())+
  theme(axis.text = element_text(size = 20)) +
  facet_grid(group~.,scales = "free_y",space = "free_y")



#Looking through GO results from each time interval in NOISeq 

#For BLAST
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

##Uncut-10m-UP -- Molecular function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_U_10m_UP

#gathering BLAST for export 
GO0008061_U_10m_UP_DEG_blastp
GO0003700_U_10m_UP_DEG_blastp
GO0004553_U_10m_UP_DEG_blastp
GO0038024_U_10m_UP_DEG_blastp
GO0098599_U_10m_UP_DEG_blastp


#GO:0008061	-- chitin binding
GO0008061_U_10m_UP_scores <- scoresInTerm(GOdata_MF_U_10m_UP,"GO:0008061", use.names = T)
GO0008061_U_10m_UP_scores_df<- as.data.frame(GO0008061_U_10m_UP_scores)
GO0008061_U_10m_UP_DEG <- subset.data.frame(GO0008061_U_10m_UP_scores_df,GO0008061_U_10m_UP_scores_df == 2 )
length(GO0008061_U_10m_UP_DEG$GO.0008061) # 2 DEG
GO0008061_U_10m_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0008061_U_10m_UP_DEG,by='row.names')


#GO:0003700-- DNA-binding tf activity 
GO0003700_U_10m_UP_scores <- scoresInTerm(GOdata_MF_U_10m_UP ,"GO:0003700", use.names = T)
GO0003700_U_10m_UP_scores_df<- as.data.frame(GO0003700_U_10m_UP_scores)
GO0003700_U_10m_UP_DEG <- subset.data.frame(GO0003700_U_10m_UP_scores_df,GO0003700_U_10m_UP_scores_df== 2 )
length(GO0003700_U_10m_UP_DEG$GO.0003700) # 3 DEG
GO0003700_U_10m_UP_DEG_blastp<- merge(HumRef2020ML_df ,GO0003700_U_10m_UP_DEG,by='row.names')

# GO:0004553 hydrolase activity, hydrolyzing O-glycosyl compounds 
GO0004553_U_10m_UP_scores <- scoresInTerm(GOdata_MF_U_10m_UP ,"GO:0004553", use.names = T)
GO0004553_U_10m_UP_scores_df<- as.data.frame(GO0004553_U_10m_UP_scores)
GO0004553_U_10m_UP_DEG <- subset.data.frame(GO0004553_U_10m_UP_scores_df,GO0004553_U_10m_UP_scores_df== 2 )
length(GO0004553_U_10m_UP_DEG$GO.0004553) # 2 DEG
GO0004553_U_10m_UP_DEG_blastp<- merge(HumRef2020ML_df ,GO0004553_U_10m_UP_DEG,by='row.names')

#GO:0038024 -- cargo receptor activity
GO0038024_U_10m_UP_scores <- scoresInTerm(GOdata_MF_U_10m_UP ,"GO:0038024", use.names = T)
GO0038024_U_10m_UP_scores_df<- as.data.frame(GO0038024_U_10m_UP_scores)
GO0038024_U_10m_UP_DEG <- subset.data.frame(GO0038024_U_10m_UP_scores_df,GO0038024_U_10m_UP_scores_df== 2 )
length(GO0038024_U_10m_UP_DEG$GO.0038024) # 2 DEG
GO0038024_U_10m_UP_DEG_blastp<- merge(HumRef2020ML_df ,GO0038024_U_10m_UP_DEG,by='row.names')

#GO:0098599--palmitoyl hydrolase activity
GO0098599_U_10m_UP_scores <- scoresInTerm(GOdata_MF_U_10m_UP ,"GO:0098599", use.names = T)
GO0098599_U_10m_UP_scores_df<- as.data.frame(GO0098599_U_10m_UP_scores)
GO0098599_U_10m_UP_DEG <- subset.data.frame(GO0098599_U_10m_UP_scores_df,GO0098599_U_10m_UP_scores_df== 2 )
length(GO0098599_U_10m_UP_DEG$GO.0098599) # 1 DEG
GO0098599_U_10m_UP_DEG_blastp<- merge(HumRef2020ML_df ,GO0098599_U_10m_UP_DEG,by='row.names')

##Uncut-10m-DOWN -- Molecular function
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_U_10m_D

#gathering BLAST for export 
GO0005509_U_10m_D_DEG_blastp
GO0016880_U_10m_D_DEG_blastp
GO0005198_U_10m_D_DEG_blastp
GO0008519_U_10m_D_DEG_blastp

# GO:0005509 -- calcium ion binding  
GO0005509_U_10m_D_scores <- scoresInTerm(GOdata_MF_U_10m_D ,"GO:0005509", use.names = T)
GO0005509_U_10m_D_scores_df <- as.data.frame(GO0005509_U_10m_D_scores)
GO0005509_U_10m_D_DEG <- subset.data.frame(GO0005509_U_10m_D_scores_df,GO0005509_U_10m_D_scores_df== 2 )
length(GO0005509_U_10m_D_DEG$GO.0005509) # 5 DEG
GO0005509_U_10m_D_DEG_blastp<- merge(HumRef2020ML_df ,GO0005509_U_10m_D_DEG,by='row.names')

#GO:0016880     acid-ammonia (or amide) ligase activity   
GO0016880_U_10m_D_scores <- scoresInTerm(GOdata_MF_U_10m_D ,"GO:0016880", use.names = T)
GO0016880_U_10m_D_scores_df <- as.data.frame(GO0016880_U_10m_D_scores)
GO0016880_U_10m_D_DEG <- subset.data.frame(GO0016880_U_10m_D_scores_df,GO0016880_U_10m_D_scores_df== 2 )
length(GO0016880_U_10m_D_DEG$GO.0016880) # 5 DEG
GO0016880_U_10m_D_DEG_blastp<- merge(HumRef2020ML_df ,GO0016880_U_10m_D_DEG,by='row.names')

#GO:0008519 ammonium transmembrane transporter activity
GO0008519_U_10m_D_scores <- scoresInTerm(GOdata_MF_U_10m_D ,"GO:0008519", use.names = T)
GO0008519_U_10m_D_scores_df <- as.data.frame(GO0008519_U_10m_D_scores)
GO0008519_U_10m_D_DEG <- subset.data.frame(GO0008519_U_10m_D_scores_df,GO0008519_U_10m_D_scores_df== 2 )
length(GO0008519_U_10m_D_DEG$GO.0008519) # 5 DEG
GO0008519_U_10m_D_DEG_blastp<- merge(HumRef2020ML_df ,GO0008519_U_10m_D_DEG,by='row.names')

#GO:0005198 -- structural molecule activity
GO0005198_U_10m_D_scores <- scoresInTerm(GOdata_MF_U_10m_D ,"GO:0005198", use.names = T)
GO0005198_U_10m_D_scores_df <- as.data.frame(GO0005198_U_10m_D_scores)
GO0005198_U_10m_D_DEG <- subset.data.frame(GO0005198_U_10m_D_scores_df,GO0005198_U_10m_D_scores_df== 2 )
length(GO0005198_U_10m_D_DEG$GO.0005198) # 2 DEG
GO0005198_U_10m_D_DEG_blastp<- merge(HumRef2020ML_df ,GO0005198_U_10m_D_DEG,by='row.names')


#10m-1h-UP -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_10m_01h_UP

#gathering BLAST for export 
GO0003700_10m_01h_UP_DEG_blastp
GO0001067_10m_01h_UP_DEG_blastp
GO0004222_10m_01h_UP_DEG_blastp

#GO:0003700-- DNA-binding tf activity 
GO0003700_10m_01h_UP_scores <- scoresInTerm(GOdata_MF_10m_01h_UP,"GO:0003700", use.names = T)
GO0003700_10m_01h_UP_scores_df<- as.data.frame(GO0003700_10m_01h_UP_scores)
GO0003700_10m_01h_UP_DEG <- subset.data.frame(GO0003700_10m_01h_UP_scores_df,GO0003700_10m_01h_UP_scores_df== 2 )
length(GO0003700_10m_01h_UP_DEG$GO.0003700) # 4 DEG
GO0003700_10m_01h_UP_DEG_blastp<- merge(HumRef2020ML_df ,GO0003700_10m_01h_UP_DEG,by='row.names')

#GO:0001067--  transcription regulatory region nucleic acid binding 
GO0001067_10m_01h_UP_scores <- scoresInTerm(GOdata_MF_10m_01h_UP,"GO:0001067", use.names = T)
GO0001067_10m_01h_UP_scores_df<- as.data.frame(GO0001067_10m_01h_UP_scores)
GO0001067_10m_01h_UP_DEG <- subset.data.frame(GO0001067_10m_01h_UP_scores_df,GO0001067_10m_01h_UP_scores_df== 2 )
length(GO0001067_10m_01h_UP_DEG$GO.0001067) # 1 DEG 
GO0001067_10m_01h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0001067_10m_01h_UP_DEG,by='row.names')


#GO:0004222-- metalloendopeptidase activity 
GO0004222_10m_01h_UP_scores <- scoresInTerm(GOdata_MF_10m_01h_UP,"GO:0004222", use.names = T)
GO0004222_10m_01h_UP_scores_df<- as.data.frame(GO0004222_10m_01h_UP_scores)
GO0004222_10m_01h_UP_DEG <- subset.data.frame(GO0004222_10m_01h_UP_scores_df,GO0004222_10m_01h_UP_scores_df== 2 )
length(GO0004222_10m_01h_UP_DEG$GO.0004222) # 2 DEG 
GO0004222_10m_01h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0004222_10m_01h_UP_DEG,by='row.names')

#10m-1h-DOWN -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_10m_01h_D
#gathering BLAST for export
GO0003700_10m_01h_D_DEG_blastp


#GO:0003700-- DNA-binding tf activity 
GO0003700_10m_01h_D_scores <- scoresInTerm(GOdata_MF_10m_01h_D,"GO:0003700", use.names = T)
GO0003700_10m_01h_D_scores_df<- as.data.frame(GO0003700_10m_01h_D_scores)
GO0003700_10m_01h_D_DEG <- subset.data.frame(GO0003700_10m_01h_D_scores_df,GO0003700_10m_01h_D_scores_df== 2 )
length(GO0003700_10m_01h_D_DEG$GO.0003700) # 2 DEG
GO0003700_10m_01h_D_DEG_blastp<- merge(HumRef2020ML_df ,GO0003700_10m_01h_D_DEG,by='row.names')


#1h-3h--UP -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_01h_03h_UP
#gathering BLAST for export
GO0005198_01h_03h_UP_DEG_blastp
GO0003924_01h_03h_UP_DEG_blastp

#GO:0005198 -- structural molecule activity  
GO0005198_01h_03h_UP_scores <- scoresInTerm(GOdata_MF_01h_03h_UP,"GO:0005198", use.names = T)
GO0005198_01h_03h_UP_scores_df<- as.data.frame(GO0005198_01h_03h_UP_scores)
GO0005198_01h_03h_UP_DEG <- subset.data.frame(GO0005198_01h_03h_UP_scores_df,GO0005198_01h_03h_UP_scores_df== 2 )
length(GO0005198_01h_03h_UP_DEG$GO.0005198) # 12 DEG 
GO0005198_01h_03h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005198_01h_03h_UP_DEG,by='row.names')


#GO:0003924 GTPase activity  
GO0003924_01h_03h_UP_scores <- scoresInTerm(GOdata_MF_01h_03h_UP,"GO:0003924", use.names = T)
GO0003924_01h_03h_UP_scores_df<- as.data.frame(GO0003924_01h_03h_UP_scores)
GO0003924_01h_03h_UP_DEG <- subset.data.frame(GO0003924_01h_03h_UP_scores_df,GO0003924_01h_03h_UP_scores_df== 2 )
length(GO0003924_01h_03h_UP_DEG$GO.0003924) # 6 DEG 
GO0003924_01h_03h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0003924_01h_03h_UP_DEG,by='row.names')


#1h-3h-D-- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_01h_03h_D
#gathering BLAST for export
GO0003723_01h_03h_D_DEG_blastp 
GO0045182_01h_03h_D_DEG_blastp
GO0003700_01h_03h_D_DEG_blastp
GO0140102_01h_03h_D_DEG_blastp
GO0015036_01h_03h_D_DEG_blastp 

#GO:0003723  -- RNA-binding 
GO0003723_01h_03h_D_scores <- scoresInTerm(GOdata_MF_01h_03h_D,"GO:0003723", use.names = T)
GO0003723_01h_03h_D_scores_df<- as.data.frame(GO0003723_01h_03h_D_scores )
GO0003723_01h_03h_D_DEG <- subset.data.frame(GO0003723_01h_03h_D_scores_df,GO0003723_01h_03h_D_scores_df== 2 )
length(GO0003723_01h_03h_D_DEG$GO.0003723) # 15 DEG 
GO0003723_01h_03h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0003723_01h_03h_D_DEG,by='row.names')

#GO:0045182 -- translation regulator activity 
GO0045182_01h_03h_D_scores <- scoresInTerm(GOdata_MF_01h_03h_D,"GO:0045182", use.names = T)
GO0045182_01h_03h_D_scores_df<- as.data.frame(GO0045182_01h_03h_D_scores )
GO0045182_01h_03h_D_DEG <- subset.data.frame(GO0045182_01h_03h_D_scores_df,GO0045182_01h_03h_D_scores_df== 2 )
length(GO0045182_01h_03h_D_DEG$GO.0045182) # 6 DEG 
GO0045182_01h_03h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0045182_01h_03h_D_DEG,by='row.names')

#GO:0003700-- DNA-binding transcription factor activity 
GO0003700_01h_03h_D_scores <- scoresInTerm(GOdata_MF_01h_03h_D,"GO:0003700", use.names = T)
GO0003700_01h_03h_D_scores_df<- as.data.frame(GO0003700_01h_03h_D_scores )
GO0003700_01h_03h_D_DEG <- subset.data.frame(GO0003700_01h_03h_D_scores_df,GO0003700_01h_03h_D_scores_df== 2 )
length(GO0003700_01h_03h_D_DEG$GO.0003700) # 9 DEG 
GO0003700_01h_03h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0003700_01h_03h_D_DEG,by='row.names')


#GO:0140102--catalytic activity, acting on a rRNA
GO0140102_01h_03h_D_scores <- scoresInTerm(GOdata_MF_01h_03h_D,"GO:0140102", use.names = T)
GO0140102_01h_03h_D_scores_df<- as.data.frame(GO0140102_01h_03h_D_scores )
GO0140102_01h_03h_D_DEG <- subset.data.frame(GO0140102_01h_03h_D_scores_df,GO0140102_01h_03h_D_scores_df== 2 )
length(GO0140102_01h_03h_D_DEG$GO.0140102) # 9 DEG 
GO0140102_01h_03h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0140102_01h_03h_D_DEG,by='row.names')

#GO:0015036 - disulfide oxidoreductase activity 
GO0015036_01h_03h_D_scores <- scoresInTerm(GOdata_MF_01h_03h_D,"GO:0015036", use.names = T)
GO0015036_01h_03h_D_scores_df<- as.data.frame(GO0015036_01h_03h_D_scores )
GO0015036_01h_03h_D_DEG <- subset.data.frame(GO0015036_01h_03h_D_scores_df,GO0015036_01h_03h_D_scores_df== 2 )
length(GO0015036_01h_03h_D_DEG$GO.0015036) # 2 DEG 
GO0015036_01h_03h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0015036_01h_03h_D_DEG,by='row.names')


#3h-6h-UP-- Molecular Function
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_03h_06h_UP
#gathering BLAST for export
GO0045182_03h_06h_UP_DEG_blastp 
GO0005509_03h_06h_UP_DEG_blastp
GO0004930_03h_06h_UP_DEG_blastp
GO0005215_03h_06h_UP_DEG_blastp
GO0005201_03h_06h_UP_DEG_blastp



#GO:0045182  -- translational regulator activity 
GO0045182_03h_06h_UP_scores <- scoresInTerm(GOdata_MF_03h_06h_UP,"GO:0045182", use.names = T)
GO0045182_03h_06h_UP_scores_df<- as.data.frame(GO0045182_03h_06h_UP_scores)
GO0045182_03h_06h_UP_DEG <- subset.data.frame(GO0045182_03h_06h_UP_scores_df,GO0045182_03h_06h_UP_scores_df== 2 )
length(GO0045182_03h_06h_UP_DEG$GO.0045182) # 18 DEG 
GO0045182_03h_06h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0045182_03h_06h_UP_DEG,by='row.names')

#GO:0004930 -- g protein coupled
GO0004930_03h_06h_UP_scores <- scoresInTerm(GOdata_MF_03h_06h_UP,"GO:0004930", use.names = T)
GO0004930_03h_06h_UP_scores_df<- as.data.frame(GO0004930_03h_06h_UP_scores)
GO0004930_03h_06h_UP_DEG <- subset.data.frame(GO0004930_03h_06h_UP_scores_df,GO0004930_03h_06h_UP_scores_df== 2 )
length(GO0004930_03h_06h_UP_DEG$GO.0004930) # 156 DEG 
GO0004930_03h_06h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0004930_03h_06h_UP_DEG,by='row.names')

# GO:0005509  -- calcium ion binding 
GO0005509_03h_06h_UP_scores <- scoresInTerm(GOdata_MF_03h_06h_UP,"GO:0005509", use.names = T)
GO0005509_03h_06h_UP_scores_df<- as.data.frame(GO0005509_03h_06h_UP_scores)
GO0005509_03h_06h_UP_DEG <- subset.data.frame(GO0005509_03h_06h_UP_scores_df,GO0005509_03h_06h_UP_scores_df== 2 )
length(GO0005509_03h_06h_UP_DEG$GO.0005509) # 81 DEG 
GO0005509_03h_06h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005509_03h_06h_UP_DEG,by='row.names')

#GO:0005215 -- transporter activity 
GO0005215_03h_06h_UP_scores <- scoresInTerm(GOdata_MF_03h_06h_UP,"GO:0005215", use.names = T)
GO0005215_03h_06h_UP_scores_df<- as.data.frame(GO0005215_03h_06h_UP_scores)
GO0005215_03h_06h_UP_DEG <- subset.data.frame(GO0005215_03h_06h_UP_scores_df,GO0005215_03h_06h_UP_scores_df== 2 )
length(GO0005215_03h_06h_UP_DEG$GO.0005215) # 106 DEG 
GO0005215_03h_06h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005215_03h_06h_UP_DEG,by='row.names')

#GO:0005201 extracellular matrix structural constituent
GO0005201_03h_06h_UP_scores <- scoresInTerm(GOdata_MF_03h_06h_UP,"GO:0005201", use.names = T)
GO0005201_03h_06h_UP_scores_df<- as.data.frame(GO0005201_03h_06h_UP_scores)
GO0005201_03h_06h_UP_DEG <- subset.data.frame(GO0005201_03h_06h_UP_scores_df,GO0005201_03h_06h_UP_scores_df== 2 )
length(GO0005201_03h_06h_UP_DEG$GO.0005201) # 10 DEG 
GO0005201_03h_06h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005201_03h_06h_UP_DEG,by='row.names')

#3h-6h-D-- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_03h_06h_D
#gathering BLAST for export
GO0003824_03h_06h_D_DEG_blastp 
GO0005198_03h_06h_D_DEG_blastp
GO0016491_03h_06h_D_DEG_blastp

#GO:0003824   -- catalytic activity 
GO0003824_03h_06h_D_scores <- scoresInTerm(GOdata_MF_03h_06h_D,"GO:0003824", use.names = T)
GO0003824_03h_06h_D_scores_df<- as.data.frame(GO0003824_03h_06h_D_scores)
GO0003824_03h_06h_D_DEG <- subset.data.frame(GO0003824_03h_06h_D_scores_df,GO0003824_03h_06h_D_scores_df== 2 )
length(GO0003824_03h_06h_D_DEG$GO.0003824) # 391 DEG 
GO0003824_03h_06h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0003824_03h_06h_D_DEG,by='row.names')

#GO:0005198 structural molecule activity
GO0005198_03h_06h_D_scores <- scoresInTerm(GOdata_MF_03h_06h_D,"GO:0005198", use.names = T)
GO0005198_03h_06h_D_scores_df<- as.data.frame(GO0005198_03h_06h_D_scores)
GO0005198_03h_06h_D_DEG <- subset.data.frame(GO0005198_03h_06h_D_scores_df,GO0005198_03h_06h_D_scores_df== 2 )
length(GO0005198_03h_06h_D_DEG$GO.0005198) # 53 DEG 
GO0005198_03h_06h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0005198_03h_06h_D_DEG,by='row.names')

#GO:0016491 -- oxidoreductase activity
GO0016491_03h_06h_D_scores <- scoresInTerm(GOdata_MF_03h_06h_D,"GO:0016491", use.names = T)
GO0016491_03h_06h_D_scores_df<- as.data.frame(GO0016491_03h_06h_D_scores)
GO0016491_03h_06h_D_DEG <- subset.data.frame(GO0016491_03h_06h_D_scores_df,GO0016491_03h_06h_D_scores_df$GO.0016491== 2 )
length(GO0016491_03h_06h_D_scores$`GO:0016491`) # 358 DEG 
GO0016491_03h_06h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0016491_03h_06h_D_DEG,by='row.names')


#6h-12h--UP -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_06h_12h_UP
#gathering BLAST for export
GO0005198_06h_12h_UP_DEG_blastp 
GO0016880_06h_12h_UP_DEG_blastp
GO0016702_06h_12h_UP_DEG_blastp

#GO:0005198 -- structural molecule activity  
GO0005198_06h_12h_UP_scores <- scoresInTerm(GOdata_06h_12h_MF_UP,"GO:0005198", use.names = T)
GO0005198_06h_12h_UP_scores_df<- as.data.frame(GO0005198_06h_12h_UP_scores)
GO0005198_06h_12h_UP_DEG <- subset.data.frame(GO0005198_06h_12h_UP_scores_df,GO0005198_06h_12h_UP_scores_df== 2 )
length(GO0005198_06h_12h_UP_DEG$GO.0005198) # 5 DEG 
GO0005198_06h_12h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005198_06h_12h_UP_DEG ,by='row.names')

#GO:0016880 -- acid-ammonia (or amide) ligase activity
GO0016880_06h_12h_UP_scores <- scoresInTerm(GOdata_06h_12h_MF_UP,"GO:0016880", use.names = T)
GO0016880_06h_12h_UP_scores_df<- as.data.frame(GO0016880_06h_12h_UP_scores)
GO0016880_06h_12h_UP_DEG <- subset.data.frame(GO0016880_06h_12h_UP_scores_df,GO0016880_06h_12h_UP_scores_df== 2 )
length(GO0016880_06h_12h_UP_DEG$GO.0016880) # 1 DEG 
GO0016880_06h_12h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0016880_06h_12h_UP_DEG ,by='row.names')

#GO:0016702 oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen 
GO0016702_06h_12h_UP_scores <- scoresInTerm(GOdata_06h_12h_MF_UP,"GO:0016702", use.names = T)
GO0016702_06h_12h_UP_scores_df<- as.data.frame(GO0016702_06h_12h_UP_scores)
GO0016702_06h_12h_UP_DEG <- subset.data.frame(GO0016702_06h_12h_UP_scores_df,GO0016702_06h_12h_UP_scores_df== 2 )
length(GO0016702_06h_12h_UP_DEG$GO.0016702) # 1 DEG 
GO0016702_06h_12h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0016702_06h_12h_UP_DEG ,by='row.names')



#6-12--Down -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_06h_12h_D
#gathering BLAST for export
MF_06h_12h_D <- list(GO0005201_06h_12h_D_DEG_blastp,
GO0016721_06h_12h_D_DEG_blastp,
GO0004869_06h_12h_D_DEG_blastp, 
GO0016859_06h_12h_D_DEG_blastp,
GO0004298_06h_12h_D_DEG_blastp)

#GO:0005201--extracellular matrix structural constituent  
GO0005201_06h_12h_D_scores <- scoresInTerm(GOdata_06h_12h_MF_D,"GO:0005201", use.names = T)
GO0005201_06h_12h_D_scores_df<- as.data.frame(GO0005201_06h_12h_D_scores)
GO0005201_06h_12h_D_DEG <- subset.data.frame(GO0005201_06h_12h_D_scores_df,GO0005201_06h_12h_D_scores_df== 2 )
length(GO0005201_06h_12h_D_DEG$GO.0005201) # 2 DEG 
GO0005201_06h_12h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0005201_06h_12h_D_DEG ,by='row.names')

#GO:0016721 oxidoreductase activity, acting on superoxide radicals as acceptor 
GO0016721_06h_12h_D_scores <- scoresInTerm(GOdata_06h_12h_MF_D,"GO:0016721", use.names = T)
GO0016721_06h_12h_D_scores_df<- as.data.frame(GO0016721_06h_12h_D_scores)
GO0016721_06h_12h_D_DEG <- subset.data.frame(GO0016721_06h_12h_D_scores_df,GO0016721_06h_12h_D_scores_df== 2 )
length(GO0016721_06h_12h_D_DEG$GO.0016721) # 1 DEG 
GO0016721_06h_12h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0016721_06h_12h_D_DEG ,by='row.names')


# GO:0004869-- cysteine-type endopeptidase inhibitor activity
GO0004869_06h_12h_D_scores <- scoresInTerm(GOdata_06h_12h_MF_D,"GO:0004869", use.names = T)
GO0004869_06h_12h_D_scores_df<- as.data.frame(GO0004869_06h_12h_D_scores)
GO0004869_06h_12h_D_DEG <- subset.data.frame(GO0004869_06h_12h_D_scores_df,GO0004869_06h_12h_D_scores_df== 2 )
length(GO0004869_06h_12h_D_DEG$GO.0004869) # 1 DEG 
GO0004869_06h_12h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0004869_06h_12h_D_DEG ,by='row.names')

#GO:0016859-- cis-trans isomerase activity    
GO0016859_06h_12h_D_scores <- scoresInTerm(GOdata_06h_12h_MF_D,"GO:0016859", use.names = T)
GO0016859_06h_12h_D_scores_df<- as.data.frame(GO0016859_06h_12h_D_scores)
GO0016859_06h_12h_D_DEG <- subset.data.frame(GO0016859_06h_12h_D_scores_df,GO0016859_06h_12h_D_scores_df== 2 )
length(GO0016859_06h_12h_D_DEG$GO.0016859) # 1 DEG 
GO0016859_06h_12h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0016859_06h_12h_D_DEG ,by='row.names')

#GO:0004298 --threonine-type endopeptidase activity
GO0004298_06h_12h_D_scores <- scoresInTerm(GOdata_06h_12h_MF_D,"GO:0004298", use.names = T)
GO0004298_06h_12h_D_scores_df<- as.data.frame(GO0004298_06h_12h_D_scores)
GO0004298_06h_12h_D_DEG <- subset.data.frame(GO0004298_06h_12h_D_scores_df,GO0004298_06h_12h_D_scores_df== 2 )
length(GO0004298_06h_12h_D_DEG$GO.0004298) # 1 DEG 
GO0004298_06h_12h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0004298_06h_12h_D_DEG ,by='row.names')

#12-24--UP -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_12h_24h_UP
#gathering BLAST for export
GO0005198_12h_24h_UP_DEG_blastp

#GO:0005198 -- structural molecule activity 
GO0005198_12h_24h_UP_scores <- scoresInTerm(GOdata_12h_24h_MF_UP,"GO:0005198", use.names = T)
GO0005198_12h_24h_UP_scores_df<- as.data.frame(GO0005198_12h_24h_UP_scores)
GO0005198_12h_24h_UP_DEG <- subset.data.frame(GO0005198_12h_24h_UP_scores_df,GO0005198_12h_24h_UP_scores_df== 2 )
length(GO0005198_12h_24h_UP_DEG$GO.0005198) # 2 DEG 
GO0005198_12h_24h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005198_12h_24h_UP_DEG ,by='row.names')

#12-24--DOWN -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_12h_24h_D
#gathering BLAST for export
GO0005200_12h_24h_D_DEG_blastp
GO0005509_12h_24h_D_DEG_blastp
GO0003924_12h_24h_D_DEG_blastp


#GO:0005200 -- constituent of cytoskeleton  
GO0005200_12h_24h_D_scores <- scoresInTerm(GOdata_12h_24h_MF_D,"GO:0005200", use.names = T)
GO0005200_12h_24h_D_scores_df<- as.data.frame(GO0005200_12h_24h_D_scores)
GO0005200_12h_24h_D_DEG <- subset.data.frame(GO0005200_12h_24h_D_scores_df,GO0005200_12h_24h_D_scores_df== 2 )
length(GO0005200_12h_24h_D_DEG$GO.0005200) # 5 DEG 
GO0005200_12h_24h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0005200_12h_24h_D_DEG ,by='row.names')

# GO:0005509 -- calcium ion binding 
GO0005509_12h_24h_D_scores <- scoresInTerm(GOdata_12h_24h_MF_D,"GO:0005509", use.names = T)
GO0005509_12h_24h_D_scores_df<- as.data.frame(GO0005509_12h_24h_D_scores)
GO0005509_12h_24h_D_DEG <- subset.data.frame(GO0005509_12h_24h_D_scores_df,GO0005509_12h_24h_D_scores_df== 2 )
length(GO0005509_12h_24h_D_DEG$GO.0005509) # 5 DEG 
GO0005509_12h_24h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0005509_12h_24h_D_DEG ,by='row.names')

#GO:0003924 -- GTPase activity 
GO0003924_12h_24h_D_scores <- scoresInTerm(GOdata_12h_24h_MF_D,"GO:0003924", use.names = T)
GO0003924_12h_24h_D_scores_df<- as.data.frame(GO0003924_12h_24h_D_scores)
GO0003924_12h_24h_D_DEG <- subset.data.frame(GO0003924_12h_24h_D_scores_df,GO0003924_12h_24h_D_scores_df== 2 )
length(GO0003924_12h_24h_D_DEG$GO.0003924) # 4 DEG 
GO0003924_12h_24h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0003924_12h_24h_D_DEG ,by='row.names')



#24-48--UP -- Molecular Function
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_24h_48h_UP
#gathering BLAST for export
GO0005200_24h_48h_UP_DEG_blastp 
GO0003924_24h_48h_UP_DEG_blastp
GO0043167_24h_48h_UP_DEG_blastp
GO0003723_24h_48h_UP_DEG_blastp

#GO:0005200 -- constituent of cytoskeleton 
GO0005200_24h_48h_UP_scores <- scoresInTerm(GOdata_24h_48h_MF_UP,"GO:0005200", use.names = T)
GO0005200_24h_48h_UP_scores_df<- as.data.frame(GO0005200_24h_48h_UP_scores)
GO0005200_24h_48h_UP_DEG <- subset.data.frame(GO0005200_24h_48h_UP_scores_df,GO0005200_24h_48h_UP_scores_df== 2 )
length(GO0005200_24h_48h_UP_DEG$GO.0005200) # 6 DEG 
GO0005200_24h_48h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0005200_24h_48h_UP_DEG ,by='row.names')

#GO:0003924 GTPase activity 
GO0003924_24h_48h_UP_scores <- scoresInTerm(GOdata_24h_48h_MF_UP,"GO:0003924", use.names = T)
GO0003924_24h_48h_UP_scores_df<- as.data.frame(GO0003924_24h_48h_UP_scores)
GO0003924_24h_48h_UP_DEG <- subset.data.frame(GO0003924_24h_48h_UP_scores_df,GO0003924_24h_48h_UP_scores_df== 2 )
length(GO0003924_24h_48h_UP_DEG$GO.0003924) # 6 DEG 
GO0003924_24h_48h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0003924_24h_48h_UP_DEG ,by='row.names')


#GO:0043167--  ion binding
GO0043167_24h_48h_UP_scores <- scoresInTerm(GOdata_24h_48h_MF_UP,"GO:0043167", use.names = T)
GO0043167_24h_48h_UP_scores_df<- as.data.frame(GO0043167_24h_48h_UP_scores)
GO0043167_24h_48h_UP_DEG <- subset.data.frame(GO0043167_24h_48h_UP_scores_df,GO0043167_24h_48h_UP_scores_df== 2 )
length(GO0043167_24h_48h_UP_DEG$GO.0043167) # 16 DEG 
GO0043167_24h_48h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0043167_24h_48h_UP_DEG ,by='row.names')

#GO:0003723 -- RNA binding
GO0003723_24h_48h_UP_scores <- scoresInTerm(GOdata_24h_48h_MF_UP,"GO:0003723", use.names = T)
GO0003723_24h_48h_UP_scores_df<- as.data.frame(GO0003723_24h_48h_UP_scores)
GO0003723_24h_48h_UP_DEG <- subset.data.frame(GO0003723_24h_48h_UP_scores_df,GO0003723_24h_48h_UP_scores_df== 2 )
length(GO0003723_24h_48h_UP_DEG$GO.0003723) # 6 DEG 
GO0003723_24h_48h_UP_DEG_blastp <- merge(HumRef2020ML_df ,GO0003723_24h_48h_UP_DEG ,by='row.names')


#24-48--DOWN -- Molecular Function 
#gathering redundancy reduced matrix for interval(Uncut-10m,etc) and direction(Up or Down)
reducedMatrix_MF_24h_48h_D
#gathering BLAST for export
GO0005509_24h_48h_D_DEG_blastp
GO0003735_24h_48h_D_DEG_blastp
GO0016859_24h_48h_D_DEG_blastp

#GO:0005509  -- Calcium ion binding  
GO0005509_24h_48h_D_scores <- scoresInTerm(GOdata_24h_48h_MF_D,"GO:0005509", use.names = T)
GO0005509_24h_48h_D_scores_df<- as.data.frame(GO0005509_24h_48h_D_scores)
GO0005509_24h_48h_D_DEG <- subset.data.frame(GO0005509_24h_48h_D_scores_df,GO0005509_24h_48h_D_scores_df== 2 )
length(GO0005509_24h_48h_D_DEG$GO.0005509) # 4 DEG 
GO0005509_24h_48h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0005509_24h_48h_D_DEG ,by='row.names')

# GO:0003735 structural constituent of ribosome 
GO0003735_24h_48h_D_scores <- scoresInTerm(GOdata_24h_48h_MF_D,"GO:0003735", use.names = T)
GO0003735_24h_48h_D_scores_df<- as.data.frame(GO0003735_24h_48h_D_scores)
GO0003735_24h_48h_D_DEG <- subset.data.frame(GO0003735_24h_48h_D_scores_df,GO0003735_24h_48h_D_scores_df== 2 )
length(GO0003735_24h_48h_D_DEG$GO.0003735) # 2 DEG 
GO0003735_24h_48h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0003735_24h_48h_D_DEG ,by='row.names')

# GO:0016859  cis-trans isomerase activity 
GO0016859_24h_48h_D_scores <- scoresInTerm(GOdata_24h_48h_MF_D,"GO:0016859", use.names = T)
GO0016859_24h_48h_D_scores_df<- as.data.frame(GO0016859_24h_48h_D_scores)
GO0016859_24h_48h_D_DEG <- subset.data.frame(GO0016859_24h_48h_D_scores_df,GO0016859_24h_48h_D_scores_df== 2 )
length(GO0016859_24h_48h_D_DEG$GO.0016859) # 1 DEG 
GO0016859_24h_48h_D_DEG_blastp <- merge(HumRef2020ML_df ,GO0016859_24h_48h_D_DEG ,by='row.names')

#Merge all of the Top GO categories in each interval with the genes


NOIS_GO <- list(MF_Uncut_10m_up = list(GO0008061_U_10m_UP_DEG_blastp,
                GO0003700_U_10m_UP_DEG_blastp,
                GO0004553_U_10m_UP_DEG_blastp,
                GO0038024_U_10m_UP_DEG_blastp,
                GO0098599_U_10m_UP_DEG_blastp),
                       MF_Uncut_10m_down = list(GO0005509_U_10m_D_DEG_blastp,
                                                GO0016880_U_10m_D_DEG_blastp,
                                                GO0005198_U_10m_D_DEG_blastp,
                                                GO0008519_U_10m_D_DEG_blastp),
                       MF_10m_01h_down= GO0003700_10m_01h_D_DEG_blastp ,
                       MF_10m_01h_up=list(GO0003700_10m_01h_UP_DEG_blastp,
                       GO0001067_10m_01h_UP_DEG_blastp,
                       GO0004222_10m_01h_UP_DEG_blastp) ,
                       MF_01h_03h_down =list(GO0003723_01h_03h_D_DEG_blastp, 
                       GO0045182_01h_03h_D_DEG_blastp,
                       GO0003700_01h_03h_D_DEG_blastp,
                       GO0140102_01h_03h_D_DEG_blastp,
                       GO0015036_01h_03h_D_DEG_blastp),
                       MF_01h_03h_up =list(GO0005198_01h_03h_UP_DEG_blastp,
                       GO0003924_01h_03h_UP_DEG_blastp),
                       MF_03h_06h_down = list(GO0003824_03h_06h_D_DEG_blastp,
                       GO0005198_03h_06h_D_DEG_blastp,
                       GO0016491_03h_06h_D_DEG_blastp),
                       MF_03h_06h_up = list(GO0045182_03h_06h_UP_DEG_blastp,
                       GO0005509_03h_06h_UP_DEG_blastp,
                       GO0004930_03h_06h_UP_DEG_blastp,
                       GO0005215_03h_06h_UP_DEG_blastp,
                       GO0005201_03h_06h_UP_DEG_blastp),
                       MF_06h_12h_down = list(GO0005201_06h_12h_D_DEG_blastp,
                       GO0016721_06h_12h_D_DEG_blastp,
                       GO0004869_06h_12h_D_DEG_blastp,
                       GO0016859_06h_12h_D_DEG_blastp,
                       GO0004298_06h_12h_D_DEG_blastp),
                       MF_06h_12h_up = list(GO0005198_06h_12h_UP_DEG_blastp,
                       GO0016880_06h_12h_UP_DEG_blastp,
                       GO0016702_06h_12h_UP_DEG_blastp),
                       MF_12h_24h_down = list(GO0005200_12h_24h_D_DEG_blastp,
                       GO0005509_12h_24h_D_DEG_blastp,
                       GO0003924_12h_24h_D_DEG_blastp),
                       MF_12h_24h_up =GO0005198_12h_24h_UP_DEG_blastp,
                       MF_24h_48h_down = list(GO0005509_24h_48h_D_DEG_blastp,
                       GO0003735_24h_48h_D_DEG_blastp,
                       GO0016859_24h_48h_D_DEG_blastp),
                       MF_24h_48h_up =list(GO0005200_24h_48h_UP_DEG_blastp,
                       GO0003924_24h_48h_UP_DEG_blastp,
                       GO0043167_24h_48h_UP_DEG_blastp,
                       GO0003723_24h_48h_UP_DEG_blastp ))

options(max.print=999999)
capture.output(NOIS_GO, file = "NOIS_GO_redunred_v4.csv")


