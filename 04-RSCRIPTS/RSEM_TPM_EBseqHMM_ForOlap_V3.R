##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
#Copyright (C) 2021 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin17.0 (64-bit)
#

### FILE NAME ###
#RNAseq_Bisect_RSEM_TPM_EBseqHMM_ForOlap.R# 
## This script uses TPM as an inout for EBseqHMM , then uses MEDIAN NORM to generate "sizes" for EBSeqHMM


## Loading EBSeqHMM

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBSeqHMM")
library(EBSeqHMM)
# V 1.28.0


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


## Preparing the data for EBSeqHMM -- need a matrix and to establish sizes 


RSEMcounts_B2_Trm_TPM_mt <- as.matrix(RSEMcounts_B2_Trm_TPM) #puts the data into matrix format

cond <- c(rep("U",3),rep("0.1",3),rep("01",3),rep("03",3),rep("06",3),rep("12",3),rep("24",3),rep("48",3))
Conditions <- factor(cond, levels=c("U","0.1","01","03","06","12","24","48"))

Sizes <- MedianNorm(RSEMcounts_B2_Trm_TPM_mt)

#The normalized matrix may be obtained by GetNormalizedMat function -- but this is just for visualization purposes

GeneNormData <- GetNormalizedMat(RSEMcounts_B2_Trm_TPM_mt, Sizes)


#Function EBSeqHMMTest may be used to estimate parameters and the posterior probability (PP) of being in each expression path.
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=RSEMcounts_B2_Trm_TPM_mt, sizeFactors=Sizes, Conditions=Conditions)

View(EBSeqHMMGeneOut)
GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
View(GeneDECalls)
length(rownames(GeneDECalls))  #1910 -- TPM
write.csv(rownames(GeneDECalls),row.names = F, file = "EBseqHMM_sig_genes_TPM_V1.28.0_V3_.csv") 

# Convert EBseq-hmm output to dataframe
GeneDECalls_df <- as.data.frame(GeneDECalls)
GeneDECalls_df2 <- data.frame(Gene= row.names(GeneDECalls_df), Most_Likely_Path = GeneDECalls_df$Most_Likely_Path,Max_PP = GeneDECalls_df$Max_PP )

#Split into a list 
GeneDECalls_df2_split<- split(GeneDECalls_df2,GeneDECalls_df2[,"Most_Likely_Path"])
length(GeneDECalls_df2_split)
View(GeneDECalls_df2_split)

# Using a cutoff of PP 0.5, looking at the top expression paths associated with the dataset 
GeneDECalls_PP0.5 <- GeneDECalls_df2[GeneDECalls_df2[,"Max_PP"] >= 0.5,]
length(rownames(GeneDECalls_PP0.5)) #122
View(GeneDECalls_PP0.5)

#Split into a list 
GeneDECalls_PP0.5_split<- split(GeneDECalls_PP0.5,GeneDECalls_PP0.5[,"Most_Likely_Path"])
View(GeneDECalls_PP0.5_split)
capture.output(GeneDECalls_PP0.5_split, file = "GeneDECalls_PP0.5_split.txt")


#BLAST  
HumRef2020ML <- read.table("../03-DATA_FILES/HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

#ML gene models and their protein blast hit 
GeneDECalls_blastp_merge <- merge(HumRef2020ML_df ,GeneDECalls_df,by='row.names')
length(rownames(GeneDECalls_blastp_merge)) #1001
GeneDECalls_blastp_merge_split<- split(GeneDECalls_blastp_merge,GeneDECalls_blastp_merge[,"Most_Likely_Path"])
View(GeneDECalls_blastp_merge_split)
capture.output(GeneDECalls_blastp_merge_split, file = "GeneDECalls_blastp_merge_split.txt")



# Using a cutoff of PP 0.5, looking at the top expression paths associated with the dataset ONLY with Blast hits 
GeneDECalls_blastp_merge_PP0.5 <- GeneDECalls_blastp_merge[GeneDECalls_blastp_merge[,"Max_PP"] >= 0.5,]
length(rownames(GeneDECalls_blastp_merge_PP0.5)) #63




#order according to Max_PP -- Which genes/paths have the highest PP?
GeneDECalls_blastp_merge_PP0.5_PP_order<- GeneDECalls_blastp_merge_PP0.5[order(GeneDECalls_blastp_merge_PP0.5$Max_PP, decreasing = T),]
write.csv(GeneDECalls_blastp_merge_PP0.5_PP_order, file = "GeneDECalls_blastp_merge_PP0.5_PP_order.csv")

# order according to Most_Likely_Path -- easier to visualize 
GeneDECalls_blastp_merge_PP0.5_MLP_order<- GeneDECalls_blastp_merge_PP0.5[order(GeneDECalls_blastp_merge_PP0.5$Most_Likely_Path),]
write.csv(GeneDECalls_blastp_merge_PP0.5_MLP_order, file = "GeneDECalls_blastp_merge_PP0.5_MLP_order.csv")



#Split into a list 
GeneDECalls_blastp_merge_PP0.5_split<- split(GeneDECalls_blastp_merge_PP0.5,GeneDECalls_blastp_merge_PP0.5[,"Most_Likely_Path"])
View(GeneDECalls_blastp_merge_PP0.5_split)
capture.output(GeneDECalls_blastp_merge_PP0.5_split, file = "GeneDECalls_blastp_merge_PP0.5_split.txt")




