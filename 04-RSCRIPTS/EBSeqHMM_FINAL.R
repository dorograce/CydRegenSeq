##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 


#R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
#Copyright (C) 2023 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin20 (64-bit)

## Loading EBSeqHMM

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EBSeqHMM")
library(EBSeqHMM)
sessionInfo()
# V 1.34.0


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
length(rownames(GeneDECalls))  #1777 
write.csv(rownames(GeneDECalls),row.names = F, file = "EBseqHMM_sig_v4.csv") 

# Convert EBseq-hmm output to dataframe
GeneDECalls_df <- as.data.frame(GeneDECalls)
GeneDECalls_df2 <- data.frame(Gene= row.names(GeneDECalls_df), Most_Likely_Path = GeneDECalls_df$Most_Likely_Path,Max_PP = GeneDECalls_df$Max_PP )

#Split into a list 
GeneDECalls_df2_split<- split(GeneDECalls_df2,GeneDECalls_df2[,"Most_Likely_Path"])
length(GeneDECalls_df2_split) #161
View(GeneDECalls_df2_split)

#The expression path with the most genes
GeneDECalls_df2_split[["Up-Down-Down-Up-Down-Down-Up"]][["Gene"]]

# Using a cutoff of PP 0.5, looking at the top expression paths associated with the dataset 
GeneDECalls_PP0.5 <- GeneDECalls_df2[GeneDECalls_df2[,"Max_PP"] >= 0.5,]
length(rownames(GeneDECalls_PP0.5)) #102
View(GeneDECalls_PP0.5)

#Split into a list 
GeneDECalls_PP0.5_split<- split(GeneDECalls_PP0.5,GeneDECalls_PP0.5[,"Most_Likely_Path"])

options(max.print=999999)
capture.output(GeneDECalls_PP0.5_split, file = "EBSeqHMM_PP0.5.csv")

#BLAST  
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

#ML gene models and their protein blast hit 
GeneDECalls_blastp_merge <- merge(HumRef2020ML_df ,GeneDECalls_df,by='row.names')
length(rownames(GeneDECalls_blastp_merge)) #925
GeneDECalls_blastp_merge_split<- split(GeneDECalls_blastp_merge,GeneDECalls_blastp_merge[,"Most_Likely_Path"])

#export list 
options(max.print=999999)
capture.output(GeneDECalls_blastp_merge, file = "EBSeqHMM_v4_DE_Blast.csv")


# Using a cutoff of PP 0.5, looking at the top expression paths associated with the dataset ONLY with Blast hits 
GeneDECalls_blastp_merge_PP0.5 <- GeneDECalls_blastp_merge[GeneDECalls_blastp_merge[,"Max_PP"] >= 0.5,]
length(rownames(GeneDECalls_blastp_merge_PP0.5)) #48


#order according to Max_PP -- Which genes/paths have the highest PP?
GeneDECalls_blastp_merge_PP0.5_PP_order<- GeneDECalls_blastp_merge_PP0.5[order(GeneDECalls_blastp_merge_PP0.5$Max_PP, decreasing = T),]
write.csv(GeneDECalls_blastp_merge_PP0.5_PP_order, file = "EBSeqHMM_blast_PP0.5_PP_order.csv")

# order according to Most_Likely_Path -- easier to visualize 
GeneDECalls_blastp_merge_PP0.5_MLP_order<- GeneDECalls_blastp_merge_PP0.5[order(GeneDECalls_blastp_merge_PP0.5$Most_Likely_Path),]
length(GeneDECalls_blastp_merge_PP0.5_MLP_order$Row.names)
write.csv(GeneDECalls_blastp_merge_PP0.5_MLP_order, file = "EBSeqHMM_blast_PP0.5_MLP_order.csv")


#Split into a list 
GeneDECalls_blastp_merge_PP0.5_split<- split(GeneDECalls_blastp_merge_PP0.5,GeneDECalls_blastp_merge_PP0.5[,"Most_Likely_Path"])
length(GeneDECalls_blastp_merge_PP0.5_split)

#export list 
options(max.print=999999)
capture.output(GeneDECalls_blastp_merge_PP0.5_split, file = "EBSeqHMM_PP0.5_blast.csv")


#Determine which genes are also in the consensus set 
# For consensus gene set annotation 
con_genenames <- read.csv(file="consensus_DE_v4.csv",1) # ML gene names of the genes identified as DE in all three methods
con_genenames_df <- data.frame(row.names= con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111...,V1 = con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111... )
EBseqDECalls_PP.05_consensus <- merge(GeneDECalls_blastp_merge_PP0.5_MLP_order,con_genenames_df,by='V1')
write.csv(EBseqDECalls_PP.05_con, file= "EBseqDECalls_PP.05_consensus.csv")




