##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

#R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
#Copyright (C) 2023 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin20 (64-bit)

#Taking the DEG results from EBseq-hmm, NOISeq and EdgeR to compare DEG identification across all three programs 
##  EBseq-hmm(TPM),  EdgeR(EC), NOISEQBIO(TPM,CPM filter)

## TPMS = Transcripts per million 

#Loading required packages 

#‘Vennerable’ version 3.0

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("RBGL","graph","reshape"))
library(RBGL)
library(graph)
library(reshape)

BiocManager::valid()

install.packages("devtools")
library(devtools)
install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library(Vennerable)

session_info()

Venn_Data_NOIS <- as.matrix(read.csv("NOISeq_sig_v4.csv", header = T))
Venn_Data_edgeR <- as.matrix(read.csv("EdgeR_sig_v4.csv",header = T))
Venn_Data_EBseqHMM <- as.matrix(read.csv("EBseqHMM_sig_v4.csv", header = T))

VNOIS_ED_EB <-list( NOIS = Venn_Data_NOIS,edgeR = Venn_Data_edgeR,EBseqHMM = Venn_Data_EBseqHMM) # creating a list of all three datasets
View(VNOIS_ED_EB)
VNOIS_ED_EB_Set <- Venn(VNOIS_ED_EB)
plot(VNOIS_ED_EB_Set,doWeights = F, type="circles")


# Extracting lists of genes from each pairing of methods 
length(VNOIS_ED_EB_Set@IntersectionSets[["100"]]) #only NOIS #399
length(VNOIS_ED_EB_Set@IntersectionSets[["110"]]) # NOIS and edgeR #17
length(VNOIS_ED_EB_Set@IntersectionSets[["010"]]) #only edgeR #60
length(VNOIS_ED_EB_Set@IntersectionSets[["001"]]) #only EBSeqHMM #1206
length(VNOIS_ED_EB_Set@IntersectionSets[["101"]]) #NOIS and EBseqHMM #300
length(VNOIS_ED_EB_Set@IntersectionSets[["011"]]) #EBseqHMM and edgeR #153
Int_NOIS.TPM_ED.EC_EB.TPM <- as.data.frame(VNOIS_ED_EB_Set@IntersectionSets[["111"]])  
length(Int_NOIS.TPM_ED.EC_EB.TPM$`VNOIS_ED_EB_Set@IntersectionSets[["111"]]`)#118
write.csv(Int_NOIS.TPM_ED.EC_EB.TPM ,row.names = F, file = "consensus_DE_v4.csv")


## Run a Blast on the 3 overlap list -- including the e value associated with it 

HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)
#ML gene models and their protein blast hit 


Int_NOIS.TPM_ED.EC_EB.TPM_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Int_NOIS.TPM_ED.EC_EB.TPM$`VNOIS_ED_EB_Set@IntersectionSets[["111"]]`)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

write.csv(Int_NOIS.TPM_ED.EC_EB.TPM_blastp ,row.names =F , file = "consensus_DE_blastp_v4.csv")



