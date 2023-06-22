##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

##R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##Copyright (C) 2021 The R Foundation for Statistical Computing
##Platform: x86_64-apple-darwin17.0 (64-bit)


#Heatmap.2 visualization of Genes -- hierarchically cluster genes in the 3 program overlap list (396), visualize genes expression clusters with a heatmap,
#visualize each cluster using a multi-line graph, and run gene ontology analysis on the 3 program overlap list and each respective cluster from the hclust)    


#Taking the DEG results from EBseq-hmm, NOISeq and EdgeR to compare DEG identification across all three programs 
## VERSION 3 --  EBseq-hmm(TPM),  EdgeR(EC), NOISEQBIO(TPM,CPM filter)



## TPMS = Transcripts per million 
## EC = Expected count 
#####################



# Need to use dplyr to calculate the medians for the replicates in each time point -- three replicates per time point 
library(dplyr)
library(tibble) # for rownames_to_column function 
library(tidyverse)

#heatmap.2 {gplots}
library(gplots)


## For making multi-line graphs 
# I need to use dplyr to calculate the medians for the replicates in each time point -- three replicates per time point 
library(ggplot2)


## Loading the gene count data 
# Uncut Replicate 1 #
unct_i_Rp1_B2_Trm <-read.delim("../03-DATA_FILES/uncut-i-rep1_S21_L002_R1_001.fastq.gztrimmed_1PRSEMout.genes.results")
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


## binding the TPM count data together for all samples 

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
                                    hab48_3_i=hab48_i_Rp3_B2_Trm$TPM
)



Int_NOIS.TPM_ED.EC_EB.TPM <- read.csv(file="Int_NOIS.TPM_ED.EC_EB.TPM_V3.csv",1) # ML gene names of the genes identified as DE in all three methods
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS <- subset.data.frame(RSEMcounts_B2_Trm_TPM,row.names(RSEMcounts_B2_Trm_TPM) %in% Int_NOIS.TPM_ED.EC_EB.TPM$VNOIS_ED_EB_Set.IntersectionSets...111...) 
# Genes in all three methods PLUS their respective TPMS at each time point and replicate 
length(row.names(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS)) #304 


# Using the pipe operator %>%
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS %>% 
  rownames_to_column('gene') %>% #switches row names to a column because the tibble format will remove row names from a standard DF 
  rowwise() %>% # calculate the median at each row within the specified columns 
  mutate(Uncut_Med = median( uncut_1_i,uncut_2_i,uncut_3_i),
         Hab00_Med = median(hab00_1,hab00_2,hab00_3),
         Hab01_Med = median(hab01_1,hab01_2,hab01_3),
         Hab03_Med = median(hab03_1,hab03_2,hab03_3),
         Hab06_Med = median(hab06_1,hab06_2,hab06_3),
         Hab12_Med = median(hab12_1,hab12_2,hab12_3),
         Hab24_Med = median(hab24_1_i,hab24_2_i,hab24_3_i),
         Hab48_Med= median(hab48_1_i,hab48_2_i,hab48_3_i))


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds, select = c(gene,Uncut_Med,Hab00_Med,Hab01_Med,Hab03_Med,Hab06_Med,Hab12_Med,Hab24_Med,Hab48_Med )) 

Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only %>% 
  column_to_rownames('gene') # switch it back to a data frame and put back in the row names



## Heatmap.2 using median values for each time point for each gene in the 396 long gene list 


# THe consensus DEG - NOISEQ(TPM), EDGE R(EC), EBSEQHMM(TPM)
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows_mt <- as.matrix(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows)

#load colorblind friendly color palette 
library(viridis)
myCol <- viridis(100,alpha = 1)
# if you get the error "Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures)" 
# try running dev.off()
dev.off()
Ed_EB_NOIS_Pr_Compl_htmp <- heatmap.2(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows_mt,
                                      col=myCol,
                                      Colv = F,
                                      main="Overlap Genes, pearson, complete",
                                      key=T, keysize=1.0,
                                      scale= "row",
                                      density.info="none",
                                      reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
                                      trace="none",
                                      cexRow=0.2,
                                      cexCol=0.8,
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="complete"))




#Now need to know which genes are in which clusters generated 
# as in Cary et al 2019 -- using cutree to cut the tree at several heights until there is more than one cluster -- resulting in 19 clusters when h =1 

Ed_EB_NOIS_Pr_Compl_htmp_hclust <- as.hclust(Ed_EB_NOIS_Pr_Compl_htmp$rowDendrogram)
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut<- as.data.frame(cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h=1))
length(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut$`cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h = 1)`) #304
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_df <- data.frame( cluster = Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut$`cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h = 1)`, gene = rownames(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble <-as_tibble(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_df)
max(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_df$cluster) #16 clusters
### pull out data frames of genes in each cluster 
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 1))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 2))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 3))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 4))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 5))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 6))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 7))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 8))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 9))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 10))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 11))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 12))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 13))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 14))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 15))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 16))



### Multi Line graph for Cluster 1 
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 1))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char)
# I want to subset the tibble I just made to ONLY
#include the columsn for each gene that show the median TPM values for each time point 

#Need to reformat the data frame containing each gene in the cluster -- short to long format
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c1 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each=41) ## 41 is the total # of genes in the cluster, so each condition will be repeated 46 x 8 conditions
Cond_c1_factor <- factor(Cond_c1, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long,Cond_c1_factor = Cond_c1_factor )
#adding condition as a factor for plot labeling 

ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long_vect,
       aes(x= Cond_c1_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #1",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi-Line graph for Cluster 2 - supplement 

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 2))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2_char)  

#Need to reformat the data frame short to long format  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)

Cond_c2 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 17 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c2_factor <- factor(Cond_c2, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2_long,Cond_c2_factor = Cond_c2_factor )

ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster2_long_vect,
       aes(x= Cond_c2_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #2",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 3 -- in supplement 

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 3))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3_char) 

#Need to reformat the data frame short to long format 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c3 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 15 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c3_factor <- factor(Cond_c3, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3_long,Cond_c3_factor = Cond_c3_factor )


ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster3_long_vect,
       aes(x= Cond_c3_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #3",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 4 - in supplement

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 4))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4_char) 

#Need to reformat the data frame containing each gene in the cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)
#The first argument after the data (conditions) tells R where to store the variable names,
#and the second (values) tells R where to store the values of each former variable. 
#The -id simply tells R to gather everything but id.


Cond_c4 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 25 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c4_factor <- factor(Cond_c4, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4_long,Cond_c4_factor = Cond_c4_factor )


ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster4_long_vect,
       aes(x= Cond_c4_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #4",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 5 - in supplement

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 5))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5_char)   

#Need to reformat the data frame containing each gene in the cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c5 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 33 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c5_factor <- factor(Cond_c5, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5_long,Cond_c5_factor = Cond_c5_factor )


ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster5_long_vect,
       aes(x= Cond_c5_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #5",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 6 - in supplement

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 6))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c6 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 6 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c6_factor <- factor(Cond_c6, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6_long,Cond_c6_factor = Cond_c6_factor )


ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster6_long_vect,
       aes(x= Cond_c6_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #6",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 7 - in supplement

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 7))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c7 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 4 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c7_factor <- factor(Cond_c7, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7_long,Cond_c7_factor = Cond_c7_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster7_long_vect,
       aes(x= Cond_c7_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #7",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 8 - in main text 

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 8))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c8 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 6 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c8_factor <- factor(Cond_c8, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8_long,Cond_c8_factor = Cond_c8_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster8_long_vect,
       aes(x= Cond_c8_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #8",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 9 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 9))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c9 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 34 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c9_factor <- factor(Cond_c9, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9_long,Cond_c9_factor = Cond_c9_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster9_long_vect,
       aes(x= Cond_c9_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #9",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 10 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 10))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c10 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 21 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c10_factor <- factor(Cond_c10, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10_long,Cond_c10_factor = Cond_c10_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster10_long_vect,
       aes(x= Cond_c10_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #10",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 11 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 11))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c11 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 21 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c11_factor <- factor(Cond_c11, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11_long,Cond_c11_factor = Cond_c11_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster11_long_vect,
       aes(x= Cond_c11_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #11",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 12 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 12))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c12 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 33 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c12_factor <- factor(Cond_c12, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12_long,Cond_c12_factor = Cond_c12_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster12_long_vect,
       aes(x= Cond_c12_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #12",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 13 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 13))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c13 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 25 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c13_factor <- factor(Cond_c13, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13_long,Cond_c13_factor = Cond_c13_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster13_long_vect,
       aes(x= Cond_c13_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #13",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 14 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 14))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c14 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 12 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c14_factor <- factor(Cond_c14, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14_long,Cond_c14_factor = Cond_c14_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster14_long_vect,
       aes(x= Cond_c14_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #14",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 15 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 15))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c15 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 4 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c15_factor <- factor(Cond_c15, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15_long,Cond_c15_factor = Cond_c15_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster15_long_vect,
       aes(x= Cond_c15_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #15",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 16 - in supplement  

Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 16))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16$gene) 

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16_char)  

#Need to reformat the data frame short to long  
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c16 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each= 7 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c16_factor <- factor(Cond_c16, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16_long,Cond_c16_factor = Cond_c16_factor )


ggplot(data = Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster16_long_vect,
       aes(x= Cond_c16_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #16",x =" Time after bisection", y = "Transcripts per million (TPM)")



## Run Gene Ontology on lists of genes in selected clusters from heatmap (ED_EB_NOIS_Pr_Cmpt_htmp)
## Mnemiopsis GO terms and TopGO code adopted from Melissa DeBiasse 
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
library("topGO")

geneID2GO <- readMappings(file = "ML2.2.aa_goterms_2.txt") #Mnemiopsis genes and their GO annotations generated from interproscan 
length(geneID2GO)
geneUniverse <- names(geneID2GO) #sets gene names from annotation file as the gene universe


## GO on consensus DEG dataset(304)
Int_NOIS.TPM_ED.EC_EB.TPM_char <- as.character(Int_NOIS.TPM_ED.EC_EB.TPM$VNOIS_ED_EB_Set.IntersectionSets...111...)
length(Int_NOIS.TPM_ED.EC_EB.TPM_char) #304
geneList_Int_NOIS.TPM_ED.EC_EB.TPM <- factor(as.integer(geneUniverse %in% Int_NOIS.TPM_ED.EC_EB.TPM_char)) #informs TOPGO where to locate the genes of interest in the gene universe
names(geneList_Int_NOIS.TPM_ED.EC_EB.TPM) <-geneUniverse #geneList lists the gene names in the gene universe that are id as genes of interest(DGE)
length(geneUniverse)
length(Filter(function(x) x == 0, geneList_Int_NOIS.TPM_ED.EC_EB.TPM)) #Genes w GO annotation not in the DEG
length(Filter(function(x) x == 1, geneList_Int_NOIS.TPM_ED.EC_EB.TPM))#DEG with GO annotation 


## Biological Process 
GOdata_BP_Int_NOIS.TPM_ED.EC_EB.TPM <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList_Int_NOIS.TPM_ED.EC_EB.TPM, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_BP_Int_NOIS.TPM_ED.EC_EB.TPM<- runTest(GOdata_BP_Int_NOIS.TPM_ED.EC_EB.TPM, algorithm="classic", statistic="fisher")
allRes_BP_Int_NOIS.TPM_ED.EC_EB.TPM <- GenTable(GOdata_BP_Int_NOIS.TPM_ED.EC_EB.TPM,classicFisher = resultFisher_BP_Int_NOIS.TPM_ED.EC_EB.TPM, orderBy = "resultFisher", ranksOf = "classicFisher")

# Cellular Component 
GOdata_CC_Int_NOIS.TPM_ED.EC_EB.TPM<- new("topGOdata", description="My project", ontology="CC", allGenes=geneList_Int_NOIS.TPM_ED.EC_EB.TPM, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_CC_Int_NOIS.TPM_ED.EC_EB.TPM <- runTest(GOdata_CC_Int_NOIS.TPM_ED.EC_EB.TPM, algorithm="classic", statistic="fisher")
allRes_CC_Int_NOIS.TPM_ED.EC_EB.TPM <- GenTable(GOdata_CC_Int_NOIS.TPM_ED.EC_EB.TPM,classicFisher = resultFisher_CC_Int_NOIS.TPM_ED.EC_EB.TPM, orderBy = "resultFisher", ranksOf = "classicFisher")


## Molecular Function 
GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM<- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_Int_NOIS.TPM_ED.EC_EB.TPM, annot = annFUN.gene2GO, gene2GO = geneID2GO) 

resultFisher_MF_Int_NOIS.TPM_ED.EC_EB.TPM <- runTest(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM, algorithm="classic", statistic="fisher")
allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM <- GenTable(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,classicFisher = resultFisher_MF_Int_NOIS.TPM_ED.EC_EB.TPM, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 100)
View(allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM)
write.table(allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM, file='allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM.tsv', quote=FALSE, sep='\t', col.names = NA)
# "allRes_MF_Int_NOIS.TPM_ED.EC_EB.EC.tsv" will be used for the GO Figure! input

#Determine which consensus DEGs contain certain GO terms
#2 == included in the DEG, 1== not included in the DEG 

#BLASTP on DEG categories
#ML gene models and their protein blast hit 
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

# DNA-binding transcription factor activity -- "GO:0003700"
GO0003700_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0003700", use.names = T)
GO0003700_scores_df<- as.data.frame(GO0003700_scores)
GO0003700_DEG <- subset.data.frame(GO0003700_scores_df,GO0003700_scores_df == 2 )
length(GO0003700_DEG$GO.0003700) # 9 DEG
GO0003700_DEG_blastp<- merge(HumRef2020ML_df ,GO0003700_DEG,by='row.names')


# Transcriptional regulator activity GO:0140110 ## GO FIGURE! merged these two together in the bubble plot
GO0140110_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0140110", use.names = T)
GO0140110_scores_df<- as.data.frame(GO0140110_scores)
GO0140110_DEG <- subset.data.frame(GO0140110_scores_df,GO0140110_scores_df == 2 )
length(GO0140110_DEG$GO.0140110) # 9 DEG
GO0140110_DEG_blastp<- merge(HumRef2020ML_df ,GO0140110_DEG,by='row.names')

# Endopeptidase GO:0004175
GO0004175_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0004175", use.names = T)
GO0004175_scores_df<- as.data.frame(GO0004175_scores)
GO0004175_DEG <- subset.data.frame(GO0004175_scores_df,GO0004175_scores_df == 2 )
length(GO0004175_DEG$GO.0004175) # 10 DEG
GO0004175_DEG_blastp<- merge(HumRef2020ML_df ,GO0004175_DEG,by='row.names')

# peptidase activity GO:0008233 ## NOT listed in GO FIGURE! output, likely merged with another term 
GO0008233_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0008233", use.names = T)
GO0008233_scores_df<- as.data.frame(GO0008233_scores)
GO0008233_DEG <- subset.data.frame(GO0008233_scores_df,GO0008233_scores_df == 2 )
length(GO0008233_DEG$GO.0008233) # 12 DEG
GO0008233_DEG_blastp<- merge(HumRef2020ML_df ,GO0008233_DEG,by='row.names')

# structural constituent of cytoskeleton -- "GO:0005200"
GO0005200_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0005200", use.names = T)
GO0005200_scores_df<- as.data.frame(GO0005200_scores)
GO0005200_DEG <- subset.data.frame(GO0005200_scores_df,GO0005200_scores_df == 2 )
length(GO0005200_DEG$GO.0005200) # 4 DEG
GO0005200_DEG_blastp<- merge(HumRef2020ML_df ,GO0005200_DEG,by='row.names')


#  scavenger receptor activity -- " GO:0005044"
GO0005044_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0005044", use.names = T)
GO0005044_scores_df<- as.data.frame(GO0005044_scores)
GO0005044_DEG <- subset.data.frame(GO0005044_scores_df,GO0005044_scores_df == 2 )
length(GO0005044_DEG$GO.0005044) # 6 DEG
GO0005044_DEG_blastp<- merge(HumRef2020ML_df ,GO0005044_DEG,by='row.names')


# metallopeptidase activity -- "GO:0008237"
GO0008237_scores <- scoresInTerm(GOdata_MF_Int_NOIS.TPM_ED.EC_EB.TPM,"GO:0008237", use.names = T)
GO0008237_scores_df<- as.data.frame(GO0008237_scores)
GO0008237_DEG <- subset.data.frame(GO0008237_scores_df,GO0008237_scores_df == 2 )
length(GO0008237_DEG$GO.0008237) # 7 DEG
GO0008237_DEG_blastp<- merge(HumRef2020ML_df ,GO0008237_DEG,by='row.names')

GO_FIGURE_TOP5 <- list(DNA_binding_transcription_factor_activity_GO0003700 = GO0003700_DEG_blastp,
                       Endopeptidase_activity = GO0004175_DEG_blastp,
                       structural_constituent_of_cytoskeleton = GO0005200_DEG_blastp,
                       scavenger_receptor_activity = GO0005044_DEG_blastp,
                       metallopeptidase_activity_GO0008237=GO0008237_DEG_blastp
                       )

capture.output(GO_FIGURE_TOP5, file = "GO_FIGURE_TOP5_blastp.txt")

#running TOPGO for Cluster #1 & #12 generated from the heatmap(cutree)


###Cluster 1 ###
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1$gene) 
length(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char)
geneList_clust1 <- factor(as.integer(geneUniverse %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char)) #tells TopGO where these interesting genes appear in the 'geneUniverse' vector
names(geneList_clust1)<-geneUniverse #the geneList object tells TopGO which genes in the gene universe are your genes of interest


##Biological Process
GOdata_BP_clust1 <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList_clust1, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_BP_clust1 <- runTest(GOdata_BP_clust1, algorithm="classic", statistic="fisher")
allRes_BP_clust1 <- GenTable(GOdata_BP_clust1, classicFisher = resultFisher_BP_clust1, orderBy = "resultFisher", ranksOf = "classicFisher")

# Cellular component 
GOdata_CC_clust1 <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList_clust1, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_CC_clust1 <- runTest(GOdata_CC_clust1, algorithm="classic", statistic="fisher")
allRes_CC_clust1 <- GenTable(GOdata_CC_clust1, classicFisher = resultFisher_CC_clust1, orderBy = "resultFisher", ranksOf = "classicFisher")


## Molecular Function 
GOdata_MF_clust1 <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_clust1, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_clust1 <- runTest(GOdata_MF_clust1, algorithm="classic", statistic="fisher")
allRes_MF_clust1 <- GenTable(GOdata_MF_clust1, classicFisher = resultFisher_MF_clust1, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 100)
write.table(allRes_MF_clust1 , file='allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM_clust1.tsv', quote=FALSE, sep='\t', col.names = NA)
# "allres_MF_clust1.tsv" will be used for the GO Figure! input

###Cluster 12 ###
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12$gene) 
length(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char)
geneList_clust12 <- factor(as.integer(geneUniverse %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char)) #tells TopGO where these interesting genes appear in the 'geneUniverse' vector
names(geneList_clust12)<-geneUniverse #the geneList object tells TopGO which genes in the gene universe are your genes of interest


##Biological Process
GOdata_BP_clust12 <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList_clust12, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_BP_clust12 <- runTest(GOdata_BP_clust12, algorithm="classic", statistic="fisher")
allRes_BP_clust12 <- GenTable(GOdata_BP_clust12, classicFisher = resultFisher_BP_clust12, orderBy = "resultFisher", ranksOf = "classicFisher")

# Cellular component 
GOdata_CC_clust12 <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList_clust12, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher_CC_clust12 <- runTest(GOdata_CC_clust12, algorithm="classic", statistic="fisher")
allRes_CC_clust12 <- GenTable(GOdata_CC_clust12, classicFisher = resultFisher_CC_clust12, orderBy = "resultFisher", ranksOf = "classicFisher")


## Molecular Function 
GOdata_MF_clust12 <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList_clust12, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
resultFisher_MF_clust12 <- runTest(GOdata_MF_clust12, algorithm="classic", statistic="fisher")
allRes_MF_clust12 <- GenTable(GOdata_MF_clust12, classicFisher = resultFisher_MF_clust12, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 100)
write.table(allRes_MF_clust12 , file='allRes_MF_Int_NOIS.TPM_ED.EC_EB.TPM_clust12.tsv', quote=FALSE, sep='\t', col.names = NA)
# "allres_MF_clust12.tsv" will be used for the GO Figure! input



#BLASTP on Clusters 
#ML gene models and their protein blast hit 
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)



## Cluster 1 
Int_NOIS.TPM_ED.EC_EB.TPM_clust1_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 2 
Int_NOIS.TPM_ED.EC_EB.TPM_clust2_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust2_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 3 
Int_NOIS.TPM_ED.EC_EB.TPM_clust3_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust3_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 4 
Int_NOIS.TPM_ED.EC_EB.TPM_clust4_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust4_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 5 
Int_NOIS.TPM_ED.EC_EB.TPM_clust5_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust5_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 6 
Int_NOIS.TPM_ED.EC_EB.TPM_clust6_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust6_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 7 
Int_NOIS.TPM_ED.EC_EB.TPM_clust7_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust7_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 8 
Int_NOIS.TPM_ED.EC_EB.TPM_clust8_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     


## Cluster 9 
Int_NOIS.TPM_ED.EC_EB.TPM_clust9_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust9_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     


## Cluster 10 
Int_NOIS.TPM_ED.EC_EB.TPM_clust10_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust10_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 11
Int_NOIS.TPM_ED.EC_EB.TPM_clust11_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust11_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     


## Cluster 12
Int_NOIS.TPM_ED.EC_EB.TPM_clust12_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust12_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 13
Int_NOIS.TPM_ED.EC_EB.TPM_clust13_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust13_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 14
Int_NOIS.TPM_ED.EC_EB.TPM_clust14_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust14_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 15
Int_NOIS.TPM_ED.EC_EB.TPM_clust15_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust15_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

## Cluster 16
Int_NOIS.TPM_ED.EC_EB.TPM_clust16_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust16_char)) # updated with the  quantification value used for each method -- important to note that BOTH EBSeqHMM and EdgeR  use 'expected count'(EC)                                     

library(xlsx)

write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust1_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust1", row.names = F)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust2_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust2", row.names= F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust3_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust3", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust4_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust4", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust5_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust5", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust6_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust6", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust7_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust7", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust8_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust8", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust9_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust9", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust10_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust10", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust11_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust11", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust12_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust12", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust13_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust13", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust14_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust14", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust15_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust15", row.names = F, append=TRUE)
write.xlsx(Int_NOIS.TPM_ED.EC_EB.TPM_clust16_blastp, file="Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx", sheetName="Clust16", row.names = F, append=TRUE)


