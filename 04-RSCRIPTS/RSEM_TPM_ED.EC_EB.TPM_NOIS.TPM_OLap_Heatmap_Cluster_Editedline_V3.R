##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

##R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##Copyright (C) 2021 The R Foundation for Statistical Computing
##Platform: x86_64-apple-darwin17.0 (64-bit)


#Heatmap.2 visualization of Genes -- hierarchically cluster genes in the 3 program overlap list (304), visualize genes expression clusters with a heatmap,
#visualize each cluster using a multi-line graph, and run gene ontology analysis on the 3 program overlap list and each respective cluster from the hclust)    

## TPMS = Transcripts per million 

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



Int_NOIS.TPM_ED.EC_EB.TPM <- read.csv(file="Int_NOIS.TPM_ED.EC_EB.TPM_V3.csv") # ML gene names of the genes identified as DE in all three methods
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



## Heatmap.2 using median values for each time point for each gene in the 304 long gene list 


# lets do the entire list of genes(204) that are DE in ALL THREE methods - NOISEQ, EDGE R, EBSEQHMM
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows_mt <- as.matrix(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows)

#load colorblind friendly color palette 
library(viridis)
myCol <- viridis(100,alpha = 1)
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


# if you get the error "Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures)" 
# try running dev.off()
# if Ed_EB_NOIS_Pr_Compl_htmp creates a plot but the object cannot be found,try expanding the plot window to increase margins 

#Now need to know which genes are in which clusters generated 
# as in Cary et al 2019 -- using cutree to cut the tree at several heights until there is more than one cluster -- resulting in 19 clusters when h =1 

Ed_EB_NOIS_Pr_Compl_htmp_hclust <- as.hclust(Ed_EB_NOIS_Pr_Compl_htmp$rowDendrogram)
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut<- as.data.frame(cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h=1))
length(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut$`cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h = 1)`) #304
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_df <- data.frame( cluster = Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut$`cutree(Ed_EB_NOIS_Pr_Compl_htmp_hclust, h = 1)`, gene = rownames(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble <-as_tibble(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_df)

### pull out data frames of genes in each cluster 
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 1))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust8 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 8))


### Multi Line graph for Cluster 1 # EDIT by taking out two outliers "ML137118a" and "ML007419a"
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1 <- as.data.frame(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_tibble %>% filter(cluster == 1))
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char <- as.character(Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1$gene) 
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char_edit1 <- gsub('ML137118a','',Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char)
Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char_edit2 <- gsub('ML007419a','',Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char_edit1)

# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1 = subset(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows, rownames(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_meds_only_rows) %in% Ed_EB_NOIS_Pr_Compl_htmp_hclust_cut_clust1_char_edit2)
# I want to subset the tibble I just made to ONLY
#include the columsn for each gene that show the median TPM values for each time point 

#Need to reformat the data frame containing each gene in the cluster -- short to long format
Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long <- Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values, -gene)


Cond_c1 <- rep(c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"),each=39) ## x is the total # of genes in the cluster, so each condition will be repeated 46 x 8 conditions
Cond_c1_factor <- factor(Cond_c1, levels = c("Uncut","10Min","01Hr","03Hr","06Hr","12Hr","24Hr","48Hr"))          


Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long_vect <- cbind(Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long,Cond_c1_factor = Cond_c1_factor )
#adding condition as a factor for plot labeling 

ggplot(data =Int_NOIS.TPM_ED.EC_EB.TPM_TPMS_med_only_cluster1_long_vect,
       aes(x= Cond_c1_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = T) + theme(text = element_text(size = 20)) + labs(title="Cluster #1",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 12 - in main text 

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
       aes(x= Cond_c12_factor, y=values, group = gene, color = gene)) + geom_line(show.legend = F) + theme(text = element_text(size = 20)) + labs(title="Cluster #12",x =" Time after bisection", y = "Transcripts per million (TPM)")

