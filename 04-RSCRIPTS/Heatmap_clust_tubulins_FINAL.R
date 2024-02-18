##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 


#Heatmap.2 visualization of Genes -- hierarchically cluster genes in the 3 program overlap list(consensus), visualize genes expression clusters with a heatmap,
#visualize each cluster using a multi-line graph, and run gene ontology analysis on the consensus list)    


#Taking the DEG results from EBseq-hmm, NOISeq and EdgeR to compare DEG identification across all three programs 
## EBseq-hmm(TPM),  EdgeR(EC), NOISEQBIO(TPM,CPM filter)


# con = consensus
## TPMS = Transcripts per million 
## EC = Expected count 
#####################


# Need to use dplyr to calculate the medians for the replicates in each time point -- three replicates per time point 
library(dplyr)
library(tibble) # for rownames_to_column function 
library(tidyverse)

#heatmap.2 {gplots}
library(gplots)
#load colorblind friendly color palette 
install.packages("viridis")
library(viridis)
myCol <- viridis(100,alpha = 1)

## For making multi-line graphs 
# I need to use dplyr to calculate the medians for the replicates in each time point -- three replicates per time point 
install.packages("ggrepel")
library(ggrepel)
library(ggplot2)


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

## binding the TPM count data together for all samples 

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
                                    hab48_3_i=hab48_i_Rp3_B2_Trm$TPM
)



con_genenames <- read.csv(file="consensus_DE_v4.csv",1) # ML gene names of the genes identified as DE in all three methods
con_geneexp <- subset.data.frame(RSEMcounts_B2_Trm_TPM,row.names(RSEMcounts_B2_Trm_TPM) %in% con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111...) 
# Genes in all three methods PLUS their respective TPMS at each time point and replicate 
length(row.names(con_geneexp)) #118


# Using the pipe operator %>%
con_med_exp <- con_geneexp  %>% 
  rownames_to_column('gene') %>% #switches row names to a column because the tibble format will remove row names from a standard DF 
  rowwise() %>% # calculate the median at each row within the specified columns 
  mutate(Uncut_Med = median(uncut_1_i,uncut_2_i,uncut_3_i),
         Hab00_Med = median(hab00_1,hab00_2,hab00_3),
         Hab01_Med = median(hab01_1,hab01_2,hab01_3),
         Hab03_Med = median(hab03_1,hab03_2,hab03_3),
         Hab06_Med = median(hab06_1,hab06_2,hab06_3),
         Hab12_Med = median(hab12_1,hab12_2,hab12_3),
         Hab24_Med = median(hab24_1_i,hab24_2_i,hab24_3_i),
         Hab48_Med= median(hab48_1_i,hab48_2_i,hab48_3_i))

con_med_exp_only = subset(con_med_exp, select = c(gene,Uncut_Med,Hab00_Med,Hab01_Med,Hab03_Med,Hab06_Med,Hab12_Med,Hab24_Med,Hab48_Med )) 

con_med_exp_only_rows <- con_med_exp_only %>% 
  column_to_rownames('gene') # switch it back to a data frame and put back in the row names



## Heatmap.2 using median values for each time point for each gene in the consensus gene list 


# THe consensus DEG - NOISEQ(TPM), EDGE R(EC), EBSEQHMM(TPM)
consensus_med_only_rows_mt <- as.matrix(con_med_exp_only_rows)


# if you get the error "Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures)" 
# try running dev.off()
dev.off()
con_htmp <- heatmap.2(consensus_med_only_rows_mt,
                                      col=myCol,
                                      Colv = T,
                                      main="Consensus Genes, pearson, complete",
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
# as in Cary et al 2019 -- using cutree to cut the tree at several heights until there is more than one cluster -- resulting in 11 clusters when h=1 

con_htmp_hclust <- as.hclust(con_htmp$rowDendrogram)
con_htmp_hclust_cut<- as.data.frame(cutree(con_htmp_hclust, h=1))
length(con_htmp_hclust_cut$`cutree(con_htmp_hclust, h = 1)`)#118
con_htmp_hclust_cut_df <- data.frame(cluster = con_htmp_hclust_cut$`cutree(con_htmp_hclust, h = 1)`,gene = rownames(con_htmp_hclust_cut))
write.csv(con_htmp_hclust_cut_df, file = "con_htmp_hclust_cut_df.csv", row.names = F)
con_htmp_hclust_cut_tibble <-as_tibble(con_htmp_hclust_cut_df)
max(con_htmp_hclust_cut_tibble$cluster) #11 clusters
### pull out data frames of genes in each cluster 
con_clust1 <- as.data.frame( con_htmp_hclust_cut_tibble %>% filter(cluster == 1))
con_clust2 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 2))
con_clust3 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 3))
con_clust4 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 4))
con_clust5 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 5))
con_clust6 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 6))
con_clust7 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 7))
con_clust8 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 8))
con_clust9 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 9))
con_clust10 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 10))
con_clust11 <- as.data.frame(con_htmp_hclust_cut_tibble %>% filter(cluster == 11))

### Multi Line graph for Cluster 1 
con_clust1_char <- as.character(con_clust1$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster1 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust1_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster1_long <- con_med_only_cluster1 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster1_long_order <- con_med_only_cluster1_long[order(con_med_only_cluster1_long$gene),]

Cond_c1 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),12 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c1_factor <- factor(Cond_c1, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster1_long_order_x <- cbind(con_med_only_cluster1_long_order, x= rep(1:8,12), Cond_c1_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster1_long_order_x_label <- con_med_only_cluster1_long_order_x                           # Modify data
con_cluster1_long_order_x_label$label <- NA
con_cluster1_long_order_x_label$label[which(con_cluster1_long_order_x_label$x == max(con_cluster1_long_order_x_label$x))] <-con_cluster1_long_order_x_label$gene[which(con_cluster1_long_order_x_label$x == max(con_cluster1_long_order_x_label$x))]

ggplot(con_cluster1_long_order_x_label,
       aes(x= Cond_c1_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #1",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 2 
con_clust2_char <- as.character(con_clust2$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster2 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust2_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster2_long <- con_med_only_cluster2 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster2_long_order <- con_med_only_cluster2_long[order(con_med_only_cluster2_long$gene),]

Cond_c2 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),19 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c2_factor <- factor(Cond_c2, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster2_long_order_x <- cbind(con_med_only_cluster2_long_order, x= rep(1:8,19), Cond_c2_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster2_long_order_x_label <- con_med_only_cluster2_long_order_x                           # Modify data
con_cluster2_long_order_x_label$label <- NA
con_cluster2_long_order_x_label$label[which(con_cluster2_long_order_x_label$x == max(con_cluster2_long_order_x_label$x))] <-con_cluster2_long_order_x_label$gene[which(con_cluster2_long_order_x_label$x == max(con_cluster2_long_order_x_label$x))]

ggplot(con_cluster2_long_order_x_label,
       aes(x= Cond_c2_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #2",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 3 
con_clust3_char <- as.character(con_clust3$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster3 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust3_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster3_long <- con_med_only_cluster3 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster3_long_order <- con_med_only_cluster3_long[order(con_med_only_cluster3_long$gene),]

Cond_c3 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),9 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c3_factor <- factor(Cond_c3, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster3_long_order_x <- cbind(con_med_only_cluster3_long_order, x= rep(1:8,9), Cond_c3_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster3_long_order_x_label <- con_med_only_cluster3_long_order_x                           # Modify data
con_cluster3_long_order_x_label$label <- NA
con_cluster3_long_order_x_label$label[which(con_cluster3_long_order_x_label$x == max(con_cluster3_long_order_x_label$x))] <-con_cluster3_long_order_x_label$gene[which(con_cluster3_long_order_x_label$x == max(con_cluster3_long_order_x_label$x))]

ggplot(con_cluster3_long_order_x_label,
       aes(x= Cond_c3_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #3",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 4 
con_clust4_char <- as.character(con_clust4$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster4 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust4_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster4_long <- con_med_only_cluster4 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster4_long_order <- con_med_only_cluster4_long[order(con_med_only_cluster4_long$gene),]

Cond_c4 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),6 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c4_factor <- factor(Cond_c4, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster4_long_order_x <- cbind(con_med_only_cluster4_long_order, x= rep(1:8,6), Cond_c4_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster4_long_order_x_label <- con_med_only_cluster4_long_order_x                           # Modify data
con_cluster4_long_order_x_label$label <- NA
con_cluster4_long_order_x_label$label[which(con_cluster4_long_order_x_label$x == max(con_cluster4_long_order_x_label$x))] <-con_cluster4_long_order_x_label$gene[which(con_cluster4_long_order_x_label$x == max(con_cluster4_long_order_x_label$x))]

ggplot(con_cluster4_long_order_x_label,
       aes(x= Cond_c4_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #4",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 5 
con_clust5_char <- as.character(con_clust5$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster5 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust5_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster5_long <- con_med_only_cluster5 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster5_long_order <- con_med_only_cluster5_long[order(con_med_only_cluster5_long$gene),]

Cond_c5 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),20 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c5_factor <- factor(Cond_c5, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster5_long_order_x <- cbind(con_med_only_cluster5_long_order, x= rep(1:8,20), Cond_c5_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster5_long_order_x_label <- con_med_only_cluster5_long_order_x                           # Modify data
con_cluster5_long_order_x_label$label <- NA
con_cluster5_long_order_x_label$label[which(con_cluster5_long_order_x_label$x == max(con_cluster5_long_order_x_label$x))] <-con_cluster5_long_order_x_label$gene[which(con_cluster5_long_order_x_label$x == max(con_cluster5_long_order_x_label$x))]

ggplot(con_cluster5_long_order_x_label,
       aes(x= Cond_c5_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #5",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 6
con_clust6_char <- as.character(con_clust6$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster6 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust6_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster6_long <- con_med_only_cluster6 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster6_long_order <- con_med_only_cluster6_long[order(con_med_only_cluster6_long$gene),]

Cond_c6 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),17 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c6_factor <- factor(Cond_c6, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster6_long_order_x <- cbind(con_med_only_cluster6_long_order, x= rep(1:8,17), Cond_c6_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster6_long_order_x_label <- con_med_only_cluster6_long_order_x                         # Modify data
con_cluster6_long_order_x_label$label <- NA
con_cluster6_long_order_x_label$label[which(con_cluster6_long_order_x_label$x == max(con_cluster6_long_order_x_label$x))] <-con_cluster6_long_order_x_label$gene[which(con_cluster6_long_order_x_label$x == max(con_cluster6_long_order_x_label$x))]

ggplot(con_cluster6_long_order_x_label,
       aes(x= Cond_c6_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #6",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 7
con_clust7_char <- as.character(con_clust7$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster7 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust7_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster7_long <- con_med_only_cluster7 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster7_long_order <- con_med_only_cluster7_long[order(con_med_only_cluster7_long$gene),]

Cond_c7 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),4 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c7_factor <- factor(Cond_c7, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster7_long_order_x <- cbind(con_med_only_cluster7_long_order, x= rep(1:8,4), Cond_c7_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster7_long_order_x_label <- con_med_only_cluster7_long_order_x                         # Modify data
con_cluster7_long_order_x_label$label <- NA
con_cluster7_long_order_x_label$label[which(con_cluster7_long_order_x_label$x == max(con_cluster7_long_order_x_label$x))] <-con_cluster7_long_order_x_label$gene[which(con_cluster7_long_order_x_label$x == max(con_cluster7_long_order_x_label$x))]

ggplot(con_cluster7_long_order_x_label,
       aes(x= Cond_c7_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #7",x =" Time after bisection", y = "Transcripts per million (TPM)")



### Multi Line graph for Cluster 8
con_clust8_char <- as.character(con_clust8$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster8 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust8_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster8_long <- con_med_only_cluster8 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster8_long_order <- con_med_only_cluster8_long[order(con_med_only_cluster8_long$gene),]

Cond_c8 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),10 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c8_factor <- factor(Cond_c8, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster8_long_order_x <- cbind(con_med_only_cluster8_long_order, x= rep(1:8,10), Cond_c8_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster8_long_order_x_label <- con_med_only_cluster8_long_order_x                         # Modify data
con_cluster8_long_order_x_label$label <- NA
con_cluster8_long_order_x_label$label[which(con_cluster8_long_order_x_label$x == max(con_cluster8_long_order_x_label$x))] <-con_cluster8_long_order_x_label$gene[which(con_cluster8_long_order_x_label$x == max(con_cluster8_long_order_x_label$x))]

ggplot(con_cluster8_long_order_x_label,
       aes(x= Cond_c8_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #8",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 9
con_clust9_char <- as.character(con_clust9$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster9 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust9_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster9_long <- con_med_only_cluster9 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster9_long_order <- con_med_only_cluster9_long[order(con_med_only_cluster9_long$gene),]

Cond_c9 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),8 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c9_factor <- factor(Cond_c9, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster9_long_order_x <- cbind(con_med_only_cluster9_long_order, x= rep(1:8,8), Cond_c9_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster9_long_order_x_label <- con_med_only_cluster9_long_order_x                         # Modify data
con_cluster9_long_order_x_label$label <- NA
con_cluster9_long_order_x_label$label[which(con_cluster9_long_order_x_label$x == max(con_cluster9_long_order_x_label$x))] <-con_cluster9_long_order_x_label$gene[which(con_cluster9_long_order_x_label$x == max(con_cluster9_long_order_x_label$x))]

ggplot(con_cluster9_long_order_x_label,
       aes(x= Cond_c9_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #9",x =" Time after bisection", y = "Transcripts per million (TPM)")

### Multi Line graph for Cluster 10
con_clust10_char <- as.character(con_clust10$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster10 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust10_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster10_long <- con_med_only_cluster10 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster10_long_order <- con_med_only_cluster10_long[order(con_med_only_cluster10_long$gene),]

Cond_c10 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),9) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c10_factor <- factor(Cond_c10, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster10_long_order_x <- cbind(con_med_only_cluster10_long_order, x= rep(1:8,9), Cond_c10_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster10_long_order_x_label <- con_med_only_cluster10_long_order_x                         # Modify data
con_cluster10_long_order_x_label$label <- NA
con_cluster10_long_order_x_label$label[which(con_cluster10_long_order_x_label$x == max(con_cluster10_long_order_x_label$x))] <-con_cluster10_long_order_x_label$gene[which(con_cluster10_long_order_x_label$x == max(con_cluster10_long_order_x_label$x))]

ggplot(con_cluster10_long_order_x_label,
       aes(x= Cond_c10_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #10",x =" Time after bisection", y = "Transcripts per million (TPM)")


### Multi Line graph for Cluster 11
con_clust11_char <- as.character(con_clust11$gene) 
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
con_med_only_cluster11 = subset(con_med_exp_only_rows, rownames(con_med_exp_only_rows) %in% con_clust11_char)

# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
con_med_only_cluster11_long <- con_med_only_cluster11 %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
con_med_only_cluster11_long_order <- con_med_only_cluster11_long[order(con_med_only_cluster11_long$gene),]

Cond_c11 <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),4) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_c11_factor <- factor(Cond_c11, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
con_med_only_cluster11_long_order_x <- cbind(con_med_only_cluster11_long_order, x= rep(1:8,4), Cond_c11_factor)


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
con_cluster11_long_order_x_label <- con_med_only_cluster11_long_order_x                         # Modify data
con_cluster11_long_order_x_label$label <- NA
con_cluster11_long_order_x_label$label[which(con_cluster11_long_order_x_label$x == max(con_cluster11_long_order_x_label$x))] <-con_cluster11_long_order_x_label$gene[which(con_cluster11_long_order_x_label$x == max(con_cluster11_long_order_x_label$x))]

ggplot(con_cluster11_long_order_x_label,
       aes(x= Cond_c11_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Cluster #11",x =" Time after bisection", y = "Transcripts per million (TPM)")



#BLASTP on Clusters 
#ML gene models and their protein blast hit 
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)
length(HumRef2020ML_df$V1 )


## Cluster 1 
con_clust1_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust1_char))                                      

## Cluster 2 
con_clust2_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust2_char))                                  

## Cluster 3 
con_clust3_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust3_char))                                  

## Cluster 4 
con_clust4_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust4_char))                                  

## Cluster 5 
con_clust5_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust5_char))                                  

## Cluster 6 
con_clust6_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust6_char))                                  

## Cluster 7 
con_clust7_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust7_char))                                  

## Cluster 8 
con_clust8_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust8_char))                                  


## Cluster 9 
con_clust9_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust9_char))                                  


## Cluster 10 
con_clust10_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust10_char))                                  

## Cluster 11
con_clust11_blastp <- subset.data.frame(HumRef2020ML_df, (rownames(HumRef2020ML_df) %in% con_clust11_char))                                  

con_clust_blastp <- list(clust1 = con_clust1_blastp,
                         clust2 = con_clust2_blastp,
                         clust3 = con_clust3_blastp,
                         clust4 = con_clust4_blastp,
                         clust5 = con_clust5_blastp,
                         clust6 = con_clust6_blastp,
                         clust7 = con_clust7_blastp,
                         clust8 = con_clust8_blastp,
                         clust9 = con_clust9_blastp,
                         clust10 = con_clust10_blastp,
                         clust11= con_clust11_blastp)


#export list 
library(data.table)
# code from https://stackoverflow.com/questions/27594541/export-a-list-into-a-csv-or-txt-file-in-r
outputfile <- "con_clust_blastp.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(nam in names(con_clust_blastp)){
  fwrite(list(nam), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- con_clust_blastp[[nam]]
  if(is.list(ele)) fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T, row.names = F) else fwrite(data.frame(matrix(ele, nrow=1)), file=outputfile, append=T) #write elements of the list
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}



#Multi line of all of the tubulins in the gene models -- supplement 
ML01482a	TUBA3D
ML021133a	TUBA1C	
ML026516a	TUBA1A
ML03349a	TUBA1C
ML056958a	TUBA1C
ML06742a	TUBA1C
ML085213a	TUBA1C
ML10374a	TUBA1A
ML10375a	TUBA1A
ML10376a	TUBA1C
ML10373a	TUBA3D
ML10377a	TUBA3D
ML10942a	TUBA1C
ML154172a	TUBA1C
ML174731a	TUBA1A
ML20831a	TUBA1A
ML234515a	TUBA1C
ML00482a	TUBA1C

# Using the pipe operator %>%
all_med_exp <- RSEMcounts_B2_Trm_TPM  %>% 
  rownames_to_column('gene') %>% #switches row names to a column because the tibble format will remove row names from a standard DF 
  rowwise() %>% # calculate the median at each row within the specified columns 
  mutate(Uncut_Med = median(uncut_1_i,uncut_2_i,uncut_3_i),
         Hab00_Med = median(hab00_1,hab00_2,hab00_3),
         Hab01_Med = median(hab01_1,hab01_2,hab01_3),
         Hab03_Med = median(hab03_1,hab03_2,hab03_3),
         Hab06_Med = median(hab06_1,hab06_2,hab06_3),
         Hab12_Med = median(hab12_1,hab12_2,hab12_3),
         Hab24_Med = median(hab24_1_i,hab24_2_i,hab24_3_i),
         Hab48_Med= median(hab48_1_i,hab48_2_i,hab48_3_i))

all_med_exp_only = subset(all_med_exp , select = c(gene,Uncut_Med,Hab00_Med,Hab01_Med,Hab03_Med,Hab06_Med,Hab12_Med,Hab24_Med,Hab48_Med )) 

all_med_exp_only_rows <- all_med_exp_only %>% 
  column_to_rownames('gene') # switch it back to a data frame and put back in the row names



tubs_char <- c("ML01482a","ML021133a","ML026516a","ML03349a","ML03349a","ML056958a","ML06742a","ML085213a","ML10374a",
               "ML10375a","ML10376a","ML10373a","ML10377a","ML10942a","ML154172a","ML174731a","ML20831a","ML234515a","ML00482a")
# Need to retrieve the  median TPMS associated with each gene in the specified cluster 
all_med_only_tubs = subset(all_med_exp_only_rows, rownames(all_med_exp_only_rows) %in% tubs_char )
length(rownames(all_med_only_tubs)) #18
# I want to subset the tibble I just made to ONLY
#include the columns for each gene that show the median TPM values for each time point 

#Need to reformat the data frame short to long format  
all_med_only_tub_long <- all_med_only_tubs %>% 
  rownames_to_column("gene") %>%
  gather(conditions, values,-gene)

#reorder by gene name 
all_med_only_tub_long_order <- all_med_only_tub_long[order(all_med_only_tub_long$gene),]

Cond_tubs <- rep(c("Uncut","10m","1h","3h","6h","12h","24h","48h"),18 ) ## each = x is the total # of genes in the cluster, so each condition will be repeated x * 8 conditions
Cond_tubs_factor <- factor(Cond_tubs, levels = c("Uncut","10m","1h","3h","6h","12h","24h","48h"))          


#add column to show x position
all_med_only_tub_long_order_x <- cbind(all_med_only_tub_long_order , x= rep(1:8,18),Cond_tubs_factor )


#Label the maximum data point so there is one label for each line -- as described https://statisticsglobe.com/add-labels-at-ends-of-lines-in-ggplot2-line-plot-r
all_med_only_tub_long_order_x_label <- all_med_only_tub_long_order_x                        # Modify data
all_med_only_tub_long_order_x_label$label <- NA
all_med_only_tub_long_order_x_label$label[which(all_med_only_tub_long_order_x_label$x == max(all_med_only_tub_long_order_x_label$x))] <-all_med_only_tub_long_order_x_label$gene[which(all_med_only_tub_long_order_x_label$x == max(all_med_only_tub_long_order_x_label$x))]

ggplot(all_med_only_tub_long_order_x_label,
       aes(x= Cond_tubs_factor, y=values, group = gene, color = gene)) + geom_line() + geom_label_repel(aes(label = label),na.rm = T,max.overlaps = 100, nudge_x = 5) + theme(text = element_text(size = 20)) + labs(title="Alpha Tubulin Genes",x =" Time after bisection", y = "Transcripts per million (TPM)")


