##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

#R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
#Copyright (C) 2023 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin20 (64-bit)


#Box and Whisker visualization of Genes -- we want to visualize genes with their TPM values(normalized from RSEM)     

## TPMS = Transcripts per million 


#####################


# Load packages
install.packages("tidyverse")
library(ggplot2)
library(dplyr)
library(RColorBrewer)

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








## binding the TPM count data togethe for all samples 

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


RSEMcounts_B2_Trm_TPM_mt <- as.matrix(RSEMcounts_B2_Trm_TPM)

options(max.print=999999)
capture.output(RSEMcounts_B2_Trm_TPM_mt, file = "RSEMcounts_B2_Trm_TPM.csv")

#Organizing labels for x axis 
mygroups = c("Uncut", "Uncut","Uncut","10m","10m","10m","1h","1h","1h","3h","3h","3h","6h","6h","6h","12h","12h","12h","24h","24h","24h","48h","48h","48h")
mygroups_order <- factor(mygroups, levels=c("Uncut", "10m", "1h", "3h","6h","12h","24h","48h"))


# Model Graph for Figure 
Model_Data <- data.frame(
        mygroups_order,
        expression = rep(c(5,10,15),8)
)

boxplot(Model_Data$expression ~ Model_Data$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main="Key",
        ylim = c(0, 20)
)




## Genes of Interest from NOISeq TOPGO analysis 
##GO Figure results 
## Transcription factors(GO.0003700) 

## ML182032a -- CREB5(ML_Fos3)
ML182032a_FOS_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML182032a",])
ML182032a_FOS_df_groups <-data.frame( group = rownames(ML182032a_FOS_df), expression = ML182032a_FOS_df$`RSEMcounts_B2_Trm_TPM_mt["ML182032a", ]`,mygroups_order = mygroups_order)



boxplot(ML182032a_FOS_df_groups$expression ~ ML182032a_FOS_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main="ML_Fos3-ML182032a"
)

## ML282527a -- ELK1(ML_Etslx4)
ML282527a_ELK_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML282527a",])
ML282527a_ELK_df_groups <- data.frame( group = rownames(ML282527a_ELK_df), expression = ML282527a_ELK_df$`RSEMcounts_B2_Trm_TPM_mt["ML282527a", ]`, mygroups_order = mygroups_order)



boxplot(ML282527a_ELK_df_groups$expression ~ ML282527a_ELK_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ML_Etslx4-ML282527a"
)



##ML1541120a- JUN(ML_Jun)
ML1541120a_JUN_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML1541120a",])
ML1541120a_JUN_df_groups <- data.frame( group = rownames(ML1541120a_JUN_df), expression = ML1541120a_JUN_df$`RSEMcounts_B2_Trm_TPM_mt["ML1541120a", ]`, mygroups_order = mygroups_order)



boxplot(ML1541120a_JUN_df_groups$expression ~ ML1541120a_JUN_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ML_Jun-ML1541120a"
)




##  ML057318a -- ATF2
ML057318a_ATF2_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML057318a",])
ML057318a_ATF2_df_groups <- data.frame( group = rownames(ML057318a_ATF2_df), expression = ML057318a_ATF2_df$`RSEMcounts_B2_Trm_TPM_mt["ML057318a", ]`, mygroups_order = mygroups_order)



boxplot(ML057318a_ATF2_df_groups$expression ~ ML057318a_ATF2_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ATF2-ML057318a"
)
## ML077623a -- CREM 
ML077623a_CREM_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML077623a",])
ML077623a_CREM_df_groups <- data.frame( group = rownames(ML077623a_CREM_df), expression = ML077623a_CREM_df$`RSEMcounts_B2_Trm_TPM_mt["ML077623a", ]` , mygroups_order = mygroups_order)



boxplot(ML077623a_CREM_df_groups$expression ~ ML077623a_CREM_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CREM-ML077623a"
)


## ML015722a -- CREB3L3
ML015722a_CREB3L3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML015722a",])
ML015722a_CREB3L3_df_groups <- data.frame( group = rownames(ML015722a_CREB3L3_df), expression = ML015722a_CREB3L3_df$`RSEMcounts_B2_Trm_TPM_mt["ML015722a", ]` , mygroups_order = mygroups_order)



boxplot(ML015722a_CREB3L3_df_groups$expression ~ ML015722a_CREB3L3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CREB3L3-ML015722a"
)
## ML08021a -- ATF6B 

ML08021a_ATF6B_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML08021a",])
ML08021a_ATF6B_df_groups <- data.frame( group = rownames(ML08021a_ATF6B_df), expression =  ML08021a_ATF6B_df$`RSEMcounts_B2_Trm_TPM_mt["ML08021a", ]`, mygroups_order = mygroups_order)



boxplot(ML08021a_ATF6B_df_groups$expression ~ ML08021a_ATF6B_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ATF6B-ML08021a"
)


## ML177213a-- GPBP1L1 

ML177213a_GPBP1L1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML177213a",])
ML177213a_GPBP1L1_df_groups <- data.frame( group = rownames(ML177213a_GPBP1L1_df), expression =  ML177213a_GPBP1L1_df$`RSEMcounts_B2_Trm_TPM_mt["ML177213a", ]`, mygroups_order = mygroups_order)



boxplot(ML177213a_GPBP1L1_df_groups$expression ~ ML177213a_GPBP1L1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "GPBP1L1-ML177213a"
)



##ML09109a-- ETV4 (ML_Etslx3)

ML09109a_ETV4_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML09109a",])
ML09109a_ETV4_df_groups <- data.frame( group = rownames(ML09109a_ETV4_df), expression = ML09109a_ETV4_df$`RSEMcounts_B2_Trm_TPM_mt["ML09109a", ]`, mygroups_order = mygroups_order)



boxplot(ML09109a_ETV4_df_groups$expression ~ ML09109a_ETV4_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ML_Etslx3-ML09109a"
)



#Peptidases GO.0004222 & GO.0003824
## ML007435a -- TLL2 (Tolloid)
ML007435a_TLL_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML007435a",])
ML007435a_TLL_df_groups <- data.frame( group = rownames(ML007435a_TLL_df), expression =ML007435a_TLL_df$`RSEMcounts_B2_Trm_TPM_mt["ML007435a", ]`, mygroups_order = mygroups_order)



boxplot(ML007435a_TLL_df_groups$expression ~ ML007435a_TLL_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TLL2-ML007435a"
)

## ML305524a -- MMP14 
ML305524a_MMP14_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML305524a",])
ML305524a_MMP14_df_groups <- data.frame( group = rownames(ML305524a_MMP14_df), expression =ML305524a_MMP14_df$`RSEMcounts_B2_Trm_TPM_mt["ML305524a", ]`, mygroups_order = mygroups_order)



boxplot(ML305524a_MMP14_df_groups$expression ~ ML305524a_MMP14_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP14-ML305524a"
)

#Peptidases in 
#ML31401a--TLL1 
ML31401a_TLL1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML31401a",])
ML31401a_TLL1_df_groups <- data.frame( group = rownames(ML31401a_TLL1_df), expression =ML31401a_TLL1_df$`RSEMcounts_B2_Trm_TPM_mt["ML31401a", ]`, mygroups_order = mygroups_order)



boxplot(ML31401a_TLL1_df_groups$expression ~ ML31401a_TLL1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TLL1-ML31401a"
)

#ML154125a -CASP3
ML154125a_CASP3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML154125a",])
ML154125a_CASP3_df_groups <- data.frame( group = rownames(ML154125a_CASP3_df), expression =ML154125a_CASP3_df$`RSEMcounts_B2_Trm_TPM_mt["ML154125a", ]`, mygroups_order = mygroups_order)



boxplot(ML154125a_CASP3_df_groups$expression ~ ML154125a_CASP3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CASP3-ML154125a"
)



## ML14875a -- MMP24 -- Middle 
ML14875a_MMP24_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML14875a",])
ML14875a_MMP24_df_groups <- data.frame( group = rownames(ML14875a_MMP24_df), expression =ML14875a_MMP24_df$`RSEMcounts_B2_Trm_TPM_mt["ML14875a", ]`, mygroups_order = mygroups_order)



boxplot(ML14875a_MMP24_df_groups$expression ~ ML14875a_MMP24_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP24-ML14875a"
)





#ML073224a --MMP28
ML073224a_MMP28_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML073224a",])
ML073224a_MMP28_df_groups <- data.frame( group = rownames(ML073224a_MMP28_df), expression =ML073224a_MMP28_df$`RSEMcounts_B2_Trm_TPM_mt["ML073224a", ]`, mygroups_order = mygroups_order)



boxplot(ML073224a_MMP28_df_groups$expression ~ ML073224a_MMP28_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP28-ML073224a"
)


#ML13379a    MMP28 
ML13379a_MMP28_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML13379a",])
ML13379a_MMP28_df_groups <- data.frame( group = rownames(ML13379a_MMP28_df), expression =ML13379a_MMP28_df$`RSEMcounts_B2_Trm_TPM_mt["ML13379a", ]`, mygroups_order = mygroups_order)



boxplot(ML13379a_MMP28_df_groups$expression ~ ML13379a_MMP28_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP28-ML13379a"
)


#ML14713a-MMP1
ML14713a_MMP1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML14713a",])
ML14713a_MMP1_df_groups <- data.frame( group = rownames(ML14713a_MMP1_df ), expression =ML14713a_MMP1_df$`RSEMcounts_B2_Trm_TPM_mt["ML14713a", ]`, mygroups_order = mygroups_order)



boxplot(ML14713a_MMP1_df_groups$expression ~ ML14713a_MMP1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP1-ML14713a"
)


#ML33825a-MMP28
ML33825a_MMP28_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML33825a",])
ML33825a_MMP28_df_groups <- data.frame( group = rownames(ML33825a_MMP28_df), expression = ML33825a_MMP28_df$`RSEMcounts_B2_Trm_TPM_mt["ML33825a", ]`, mygroups_order = mygroups_order)



boxplot(ML33825a_MMP28_df_groups$expression ~ ML33825a_MMP28_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP28-ML33825a"
)



# Tubulins - strucural consituent of the cytoskeleton(GO:0005200)

## ML026516a -- TUBA1A
ML026516a_TUBA1A_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML026516a",])
ML026516a_TUBA1A_df_groups <- data.frame( group = rownames(ML026516a_TUBA1A_df), expression = ML026516a_TUBA1A_df$`RSEMcounts_B2_Trm_TPM_mt["ML026516a", ]`, mygroups_order = mygroups_order)



boxplot(ML026516a_TUBA1A_df_groups$expression ~ ML026516a_TUBA1A_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA1A-ML026516a"
)





## ML06742a -- TUBA1C
ML06742a_TUBA1C_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML06742a",])
ML06742a_TUBA1C_df_groups <- data.frame( group = rownames(ML06742a_TUBA1C_df), expression = ML06742a_TUBA1C_df$`RSEMcounts_B2_Trm_TPM_mt["ML06742a", ]`, mygroups_order = mygroups_order)



boxplot(ML06742a_TUBA1C_df_groups$expression ~ ML06742a_TUBA1C_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA1C-ML06742a"
)




#ML056958a TUBA1C 
ML056958a_TUBA1C_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML056958a",])
ML056958a_TUBA1C_df_groups <- data.frame( group = rownames(ML056958a_TUBA1C_df), expression =ML056958a_TUBA1C_df$`RSEMcounts_B2_Trm_TPM_mt["ML056958a", ]` , mygroups_order = mygroups_order)



boxplot(ML056958a_TUBA1C_df_groups$expression ~ ML056958a_TUBA1C_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA1C-ML056958a"
)





## ML01482a -- TUBA3D 
ML01482a_TUBA3D_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML01482a",])
ML01482a_TUBA3D_df_groups <- data.frame( group = rownames(ML01482a_TUBA3D_df), expression =ML01482a_TUBA3D_df$`RSEMcounts_B2_Trm_TPM_mt["ML01482a", ]` , mygroups_order = mygroups_order)



boxplot(ML01482a_TUBA3D_df_groups$expression ~ ML01482a_TUBA3D_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA3D-ML01482a"
)






#ML24281a-TUBB2B 
ML24281a_TUBB2B_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML24281a",])
ML24281a_TUBB2B_df_groups <- data.frame( group = rownames(ML24281a_TUBB2B_df ), expression =ML24281a_TUBB2B_df$`RSEMcounts_B2_Trm_TPM_mt["ML24281a", ]`, mygroups_order = mygroups_order)



boxplot(ML24281a_TUBB2B_df_groups$expression ~ ML24281a_TUBB2B_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBB2B-ML10377a"
)




#Supplement


#6h-12h UP ribosomal associated genes(GO.0005198) structural molecular activity 
#ML00985a  RPL34 
ML00985a_RPL34_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML00985a",])
ML00985a_RPL34_df_groups <- data.frame( group = rownames(ML00985a_RPL34_df), expression = ML00985a_RPL34_df$`RSEMcounts_B2_Trm_TPM_mt["ML00985a", ]`, mygroups_order = mygroups_order)



boxplot(ML00985a_RPL34_df_groups$expression ~ML00985a_RPL34_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "RPL34-ML00985a "
)


#ML22072a  RPS17 
ML22072a_RPS17_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML22072a",])
ML22072a_RPS17_df_groups <- data.frame( group = rownames(ML22072a_RPS17_df), expression = ML22072a_RPS17_df$`RSEMcounts_B2_Trm_TPM_mt["ML22072a", ]`, mygroups_order = mygroups_order)



boxplot(ML22072a_RPS17_df_groups$expression ~ML22072a_RPS17_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "RPS17-ML22072a "
)

#ML358815a  RPL31
ML358815a_RPL31_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML358815a",])
ML358815a_RPL31_df_groups <- data.frame( group = rownames(ML358815a_RPL31_df), expression = ML358815a_RPL31_df$`RSEMcounts_B2_Trm_TPM_mt["ML358815a", ]`, mygroups_order = mygroups_order)



boxplot(ML358815a_RPL31_df_groups$expression ~ML358815a_RPL31_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "RPL31-ML358815a "
)


#ML03874a  UBA52
ML03874a_UBA52_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML03874a",])
ML03874a_UBA52_df_groups <- data.frame( group = rownames(ML03874a_UBA52_df ), expression = ML03874a_UBA52_df$`RSEMcounts_B2_Trm_TPM_mt["ML03874a", ]`, mygroups_order = mygroups_order)



boxplot(ML03874a_UBA52_df_groups$expression ~ML03874a_UBA52_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "UBA52-ML03874a "
)



# 3h-6h down in catalytic activity (GO.0003824)

#ML218826a     GLUL
ML218826a_GLUL_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML218826a",])
ML218826a_GLUL_df_groups <- data.frame( group = rownames(ML218826a_GLUL_df ), expression =ML218826a_GLUL_df$`RSEMcounts_B2_Trm_TPM_mt["ML218826a", ]`, mygroups_order = mygroups_order)



boxplot(ML218826a_GLUL_df_groups$expression ~ML218826a_GLUL_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "GLUL-ML218826a"
)

#ML05027a GLUL
ML05027a_GLUL_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML05027a",])
ML05027a_GLUL_df_groups <- data.frame( group = rownames(ML05027a_GLUL_df ), expression =ML05027a_GLUL_df$`RSEMcounts_B2_Trm_TPM_mt["ML05027a", ]`, mygroups_order = mygroups_order)



boxplot(ML05027a_GLUL_df_groups$expression ~ML05027a_GLUL_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "GLUL-ML05027a "
)
#ML17724a    GLUD2
ML17724a_GLUD2_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML17724a",])
ML17724a_GLUD2_df_groups <- data.frame( group = rownames(ML17724a_GLUD2_df ), expression = ML17724a_GLUD2_df$`RSEMcounts_B2_Trm_TPM_mt["ML17724a", ]`, mygroups_order = mygroups_order)



boxplot(ML17724a_GLUD2_df_groups$expression ~ML17724a_GLUD2_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "GLUD2-ML17724a "
)

#ML05654a GLUL
ML05654a_GLUL_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML05654a",])
ML05654a_GLUL_df_groups <- data.frame( group = rownames(ML05654a_GLUL_df ), expression = ML05654a_GLUL_df$`RSEMcounts_B2_Trm_TPM_mt["ML05654a", ]`, mygroups_order = mygroups_order)



boxplot(ML05654a_GLUL_df_groups$expression ~ML05654a_GLUL_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "GLUL-ML05654a "
)




#ML310320a MAPKAPK2
ML310320a_MAPKAPK2_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML310320a",])
ML310320a_MAPKAPK2_df_groups <- data.frame( group = rownames(ML310320a_MAPKAPK2_df ), expression = ML310320a_MAPKAPK2_df$`RSEMcounts_B2_Trm_TPM_mt["ML310320a", ]`, mygroups_order = mygroups_order)



boxplot(ML310320a_MAPKAPK2_df_groups$expression ~ML310320a_MAPKAPK2_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MAPKAPK2-ML310320a "
)
#ML07243a   MAP3K3 
ML07243a_MAP3K3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML07243a",])
ML07243a_MAP3K3_df_groups <- data.frame( group = rownames(ML07243a_MAP3K3_df), expression = ML07243a_MAP3K3_df$`RSEMcounts_B2_Trm_TPM_mt["ML07243a", ]`, mygroups_order = mygroups_order)



boxplot(ML07243a_MAP3K3_df_groups$expression ~ML07243a_MAP3K3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MAP3K3-ML07243a "
)



#ML07242a   MAP3K5
ML07242a_MAP3K5_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML07242a",])
ML07242a_MAP3K5_df_groups <- data.frame( group = rownames(ML07242a_MAP3K5_df), expression = ML07242a_MAP3K5_df$`RSEMcounts_B2_Trm_TPM_mt["ML07242a", ]`, mygroups_order = mygroups_order)



boxplot(ML07242a_MAP3K5_df_groups$expression ~ML07242a_MAP3K5_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = " MAP3K5-ML07242a "
)





#ML08261a--MAPK9
ML08261a_MAPK9_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML08261a",])
ML08261a_MAPK9_df_groups <- data.frame( group = rownames(ML08261a_MAPK9_df), expression = ML08261a_MAPK9_df$`RSEMcounts_B2_Trm_TPM_mt["ML08261a", ]`, mygroups_order = mygroups_order)



boxplot(ML08261a_MAPK9_df_groups$expression ~ML08261a_MAPK9_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MAPK9-ML08261a "
)



# Uncut-10m up Chitin binding (GO.0008061)
#ML03134a CHI3L1
ML03134a_CHI3L1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML03134a",])
ML03134a_CHI3L1_df_groups <- data.frame( group = rownames(ML03134a_CHI3L1_df), expression = ML03134a_CHI3L1_df$`RSEMcounts_B2_Trm_TPM_mt["ML03134a", ]`, mygroups_order = mygroups_order)



boxplot(ML03134a_CHI3L1_df_groups$expression ~ML03134a_CHI3L1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CHI3L1-ML03134a "
)

#ML368913a  CHIT1
ML368913a_CHIT1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML368913a",])
ML368913a_CHIT1_df_groups <- data.frame( group = rownames(ML368913a_CHIT1_df), expression = ML368913a_CHIT1_df$`RSEMcounts_B2_Trm_TPM_mt["ML368913a", ]`, mygroups_order = mygroups_order)



boxplot(ML368913a_CHIT1_df_groups$expression ~ ML368913a_CHIT1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CHIT1-ML368913a "
)

#Uncut-10m Down - Ca ion binding 
#ML085730b- MlPhotoprotein2
ML085730b_MlPhotoprotein2_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML085730b",])
ML085730b_MlPhotoprotein2_df_groups <- data.frame( group = rownames(ML085730b_MlPhotoprotein2_df), expression = ML085730b_MlPhotoprotein2_df$`RSEMcounts_B2_Trm_TPM_mt["ML085730b", ]`, mygroups_order = mygroups_order)



boxplot(ML085730b_MlPhotoprotein2_df_groups$expression ~ ML085730b_MlPhotoprotein2_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MlPhotoprotein2-ML085730b "
)

#ML085733a MlPhotoprotein5
ML085733a_MlPhotoprotein5_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML085733a",])
ML085733a_MlPhotoprotein5_df_groups <- data.frame( group = rownames(ML085733a_MlPhotoprotein5_df), expression = ML085733a_MlPhotoprotein5_df$`RSEMcounts_B2_Trm_TPM_mt["ML085733a", ]`, mygroups_order = mygroups_order)



boxplot(ML085733a_MlPhotoprotein5_df_groups$expression ~ ML085733a_MlPhotoprotein5_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MlPhotoprotein5-ML085733a "
)

#G-coupled protein receptor - 3-6h UP
#ML01515a  ARRDC3 
ML01515a_ARRDC3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML01515a",])
ML01515a_ARRDC3_df_groups <- data.frame( group = rownames(ML01515a_ARRDC3_df), expression = ML01515a_ARRDC3_df$`RSEMcounts_B2_Trm_TPM_mt["ML01515a", ]`, mygroups_order = mygroups_order)



boxplot(ML01515a_ARRDC3_df_groups$expression ~ ML01515a_ARRDC3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ARRDC3-ML01515a "
)

#ML01516a  ADGRG7 
ML01516a_ADGRG7_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML01516a",])
ML01516a_ADGRG7_df_groups <- data.frame( group = rownames(ML01516a_ADGRG7_df), expression = ML01516a_ADGRG7_df$`RSEMcounts_B2_Trm_TPM_mt["ML01516a", ]`, mygroups_order = mygroups_order)



boxplot(ML01516a_ADGRG7_df_groups$expression ~ ML01516a_ADGRG7_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ADGRG7-ML01516a "
)



# WNT Pathway components  
##ML11223a- TCF
ML11223a_TCF_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML11223a",])
ML11223a_TCF_df_groups <- data.frame( group = rownames(ML11223a_TCF_df), expression = ML11223a_TCF_df$`RSEMcounts_B2_Trm_TPM_mt["ML11223a", ]`, mygroups_order = mygroups_order)



boxplot(ML11223a_TCF_df_groups$expression ~ ML11223a_TCF_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TCF-ML11223a"
)


##ML073715a - beta catenin 
ML073715a_bcat_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML073715a",])
ML073715a_bcat_df_groups <- data.frame( group = rownames(ML073715a_bcat_df), expression = ML073715a_bcat_df$`RSEMcounts_B2_Trm_TPM_mt["ML073715a", ]`  , mygroups_order = mygroups_order)


boxplot(ML073715a_bcat_df_groups$expression ~ ML073715a_bcat_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "Beta-catenin-ML073715a"
)












#ML01051a - WntX
ML01051a_WntX_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML01051a",])
ML01051a_WntX_df_groups <- data.frame( group = rownames(ML01051a_WntX_df), expression = ML01051a_WntX_df$`RSEMcounts_B2_Trm_TPM_mt["ML01051a", ]`  , mygroups_order = mygroups_order)


boxplot(ML01051a_WntX_df_groups$expression ~ ML01051a_WntX_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "WntX-ML01051a"
)

