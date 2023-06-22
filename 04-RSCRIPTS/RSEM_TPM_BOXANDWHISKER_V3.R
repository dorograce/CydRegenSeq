##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

##R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##Copyright (C) 2021 The R Foundation for Statistical Computing
##Platform: x86_64-apple-darwin17.0 (64-bit)

#Box and Whisker visualization of Genes -- we want to visualize genes with their TPM values(normalized from RSEM)     

## TPMS = Transcripts per million 


#####################


# Load packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)

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






## binding the TPM count data togethe for all samples 

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


RSEMcounts_B2_Trm_TPM_mt <- as.matrix(RSEMcounts_B2_Trm_TPM)

#Organizing labels for x axis 
mygroups = c("Uncut", "Uncut","Uncut","00","00","00","01","01","01","03","03","03","06","06","06","12","12","12","24","24","24","48","48","48")
mygroups_order <- factor(mygroups, levels=c("Uncut", "00", "01", "03","06","12","24","48"))


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




## Genes of Interest 

# DNA-binding transcription factor activity -- "GO:0003700"
## ML182032a -- FOS
ML182032a_FOS_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML182032a",])
ML182032a_FOS_df_groups <-data.frame( group = rownames(ML182032a_FOS_df), expression = ML182032a_FOS_df$`RSEMcounts_B2_Trm_TPM_mt["ML182032a", ]`,mygroups_order = mygroups_order)



boxplot(ML182032a_FOS_df_groups$expression ~ ML182032a_FOS_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main="FOS"
)

## ML282527a -- ELk1
ML282527a_ELK_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML282527a",])
ML282527a_ELK_df_groups <- data.frame( group = rownames(ML282527a_ELK_df), expression = ML282527a_ELK_df$`RSEMcounts_B2_Trm_TPM_mt["ML282527a", ]`, mygroups_order = mygroups_order)



boxplot(ML282527a_ELK_df_groups$expression ~ ML282527a_ELK_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "Elk"
)


## ML09109a-- ETV4
ML09109a_ETV_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML09109a",])
ML09109a_ETV_df_groups <- data.frame( group = rownames(ML09109a_ETV_df), expression = ML09109a_ETV_df$`RSEMcounts_B2_Trm_TPM_mt["ML09109a", ]` , mygroups_order = mygroups_order)



boxplot(ML09109a_ETV_df_groups$expression ~ ML09109a_ETV_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ETV"
)

##ML1541120a- JUN
ML1541120a_JUN_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML1541120a",])
ML1541120a_JUN_df_groups <- data.frame( group = rownames(ML1541120a_JUN_df), expression = ML1541120a_JUN_df$`RSEMcounts_B2_Trm_TPM_mt["ML1541120a", ]`, mygroups_order = mygroups_order)



boxplot(ML1541120a_JUN_df_groups$expression ~ ML1541120a_JUN_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "JUN"
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
        main = "ATF2"
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
        main = "CREM"
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
        main = "CREB3L3"
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
        main = "ATF6B"
)

##ML23712a  -- FOXJ1 

ML23712a_FOXJ1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML23712a",])
ML23712a_FOXJ1_df_groups <- data.frame( group = rownames(ML23712a_FOXJ1_df), expression =  ML23712a_FOXJ1_df$`RSEMcounts_B2_Trm_TPM_mt["ML23712a", ]`, mygroups_order = mygroups_order)



boxplot(ML23712a_FOXJ1_df_groups$expression ~ ML23712a_FOXJ1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "FOXJ1"
)



# structural constituent of cytoskeleton -- "GO:0005200"
## ML026516a -- TUBA1A
ML026516a_TUBA1A_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML026516a",])
ML026516a_TUBA1A_df_groups <- data.frame( group = rownames(ML026516a_TUBA1A_df), expression = ML026516a_TUBA1A_df$`RSEMcounts_B2_Trm_TPM_mt["ML026516a", ]`, mygroups_order = mygroups_order)



boxplot(ML026516a_TUBA1A_df_groups$expression ~ ML026516a_TUBA1A_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA1A"
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
        main = "TUBA1C"
)


## ML01482a -- TUBA3D - 1
ML01482a_TUBA3D_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML01482a",])
ML01482a_TUBA3D_df_groups <- data.frame( group = rownames(ML01482a_TUBA3D_df), expression =ML01482a_TUBA3D_df$`RSEMcounts_B2_Trm_TPM_mt["ML01482a", ]` , mygroups_order = mygroups_order)



boxplot(ML01482a_TUBA3D_df_groups$expression ~ ML01482a_TUBA3D_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA3D-1"
)
## ML10373a -- TUBA3D - 2
ML10373a_TUBA3D_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML10373a",])
ML10373a_TUBA3D_df_groups <- data.frame( group = rownames(ML10373a_TUBA3D_df), expression = ML10373a_TUBA3D_df$`RSEMcounts_B2_Trm_TPM_mt["ML10373a", ]` , mygroups_order = mygroups_order)



boxplot(ML10373a_TUBA3D_df_groups$expression ~ ML10373a_TUBA3D_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TUBA3D-2"
)
##  scavenger receptor activity -- " GO:0005044"
## ML11643a -- PRSS22
ML11643a_PRSS22_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML11643a",])
ML11643a_PRSS22_df_groups <- data.frame( group = rownames(ML11643a_PRSS22_df), expression = ML11643a_PRSS22_df$`RSEMcounts_B2_Trm_TPM_mt["ML11643a", ]`, mygroups_order = mygroups_order)



boxplot(ML11643a_PRSS22_df_groups$expression ~ ML11643a_PRSS22_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "PRSS22"
)

## ML005015a -- DMBT1 - 1
ML005015a_DMBT1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML005015a",])
ML005015a_DMBT1_df_groups <- data.frame( group = rownames(ML005015a_DMBT1_df), expression = ML005015a_DMBT1_df$`RSEMcounts_B2_Trm_TPM_mt["ML005015a", ]`, mygroups_order = mygroups_order)



boxplot(ML005015a_DMBT1_df_groups$expression ~ ML005015a_DMBT1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "DMBT1-1"
)

## ML00576a -- TMPRSS3
ML00576a_TMPRSS3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML00576a",])
ML00576a_TMPRSS3_df_groups <- data.frame( group = rownames(ML00576a_TMPRSS3_df), expression = ML00576a_TMPRSS3_df$`RSEMcounts_B2_Trm_TPM_mt["ML00576a", ]`, mygroups_order = mygroups_order)



boxplot(ML00576a_TMPRSS3_df_groups$expression ~ ML00576a_TMPRSS3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TMPRSS3"
)
## ML006919a -- DMBT1 - 2
ML006919a_DMBT1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML006919a",])
ML006919a_DMBT1_df_groups <- data.frame( group = rownames(ML006919a_DMBT1_df), expression = ML006919a_DMBT1_df$`RSEMcounts_B2_Trm_TPM_mt["ML006919a", ]`, mygroups_order = mygroups_order)



boxplot(ML006919a_DMBT1_df_groups$expression ~ ML006919a_DMBT1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "DMBT1-2"
)

## ML20572a -- CR1L
ML20572a_CR1L_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML20572a",])
ML20572a_CR1L_df_groups <- data.frame( group = rownames(ML20572a_CR1L_df), expression = ML20572a_CR1L_df$`RSEMcounts_B2_Trm_TPM_mt["ML20572a", ]`, mygroups_order = mygroups_order)



boxplot(ML20572a_CR1L_df_groups$expression ~ ML20572a_CR1L_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CR1L"
)


## ML02243a -- DMBT1 - 3
ML02243a_DMBT1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML02243a",])
ML02243a_DMBT1_df_groups <- data.frame( group = rownames(ML02243a_DMBT1_df), expression =ML02243a_DMBT1_df$`RSEMcounts_B2_Trm_TPM_mt["ML02243a", ]`, mygroups_order = mygroups_order)



boxplot(ML02243a_DMBT1_df_groups$expression ~ ML02243a_DMBT1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "DMBT1-3"
)


# metallopeptidase activity -- "GO:0008237"
## ML007435a -- TLL2 (Tolloid)
ML007435a_TLL_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML007435a",])
ML007435a_TLL_df_groups <- data.frame( group = rownames(ML007435a_TLL_df), expression =ML007435a_TLL_df$`RSEMcounts_B2_Trm_TPM_mt["ML007435a", ]`, mygroups_order = mygroups_order)



boxplot(ML007435a_TLL_df_groups$expression ~ ML007435a_TLL_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TLL"
)
## ML030513a -- ERAP2
ML030513a_ERAP2_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML030513a",])
ML030513a_ERAP2_df_groups <- data.frame( group = rownames(ML030513a_ERAP2_df), expression =ML030513a_ERAP2_df$`RSEMcounts_B2_Trm_TPM_mt["ML030513a", ]`, mygroups_order = mygroups_order)



boxplot(ML030513a_ERAP2_df_groups$expression ~ ML030513a_ERAP2_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "ERAP2"
)

## ML14875a -- MMP24
ML14875a_MMP24_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML14875a",])
ML14875a_MMP24_df_groups <- data.frame( group = rownames(ML14875a_MMP24_df), expression =ML14875a_MMP24_df$`RSEMcounts_B2_Trm_TPM_mt["ML14875a", ]`, mygroups_order = mygroups_order)



boxplot(ML14875a_MMP24_df_groups$expression ~ ML14875a_MMP24_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP24"
)

## ML42444a -- MMP28 - 1
ML42444a_MMP28_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML42444a",])
ML42444a_MMP28_df_groups <- data.frame( group = rownames(ML42444a_MMP28_df), expression =ML42444a_MMP28_df$`RSEMcounts_B2_Trm_TPM_mt["ML42444a", ]`, mygroups_order = mygroups_order)



boxplot(ML42444a_MMP28_df_groups$expression ~ ML42444a_MMP28_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP28-1"
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
        main = "MMP14"
)
##  ML064958a --  MMP28 - 2 
ML064958a_MMP28_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML064958a",])
ML064958a_MMP28_df_groups <- data.frame( group = rownames(ML064958a_MMP28_df), expression =ML064958a_MMP28_df$`RSEMcounts_B2_Trm_TPM_mt["ML064958a", ]`, mygroups_order = mygroups_order)



boxplot(ML064958a_MMP28_df_groups$expression ~ ML064958a_MMP28_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MMP28-2"
)


##Endopeptidase GO:0004175
#ML15195a - CTSE 
ML15195a_CTSE_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML15195a",])
ML15195a_CTSE_df_groups <-data.frame( group = rownames(ML15195a_CTSE_df), expression = ML15195a_CTSE_df$`RSEMcounts_B2_Trm_TPM_mt["ML15195a", ]`,mygroups_order = mygroups_order)

boxplot(ML15195a_CTSE_df_groups$expression ~ ML15195a_CTSE_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CTSE"
)



#ML154125a_CASP3
ML154125a_CASP3_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML154125a",])
ML154125a_CASP3_df_groups <-data.frame( group = rownames(ML154125a_CASP3_df), expression = ML154125a_CASP3_df$`RSEMcounts_B2_Trm_TPM_mt["ML154125a", ]`,mygroups_order = mygroups_order)

boxplot(ML154125a_CASP3_df_groups$expression ~ ML154125a_CASP3_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "CASP3"
)


#ML207912a_PCSK1
ML207912a_PCSK1_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML207912a",])
ML207912a_PCSK1_df_groups <-data.frame( group = rownames(ML207912a_PCSK1_df), expression = ML207912a_PCSK1_df$`RSEMcounts_B2_Trm_TPM_mt["ML207912a", ]`,mygroups_order = mygroups_order)

boxplot(ML207912a_PCSK1_df_groups$expression ~ ML207912a_PCSK1_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "PCSK1"
)

# WNT Pathway 
##ML11223a- TCF
ML11223a_TCF_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML11223a",])
ML11223a_TCF_df_groups <- data.frame( group = rownames(ML11223a_TCF_df), expression = ML11223a_TCF_df$`RSEMcounts_B2_Trm_TPM_mt["ML11223a", ]`, mygroups_order = mygroups_order)



boxplot(ML11223a_TCF_df_groups$expression ~ ML11223a_TCF_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TCF"
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
        main = "Beta-catenin"
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
        main = "WntX"
)

## MUCIN 
##ML18936a - VWF # NOISeq and EDGER
ML18936a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML18936a",])
ML18936a_df_groups <- data.frame( group = rownames(ML18936a_df), expression = ML18936a_df$`RSEMcounts_B2_Trm_TPM_mt["ML18936a", ]` , mygroups_order = mygroups_order)



boxplot( ML18936a_df_groups$expression ~ ML18936a_df_groups$mygroups_order,
         col= brewer.pal(8,"Set1"),
         xlab = "Hours after bisection (HAB)",
         ylab = "Transcripts per million (TPM)",
         range = 0,
         outliers = T,
         main = "VWF"
)

##ML019233a- MUC19 - consensus 
ML019233a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML019233a",])
ML019233a_df_groups <- data.frame( group = rownames(ML019233a_df), expression = ML019233a_df$`RSEMcounts_B2_Trm_TPM_mt["ML019233a", ]`, mygroups_order = mygroups_order)



boxplot(ML019233a_df_groups$expression ~ ML019233a_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MUC19"
)


## ML173711a - MUC5B #EdgeR & EBseq only
ML173711a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML173711a",])
ML173711a_df_groups <- data.frame( group = rownames(ML173711a_df), expression = ML173711a_df$`RSEMcounts_B2_Trm_TPM_mt["ML173711a", ]`  , mygroups_order = mygroups_order)


boxplot(ML173711a_df_groups$expression ~ ML173711a_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "MUC5B"
)

## SOX
##  ML047927a - SOX1 ## NOT INCLUDED -- not found in any of the DE lists 
ML047927a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML047927a",])
ML047927a_df_groups <- data.frame( group = rownames(ML047927a_df), expression = ML047927a_df$`RSEMcounts_B2_Trm_TPM_mt["ML047927a", ]`  , mygroups_order = mygroups_order)


boxplot(ML047927a_df_groups$expression ~ ML047927a_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "SOX1"
)




##  ML234028a - SOX2
ML234028a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML234028a",])
ML234028a_df_groups <- data.frame( group = rownames(ML234028a_df), expression = ML234028a_df$`RSEMcounts_B2_Trm_TPM_mt["ML234028a", ]`  , mygroups_order = mygroups_order)


boxplot(ML234028a_df_groups$expression ~ ML234028a_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "SOX2"
)


##  ML102235a - MlTGFbA
ML102235a_df <- as.data.frame(RSEMcounts_B2_Trm_TPM_mt["ML102235a",])
ML102235a_df_groups <- data.frame( group = rownames(ML102235a_df), expression = ML102235a_df$`RSEMcounts_B2_Trm_TPM_mt["ML102235a", ]`  , mygroups_order = mygroups_order)


boxplot(ML102235a_df_groups$expression ~ ML102235a_df_groups$mygroups_order,
        col= brewer.pal(8,"Set1"),
        xlab = "Hours after bisection (HAB)",
        ylab = "Transcripts per million (TPM)",
        range = 0,
        outliers = T,
        main = "TGFbA"
)





