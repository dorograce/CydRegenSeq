##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

#R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
#Copyright (C) 2023 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin20 (64-bit)


## EC = expected count 

## Using gene counts (EC) from RSEM to run differential gene expression (DGE) with EdgeR 


## Loading EdgeR ## - expects negative binomial
###Pairwise comparisons using Edge R ***Make sure to use unfiltered EC from RSEM  

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)
sessionInfo()
## version : 3.42.4 

# Loading gene quantification results # 
## Expected count is needed for edgeR

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

## preparing the data 

RSEMcounts_B2_Trm_EC <- data.frame(row.names = hab0_Rp1_B2_Trm$gene_id,
                                    uncut_1=unct_i_Rp1_B2_Trm$expected_count,
                                    uncut_2=unct_i_Rp2_B2_Trm$expected_count,
                                    uncut_3=unct_i_Rp3_B2_Trm$expected_count,
                                    hab00_1=hab0_Rp1_B2_Trm$expected_count,
                                    hab00_2=hab0_Rp2_B2_Trm$expected_count,
                                    hab00_3=hab0_Rp3_B2_Trm$expected_count,
                                    hab01_1=hab1_Rp1_B2_Trm$expected_count,
                                    hab01_2=hab1_Rp2_B2_Trm$expected_count,
                                    hab01_3=hab1_Rp3_B2_Trm$expected_count,
                                    hab03_1=hab3_Rp1_B2_Trm$expected_count,
                                    hab03_2=hab3_Rp2_B2_Trm$expected_count,
                                    hab03_3=hab3_Rp3_B2_Trm$expected_count,
                                    hab06_1=hab6_Rp1_B2_Trm$expected_count,
                                    hab06_2=hab6_Rp2_B2_Trm$expected_count,
                                    hab06_3=hab6_Rp3_B2_Trm$expected_count,
                                    hab12_1=hab12_Rp1_B2_Trm$expected_count,
                                    hab12_2=hab12_Rp2_B2_Trm$expected_count,
                                    hab12_3=hab12_Rp3_B2_Trm$expected_count,
                                    hab24_1=hab24_i_Rp1_B2_Trm$expected_count,
                                    hab24_2=hab24_i_Rp2_B2_Trm$expected_count,
                                    hab24_3=hab24_i_Rp3_B2_Trm$expected_count,
                                    hab48_1=hab48_i_Rp1_B2_Trm$expected_count,
                                    hab48_2=hab48_i_Rp2_B2_Trm$expected_count,
                                    hab48_3=hab48_i_Rp3_B2_Trm$expected_count)

length(RSEMcounts_B2_Trm_EC$uncut_1) #16548 


## Running edgeR
## Make sure that the columns are labeled such that the ladder time point is positioned first i.e. hab00-000 is listed at the top of myedgeR_decide_U00 and the same for hab01-hab00 in  myedgeR_decide_00_01 and so on

#comparing uncut with 0 hrs
RSEMcounts_B2_Trm_U00<- subset(RSEMcounts_B2_Trm_EC, select = c(uncut_1,uncut_2,uncut_3,hab00_1,hab00_2,hab00_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compareU00 = colnames(RSEMcounts_B2_Trm_U00)[c(grep("uncut", colnames(RSEMcounts_B2_Trm_U00)),grep("hab00",colnames(RSEMcounts_B2_Trm_U00)))]
mygroupsU00 = rep(c("000","000","000","hab00","hab00","hab00"))
myedgeR_U00 = DGEList(counts=RSEMcounts_B2_Trm_U00[,cond2compareU00], group=mygroupsU00)
myedgeR_U00 = calcNormFactors(myedgeR_U00)  
myedgeR_U00 = estimateCommonDisp(myedgeR_U00)
myedgeR_U00 = estimateTagwiseDisp(myedgeR_U00, trend = "movingave")
myedgeR_U00 = exactTest(myedgeR_U00)
pvaluesedgeR_U00 = as.data.frame(myedgeR_U00$table$PValue)
logFCedgeR_U00 = as.data.frame(myedgeR_U00$table$logFC)
myedgeR_decide_U00 = decideTestsDGE(myedgeR_U00, adjust.method = "BH", p.value = 0.05)
myedgeR_decide_U00_df <- as.data.frame(myedgeR_decide_U00)
myedgeR_decide_U00_pval_logFC <- cbind(myedgeR_decide_U00_df,pvaluesedgeR_U00,logFCedgeR_U00 )

## Getting a list of only significant genes between Uncut and 00hab
EdgeR_sig_U_00 <- myedgeR_decide_U00[myedgeR_decide_U00[, "hab00-000"] != 0,] #which genes do not equal 0, therefore are significantly down OR upregulated
EdgeR_names_sig_U_00 <- as.data.frame(EdgeR_sig_U_00)
length(rownames(EdgeR_names_sig_U_00))# 30 v3.42.4 
EdgeR_sig_U_00_pval_logFC <- myedgeR_decide_U00_pval_logFC[myedgeR_decide_U00_pval_logFC[, "hab00-000"] != 0,]

########################
# Comparing 0 with 1  
RSEMcounts_B2_Trm_0001 <- subset(RSEMcounts_B2_Trm_EC, select = c(hab00_1,hab00_2,hab00_3,hab01_1,hab01_2,hab01_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare_00_01 = colnames(RSEMcounts_B2_Trm_0001)[c(grep("hab00", colnames(RSEMcounts_B2_Trm_0001)),grep("hab01",colnames(RSEMcounts_B2_Trm_0001)))]
mygroups_00_01 = rep(c("hab00","hab00","hab00","hab01","hab01","hab01"))
myedgeR_00_01 <- DGEList(counts=RSEMcounts_B2_Trm_0001[,cond2compare_00_01], group=mygroups_00_01)
myedgeR_00_01 = calcNormFactors(myedgeR_00_01)  
myedgeR_00_01 = estimateCommonDisp(myedgeR_00_01)
myedgeR_00_01 = estimateTagwiseDisp(myedgeR_00_01, trend = "movingave")
myedgeR_00_01 = exactTest(myedgeR_00_01)
pvaluesedgeR_0001 = as.data.frame(myedgeR_00_01$table$PValue)
logFCedgeR_0001 = as.data.frame(myedgeR_00_01$table$logFC)
myedgeR_decide_00_01 = decideTestsDGE(myedgeR_00_01, adjust.method = "BH", p.value = 0.05)
myedgeR_decide_00_01_df <- as.data.frame(myedgeR_decide_00_01)
myedgeR_decide_00_01_pval_logFC <- cbind(myedgeR_decide_00_01_df,pvaluesedgeR_0001, logFCedgeR_0001)



## Getting a list of only significant genes between 00hab and 01hab
EdgeR_sig_00_01 <- myedgeR_decide_00_01[myedgeR_decide_00_01[, "hab01-hab00"] != 0,]
EdgeR_names_sig_00_01 <- as.data.frame(EdgeR_sig_00_01)
EdgeR_sig_00_01_pval_logFC <- myedgeR_decide_00_01_pval_logFC[myedgeR_decide_00_01_pval_logFC[, "hab01-hab00"] != 0,]
EdgeR_names_sig_00_01_pval_logFC <- as.data.frame(EdgeR_sig_00_01_pval_logFC)
length(EdgeR_names_sig_00_01_pval_logFC$`hab01-hab00`) #28 v3.42.4 

#comparing 1 with 3 hrs
RSEMcounts_B2_Trm_0103<- subset(RSEMcounts_B2_Trm_EC, select = c(hab01_1,hab01_2,hab01_3,hab03_1,hab03_2,hab03_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare01_03 = colnames(RSEMcounts_B2_Trm_0103)[c(grep("hab01", colnames(RSEMcounts_B2_Trm_0103)),grep("hab03",colnames(RSEMcounts_B2_Trm_0103)))]
cond2compare01_03
mygroups_01_03 = rep(c("hab01","hab01","hab01","hab03","hab03","hab03"))
myedgeR_01_03 <- DGEList(counts=RSEMcounts_B2_Trm_0103[,cond2compare01_03], group=mygroups_01_03)
myedgeR_01_03 = calcNormFactors(myedgeR_01_03)  
myedgeR_01_03 = estimateCommonDisp(myedgeR_01_03)
myedgeR_01_03 = estimateTagwiseDisp(myedgeR_01_03, trend = "movingave")
myedgeR_01_03 = exactTest(myedgeR_01_03)
pvaluesedgeR_0103 = as.data.frame(myedgeR_01_03$table$PValue)
logFCedgeR_0103 = as.data.frame(myedgeR_01_03$table$logFC)
myedgeR_decide_01_03 = decideTestsDGE(myedgeR_01_03,adjust.method= "BH",p.value = 0.05)

myedgeR_decide_01_03_df <- as.data.frame(myedgeR_decide_01_03)
myedgeR_decide_01_03_pval_logFC <- cbind(myedgeR_decide_01_03_df ,pvaluesedgeR_0103, logFCedgeR_0103 )


## Getting a list of only significant genes between 01hab and 03hab

EdgeR_sig_01_03 <- myedgeR_decide_01_03[myedgeR_decide_01_03[, "hab03-hab01"] != 0,]
EdgeR_names_sig_01_03 <- as.data.frame(EdgeR_sig_01_03)
EdgeR_sig_01_03_pval_logFC <- myedgeR_decide_01_03_pval_logFC[myedgeR_decide_01_03_pval_logFC[, "hab03-hab01"] != 0,]
length(EdgeR_sig_01_03_pval_logFC$`hab03-hab01`) # 64  v3.42.4 


#comparing 3 with 6 hrs
RSEMcounts_B2_Trm_0306 <- subset(RSEMcounts_B2_Trm_EC, select = c(hab03_1,hab03_2,hab03_3,hab06_1,hab06_2,hab06_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare03_06 = colnames(RSEMcounts_B2_Trm_0306)[c(grep("hab03", colnames(RSEMcounts_B2_Trm_0306 )),grep("hab06",colnames(RSEMcounts_B2_Trm_0306 )))]
cond2compare03_06
mygroups_03_06 = rep(c("hab03","hab03","hab03","hab06","hab06","hab06"))

myedgeR_03_06 <- DGEList(counts=RSEMcounts_B2_Trm_0306[,cond2compare03_06], group=mygroups_03_06)
myedgeR_03_06 = calcNormFactors(myedgeR_03_06)  
myedgeR_03_06 = estimateCommonDisp(myedgeR_03_06)
myedgeR_03_06 = estimateTagwiseDisp(myedgeR_03_06, trend = "movingave")
myedgeR_03_06 = exactTest(myedgeR_03_06)
pvaluesedgeR_0306 = as.data.frame(myedgeR_03_06$table$PValue)
logFCedgeR_0306 = as.data.frame(myedgeR_03_06$table$logFC)

myedgeR_decide_03_06 = decideTestsDGE(myedgeR_03_06, adjust.method = "BH", p.value = 0.05)
myedgeR_decide_03_06_df <- as.data.frame(myedgeR_decide_03_06)
myedgeR_decide_03_06_pval_logFC <- cbind(myedgeR_decide_03_06_df,pvaluesedgeR_0306,logFCedgeR_0306)

## Getting a list of only significant genes between 03hab and 06hab

EdgeR_sig_03_06 <- myedgeR_decide_03_06[myedgeR_decide_03_06[, "hab06-hab03"] != 0,]
EdgeR_names_sig_03_06 <- as.data.frame(EdgeR_sig_03_06)
EdgeR_sig_03_06_pval_logFC <- myedgeR_decide_03_06_pval_logFC[myedgeR_decide_03_06_pval_logFC[, "hab06-hab03"] != 0,]
length(EdgeR_sig_03_06_pval_logFC$`hab06-hab03`) #244 v3.42.4 


#comparing 6 with 12 hrs
RSEMcounts_B2_Trm_0612 <- subset(RSEMcounts_B2_Trm_EC, select = c(hab06_1,hab06_2,hab06_3 ,hab12_1,hab12_2,hab12_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare06_12 = colnames(RSEMcounts_B2_Trm_0612)[c(grep("hab06", colnames(RSEMcounts_B2_Trm_0612)),grep("hab12",colnames(RSEMcounts_B2_Trm_0612)))]
cond2compare06_12
mygroups_06_12 = rep(c("hab06","hab06","hab06","hab12","hab12","hab12"))
myedgeR_06_12 <- DGEList(counts=RSEMcounts_B2_Trm_0612[,cond2compare06_12], group=mygroups_06_12)
myedgeR_06_12 = calcNormFactors(myedgeR_06_12)  
myedgeR_06_12 = estimateCommonDisp(myedgeR_06_12)
myedgeR_06_12 = estimateTagwiseDisp(myedgeR_06_12, trend = "movingave")
myedgeR_06_12 = exactTest(myedgeR_06_12)
pvaluesedgeR_0612 = as.data.frame(myedgeR_06_12$table$PValue)
logFCedgeR_0612 = as.data.frame(myedgeR_06_12$table$logFC)
myedgeR_decide_06_12 = decideTestsDGE(myedgeR_06_12, adjust.method = "BH", p.value = 0.05)
myedgeR_decide_06_12_df <- as.data.frame(myedgeR_decide_06_12)
myedgeR_decide_06_12_pval_logFC <- cbind(myedgeR_decide_06_12_df,pvaluesedgeR_0612,logFCedgeR_0612 )


## Getting a list of only significant genes between 06hab and 12hab

EdgeR_sig_06_12 <- myedgeR_decide_06_12[myedgeR_decide_06_12[, "hab12-hab06"] != 0,]
EdgeR_names_sig_06_12 <- as.data.frame(EdgeR_sig_06_12)
EdgeR_sig_06_12_pval_logFC <- myedgeR_decide_06_12_pval_logFC[myedgeR_decide_06_12_pval_logFC[, "hab12-hab06"] != 0,]
length(EdgeR_sig_06_12_pval_logFC$`hab12-hab06`) #8 v3.42.4



#comparing 12 with 24 hrs
RSEMcounts_B2_Trm_1224<- subset(RSEMcounts_B2_Trm_EC, select = c(hab12_1,hab12_2,hab12_3,hab24_1,hab24_2,hab24_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare12_24 = colnames(RSEMcounts_B2_Trm_1224)[c(grep("hab12", colnames(RSEMcounts_B2_Trm_1224)),grep("hab24",colnames(RSEMcounts_B2_Trm_1224)))]
cond2compare12_24
mygroups_12_24 = rep(c("hab12","hab12","hab12","hab24","hab24","hab24"))
myedgeR_12_24 <- DGEList(counts=RSEMcounts_B2_Trm_1224[,cond2compare12_24], group=mygroups_12_24)
myedgeR_12_24 = calcNormFactors(myedgeR_12_24)  
myedgeR_12_24 = estimateCommonDisp(myedgeR_12_24)
myedgeR_12_24 = estimateTagwiseDisp(myedgeR_12_24, trend = "movingave")
myedgeR_12_24 = exactTest(myedgeR_12_24)
myedgeR_decide_12_24 = decideTestsDGE(myedgeR_12_24, adjust.method = "BH", p.value = 0.05)
myedgeR_decide_12_24_df <- as.data.frame(myedgeR_decide_12_24)


## Getting a list of only significant genes between 12hab and 24hab

EdgeR_sig_12_24 <- myedgeR_decide_12_24[myedgeR_decide_12_24[, "hab24-hab12"] != 0,]
EdgeR_names_sig_12_24 <- as.data.frame(EdgeR_sig_12_24)
#0 sig genes!! v3.42.4


#comparing 24 with 48 hrs
RSEMcounts_B2_Trm_2448 <- subset(RSEMcounts_B2_Trm_EC, select = c(hab24_1,hab24_2,hab24_3,hab48_1,hab48_2,hab48_3))#only columns of interest from EDGER filtered data ( non normalized) but no technical replicates
cond2compare24_48 = colnames(RSEMcounts_B2_Trm_2448)[c(grep("hab24", colnames(RSEMcounts_B2_Trm_2448)),grep("hab48",colnames(RSEMcounts_B2_Trm_2448)))]
cond2compare24_48
mygroups_24_48 = rep(c("hab24","hab24","hab24","hab48","hab48","hab48"))
myedgeR_24_48 <- DGEList(counts=RSEMcounts_B2_Trm_2448[,cond2compare24_48], group=mygroups_24_48)
myedgeR_24_48 = calcNormFactors(myedgeR_24_48)  
myedgeR_24_48 = estimateCommonDisp(myedgeR_24_48)
myedgeR_24_48 = estimateTagwiseDisp(myedgeR_24_48, trend = "movingave")
myedgeR_24_48 = exactTest(myedgeR_24_48)
pvaluesedgeR_2448 = as.data.frame(myedgeR_24_48$table$PValue)
logFCedgeR_2448 = as.data.frame(myedgeR_24_48$table$logFC)

myedgeR_decide_24_48 = decideTestsDGE(myedgeR_24_48, adjust.method = "BH", p.value = 0.05)

myedgeR_decide_24_48_df <- as.data.frame(myedgeR_decide_24_48)
myedgeR_decide_24_48_pval_logFC <- cbind(myedgeR_decide_24_48_df,pvaluesedgeR_2448,logFCedgeR_2448 )

## Getting a list of only significant genes between 24hab and 48hab

EdgeR_sig_24_48 <- myedgeR_decide_24_48[myedgeR_decide_24_48[, "hab48-hab24"] != 0,]
EdgeR_names_sig_24_48 <- as.data.frame(EdgeR_sig_24_48)
EdgeR_sig_24_48_pval_logFC <- myedgeR_decide_24_48_pval_logFC[myedgeR_decide_24_48_pval_logFC[, "hab48-hab24"] != 0,]
length(EdgeR_sig_24_48_pval_logFC$`hab48-hab24`) #9 v3.42.4



EdgeR_sig_genes_EC_ALL <- list(rownames(EdgeR_names_sig_24_48),
                               rownames(EdgeR_names_sig_U_00),
                               rownames(EdgeR_names_sig_06_12),
                               rownames(EdgeR_names_sig_00_01),
                               rownames(EdgeR_names_sig_03_06),
                               rownames(EdgeR_names_sig_01_03))

## The list of significantly DE across all pairwise time intervals contains replicates, as some genes are DE in more than one time interval
EdgeR_sig_genes_EC_ALL <- lapply(EdgeR_sig_genes_EC_ALL, function(x) as.matrix(x))
EdgeR_sig_genes_EC_ALL <- do.call(rbind,EdgeR_sig_genes_EC_ALL)
length(EdgeR_sig_genes_EC_ALL) #383 v3.42.4
EdgeR_sig_genes_EC_nodups <- EdgeR_sig_genes_EC_ALL[!duplicated(EdgeR_sig_genes_EC_ALL)]
length(EdgeR_sig_genes_EC_nodups )#348 V 3.42.4 
write.csv(EdgeR_sig_genes_EC_nodups,row.names = F,file ="EdgeR_sig_v4.csv")


#Separate up from down regulation in each pairwise grouping & BLAST

## For BLAST
#ML gene models and their protein blast hit 
HumRef2020ML <- read.table("HumRef2020ML_eval.txt",header = F)
HumRef2020ML_df <- as.data.frame(HumRef2020ML, row.names = HumRef2020ML$V1)

# For consensus gene set annotation 
con_genenames <- read.csv(file="consensus_DE_v4.csv",1) # ML gene names of the genes identified as DE in all three methods
con_genenames_df <- data.frame(row.names= con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111...,V1 = con_genenames$VNOIS_ED_EB_Set.IntersectionSets...111... )


#uncut - 10 mins
#down
EdgeR_sig_U_00_pval_logFC
EdgeR_sig_U_00_down_pval_logFC <- EdgeR_sig_U_00_pval_logFC[EdgeR_sig_U_00_pval_logFC[, "hab00-000"] == -1,]
length(EdgeR_sig_U_00_down_pval_logFC$`uncut-hab00`) #0

#up
EdgeR_sig_U_00_up_pval_logFC <- EdgeR_sig_U_00_pval_logFC[EdgeR_sig_U_00_pval_logFC[, "hab00-000"] == 1,]
length(EdgeR_sig_U_00_up_pval_logFC$`hab00-000`) #30
#BLAST 
EdgeR_sig_U_00_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_U_00_up_pval_logFC,by='row.names')
#order P-values 
EdgeR_sig_U_00_up_pval_logFC_blastp_merge_order<- EdgeR_sig_U_00_up_pval_logFC_blastp_merge[order(EdgeR_sig_U_00_up_pval_logFC_blastp_merge$`myedgeR_U00$table$PValue`, decreasing = F),]
length(EdgeR_sig_U_00_up_pval_logFC_blastp_merge_order$Row.names) #16

#Consensus UP order
EdgeR_sig_U_00_up_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_U_00_up_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_U_00_up_pval_logFC_blastp_con_merge_order<- EdgeR_sig_U_00_up_pval_logFC_blastp_con_merge [order(EdgeR_sig_U_00_up_pval_logFC_blastp_con_merge$`myedgeR_U00$table$PValue`, decreasing = F),]


# 10 mins - 1 hour
#down
EdgeR_sig_00_01_pval_logFC
EdgeR_sig_00_01_down_pval_logFC <- EdgeR_sig_00_01_pval_logFC[EdgeR_sig_00_01_pval_logFC[, "hab01-hab00"] == -1,]
length(EdgeR_sig_00_01_down_pval_logFC$`hab01-hab00`) #3
#BLAST
EdgeR_sig_00_01_down_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_00_01_down_pval_logFC,by='row.names')
#order P-values 
EdgeR_sig_00_01_down_pval_logFC_blastp_merge_order<- EdgeR_sig_00_01_down_pval_logFC_blastp_merge[order(EdgeR_sig_00_01_down_pval_logFC_blastp_merge$`myedgeR_00_01$table$PValue`, decreasing = F),]
length(EdgeR_sig_00_01_down_pval_logFC_blastp_merge_order$Row.names)#2

#Consensus DOWN order
EdgeR_sig_00_01_down_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_00_01_down_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_00_01_down_pval_logFC_blastp_con_merge_order <- EdgeR_sig_00_01_down_pval_logFC_blastp_con_merge[order(EdgeR_sig_00_01_down_pval_logFC_blastp_con_merge$`myedgeR_00_01$table$PValue`, decreasing = F),]


#up
EdgeR_sig_00_01_up_pval_logFC <- EdgeR_sig_00_01_pval_logFC[EdgeR_sig_00_01_pval_logFC[, "hab01-hab00"] == 1,]
length(EdgeR_sig_00_01_up_pval_logFC$`hab01-hab00`) #25

#BLAST
EdgeR_sig_00_01_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_00_01_up_pval_logFC,by='row.names')
#order P-values 
EdgeR_sig_00_01_up_pval_logFC_blastp_merge_order<- EdgeR_sig_00_01_up_pval_logFC_blastp_merge[order(EdgeR_sig_00_01_up_pval_logFC_blastp_merge$`myedgeR_00_01$table$PValue`, decreasing = F),]
length(EdgeR_sig_00_01_up_pval_logFC_blastp_merge_order$Row.names)#17

#Consensus UP order
EdgeR_sig_00_01_up_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_00_01_up_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_00_01_up_pval_logFC_blastp_con_merge_order <- EdgeR_sig_00_01_up_pval_logFC_blastp_con_merge[order(EdgeR_sig_00_01_up_pval_logFC_blastp_con_merge$`myedgeR_00_01$table$PValue`, decreasing = F),]


# 1 hour - 3 hours
#down
EdgeR_sig_01_03_pval_logFC
EdgeR_sig_01_03_down_pval_logFC <- EdgeR_sig_01_03_pval_logFC[EdgeR_sig_01_03_pval_logFC[, "hab03-hab01"] == -1,]
length(EdgeR_sig_01_03_down_pval_logFC$`hab03-hab01`) #39
#BLAST
EdgeR_sig_01_03_down_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_01_03_down_pval_logFC,by='row.names') #25
#order P-values 
EdgeR_sig_01_03_down_pval_logFC_blastp_merge_order<- EdgeR_sig_01_03_down_pval_logFC_blastp_merge[order(EdgeR_sig_01_03_down_pval_logFC_blastp_merge$`myedgeR_01_03$table$PValue`, decreasing = F),]
length(EdgeR_sig_01_03_down_pval_logFC_blastp_merge_order$Row.names)#25

#Consensus DOWN order
EdgeR_sig_01_03_down_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_01_03_down_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_01_03_down_pval_logFC_blastp_con_merge_order <- EdgeR_sig_01_03_down_pval_logFC_blastp_con_merge[order(EdgeR_sig_01_03_down_pval_logFC_blastp_con_merge$`myedgeR_01_03$table$PValue`, decreasing = F),]


#up 
EdgeR_sig_01_03_up_pval_logFC <- EdgeR_sig_01_03_pval_logFC[EdgeR_sig_01_03_pval_logFC[, "hab03-hab01"] == 1,]
length(EdgeR_sig_01_03_up_pval_logFC$`hab03-hab01`) #25
#BLAST
EdgeR_sig_01_03_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_01_03_up_pval_logFC,by='row.names')
#order P-values 
EdgeR_sig_01_03_up_pval_logFC_blastp_merge_order<- EdgeR_sig_01_03_up_pval_logFC_blastp_merge[order(EdgeR_sig_01_03_up_pval_logFC_blastp_merge$`myedgeR_01_03$table$PValue`, decreasing = F),]
length(EdgeR_sig_01_03_up_pval_logFC_blastp_merge_order$Row.names)#14

#Consensus up order
EdgeR_sig_01_03_up_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_01_03_up_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_01_03_up_pval_logFC_blastp_con_merge_order <- EdgeR_sig_01_03_up_pval_logFC_blastp_con_merge[order(EdgeR_sig_01_03_up_pval_logFC_blastp_con_merge$`myedgeR_01_03$table$PValue`, decreasing = F),]




# 3 hour - 6 hours

#down
EdgeR_sig_03_06_pval_logFC
EdgeR_sig_03_06_down_pval_logFC <- EdgeR_sig_03_06_pval_logFC[EdgeR_sig_03_06_pval_logFC[, "hab06-hab03"] == -1,]
length(EdgeR_sig_03_06_down_pval_logFC$`hab06-hab03`) #115
#BLAST
EdgeR_sig_03_06_down_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_03_06_down_pval_logFC,by='row.names') 
#order P-values 
EdgeR_sig_03_06_down_pval_logFC_blastp_merge_order<- EdgeR_sig_03_06_down_pval_logFC_blastp_merge[order(EdgeR_sig_03_06_down_pval_logFC_blastp_merge$`myedgeR_03_06$table$PValue`, decreasing = F),]
length(EdgeR_sig_03_06_down_pval_logFC_blastp_merge_order$Row.names)#70

#Consensus down order
EdgeR_sig_03_06_down_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_03_06_down_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_03_06_down_pval_logFC_blastp_con_merge_order <- EdgeR_sig_03_06_down_pval_logFC_blastp_con_merge[order(EdgeR_sig_03_06_down_pval_logFC_blastp_con_merge$`myedgeR_03_06$table$PValue`, decreasing = F),]


#up
EdgeR_sig_03_06_up_pval_logFC <- EdgeR_sig_03_06_pval_logFC[EdgeR_sig_03_06_pval_logFC[, "hab06-hab03"] == 1,]
length(EdgeR_sig_03_06_up_pval_logFC$`hab06-hab03`) #129 
#BLAST 
EdgeR_sig_03_06_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_03_06_up_pval_logFC,by='row.names') 
#order P-values 
EdgeR_sig_03_06_up_pval_logFC_blastp_merge_order<- EdgeR_sig_03_06_up_pval_logFC_blastp_merge[order(EdgeR_sig_03_06_up_pval_logFC_blastp_merge$`myedgeR_03_06$table$PValue`, decreasing = F),]
length(EdgeR_sig_03_06_up_pval_logFC_blastp_merge_order$Row.names) #88


#Consensus up order
EdgeR_sig_03_06_up_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_03_06_up_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_03_06_up_pval_logFC_blastp_con_merge_order <- EdgeR_sig_03_06_up_pval_logFC_blastp_con_merge[order(EdgeR_sig_03_06_up_pval_logFC_blastp_con_merge$`myedgeR_03_06$table$PValue`, decreasing = F),]



# 6 hour - 12 hours
#down
EdgeR_sig_06_12_pval_logFC
EdgeR_sig_06_12_down_pval_logFC <- EdgeR_sig_06_12_pval_logFC[EdgeR_sig_06_12_pval_logFC[, "hab12-hab06"] == -1,]
length(EdgeR_sig_06_12_down_pval_logFC$`hab12-hab06`) #3

# BLAST 
EdgeR_sig_06_12_down_pval_logFC_blastp_merge<- merge(HumRef2020ML_df ,EdgeR_sig_06_12_down_pval_logFC,by='row.names')
length(EdgeR_sig_06_12_down_pval_logFC_blastp_merge$Row.names) #2


#Consensus down order
EdgeR_sig_06_12_down_pval_logFC_blastp_con_merge <- merge(EdgeR_sig_06_12_down_pval_logFC_blastp_merge,con_genenames_df,by='V1')
EdgeR_sig_06_12_down_pval_logFC_blastp_con_merge_order <- EdgeR_sig_06_12_down_pval_logFC_blastp_con_merge[order(EdgeR_sig_06_12_down_pval_logFC_blastp_con_merge$`myedgeR_06_12$table$PValue`, decreasing = F),]



#up
EdgeR_sig_06_12_up_pval_logFC <- EdgeR_sig_06_12_pval_logFC[EdgeR_sig_06_12_pval_logFC[, "hab12-hab06"] == 1,]
length(EdgeR_sig_06_12_up_pval_logFC$`hab12-hab06`) #5
#BLAST
EdgeR_sig_06_12_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_06_12_up_pval_logFC,by='row.names') 

#order P-values 
EdgeR_sig_06_12_up_pval_logFC_blastp_merge_order<- EdgeR_sig_06_12_up_pval_logFC_blastp_merge[order(EdgeR_sig_06_12_up_pval_logFC_blastp_merge$`myedgeR_06_12$table$PValue`, decreasing = F),]
length(EdgeR_sig_06_12_up_pval_logFC_blastp_merge_order$Row.names)#3

#Consensus up order
EdgeR_sig_06_12_up_pval_logFC_blastp_con_merge  <- merge(EdgeR_sig_06_12_up_pval_logFC_blastp_merge ,con_genenames_df,by='V1')
EdgeR_sig_06_12_up_pval_logFC_blastp_con_merge_order <- EdgeR_sig_06_12_up_pval_logFC_blastp_con_merge[order(EdgeR_sig_06_12_up_pval_logFC_blastp_con_merge$`myedgeR_06_12$table$PValue`, decreasing = F),]





# 12 hour - 24 hour -- no sig changes! 

# 24 hour - 48 hour
#down
EdgeR_sig_24_48_pval_logFC
EdgeR_sig_24_48_down_pval_logFC <- EdgeR_sig_24_48_pval_logFC[EdgeR_sig_24_48_pval_logFC[, "hab48-hab24"] == -1,]
length(EdgeR_sig_24_48_down_pval_logFC$`hab48-hab24`) #2
#BLAST
EdgeR_sig_24_48_down_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_24_48_down_pval_logFC,by='row.names')#1

#Consensus down order
EdgeR_sig_24_48_down_pval_logFC_blastp_con_merge  <- merge(EdgeR_sig_24_48_down_pval_logFC_blastp_merge ,con_genenames_df,by='V1')
EdgeR_sig_24_48_down_pval_logFC_blastp_con_merge_order <- EdgeR_sig_24_48_down_pval_logFC_blastp_con_merge[order(EdgeR_sig_24_48_down_pval_logFC_blastp_con_merge$`myedgeR_24_48$table$PValue`, decreasing = F),]



#up
EdgeR_sig_24_48_up_pval_logFC <- EdgeR_sig_24_48_pval_logFC[EdgeR_sig_24_48_pval_logFC[, "hab48-hab24"] == 1,]
length(EdgeR_sig_24_48_up_pval_logFC$`hab48-hab24`) #7
#BLAST
EdgeR_sig_24_48_up_pval_logFC_blastp_merge <- merge(HumRef2020ML_df ,EdgeR_sig_24_48_up_pval_logFC,by='row.names')
length(EdgeR_sig_24_48_up_pval_logFC_blastp_merge$Row.names)#3

#order P-values 
EdgeR_sig_24_48_up_pval_logFC_blastp_merge_order<- EdgeR_sig_24_48_up_pval_logFC_blastp_merge[order(EdgeR_sig_24_48_up_pval_logFC_blastp_merge$`myedgeR_24_48$table$PValue`, decreasing = F),]

#Consensus up order
EdgeR_sig_24_48_up_pval_logFC_blastp_con_merge  <- merge(EdgeR_sig_24_48_up_pval_logFC_blastp_merge ,con_genenames_df,by='V1')
EdgeR_sig_24_48_up_pval_logFC_blastp_con_merge_order <- EdgeR_sig_24_48_up_pval_logFC_blastp_con_merge[order(EdgeR_sig_24_48_up_pval_logFC_blastp_con_merge$`myedgeR_24_48$table$PValue`, decreasing = F),]





#Merge all of the sig genes in each interval with their corresponding p-values and LogFC into a list 
EdgeR_sig_pval_logFC <- list(sig_U_00_up = EdgeR_sig_U_00_up_pval_logFC_blastp_merge_order,
                             sig_00_01_down= EdgeR_sig_00_01_down_pval_logFC_blastp_merge_order,
                             sig_00_01_up= EdgeR_sig_00_01_up_pval_logFC_blastp_merge_order,
                             sig_01_03_down = EdgeR_sig_01_03_down_pval_logFC_blastp_merge_order,
                             sig_01_03_up = EdgeR_sig_01_03_up_pval_logFC_blastp_merge_order,
                             sig_03_06_down = EdgeR_sig_03_06_down_pval_logFC_blastp_merge_order,
                             sig_03_06_up = EdgeR_sig_03_06_up_pval_logFC_blastp_merge_order,
                             sig_06_12_down = EdgeR_sig_06_12_down_pval_logFC_blastp_merge,
                             sig_06_12_up = EdgeR_sig_06_12_up_pval_logFC_blastp_merge_order,
                             sig_24_48_down = EdgeR_sig_24_48_down_pval_logFC_blastp_merge,
                             sig_24_48_up = EdgeR_sig_24_48_up_pval_logFC_blastp_merge_order)
options(max.print=999999)
capture.output(EdgeR_sig_pval_logFC, file = "EdgeR_sig_pval_logFC_v4.csv")


#Merge all of the sig genes in each interval that are also found in the consensus with their corresponding p-values and LogFC into a list 
EdgeR_sig_pval_logFC_consensus <- list(sig_U_00_up = EdgeR_sig_U_00_up_pval_logFC_blastp_con_merge_order,
                             sig_00_01_down= EdgeR_sig_00_01_down_pval_logFC_blastp_con_merge_order,
                             sig_00_01_up= EdgeR_sig_00_01_up_pval_logFC_blastp_con_merge_order,
                             sig_01_03_down = EdgeR_sig_01_03_down_pval_logFC_blastp_con_merge_order,
                             sig_01_03_up = EdgeR_sig_01_03_up_pval_logFC_blastp_con_merge_order,
                             sig_03_06_down = EdgeR_sig_03_06_down_pval_logFC_blastp_con_merge_order,
                             sig_03_06_up = EdgeR_sig_03_06_up_pval_logFC_blastp_con_merge_order,
                             sig_06_12_down = EdgeR_sig_06_12_down_pval_logFC_blastp_con_merge_order,
                             sig_06_12_up = EdgeR_sig_06_12_up_pval_logFC_blastp_con_merge_order,
                             sig_24_48_down = EdgeR_sig_24_48_down_pval_logFC_blastp_con_merge_order,
                             sig_24_48_up = EdgeR_sig_24_48_up_pval_logFC_blastp_con_merge_order)
options(max.print=999999)
capture.output(EdgeR_sig_pval_logFC_consensus, file = "EdgeR_sig_pval_logFC_consensus_v4.csv")





# Getting a dataframe of # of DEG(with or without BLAST) for each time interval 
EdgeR_sig_numbers <- data.frame(sig_U_00_up = length(EdgeR_sig_U_00_up_pval_logFC$`hab00-000`),
                                sig_U_00_down = length(EdgeR_sig_U_00_down_pval_logFC$`uncut-hab00`),
                                sig_00_01_up= length(EdgeR_sig_00_01_up_pval_logFC$`hab01-hab00`),
                                sig_00_01_down= length(EdgeR_sig_00_01_down_pval_logFC$`hab01-hab00`),
                                sig_01_03_up = length(EdgeR_sig_01_03_up_pval_logFC$`hab03-hab01`),
                                sig_01_03_down = length(EdgeR_sig_01_03_down_pval_logFC$`hab03-hab01`),
                                sig_03_06_up = length(EdgeR_sig_03_06_up_pval_logFC$`hab06-hab03`),
                                sig_03_06_down = length(EdgeR_sig_03_06_down_pval_logFC$`hab06-hab03`),
                                sig_06_12_up = length(EdgeR_sig_06_12_up_pval_logFC$`hab12-hab06`),
                                sig_06_12_down = length(EdgeR_sig_06_12_down_pval_logFC$`hab12-hab06`),
                                sig_24_48_up = length(EdgeR_sig_24_48_up_pval_logFC$`hab48-hab24`),
                                sig_24_48_down = length(EdgeR_sig_24_48_down_pval_logFC$`hab48-hab24`)
                                 )

