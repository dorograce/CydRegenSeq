##R scripts for Ctenophore Cydippid regeneration RNA-seq differential gene expression analysis

NOTE: run scripts in the following order. Some data files needed by scripts are generated by earliers scripts.

#### RSEM_TPM_NOISeq_CPMfilte_ForOlap_V3.R

Pair-wise differential gene expression with NOISeq v2.38.0.  Must import transcript quantification files under folder RSEM.genes.out.FilesTransfer.  Generates "NOISeqbio_sig_genes_TPM_nodups_v2.38.0_2_V3.csv" file -- contains all DEG across all intervals and no gene duplicates.  Must import Blast annotations with "HumRef2020ML_eval.txt". Generates "NOIS_DE_log2FC_V3.csv" -- DEG with corresponding log2FC value in each interval with a BLAST hit 

#### RSEM_EC_EdgeR_forOlap_V3.R

Pair-wise differential gene expression with edgeR v3.36.0.  Must import transcript quantification files under folder RSEM.genes.out.FilesTransfer.  Generates "EdgeR_sig_genes_EC_nodups_V3.36.0_V3.csv" -- contains all DEG across all intervals and no gene duplicates.  Must import Blast annotations with "HumRef2020ML_eval.txt".  Generates "EdgeR_sig_pval_logFC_V3.tsv" -- DEG with corresponding logFC and p-value value in each interval with a BLAST hit.

#### RSEM_TPM_EBseqHMM_ForOlap_V3.R

Time-course differential gene expression with EBseqHMM v1.28.0.  Must import transcript quantification files under folder RSEM.genes.out.FilesTransfer.  Generates "EBseqHMM_sig_genes_TPM_V1.28.0_V3_.csv" file -- contains all DEG across entire time course.  Must import Blast annotations with "HumRef2020ML_eval.txt".  Generates "GeneDECalls_blastp_merge_PP0.5_split.txt" -- DEG in each path with a posterior probaility cutoff of 0.5 with corresponding best BLAST hit.

#### RSEM_NOISeq_EBseqHMM_edgeR_venn_Blastp_V3.R

Venn diagram of DEG found in EBseqHMM, edgeR, and NOISeq 
Must import DEG lists for NOISeq -- "NOISeqbio_sig_genes_TPM_nodups_v2.38.0_2_V3.csv", edgeR -- "EdgeR_sig_genes_EC_nodups_V3.36.0_V3.csv", EBseqhmm -- "EBseqHMM_sig_genes_TPM_V1.28.0_V3_.csv".  Must import Blast annotations with "HumRef2020ML_eval.txt". Generates consensus DEG list with best BLAST hit "Int_NOIS.TPM_ED.EC_EB.TPM_V3_blastp.csv".

#### RSEM_TPM_ED.EC_EB.TPM_NOIS.TPM_OLap_Heatmap_Cluster_TOPGO_BLASTP_V3

Heatmap of gene expression in the conensus DEG.  Cluster analysis from the generated heatmap.  Must import "Int_NOIS.TPM_ED.EC_EB.TPM_V3.csv".  Multi line graphs of Cluster #1-1.  Gene Ontolgy analysis of conensus DEG.  Must import "ML2.2.aa_goterms_2.txt" -- #Mnemiopsis genes and their GO annotations generated from interproscan in Babonis et al 2018.  Must import Blast annotations with "HumRef2020ML_eval.txt.  Generates "GO_FIGURE_TOP5_blastp.txt" -- the top 5 GO categories and the genes in each with corresponding best BLAST hit.  Must import Blast annotations with "HumRef2020ML_eval.txt.  Generates "Int_NOIS.TPM_ED.EC_EB.TPM_V3_hclust_blastp.xlsx" -- Genes in each cluster with corresponding best BLAST hit.

#### RSEM_TPM_ED.EC_EB.TPM_NOIS.TPM_OLap_Heatmap_Cluster_Editedline_V3.R

Heatmap of gene expression in the conensus DEG.  Cluster analysis from the generated heatmap.  Must import "Int_NOIS.TPM_ED.EC_EB.TPM_V3.csv".  Multi line graphs of Cluster #1 and Cluster #12 -- edited for outliers.

#### RSEM_TPM_BOXANDWHISKER_V3.R
Box and whisker plots showing temporal gene expression of individually selected genes. Must import transcript quantification files under folder RSEM.genes.out.FilesTransfer.

#### ED.EC_NOIS.TPM_DGEstats_Bargraph.R

mirror bargraphs showing # of DEG in each time interval for both pairwise methods (edgeR and NOISeq)
