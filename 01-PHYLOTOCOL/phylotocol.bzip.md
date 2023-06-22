# <TITLE>
 Principal Investigator: Dorothy Mitchell, Allison Edgar, Joseph Ryan, Mark Martindale
 Draft or Version Number: v.0.01
 Date: 6-12-2021

## SUMMARY OF CHANGES FROM PREVIOUS VERSION
\<affected section> \<summary of revisions made> \<rationale>

## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _OBJECTIVE_  

Goal is to classify Mnemiopsis bzip domain containing proteins

## 2 STUDY DESIGN AND ENDPOINTS  

#### 2.1 <step 1>

use bzip_1 domain from pfam to pull out bzip proteins from ML2.2, Drosophila melanogaster, Nematostella vectensis, and Homo sapiens using hmm2aln.pl

#### 2.2 <step 2>

produce an alignment of residuse that match the hmm

```hmm2aln.pl -hmm=PF00170.hmm --name=bzip_1 --fasta_dir=../01-AA_SEQS --threads=45 --no_clean --nofillcnf=nofill.bzip_1.conf

#### 2.3 <step 2>

generate a maximum likelihood tree

```iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m TEST -pre [output prefix] > iq.out 2> iq.err```


## 3 WORK COMPLETED SO FAR WITH DATES  

05-12-2021	nothing
\<day month year> \<steps above that have been completed>

## 4 LITERATURE REFERENCED  

\<literature>

## APPENDIX

Version&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;Date&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Significant Revisions  
1.1  
1.2  
1.3  
1.4  

