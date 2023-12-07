# Commands used in Mitchell et al. - Mnemiopsis cydippid  regeneration 
#### 1\. Run trimmomatic on all FASTQ files ($f = the prefix of the file followed by "_R1_001.fastq.gz")

`trimmomatic PE -basein ${f} -baseout ${f}trimmed LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

#### 2\. Run RSEM with bowtie2 on trimmed reads to estimate counts. ($f = the prefix of the file followed by "_R1_001.fastq.gz")

requires ML2.2.nt which can be downloaded from https://research.nhgri.nih.gov/mnemiopsis/download/transcriptome/ML2.2.nt.gz

`rsem-calculate-expression --bowtie2 --paired-end ${f} *_R1_001.fastq.gztrimmed_2P ML2.2.nt ${f}RSEMout`

#### 3\. BLASTP Mnemiopsis protein models against the human protein database 

requires ML2.2.aa which can be downloaded from https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

requires HumRefSelect2020 which was created using reduce_refseq (https://github.com/josephryan/reduce_refseq) and is available in the 03-DATA_FILES directory of this repo

`blastp -db HumRefSelect2020.fa -query ML2.2.aa.fa -out ML2.2_v_HumRefSelect2020.blastp -evalue .001 -outfmt 6`

#### 4\. Parse the BLAST output and create a tab-delimited file that includes each ML2.2 protein, its top human hit, and e-value of the hit.

requires parse_blast.pl which is available in the 02-SCRIPTS directory of this repo

`perl parse_blast.pl ML2.2_v_HumRefSelect2020.blastp > parse_blast.out`

#### 5\. Use hmm2aln.pl to identify protein sequences with BZIP and ETS domains from human, Drosophila, Nematostella, and Mnemiopsis.

requires: HMMer (http://hmmer.org/) and hmm2aln.pl (https://github.com/josephryan/hmm2aln.pl) and the PFAM domains for each of these domains (https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00170?annotation=hmm and https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00178?annotation=hmm) and the nofill.bzip.conf file in the 03-DATA_FILES` directory of this repo

`hmm2aln.pl --hmm=PF00170.hmm --name=BZIP --fasta_dir=HsDmNvMl_seqs --threads=40 --no_clean --nofillcnf=nofill.bzip.conf > bzip.aln.fa

`hmm2aln.pl --hmm=PF00178.hmm --name=ETS --fasta_dir=HsDmNvMl_seqs --threads=40 --no_clean > ets.aln.fa

#### 6\. Run iqtree 

requires IQ-tree (http://www.iqtree.org/)

`iqtree-omp -s bzip.aln.fa -nt AUTO -bb 1000 -m TEST -pre bzip_1_4species > iq.bzip.out`

`iqtree-omp -s est.aln.fa -nt AUTO -bb 1000 -m TEST -pre ets_4species > iq.out`

#### 7\. Run differential gene expression R script in 04-RSCRIPTs directory in this repo

