# Kmer - Based Genome Survey
This chapter allows you to obtain information about the short-reads (obtained through sequencing) quality.


## Reads quality check:
Controlling the short-reads quality

```bash
#|assembly|
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```


## Trimmomatic
Trimming ("cutting") the reads, removing non-optimal portions according to parameters chosen by us.

```bash
#|assembly|
trimmomatic PE -threads 8 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log 
```


## KAT
It analyzes the k-mers obtained after the trimming.
- '-t' : inform the program about how many thread (parallelisation) it should use (4).
- '-m' : is the dimension of the k-mer we will use (27).
- '-o' : determines the prefix of the output.
- finally, append the two input files.

```bash
#|kat|
kat hist -t 6 -m 27 -o Anoste SRR11672503_1_paired_fastqc.html SRR11672503_2_paired_fastqc.html
```
