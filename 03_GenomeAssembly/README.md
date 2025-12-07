g# Genome Assembly
A genome assembly is a computational representation of a genome sequence. Because we are not able to sequence along the complete length of a chromosome, each chromosome assembly is made up of short stretches of sequenced DNA pasted together.

* Contig assembly + error correction (depending on the type of reads)
* Scaffolding
* Gap filling

 ### Some defitinitions

**Contig** : “a contiguous sequence generated from determining the non-redundant path along an order set of component sequences. A contig should contain no gaps but often the terms contig and scaffold are used interchangeably”

**Scaffold** : “an ordered and oriented set of contigs. A scaffold will contain gaps, but there is typically some evidence to support the contig order, orientation and gap size estimates.”

**Assembly** : “a set of chromosomes, unlocalized and unplaced (random) sequences and alternate loci used to represent an organism’s genome. Most current assemblies are a haploid representation of an organism’s genome, although some loci may be represented more than once (see Alternate locus, above). This representation may be obtained from a single individual (e.g. chimp or mouse) or multiple individuals (e.g. human reference assembly). Except in the case of organisms which have been bred to homozygosity, the haploid assembly does not typically represent a single haplotype, but rather a mixture of haplotypes. As sequencing technology evolves, it is anticipated that diploid sequences representing an individual’s genome will become available.”

**Diploid assembly** : “A genome assembly for which a Chromosome Assembly is available for both sets of an individual’s chromosomes, as defined by the NCBI Assembly model. It is anticipated that a diploid genome assembly is representing the genome of an individual. Therefore it is not anticipated that alternate loci will be defined for this assembly, although it is possible that unlocalized or unplaced sequences could be part of the assembly.”

**Primary Assembly** : “Relevant for haploid assemblies only. The primary assemblies represents the collection of assembled chromosomes, unlocalized and unplaced sequences that, when combined, should represent a non-redundant haploid genome. This excludes any of the alternate locus groups.”

**Unlocalized sequence/scaffold** : “A sequence found in an assembly that is associated with a specific chromosome but cannot be ordered or oriented on that chromosome.”

**Unplaced Sequence/scaffold** : “A sequence found in an assembly that is not associated with any chromosome.”

Source: [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc/help/definitions/#:~:text=The%20primary%20assemblies%20represents%20the,a%20non%2Dredundant%20haploid%20genome.)

### What we're going to do:

Once the scaffolding step is completed, I can reuse the reads (including those that previously failed to align) to close the remaining gaps (reads that originally mapped across the boundaries of two contigs were not used during the polishing stage).
In the end, we obtain a Genome Draft: with chromosomal information available, we can state that an entire scaffold corresponds to a chromosome.
This allows us to assemble the full genome of a species as accurately and contiguously as possible, enabling further study.

## Methods for assessing assembly quality:

### N50: the length of the shortest contig whose cumulative length accounts for at least 50% of the entire assembly.
The larger this value is, the more contiguous the assembly.
```bash
#|assembly|
assembly-stats Anoste_raw.fasta > Anoste_raw.stats
```
### BUSCO: a tool that performs sequence searches to compare our assembly against a curated dataset of genes known to be present in a given group of organisms.
This dataset contains a set of genes shared by all organisms studied within that group. Naturally, the closer the taxa are evolutionarily (i.e., the narrower the group), the larger the pool of shared genes.
The goal is to compare the assembled genome with this dataset to check whether the expected genes are present. If many are missing, it likely indicates problems in the assembly process.
```bash
#|sequence|
export NUMEXPR_MAX_THREADS=80
busco -m geno -l $BUSCO/culicidae_odb12 -c 6 -o Anoste_raw_busco -i Anoste_raw.fasta 
```

The first `export` is needed to set the `NUMEXPR_MAX_THREADS` equals to the number of the cores that we will be using during the real busco command. The option are:

- '-m' #define if we are using genomes or proteomes
- '-l' #redirect to the link of the library of reference
- '-c' #the numebr of cores for parallelise the process
- '-o' #name of the output **folder**
- '-i' #the input. In out case the raw assembly

### spectra-cn (KAT): remaps k-mers onto the assembly to visualize how they have been used.
It provides a representation of the frequency of k-mers included or not included in the assembly.

## Polishing
We will now use the short reads to refine the assembly → genome polishing.

### Mapping of the reads
Mapping: the process by which we realign reads to a specific assembly.

>We're creating the file mapping.sh and writing the following code inside of it:

```bash
#|assembly|
#Short reads
minimap2 -ax sr –MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq Anoste_raw.fasta SRR11672503_2_paired.fastq > Anoste_raw.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw_sr.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr_sorted.bam
rm Anoste_pol_sr.bam

#Long Reads
minimap2 -ax sr –MD -t 6 Anoste_raw.fasta SRR11672503.fastq.gz > Anoste_raw_lr.sam 
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_pol_lr.bam
```

I'm then saving this file and after that:

```bash
-	bash mapping.sh
```

### Assembly polishing:
#### Hypo:
Hypo is used to improve (polish) a draft assembly (contigs) by leveraging raw reads (short reads and/or long reads) mapped onto the assembly.
Its goal is to correct errors, indels, and incorrect bases to obtain a more accurate version of the genome.

*	R1 = SRR11672503_1_paired.fastq
*	R2 = SRR11672503_2_paired.fastq

```bash
#|assembly|
echo e- "$R1\n$R2" > Sr.path
hypo -d Anoste_raw.fasta -r @Sr.path -s 227054799 -c 136 -b Anoste_raw_sr_sorted.bam -B Anoste_raw_lr_sorted.bam -t 6 
```
### Assembly quality check:
Now we need to perform assembly quality control again using N50 and BUSCO.
assembly-stats is a program (often installed via conda/bioconda) that reads a FASTA assembly file and calculates:

* number of contigs
* total assembly size
* maximum and minimum contig length
* N50, N75, N90
* L50, L75, L90
* GC%
* other useful metrics for quickly assessing assembly quality

#### N50 
```bash
#|assembly|
assembly-stats Anoste_pol.fasta > Anoste_pol.stats
```
#### Busco
```bash
#|sequence|
busco -m geno -l $BUSCO/culicidae_odb12 -c 8 -o Anoste_pol_busco -i Anoste_pol.fasta
```
Looking at the BUSCO results: my assembly genomically contains 98.8% of the 7,027 coding genes found across the 13 culicid genomes included in the Culicidae dataset used for the analysis.
When we broaden the clade analyzed with BUSCO, the number of orthologs shared by all species decreases (the genomic core defining that clade becomes smaller as more, and more diverse, species are included).

* This is why the analysis using the Culicidae database yields more orthologs than the analysis using the Diptera database.

#### spectra-cn (KAT)
```bash
#|kat|
kat comp -t 8 -o Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
```
KAT (kat comp): compares the k-mer distribution between the raw reads (input) and the assembly (FASTA).
It helps us understand which read sequences were included in the assembly, which parts were missed (k-mers present in the reads but absent from the assembly → missing), whether there are duplications (k-mers over-represented in the assembly → duplicated), and whether the assembly error rate is high or low.

* KAT produces a k-mer comparison plot (usually a .png) and accompanying statistics files.
##### K-mer comparison plot (obtained):

While GenomeScope shows the distribution of k-mer variability, here we visualize the multiplicity of k-mers with coverage over the assembly. The red color indicates k-mers that are duplicated only once.
This means that our distinct k-mers were used only once across the entire assembly (there is no genome duplication complexity).

Sometimes, there may be regions of the genome with high multiplicity, but the k-mers are duplicated (shown in different colors on the plot).
We also see 0x coverage (black), which represents k-mers that were not used to produce the assembly. These could be sequencing errors, k-mer calculation errors, or k-mers from contaminants (from other assemblies).
It could also represent alternative alleles (in heterozygous regions, only one of the two alleles is included in the assembly).

# Decontamination:

By sequencing and extracting DNA, we also obtained environmental bacterial DNA, DNA from the organism itself (including the gut), and viral DNA.
Tools are needed to remove reads that do not belong to the taxon of interest. This is done partially using taxonomic annotation (e.g., BLAST → assigning a taxonomic identifier to each contig).
Further cleaning is performed based on assembly characteristics (e.g., differing GC content depending on the taxon or DNA type, such as mitochondrial DNA).

### 1-Mapping of the contaminants:

I'm creating the file mapping_contaminants.sh, that'll contain my pipeline:

```bash
#|assembly|
minimap2 -ax sr –MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq Anoste_raw.fasta SRR11672503_2_paired.fastq | samtools view -Sb - > Anoste_pol_sr.bam
samtools sort -@6 -o Anoste_pol_sr_sorted.bam Anoste_pol_sr.bam
samtools index Anoste_pol_sr_sorted.bam
rm Anoste_pol_sr.bam
```
### 2-Taxonomic annotation of the contigs:

We need to use BLASTn: a BLAST process done on nucleotides.

This is the structure:

```bash
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
```

#### Building BlobDB:

For this, we'll need the following files:

* Anoste_pol.fasta
* Anoste_pol_sorted.bam
* Anoste_pol_sorted.bam.bai
* Anoste_blast.tsv

```bash
blobtools create -i <ASSEMBLY> -b <MAPPING_FILE> -t <BLASTN_FILE> -o <OUTPUT_PREFIX>
```

Now we're launching the following commands, using the files created with the last one:

```bash
* blobtools view -i <JSON_FILE> -o <OUTPUT_PREFIX>
* blobtools plot -i <JSON_FILE> -o <OUTPUT_PREFIX>
```
#### BlobDB results:

It generates two PNG files:

* Covariance
* Plot → showing identified phyla along with GC content and coverage

It is also possible to create similar plots at the genus or species level to achieve higher taxonomic resolution.
Caution is needed when filtering at very specific levels, especially for poorly studied groups, as this may result in over-filtering.

BlobDB works with phyla, showing the percentage of mapped reads:

* 91% of our contigs are identified as Arthropoda
* 6% of reads are Pseudomonas (gut)
* Other stuff in a small %

#### Removing the non-Arthropoda portions:

We need to remove the portions of the assembly that are non-Arthropoda (since we're studying a mosquito).

* First, check the table.txt file just created (less Anoste.Anoste_blob.blobDB.table.txt)

I want to keep only the contigs that have “Arthropoda” in the third-to-last column.

Save the list of Arthropoda contigs to a new file:

```bash
grep "Arthropoda" Anoste.Anoste_blob.blobDB.table.txt | wc -l > contig_arthropoda.tsv
```

#### Decontaminating:

```bash
awk ' { if ((NR>1)&&($0 ~/^>/)) { printf("\n%s", $0) } else if (NR==1) { printf("\t%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta | grep -w -v -Ff <(cut -f1 contig_arthropoda.tsv) - | tr "\t" "\n" > Anoste_decontaminated.fasta
```
The first part ensures the FASTA is in oneline form and reformat the file obtaing ">header\tsequence". In this way it is easier to process the file with grep. Then filters out everything is contained inside the patterns file, outputting a fasta file oneline with only desired contigs.

Using the list of “good” contigs contained in the file contig_arthropoda.tsv, a new FASTA file (Anoste_decontaminated.fasta) is generated, containing only the contigs to keep (those not removed by the filtering step).

Using awk, I obtained a file where each line is concatenated directly to its header (one-line format) [awk is a program used to interact with columns in a file].

Single quotes after awk indicate what actions to perform; curly braces {} define a sequence of commands (e.g., “do this, then do that, etc.”).

The file or standard output to process is specified afterward.

* The - symbol in a command represents the standard output of the previous command.
* -Ff allows splitting each line of a file into fields based on a delimiter.
* cut -F1 extracts the first field (i.e., the first column).
* < means: treat what follows as input for a temporary file → it’s a command substitution.
* -w ensures that only the exact pattern is matched; without it, searching for ctg1 could also match ctg10, ctg11, etc.
* \n represents a new line.

After this, the following grep keeps only the contigs that are not contaminants (thus retaining Arthropoda contigs).

### 3-Scaffolding

After polishing and decontamination, scaffolding is performed: given prior information (a reference genome or experiments indicating contig proximity), we assemble contigs as accurately as possible, allowing gaps represented by “N”.

The program RagTag matches our contigs to a reference genome (e.g., contig27 is likely on chromosome 5).
It is a collection of software tools for scaffolding and improving modern genome assemblies. Its tasks include:

* Homology-based misassembly correction
* Homology-based assembly scaffolding and patching
* Scaffold merging

```bash
ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
```
-C moves the contigs without a location in the 0-Chromosome, that contains "junk": we agree to eliminate it afterwards, losing a part of the information.
