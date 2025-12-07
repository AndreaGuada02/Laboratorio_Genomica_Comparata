# Genome Assembly
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
conda activate assembly
assembly-stats Anoste_raw.fasta > Anoste_raw.stats
```
### BUSCO: a tool that performs sequence searches to compare our assembly against a curated dataset of genes known to be present in a given group of organisms.
This dataset contains a set of genes shared by all organisms studied within that group. Naturally, the closer the taxa are evolutionarily (i.e., the narrower the group), the larger the pool of shared genes.
The goal is to compare the assembled genome with this dataset to check whether the expected genes are present. If many are missing, it likely indicates problems in the assembly process.
```bash
conda activate sequence
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

# Polishing
Operations in order to polish the genome.

### Mapping of the reads
Mapping the short-reads on the whole genome.

>Il codice in questione è stato scritto all'interno del file mapping.sh
```bash
#|assembly|
minimap2 -ax --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired_fastq SRR11672503_1_paired_fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw_sr.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr_sorted.bam
rm Anoste_raw_sr.bam
```

Mappatura delle long reads sul genoma assemblato

```bash
#|assembly|
minimap2 -ax --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired_fastq SRR11672503_1_paired_fastq > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_raw_lr.bam
```


### Pulizia dell'assemblaggio
Miglioramento della sequenza del genoma frazie al confronto con le reads mappate precedentemente

```bash
#|assembly|
echo e- "$R1\n$R2" > Sr.path
hypo -d Anoste_raw.fasta -r @Sr.path -s 227054799 -c 136 -b Anoste_raw_sr_sorted.bam -B Anoste_raw_lr_sorted.bam -t 6 
```


### Controllo qualità del genoma pulito
Verifica delle statistiche inerenti al genoma di i passaggi di pulizia. il controllo qualitativo è svolto con gli stessi metodi adoperati in precedenza con l'assemblaggio "raw" del genoma.

##### N50 
```bash
#|assembly|
assembly-stats Anoste_pol.fasta > Anoste_pol.stats
```
##### Busco
```bash
#|sequence|
busco -m geno -l $BUSCO/culicidae_odb12 -c 8 -o Anoste_pol_busco -i Anoste_pol.fasta
```

##### Spectra-cn (KAT)
```bash
#|kat|
kat comp -t 8 -o Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
```
