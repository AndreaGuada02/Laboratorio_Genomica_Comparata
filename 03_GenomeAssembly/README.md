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
