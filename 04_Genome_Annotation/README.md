# Genome Annotation:

Once you have a de-novo assembled genome, what you usually want to know is what genomic features are on that genome. For this, we need an annotation pipeline.
* 1.	Annotation of repetitive elements in the genome. Repetitive regions and transposons must be excluded from gene annotation because they have completly different structures.
* 2.	Alignment of high quality external evidences (Transcriptomes and proteins).
* 3.	Extraction of an initial set of highly confident gene models.
* 4.	Training of ab-initio gene predictors softwares.
* 5.	Re - annotation of the genome through trained models.
Steps 3, 4 and 5 are usually repeated multiple times to improve gene prediction at each iteration
* 6.	Collection of all gene models and creation of a consensus set based on support given by external evidences (usually gene models low or not supported by external evidences are discared).
* 7.	Evaluation of our final gene set  Busco.

Annotation can be performed from the beginning using programs that generate HMMs by leveraging protein features in certain species to create predictive models of how specific proteins might be structured.

* A group of proteins from closely related species is modeled this way: if these proteins are also found in the genome of interest, they are input into software that predicts how the proteins are structured.
* Each species can have particular characteristics in protein structure (e.g., codon bias: some codons may be preferred over others for building an amino acid).
* These models predict how the protein should appear within the genome. Even proteins not found by BLAST could still be identified through this predictive step → this is the training phase.

Predictive models are applied iteratively:

* After generating the models, they are extracted to reinforce them.
* Protein searches are then repeated on the genome.
* Models are re-extracted and the process is repeated 2–4 times (avoiding over-parameterization). The iteration stops when the proteome size stabilizes or decreases.
* Models identify proteins with a statistical score based on similarity. Proteins can be filtered based on this score.

Finally, the results are evaluated. BUSCO can also be applied to proteomes: the proteins extracted during annotation are analyzed, and the performance of the models used can also be assessed.

### Maker:

MAKER works with a configuration file → it takes input from a file that must be prepared beforehand. Once MAKER produces a result, it creates a folder. We then need to compress it and generate a summary of the results:

```bash
maker -CTL
```
Now we obtain the three files to modify, in order to use Maker

```bash
maker -base folder_name
```

After completion, the FASTA and GFF files must be merged:

```bash
fasta_merge -d <DATASTORE INDEX FILE>
gff3_merge -d <DATASTORE INDEX FILE>
```

We can also summarize results of Evidence-based gene annotation and RepeatMasker using the AGAT package, a very usefull set of perl scripts to manage gff3 files, print the help and run the script:

### AGAT:

We can also summarize results of Evidence-based gene annotation and RepeatMasker using the AGAT package, a very usefull set of perl scripts to manage gff3 files, print the help and run the script:

```bash
agat_sp_statistics.pl --gff file.gff -o <output_file>
agat_sq_repeats_analyzer.pl -i <input_file> -o <output_file>
```
Summary statistics of gene models (1) and of repeats (2)

### SNAP: a gene prediction tool that uses previously obtained data.

SNAP is designed to identify protein-coding genes in genomic DNA sequences, particularly in newly sequenced genomes. It uses a probabilistic approach based on HMMs to model various gene features and can be trained with a set of known gene models to improve accuracy for a specific organism.
* It is often used alongside other gene prediction tools, and results from multiple tools can be combined to create a consensus annotation.

### Augustus: (can also be run via BUSCO)

Augustus training is more complex and computationally intensive but is one of the highest-performing predictors.

Since this process is computationally heavy and the trained models must be in a specific path where Augustus will search for them, I have already trained Augustus for you. The model is called Aste (just use Aste when required).

### Evaluating genome annotation:

* We have multiple options to evaluate the predicted genome annotation:
* Compare gene statistics with bibliographic knowledge (e.g., mean gene length, mean exon length, mean number of exons and introns per gene, etc.).
* Run BUSCO on the predicted proteins.
* Summarize AED values.
* Align the de novo assembled transcriptome to the genome-based predicted annotation.
