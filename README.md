## Description

This repository is part of the CELEBRIMBOR pangenome analysis pipeline, it provides rust code that labels genes as core, rare, or neither depending on the number of observations of the gene over all genome samples. The code tries to account for incomplete genome samples by using the genome completeness score from software CheckM. 

The following people have contributed to writing the rust code and fitting it into the CELEBRIMBOR pipeline:
- Joel Hellewell
- John Lees
- Sam Horsfield
- Johanna Von Wachsmann

## Example
You can run the code on on checkM output called `genome_metadata.tsv` and a presence-absence matrix (generated earlier in the CELEBRIMBOR snakemake pipeline) `gene_presence_absence.Rtab`. The `completeness-column 7` argument specifies the column in `genome_metadata.tsv` that contains the completeness score for each genome sample.

`target/release/cgt_rust genome_metadata.tsv gene_presence_absence.Rtab --completeness-column 7`

