## Description

This repository is part of the CELEBRIMBOR pangenome analysis pipeline, it provides rust code that labels genes as core, rare, or neither depending on the number of observations of the gene over all genome samples. The code tries to account for incomplete genome samples by using the genome completeness score from software CheckM. 

The following people have contributed to writing the rust code and fitting it into the CELEBRIMBOR pipeline:
- Joel Hellewell
- John Lees
- Sam Horsfield
- Johanna Von Wachsmann

## Example
You can run the code on on checkM output called `genome_metadata.tsv` and a presence-absence matrix (generated earlier in the CELEBRIMBOR snakemake pipeline) `gene_presence_absence.Rtab`. The `completeness-column 7` argument specifies the column in `genome_metadata.tsv` that contains the completeness score for each genome sample.

First build the crate using `cargo build --release` in this directory. Then you can run the program on the example data provided with the following command:

`target/release/cgt_bacpop example_data/genome_metadata.tsv example_data/gene_presence_absence.Rtab --completeness-column 7`

## Description

CGT tries to account for the fact that genes might be missed (assumed completely at random) during the sequencing process. This leads to underestimation of the true gene frequency in pangenome data, which is bad if we are trying to label the core genome.
To do this it works outs by simulation a new, lower threshold for assigning genes as part of the core genome. We assume that the true gene frequencies for all genes in the genome belong to a Beta distribution (the default priors for this distribution give it a U-shape so genes are assumed to have mostly low or high true frequency).
The probability of observing gene A, with underlying true frequency a, in a genome sample is Bernoulli(ac) where c is the completeness score of the genome sample.
Across N genome samples in a pangenome, the sum of these independent Bernoulli trials with varying probabilities of success (due to different genome completeness scores and different true frequencies) is a Poisson binomial distribution.
We use many samples from the Poisson binomial distribution to calculate the number of observations X for which Probability(gene A observed less than X | gene A has true frequency >= 95%) = 5%
So there is 5% (or less) chance that we have missed core genes in our labelling process by saying that any gene with X or more observations out of N genomes is core. However, it is also possible that genes with less than 95% true frequency appear X or more times and are wrongly labelled as core. We also calculate Prob(gene A observed more than X | gene A has true frequency < 95%) and we report this probability.
Any threshold number of observations, X, for choosing when to label genes as core will have a corresponding probability of missing core genes or wrongly assigning non-core genes as core. We aim to provide the information to the user so that they can choose the balance of these probabilities that they are comfortable with, some people may prefer to miss fewer core genes even though this incurs a greater possibility of wrongly assigning non-core genes as core.

All of the values such as the Beta distribution priors, the 95% true frequency threshold, and the 5% probability of missing core genes can be changed by the user
