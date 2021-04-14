**How are features detected?**
The submitted plasmid sequence is queried against several databases using [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [DIAMOND](https://diamond.readthedocs.io), and [Infernal](http://eddylab.org/infernal/). Features are predicted if there is a match with ≥95% sequence identity to the database. Matches with significant overlap are filtered so that only the best consensus set of overall features is shown.

**What databases are used?**
pLannotate finds nucleotide matches to a set of features extracted from [AddGene](https://www.addgene.org/) plasmids and identifies translated nucleotide matches to entries in [fpbase](https://www.fpbase.org) and the [Swiss-Prot database](https://www.uniprot.org/statistics/Swiss-Prot). Only the subset of Swiss-Prot with an [annotation score](https://www.uniprot.org/help/annotation_score) of 3 or higher is used. [Rfam](https://rfam.xfam.org/) is used to detect non-coding RNAs.

**What is the "%" assigned to annotations?**
This is the percentage of the database feature that is matched in the plasmid. It is calculated as the percent match length multiplied by the percent identity within the matching region.

**What do unfilled features on the plot mean?**
These are *incomplete features*; the sequence match in the plasmid covers less than 95% of the full length of the feature in the database.

**Why should I care about incomplete features?**
These elements may be leftover fragments from earlier cloning steps used to create a plasmid. If they include only a small fraction of the feature, they likely do not still have the annotated function. However, even small feature fragments may affect plasmid function if they result in cryptic gene expression or are inadvertently combined with other elements during later cloning steps.

**What is "Linear plasmid annotation"?**
pLannotate can also annotate engineered DNA sequences that are not circular, such as [pJAZZ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847241/) linear plasmids, or [Sanger sequencing](https://en.wikipedia.org/wiki/Sanger_sequencing) fragments. pLannotate still displays the plot as a circular plasmid map, though the break between the ends is denoted with a heavy black line. 