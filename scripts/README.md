**scripts/convergent_table.R** is used to produce Table 1 of the paper.

**scripts/circos.R** is used to produce the circos diagrams which illustrate the important J parameters fitted for different organisms, 
as well as the proportion of times a given nucleic acid was represented at a given site (i.e. sequence logo)
The 36 circular ring boxes represent the 4 nucleotides per 9-site donor sequence. 
Warm colors were used for the three exonic sites whereas cold colors for intronic ones. 
The area of each box is proportional to the nucleotide-site obserevd probability *f<sub>i</sub>(s<sub>i</sub>)*. 
Positive and negative couplings are depicted with green and red lines respectively

**scripts/script_freqs.R** and its GT counterpart were used to generate *f<sub>i</sub>* and *f<sub>ij</sub>* frequencies from the organism/sample's corresponding sequences, from the full sample or only its GT-bearing sequences respectively.

**scripts/script_freqs.R** contains utility functions for calculating the energy of sequences using provided *h<sub>i</sub>* and *h¿J<sub>ij</sub>* parameters.

**scripts/script_corr.R** was used for evaluating two-site correlations *f<sub>ij</sub> - f<sub>i</sub> f<sub>j</sub>* for significance.

**scripts/script_uniques.R** was used for computing number of sequences (total, GT, unique, unique GT).

**scripts/convergent_table.R** was used for computing interactions between/within the exon/intron regions, and outputting these interactions in LaTeX-compatible format.

**scripts/convergence_check.R** was used to check whether sufficient iterations of a given run had elapsed as to ensure convergence.


All scripts should be run using the base directory as their working directory, but parameters within the files can be changed in order to use them elsewhere.