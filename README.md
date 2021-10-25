# QuasiSpecies Analyser (QSA)

Quasispecies viruses are viral populations composed of an ensemble of variant sequences.

The program provides the possibility of carrying out analyses on quasispecies viruses data from BAM files.

# Features
This tool is far from finished, but it can already be used to perform some initial analyses.

## Efficiency (Normalized Entropy)
The efficiency of each position in a sequence can provide useful data on the mutability of the sequence. It is then possible to understand if some specific regions of the analysed sequence are less susceptible to selective pressure.

Shannon's entropy is used to calculate entropy values at each position in the sequence. Afterwards, these values are divided by the maximum value of Shannon's entropy to find the efficiency at each position.

## α-diversity
A sample's α-diversity is calculated as the sum of the entropy at each position in the sequence normalized by the sequence's length. This normalisation is necessary to obtain comparable values.

## β-diversity
The β-diversity between two samples is defined as the difference in α-diversity.

The matrix obtained by calculating the difference between each sample's α-diversity is then used to create an undirected weighted graph. The resulting network can be used to visualise the difference in mutability between the samples.
