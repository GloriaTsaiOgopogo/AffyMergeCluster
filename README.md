# AffyMergeCluster
A Procedure for Cross-experiment Clustering of Affymetrix Oligonucleotide Microarray

# Usage:
# step 1: matrix building

source("MatrixConstructor.R")

cels <- read.cel.files(celPath)
\#\# celPath: path of a directory contians .CEL.gz files

Probecalls <- binary.expr(cels)

gpl <- cels@annotation

PDOmatrix <- probesetID2PDO(Probecalls, gpl)

clean_matrix <- matrix.shrinker(PDOmatrix)

write.phy(clean_matrix, taxnames, "./tree.phy")

# step 2: clustering
# compling

gcc bootdist.c -o bootdist [-lm]

gcc neighbor.c -o neighbor [-lm]

gcc consense.c -o .\consense

# running

seqboot.exe [infile] [outfile] [reps]

neighbor_test.exe [infile] [outfile] [reps]

consense.exe [infile] [outfile]

\# reps is an integer(recommend bigger than 100) which indicates bootstrap replicates
