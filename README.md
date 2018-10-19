# AffyMergeCluster
A Procedure for Cross-experiment Clustering of Affymetrix Oligonucleotide Microarray

# Usage:
source("MatrixConstructor.R")

cels <- read.cel.files(celPath)
\#\# celPath: a directory contians .CEL.gz files

Probecalls <- binary.expr(cels)

gpl <- cels@annotation

PDOmatrix <- probesetID2PDO(Probecalls, gpl)

clean_matrix <- matrix.shrinker(PDOmatrix)

write.phy(clean_matrix, taxnames, "./tree.phy")
