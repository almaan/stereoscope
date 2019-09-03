#!/usr/bin/bash

# Script to generate joint single cell count and meta files
# from 10X public data found at : https://support.10xgenomics.com/single-cell-gene-expression/datasets/
# in the subsection "Single Cell 3' Paper: Zheng et al. 2017"

# Create necessary folders

for ndir in raw tars all; do
    if [ ! -d $ndir ]; then
        mkdir $ndir
    fi
done

# Extract filtered tar.gz matrices located in tars folder

for ii in *tar.gz; do
    name=$( echo $ii | cut -f 1-2 -d _ )
    mkdir raw/$name
    tar -xf $ii -C raw/$name
    mv $ii tars/
done

# Assemble separate matrices for each one of the cell types

while read -r type; do
    ./assemble_matrix.py -mp raw/$type/filtered_matrices_mex/hg19 -n 500 -t $type -o all
done < types.dat

# Concatenate generated matrices

./concat_matrices.py -cp all/*cnt*tsv -mp all/*mta*tsv -o all
