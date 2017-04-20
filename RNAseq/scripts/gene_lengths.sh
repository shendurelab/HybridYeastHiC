#!/bin/bash
# gene_lengths.sh
# extracts gene lengths from GFF file
# Seungsoo Kim

dir='../nobackup/references'
sed 's/;.*//' $dir/sacCer3_genes.gff | sed 's/ID=//' | awk '{OFS="\t"; print $9, $5-$4}' > $dir/sacCer3_genes_length.txt
