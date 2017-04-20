#!/bin/bash
# revise_Suva.sh
# revise S. uvarum genome, correcting rearrangement on chromosome 3 and adding new sequence to chromosome 3, and swapping chromosomes 10 and 12
# Seungsoo Kim
#
# requires: 
#   Suva.fa
#   Suva.gff
#   Suva3_revision.bed
#   Suva3_new.bed
#
# output:
#   Suva_revised.fa
#   Suva_revised.gff

# reference directory
refs="../references"

# load bedtools module
module load bedtools/2.24.0

# extract sequence of segments of chromosome 3 from genome fasta
bedtools getfasta -fi $refs/Suva.fa -bed $refs/Suva3_revision.bed -fo $refs/Suva3_revision.fa

# create fasta of revised genome
cat <(head -n 5 $refs/Suva.fa) \
    <(cat <(sed -n '2p' $refs/Suva3_revision.fa) <(sed -n '2p' $refs/Suva3_new.fa) <(tail -n 1 $refs/Suva3_revision.fa)) \
    <(sed -n '7,19p' $refs/Suva.fa) \
    <(sed -n '24p' $refs/Suva.fa) \
    <(sed -n '21,23p' $refs/Suva.fa) \
    <(sed -n '20p' $refs/Suva.fa) \
    <(tail -n +25 $refs/Suva.fa) > $refs/Suva_revised.fa

# revise gff file
awk -F'\t' '{OFS="\t"; 
             if($1=="3") 
               {if ($4 >= 219499) 
                 print $1, $2, $3, $4-219499, $5-219499, $6, $7, $8, $9; 
               else print $1, $2, $3, $4+75853, $5+75853, $6, $7, $8, $9} 
             else print $1, $2, $3, $4, $5, $6, $7, $8, $9}' $refs/Suva.gff | \
  awk '{OFS="\t";
       if ($1=="10") 
         print "12", $2, $3, $4, $5, $6, $7, $8, $9; 
       else if ($1=="12") 
         print "10", $2, $3, $4, $5, $6, $7, $8, $9; 
       else print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
  sort -k1,1n -k4,4n > $refs/Suva_revised.gff
