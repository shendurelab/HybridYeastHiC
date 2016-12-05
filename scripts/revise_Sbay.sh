#!/bin/bash
# revise_Sbay.sh
# revise S. bayanus genome, correcting rearrangement on chromosome 3 and adding new sequence to chromosome 3, and swapping chromosomes 10 and 12
# Seungsoo Kim
#
# requires: 
#   Sbay.fa
#   Sbay.gff
#   Sbay3_revision.bed
#   Sbay3_new.bed
#
# output:
#   Sbay_revised.fa
#   Sbay_revised.gff

# reference directory
refs="../references"

# load bedtools module
module load bedtools/2.24.0

# extract sequence of segments of chromosome 3 from genome fasta
bedtools getfasta -fi $refs/Sbay.fa -bed $refs/Sbay3_revision.bed -fo $refs/Sbay3_revision.fa

# create fasta of revised genome
cat <(head -n 5 $refs/Sbay.fa) \
    <(cat <(sed -n '2p' $refs/Sbay3_revision.fa) <(sed -n '2p' $refs/Sbay3_new.fa) <(tail -n 1 $refs/Sbay3_revision.fa)) \
    <(sed -n '7,19p' $refs/Sbay.fa) \
    <(sed -n '24p' $refs/Sbay.fa) \
    <(sed -n '21,23p' $refs/Sbay.fa) \
    <(sed -n '20p' $refs/Sbay.fa) \
    <(tail -n +25 $refs/Sbay.fa) > $refs/Sbay_revised.fa

# revise gff file
awk -F'\t' '{OFS="\t"; 
             if($1=="3") 
               {if ($4 >= 219499) 
                 print $1, $2, $3, $4-219499, $5-219499, $6, $7, $8, $9; 
               else print $1, $2, $3, $4+75853, $5+75853, $6, $7, $8, $9} 
             else print $1, $2, $3, $4, $5, $6, $7, $8, $9}' $refs/Sbay.gff | \
  awk '{OFS="\t";
       if ($1=="10") 
         print "12", $2, $3, $4, $5, $6, $7, $8, $9; 
       else if ($1=="12") 
         print "10", $2, $3, $4, $5, $6, $7, $8, $9; 
       else print $1, $2, $3, $4, $5, $6, $7, $8, $9}' | \
  sort -k1,1n -k4,4n > $refs/Sbay_revised.gff
