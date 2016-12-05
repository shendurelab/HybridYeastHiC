#!/bin/bash
# revise_Spar.sh
# revise S. paradoxus genome
# Seungsoo Kim
#
# requires:
#   Spar.fa
#   Spar.gff
#   Spar_revision.fa
#
# output:
#   Spar_revised.fa
#   Spar_revised.gff

refs="../references"

# load bedtools module
module load bedtools/2.24.0

# extract chromosome 4 sequences
bedtools getfasta -fi $refs/Spar.fa -bed $refs/Spar_revision.bed -fo $refs/Spar_revision.fa

# create revised fasta file
cat <(head -n 7 $refs/Spar.fa) <(grep -v '>' $refs/Spar_revision.fa | head -n 3) <(sed -n '8p' $refs/Spar_revision.fa | tr 'ATCG' 'TAGC' | rev) <(tail -n 1 $refs/Spar_revision.fa) <(tail -n +9 $refs/Spar.fa) > $refs/Spar_revised.fa

# revise gff file
awk -F"\t" '{OFS="\t"; 
            if($1=="4" && $4 >= 943469 && $4 < 1027717)
              {if ($7=="+") 
                 print $1, $2, $3, 2136497-$5, 2136497-$4, $6, "-", $8, $9; 
               else 
                 print $1, $2, $3, 2136497-$5, 2136497-$4, $6, "+", $8, $9} 
            else if ($1=="4" && $4 >= 1027717 && $4 < 1029252) 
              print $1, $2, $3, $4+79528, $5+79528, $6, $7, $8, $9; 
            else if ($1=="4" && $4 >= 1029252 && $4 < 1193028) 
              print $1, $2, $3, $4-85783, $5-85783, $6, $7, $8, $9; 
            else 
              print $1, $2, $3, $4, $5, $6, $7, $8, $9}' $refs/Spar.gff | sort -k1,1n -k4,4n > $refs/Spar_revised.gff
