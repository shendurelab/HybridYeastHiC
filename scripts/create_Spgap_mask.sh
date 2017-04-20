#!/bin/bash
# create_ScSp_mask.sh
# Seungsoo Kim

# load modules
module load samtools/1.2 bedtools/2.24.0

dir='../nobackup/aligned'
Sparlib='ScSp/YDG613_exponential_Sau3AI'
Scerlib='ScSp/YMD1797_raffinose_Sau3AI'

# find all regions where more than 1 S. paradoxus read maps with MAPQ >= 30 to S. cerevisiae (when mapped to combined reference)
cat <(samtools view $dir/$Sparlib.1.bam) <(samtools view $dir/$Sparlib.2.bam) |\
  grep Scer | grep -v MT | awk '$5 >= 30 {OFS="\t"; print $3, $4, $4+80}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | awk '$4 > 1' > ../masks/Spgap.bed
# find all regions where more than 1 S. cerevisiae read maps with MAPQ >= 30 to S. paradoxus (when mapped to combined reference)
cat <(samtools view $dir/$Scerlib.1.bam) <(samtools view $dir/$Scerlib.2.bam) |\
  grep Spar | grep -v MT | awk '$5 >= 30 {OFS="\t"; print $3, $4, $4+80}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | awk '$4 > 1' >> ../masks/Spgap.bed
