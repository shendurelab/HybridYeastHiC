#!/bin/bash
# create_ScSp_mask.sh
# Seungsoo Kim

# load modules
module load samtools/1.2 bedtools/2.24.0

# find all regions where more than 1 S. paradoxus read maps with MAPQ >= 30 to S. cerevisiae (when mapped to combined reference)
cat <(samtools view ../nobackup/aligned/ScSp_2micron_revised/YDG613_asynch.1.bam) <(samtools view ../nobackup/aligned/ScSp_2micron_revised/YDG613_asynch.2.bam) |\
  grep Scer | grep -v MT | awk '$5 >= 30 {OFS="\t"; print $3, $4, $4+80}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | awk '$4 > 1' > ../masks/Spgap.bed
# find all regions where more than 1 S. cerevisiae read maps with MAPQ >= 30 to S. paradoxus (when mapped to combined reference)
cat <(samtools view ../nobackup/aligned/ScSp_2micron_revised/YMD1797_raff.1.bam) <(samtools view ../nobackup/aligned/ScSp_2micron_revised/YMD1797_raff.2.bam) |\
  grep Spar | grep -v MT | awk '$5 >= 30 {OFS="\t"; print $3, $4, $4+80}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | awk '$4 > 1' >> ../masks/Spgap.bed
