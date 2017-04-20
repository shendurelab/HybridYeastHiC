#!/bin/bash
# centromere_annotations.sh
# Seungsoo Kim

refs="../references"
bt2="../nobackup/bowtie2"

# load modules
module load samtools/1.2 bedtools/2.24.0 bowtie2/2.2.3

# extract centromeres from S. cerevisiae
awk '$3=="centromere" {OFS="\t"; print "Scer_" $1, $4, $5}' $refs/Scer.gff > $refs/Scer_cen.gff

# extract CEN sequences
bedtools getfasta -fi $refs/Scer.fa -bed $refs/Scer_cen.gff -fo $refs/Scer_cen.fa

# extract centromeres from Suva.gff
awk -F"\t" '/CEN/ {OFS="\t"; print "Suva_" $1, $4, $5}' $refs/Suva.gff > $refs/Suva_cen.gff

# extract centromeres from Suva_revised.gff
awk -F"\t" '/CEN/ {OFS="\t"; print "Suva_" $1, $4, $5}' $refs/Suva_revised.gff > $refs/Suva_revised_cen.gff

# find S. paradoxus centromeres
sed 's/[ACTGN]\{2\}TCAC[AG]TG[ACTGN]\{95,100\}CCGAA[ACTGN]\{6\}/	/g' $refs/Spar_cen_regions.fa | \
  sed 's/[ACTGN]\{6\}TTCGG[ACTGN]\{95,100\}CA[TC]GTGA[ATCGN]\{2\}/	/g'  | \
  awk '{OFS="\t"; if (substr($0,0,1)==">") 
                    print $0; 
                  else print length($1), length($2)}' | \
  sed 's/:/	/' | sed 's/-/	/' | tr '\n' '\t' | sed 's/>/X/g' | tr 'X' '\n' | tail -n +2 | \
  awk '{OFS="\t"; if ($5==0) 
                    print $1, int(($2+$3)/2-60), int(($2+$3)/2+60); 
                  else print $1, $2+$4, $3-$5}' > $refs/Spar_cen.gff
cp $refs/Spar_cen.gff $refs/Spar_revised_cen.gff
bedtools getfasta -fi $refs/Spar.fa -bed $refs/Spar_cen.gff -fo $refs/Spar_cen.fa

# find Y12, DBVPG6044 centromeres
bowtie2 -f --very-sensitive -x $bt2/Y12 -U $refs/Scer_cen.fa | \
  grep Y12 | grep -v '^@' | cut -f3,4,6 | \
  sed 's/M/+/g' | sed 's/D/+/g' | sed 's/[0-9]*I//' | sed 's/$/0/' | \
  while read a b c; do echo $a $b $(($c)); done | \
  awk '{OFS="\t"; print $1, $2, $2+$3}' > $refs/Y12_cen.gff
bowtie2 -f --very-sensitive -x $bt2/DBVPG6044 -U $refs/Scer_cen.fa | \
  grep DBVPG6044 | grep -v '^@' | cut -f3,4,6 | \
  sed 's/M/+/g' | sed 's/D/+/g' | sed 's/[0-9]*I//' | sed 's/$/0/' | \
  while read a b c; do echo $a $b $(($c)); done | \
  awk '{OFS="\t"; print $1, $2, $2+$3}' > $refs/DBVPG6044_cen.gff

# concatenate for hybrid genomes
cat $refs/Scer_cen.gff $refs/Suva_revised_cen.gff > $refs/ScSu_2micron_revised_cen.gff
cat $refs/Scer_cen.gff $refs/Spar_revised_cen.gff > $refs/ScSp_2micron_revised_cen.gff
cat $refs/Spar_cen.gff $refs/Suva_revised_cen.gff > $refs/SpSu_2micron_revised_cen.gff
cat $refs/Y12_cen.gff $refs/DBVPG6044_cen.gff > $refs/Y12xDBVPG6044_cen.gff
