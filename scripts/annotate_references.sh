#!/bin/bash
# annotate_references.sh
# Seungsoo Kim

# load modules
module load bowtie2/2.2.3

# PARAMETERS
# directory for output
out="../nobackup"
# directory for reference files
dir="../nobackup/annotations"
if [ ! -d "$dir" ]; then
  mkdir "$dir"
fi
# subdirectory for bowtie2 indices
bt2="../nobackup/bowtie2"
if [ ! -d "$bt2" ]; then
  mkdir "$bt2"
fi
# file with restriction enzymes - tab delimited file with following columns:
# 1) restriction enzyme name, 2) restriction site, 3) number of nt offset for cut
restriction_enzymes="../restr_enz.txt"
# file with references - tab-delimited file with following columns:
# 1) reference name, 2) reference fasta file path, 3) reference centromere annotation file path, 4) list of binsizes
refs="../refs.txt"

while read ref ref_file ref_cens binsizes
do
  for binsize in $(echo $binsizes | tr "," "\n")
  do
    #bed file of bins
    paste <(awk -v b=$binsize 'BEGIN{OFS="\t"; i=0}
                                    {for (j = 0; j < $2; j = j + b){
                                      if ($2 < j + b)
                                        print $1, j, $2, i;
                                      else print $1, j, j + b, i;
                                      i = i + 1}}' $dir/$ref.chrom_lengths) > $dir/$ref.$binsize.bed

    #annotations by chromosome - chromosome name, bin start, bin end, bin cen, chromosome number
    paste <(awk -v b=$binsize 'BEGIN{OFS="\t"; start = 0; end = 0}
                                    {end = end + int($2/b) + 1; print $1, start, end; start = start + int($2/b) + 1}' $dir/$ref.chrom_lengths) \
          <(awk -v b=$binsize '{OFS="\t"; print int(($2+$3)/(2*b))}' ../${ref_cens}) | \
      awk '{OFS="\t"; print $1, $1, $2, $3, $2+$4}' | \
      sed 's/_/        /' | \
      awk '{OFS="\t"; print $3, $4, $5, $6, $2}' > $dir/$ref.$binsize.chr_annotations

    #annotations by bin - bin number, chromosome name, distance from centromere, chromosome arm length
    awk '{OFS="\t"; for (i=$2; i<$3; i++) 
                      if (i<$4) print i, $1, $4-i, $4-$2; 
                      else if (i==$4) print i, $1, 0, 0;
                      else print i, $1, i-$4, $3-$4-1}' $dir/$ref.$binsize.chr_annotations | \
      awk '{OFS="\t"; print $2, $0}' | \
      sed 's/_/   /' | \
      awk '{OFS="\t"; print $3, $4, $5, $6, $7, $2}' > $dir/$ref.$binsize.bin_annotations

    #create subtelomere bed file
    awk -v b=$binsize '{OFS="\t"; print $1, 0, b; print $1, $2-b, $2}' $dir/$ref.chrom_lengths | awk '{OFS="\t"; print $1, $2, $3, NR-1}' > $dir/$ref.$binsize.subtel.bed

    #create pericentromeric bed file
    awk -v b=$binsize '{OFS="\t"; print $1, int(($2+$3-b)/2), int(($2+$3+b)/2), NR-1}' ../${ref_cens} > $dir/$ref.$binsize.pericen.bed
  done
done < $refs
