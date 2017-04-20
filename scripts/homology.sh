#!/bin/bash
# homology.sh
# Seungsoo Kim

refs="../references"
anns="../nobackup/annotations"
bt2="../nobackup/bowtie2"
out="../nobackup"

# load modules
module load bowtie2/2.2.3 samtools/1.2

bsize=32000

# INTERSPECIFIC HYBRIDS
# for interspecies homology, find homologous genes and look at end positions

# find Scer genes
echo "Finding S. cerevisiae genes"
awk '$3=="CDS"' $refs/Scer.gff | grep -v "MT" | sed 's/Parent.*Name=//' | sed 's/_CDS.*//' | awk '{OFS="\t"; print $9, $1, $4, $5, $7}' | sort -k1,1 > $anns/Scer_genes.txt
# find Spar unique genes
echo "Finding S. paradoxus unique genes"
awk '$3=="CDS"' $refs/Spar_revised.gff | sed 's/;ncbi.*//' | sed 's/ID.*SGD=//' | sed 's/;.*BLAST=/	/' | awk -F"\t" '$9==$10 {OFS="\t"; print $9, $1, $4, $5, $7}' | grep Y | sort -k1,1 > $anns/Spar_genes.txt
grep -f <(cut -f1 $anns/Spar_genes.txt | uniq -c | awk '$1==1 {print $2}') $anns/Spar_genes.txt > $anns/Spar_genes_unique.txt
awk '$3=="CDS"' $refs/Spar.gff | sed 's/;ncbi.*//' | sed 's/ID.*SGD=//' | sed 's/;.*BLAST=/	/' | awk -F"\t" '$9==$10 {OFS="\t"; print $9, $1, $4, $5, $7}' | grep Y | sort -k1,1 > $anns/Spar_original_genes.txt
grep -f <(cut -f1 $anns/Spar_original_genes.txt | uniq -c | awk '$1==1 {print $2}') $anns/Spar_original_genes.txt > $anns/Spar_original_genes_unique.txt
# find Suva unique genes
echo "Finding S. bayanus unique genes"
awk '$3=="CDS"' $refs/Suva_revised.gff | sed 's/;ncbi.*//' | sed 's/ID.*SGD=//' | sed 's/;.*BLAST=/	/' | awk -F"\t" '$9==$10 {OFS="\t"; print $9, $1, $4, $5, $7}' | grep Y | sort -k1,1 > $anns/Suva_genes.txt
grep -f <(cut -f1 $anns/Suva_genes.txt | uniq -c | awk '$1==1 {print $2}') $anns/Suva_genes.txt > $anns/Suva_genes_unique.txt
awk '$3=="CDS"' $refs/Suva.gff | sed 's/;ncbi.*//' | sed 's/ID.*SGD=//' | sed 's/;.*BLAST=/	/' | awk -F"\t" '$9==$10 {OFS="\t"; print $9, $1, $4, $5, $7}' | grep Y | sort -k1,1 > $anns/Suva_original_genes.txt
grep -f <(cut -f1 $anns/Suva_original_genes.txt | uniq -c | awk '$1==1 {print $2}') $anns/Suva_original_genes.txt > $anns/Suva_original_genes_unique.txt

# gene homology for original references
join $anns/Scer_genes.txt $anns/Spar_original_genes_unique.txt | awk -v a=Scer -v b=Spar '{OFS="\t"; print $1, a "_" $2, $3, $4, $5, b "_" $6, $7, $8, $9}' > $anns/Scer_Spar_original_gene_homology.txt
join $anns/Scer_genes.txt $anns/Suva_original_genes_unique.txt | awk -v a=Scer -v b=Suva '{OFS="\t"; print $1, a "_" $2, $3, $4, $5, b "_" $6, $7, $8, $9}' > $anns/Scer_Suva_original_gene_homology.txt

abbr=( 'Sc' 'Sp' 'Su' )
spec=( 'Scer' 'Spar_original' 'Suva_original' )
spec2=( 'Scer' 'Spar' 'Suva' )
genes=( 'Scer_genes.txt' 'Spar_genes_unique.txt' 'Suva_genes_unique.txt' )

# for all pairs (i,j) where i < j and 0 <= i,j <= 2
for i in {0..2}
do
  # counter
  j=$(($i+1))

  while [ $j -lt 3 ]; do
    echo "Finding homology for" ${spec2[$i]} "and" ${spec2[$j]}

    # match up genes
    echo "Matching genes"
    join $anns/${genes[$i]} $anns/${genes[$j]} | awk -v a=${spec2[$i]} -v b=${spec2[$j]} '{OFS="\t"; print $1, a "_" $2, $3, $4, $5, b "_" $6, $7, $8, $9}' > $anns/${spec2[$i]}_${spec2[$j]}_gene_homology.txt

    # map ends of homologous genes to bins, then tally how many are in each bin
    echo "Counting gene ends"
    awk '{OFS="\t"; if ($5==$9) {print $2, $3, $6, $7; print $2, $4, $6, $8} else {print $2, $3, $6, $8; print $2, $4, $6, $7}}' $anns/${spec2[$i]}_${spec2[$j]}_gene_homology.txt | \
      awk -v b=$bsize '{OFS="\t"; print $1, int($2/b), $3, int($4/b)}' | sort -k1,1 -k2,2n -k3,3 -k4,4n | uniq -c | sort -k2,2 -k3,3n -k1,1nr > $anns/${abbr[$i]}${abbr[$j]}.$bsize.homology

    # find all and best homologous bins, excluding singletons and finding neighbors
    echo "Finding all and best homologous bins"
    ./homology_bins.sh $bsize ${abbr[$i]}${abbr[$j]} 1
    
    # increment counter
    let j+=1
  done
done

# S. CEREVISIAE INTRASPECIFIC HYBRIDS
# for within species homology, map simulated reads to each other

ref="Y12xDBVPG6044"
ref1="Y12"
ref2="DBVPG6044"
bsize=32000
win=10
rlen=150

# compile genome2reads
if [ ! -x genome2reads ]; then
  g++ genome2reads.cpp -o genome2reads -O3
fi

if [ ! -d $out/simreads ]; then
  mkdir $out/simreads
fi

# create simulated reads from Y12 genome
./genome2reads $refs/$ref1.fa $win $rlen 150 > $out/simreads/${ref1}_w${win}_r${rlen}.fa

# map simulated Y12 reads to DBVPG6044 reference
bowtie2 -f --very-sensitive -p 8 --reorder \
  -x $bt2/$ref2 \
  -U $out/simreads/${ref1}_w${win}_r${rlen}.fa \
  2> $out/simreads/${ref1}_w${win}_r${rlen}.${ref2}.out | \
  samtools view -b - > $out/simreads/${ref1}_w${win}_r${rlen}.$ref2.bam

# convert to bin units, tally counts of read mappings
samtools view $out/simreads/${ref1}_w${win}_r${rlen}.$ref2.bam | \
  awk '$6 >= 30' | \
  sed 's/-/      /' | \
  awk -v b=$bsize '{OFS="\t"; print $1, int($2/b), $4, int($5/b)}' | \
  sort -k1,1 -k2,2n -k3,3 -k4,4n | uniq -c | sort -k2,2 -k3,3n -k1,1nr > $anns/$ref.$bsize.homology

# find all and best homologous bins
./homology_bins.sh $bsize $ref 15
