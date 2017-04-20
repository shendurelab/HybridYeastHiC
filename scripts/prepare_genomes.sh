# prepare_genomes.sh
# Seungsoo Kim

refs="../references"
anns="../nobackup/annotations"
bt2="../nobackup/bowtie2"

# load modules
module load samtools/1.2 bedtools/2.24.0 bowtie2/2.2.3

# rename/reorder chromosomes
sed 's/r0/r/' $refs/saccer3.fsa | sed 's/chr/Scer_/' > $refs/Scer.fa
sed 's/ .*//' $refs/SbayDRS.ultrascaf > $refs/Sbay.fa
sed 's/ .*//' $refs/Spar.ultrascaf > $refs/Spar.fa
sed 's/Sbay/Suva/' $refs/Sbay.fa > $refs/Suva.fa

Y12="$refs/Y12_FALCON_quiver.fa"
paste <(cat <(seq 1 16 | sed 's/^/>Y12_/')) <(cat <(sed '34q;d' $Y12) \
    <(sed '12q;d' $Y12) \
    <(sed '30q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '2q;d' $Y12) \
    <(sed '20q;d' $Y12) \
    <(sed '32q;d' $Y12) \
    <(sed '4q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '24q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '28q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '14q;d' $Y12) \
    <(sed '18q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(paste <(sed '26q;d' $Y12) <(sed '22q;d' $Y12) | sed 's/   //') \
    <(sed '10q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(paste <(sed '36q;d' $Y12) <(sed '16q;d' $Y12) | sed 's/   //') \
    <(sed '6q;d' $Y12) \
    <(sed '8q;d' $Y12 | tr 'AaTtCcGg' 'TtAaGgCc' | rev)) | tr '	' '\n' > $refs/Y12.fa

sed 's/quiver/quiver	/' $refs/DBVPG6044_FALCON_quiver.fa | tr '\n' 'X' | sed "s/>/\n>/g" | sort -k1,1n | tr '\t' '\n' | sed 's/X//g' | tail -n +2 > $refs/DBVPG6044_FALCON_quiver_sorted.fa
DBVPG6044="$refs/DBVPG6044_FALCON_quiver_sorted.fa"
paste <(cat <(seq 1 16 | sed 's/^/>DBVPG6044_/')) <(cat <(sed '34q;d' $DBVPG6044) \
    <(sed '12q;d' $DBVPG6044) \
    <(sed '30q;d' $DBVPG6044) \
    <(sed '2q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '22q;d' $DBVPG6044) \
    <(sed '32q;d' $DBVPG6044) \
    <(sed '4q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '24q;d' $DBVPG6044) \
    <(sed '28q;d' $DBVPG6044) \
    <(sed '16q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '18q;d' $DBVPG6044) \
    <(paste <(sed '26q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) <(sed '20q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) | sed 's/	//') \
    <(sed '10q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '14q;d' $DBVPG6044 | tr 'AaTtCcGg' 'TtAaGgCc' | rev) \
    <(sed '6q;d' $DBVPG6044) \
    <(sed '8q;d' $DBVPG6044)) | tr '	' '\n' > $refs/DBVPG6044.fa

# rename Sbay to Suva in S. uvarum gff file
sed 's/Sbay/Suva/g' $refs/Sbay.gff > $refs/Suva.gff

# revise S. paradoxus genome
./revise_Spar.sh

# revise S. uvarum genome
./revise_Suva.sh

# rename chromosomes in S. cerevisiae gff file
cat $refs/saccharomyces_cerevisiae_R64-2-1_20150113.gff | \
  sed 's/chrXVI/16/' | \
  sed 's/chrXV/15/' | \
  sed 's/chrXIV/14/' | \
  sed 's/chrXIII/13/' | \
  sed 's/chrXII/12/' | \
  sed 's/chrXI/11/' | \
  sed 's/chrX/10/' | \
  sed 's/chrIX/9/' | \
  sed 's/chrVIII/8/' | \
  sed 's/chrVII/7/' | \
  sed 's/chrVI/6/' | \
  sed 's/chrV/5/' | \
  sed 's/chrIV/4/' | \
  sed 's/chrIII/3/' | \
  sed 's/chrII/2/' | \
  sed 's/chrI/1/' | \
  sed 's/chrmt/MT/' > $refs/Scer.gff

# index references
samtools faidx $refs/Scer.fa
samtools faidx $refs/Spar_revised.fa
samtools faidx $refs/Suva_revised.fa
samtools faidx $refs/Y12.fa
samtools faidx $refs/DBVPG6044.fa

# create hybrid references
cat $refs/Scer.fa <(echo '\n') $refs/Suva_revised.fa $refs/2micron.fa > $refs/ScSu_2micron_revised.fa
cat $refs/Scer.fa <(echo '\n') $refs/Spar_revised.fa $refs/2micron.fa > $refs/ScSp_2micron_revised.fa
cat $refs/Spar_revised.fa $refs/Suva_revised.fa $refs/2micron.fa > $refs/SpSu_2micron_revised.fa
cat $refs/Y12.fa $refs/DBVPG6044.fa > $refs/Y12xDBVPG6044.fa
