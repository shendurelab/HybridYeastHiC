#!/bin/bash

module load bowtie2/2.2.3 samtools/1.2

refdir="../references"
bt2="../nobackup/bowtie2"
sim="../nobackup/simreads"
out="../nobackup/mappability"
refs=( 'Y12xDBVPG6044_quiver' 'ScSu_2micron_revised' 'ScSp_2micron_revised' 'SpSu_2micron_revised')
win=10
rlens=( 150 80 80 80)
bsizes=( 32000 )
qthrs=( 6 20 30 )

# for each reference/read-length combination
for i in {0..3}
do
  ref=${refs[$i]}
  rlen=${rlens[$i]}

  # simulate reads
  ./genome2reads $refdir/$ref.fa ${win} ${rlen} > $sim/${ref}_w${win}_r${rlen}.fa

  # align reads
  bowtie2 -f --very-sensitive -p 8 --reorder -x $bt2/$ref -U $sim/${ref}_w${win}_r${rlen}.fa 2> $sim/${ref}_w${win}_r${rlen}.out | samtools view -b - \
    > $sim/${ref}_w${win}_r${rlen}.bam

  # for each MAPQ threshold:
  for qthr in "${qthrs[@]}"
  do
    # call reads as mapped or unmapped (incorrectly mapped)
    samtools view $sim/${ref}_w${win}_r${rlen}.bam | grep -v ^@ | sed 's/-/	/' | \
      awk -v q=$qthr '{OFS="\t"; if (($6 < q) || ($1 != $4) || ($2 != $5)) 
                                   print $1, $2, "unmapped"; 
                                 else print $1, $2, "mapped"}' > $out/${ref}_w${win}_r${rlen}_q${qthr}.mapping

    awk '{if ($3 == "unmapped") a = a + 1; else b = b + 1} END {print b/(a+b)}' $out/${ref}_w${win}_r${rlen}_q${qthr}.mapping > $out/${ref}_w${win}_r${rlen}_q${qthr}.map_rate
:<<END
    # tally percentage of reads mapped
    for bsize in "${bsizes[@]}"
    do
      awk -v b=$bsize 'BEGIN{OFS="\t"; chr=""; bin=""; unmapped=0; mapped=0}
                       {if ((NR!=1) && (chr != $1 || bin != int($2/b))) 
                          {print chr, bin, unmapped/(unmapped+mapped); unmapped=0; mapped=0}
                        else if ($3=="unmapped") 
                          unmapped = unmapped + 1;
                        else 
                          mapped = mapped + 1}
                       END{print chr, bin, unmapped/(unmapped+mapped)}' $out/${ref}_w${win}_r${rlen}_q${qthr}.mapping | \
        > $out/${ref}_w${win}_r${rlen}_q${qthr}.$bsize.mappability
    done
END
  done
done
