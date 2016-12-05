#!/bin/bash
# pipeline.sh
# runs all scripts needed to analyze Hi-C data
# Seungsoo Kim

# output directory
out='../nobackup'
if [ ! -d $out ]; then
  mkdir $out
fi
if [ ! -d $out/sge_out ]; then
  mkdir $out/sge_out
fi
if [ ! -d $out/sge_error ]; then
  mkdir $out/sge_error
fi


if [[ "$1" =~ ^((-{1,2})([Hh]$|[Hh][Ee][Ll][Pp])|)$ ]]; then
  echo "USAGE: ./pipeline.sh [-r/--refs REFERENCE] [-t/--trim RAW_FILES] [-p/--process SAMPLES] [-m/--matrices MATRICES]"; exit 1
else
  while [[ $# -gt 0 ]]; do
    opt="$1"
    shift;
    current_arg="$1"
    if [[ "$current_arg" =~ ^-{1,2}.* ]]; then
      echo "WARNING: You may have left an argument blank. Double check your command." 
    fi
    case "$opt" in
      "-r"|"--refs"      ) refs="$1"; shift;;
      "-t"|"--trim"      ) trim="$1"; shift;;
      "-p"|"--process"   ) samps="$1"; shift;;
      "-m"|"--matrices"  ) matrices="$1"; shift;;
      *                  ) echo "ERROR: Invalid option: \""$opt"\"" >&2
                            exit 1;;
    esac
  done
fi

# ARGUMENTS
# reference file - tab-delimited file with columns for:
# 1) reference name, 2) path to the FASTA file, 3) path to centromere annotation file, 4) comma-separated list of bin sizes

# raw reads file list - tab-delimited file with columns for:
# 1) sample name, 2) restriction enzyme name, 3) read 1 FASTQ file, 4) read 2 FASTQ file

# sample file - tab-delimited file with columns for:
# 1) sample name, 2) sample type, 3) reference name, 4) ref 1 name (needed for dual only), 5) ref 2 name (needed for dual only), 6) restriction enzyme name

# matrix file - tab-delimited file with columns for:
# 1) sample name, 2) reference name, 3) bin size, 4) mask bed file (exclude .bed here), 5) minimum average across row

# make output folders
for d in adapt_trim junc_trim aligned assigned matrix
do
  if [ ! -d "$out/$d" ]; then
    mkdir "$out/$d"
  fi
done

# build reference files
if [ "$refs" != "" ]; then
  echo "Building references ..."
  ./build_references.sh $refs
fi

# trim reads - adapters and to restriction junctions
if [ "$trim" != "" ]; then
  echo "Trimming reads ..."
  while read samp renz r1 r2
  do
    qsub -N trim_${samp} -o $out/sge_out/ -e $out/sge_error/ -pe serial 1-8 -l mfree=4G -cwd -b y ./trim_reads.sh $samp $renz $r1 $r2
  done < $trim
fi

# process reads - trim, map, deduplicate, and assign to restriction fragments, then identify valid ligations
if [ "$samps" != "" ]; then
  echo "Processing reads ..."
  while read samp type ref ref1 ref2 renz r1 r2
  do
    qsub -N process_${samp}_${ref} -o $out/sge_out/ -e $out/sge_error/ -pe serial 1-8 -l mfree=4G -cwd -b y ./process_reads.sh $samp $type $ref $ref1 $ref2 $renz
  done < $samps
fi
# compile make_matrix.cpp
if [ ! -x make_matrix ]; then
  g++ make_matrix.cpp -o make_matrix -O3
fi

# make matrices
if [ "$matrices" != "" ]; then
  echo "Building matrices ..."
  while read samp ref bsize mask minavg
  do
    qsub -N matrix_${samp} -o $out/sge_out/ -e $out/sge_error/ -hold_jid process_${samp}_${ref} -cwd -b y ./make_matrix.sh $samp $ref $bsize $mask $minavg
  done < $matrices
fi
