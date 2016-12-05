#!/bin/bash
# masks.sh
# Seungsoo Kim

echo "Sbrep"
./create_Sbrep_mask.sh

echo "Spgap"
./create_Spgap_mask.sh

echo "Y12xDBVPG6044"
./create_Y12xDBVPG6044_mask.sh

echo "adding rDNA"
cat ../masks/Sbay_repeats.bed ../masks/rDNA.bed > ../masks/Sbrep-rDNA.bed
cat ../masks/Spgap.bed ../masks/rDNA.bed > ../masks/Spgap-rDNA.bed

cat ../masks/Sbrep-rDNA.bed ../masks/Scer_ymr290-291.bed > ../masks/Sbrep-rDNA-ScYMR290-291.bed
