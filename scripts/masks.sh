#!/bin/bash
# masks.sh
# Seungsoo Kim

echo "Surep"
./create_Surep_mask.sh

echo "Spgap"
./create_Spgap_mask.sh

echo "Y12xDBVPG6044"
./create_Y12xDBVPG6044_mask.sh

echo "adding rDNA"
cat ../masks/Suva_repeats.bed ../masks/rDNA.bed > ../masks/Surep-rDNA.bed
cat ../masks/Spgap.bed ../masks/rDNA.bed > ../masks/Spgap-rDNA.bed

cat ../masks/Surep-rDNA.bed ../masks/Scer_ymr290-291.bed > ../masks/Surep-rDNA-ScYMR290-291.bed
