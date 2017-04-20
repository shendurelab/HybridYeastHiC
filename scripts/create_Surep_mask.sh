#!/bin/bash

refs="../references"
grep repeat $refs/Suva_revised.gff | awk '{OFS="\t"; print "Suva_" $1, $4, $5}' > ../masks/Suva_repeats.bed
