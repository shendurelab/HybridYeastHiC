#!/bin/bash

refs="../references"
grep repeat $refs/Sbay_revised.gff | awk '{OFS="\t"; print "Sbay_" $1, $4, $5}' > ../masks/Sbay_repeats.bed
