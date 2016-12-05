#!/bin/bash
# run.sh
#
# Shell script to run structure simulation for a single iteration and process the output 
# PDB to count all pairs of beads within 45 nm as a "Hi-C interaction"

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load boost/1.52.0 python/2.7.2 imp/2.0.1 numpy/1.8.0

# run structure simulation
python diploid_landmark.py ${SGE_TASK_ID} $1 $2 $3 $4 > $1/${SGE_TASK_ID}.out

# process PDB output
./pdb2interactions 45 $1/${SGE_TASK_ID}.pdb > $1/${SGE_TASK_ID}.c45.interactions.txt
