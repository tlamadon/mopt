#!/bin/bash

echo "starting qsub script file"
source ~/.bash_profile
date 

RSNOW=`Rscript -e "cat(.libPaths()[[1]])"`

echo $RSNOW
# here's the SGE directives
# ------------------------------------------
#$ -q batch.q   # <- the name of the Q you want to submit to
#$ -pe mpich 100 #  <- load the openmpi parallel env w/ 3 slots
#$ -S /bin/bash   # <- run the job under bash
#$ -N mopt-example # <- name of the job in the qstat output
#$ -o std.out # <- name of the output file.
#$ -e err.out # <- name of the stderr file.
#$ -cwd

module load openmpi
module load open64
module load gcc
module load r/2.15.2

echo "loaded modules"
module list

echo "calling mpirun now"
mpirun -np 3 $RSNOW/snow/RMPISNOW -q < example.bgp.mpi.r > example.Rout


