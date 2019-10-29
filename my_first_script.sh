#!/bin/bash
#
#SBATCH -p ncf_holy # partition (queue)
#SBATCH --mem 100 # memory 
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o myscript_%j_output.out # STDOUT
#SBATCH --mail-type=END # notifications for job done

for i in {1..100000}; do
    echo $RANDOM >> SomeRandomNumbers.txt
done

sort SomeRandomNumbers.txt
