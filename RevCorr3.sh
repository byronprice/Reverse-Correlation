#!/bin/bash -l

# 24 hour time limit
#$ -l h_rt=50:00:00
#$ -N RevCorr3
#$ -j y
#$ -o RevCorr3Output.txt
#$ -l mem_per_core=16G
#$ -pe omp 16

module load matlab/2017a

matlab -nodisplay -r "RevCorrMov3(43652,20170421,'pink');exit"
