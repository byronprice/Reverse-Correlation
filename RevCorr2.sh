#!/bin/bash -l

# 72 hour time limit
#$ -l h_rt=72:00:00
#$ -N RevCorr2
#$ -j y
#$ -o RevCorr2Outputs.txt
#$ -l mem_per_core=16G
#$ -pe omp 16

module load matlab/2017a

matlab -nodisplay -r "RevCorrMovies(43652,20170421,'pink');exit"
