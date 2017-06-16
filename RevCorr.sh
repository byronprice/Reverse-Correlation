#!/bin/bash -l

# 24 hour time limit
#$ -l h_rt=12:00:00
#$ -N RevCorr
#$ -j y
#$ -o RevCorrOutput.txt
#$ -l mem_per_core=16G
#$ -pe omp 8
#$ -m ea
#$ -M byron.h.price@gmail.com


module load matlab/2017a

matlab -nodisplay -r "BatchRunSCC;exit"
