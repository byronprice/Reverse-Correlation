#!/bin/bash -l

# 300 hour time limit
#$ -l h_rt=300:00:00
#$ -N BayesReceptiveFields
#$ -j y
#$ -o BayesReceptiveFields.txt
#$ -l mem_per_core=8G
#$ -pe omp 8
#$ -m ea
#$ -M byron.h.price@gmail.com

module load matlab/2017a

matlab -nodisplay -r "for ii=0:2;[~] = BayesianFitSRFQuick(81071,20170620+ii,'pink');end;exit"