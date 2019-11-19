#!/bin/sh 
#BSUB -J xclone_input
#BSUB -n 30
#BSUB -q long
#BSUB -e xclone_input.log 
#BSUB -o xclone_input.txt 

module load R/3.6.0

Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/studio_xclone_input.R
