#!/bin/sh 
#BSUB -J xclone_input
#BSUB -n 20
#BSUB -q long
#BSUB -e xclone_input.log 
#BSUB -o xclone_input.txt 

module load R/3.6.0 bedtools/2.24.0 bedops/2.4.14

single_cells_cn="/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/single_cells_cn.txt"
unit="/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/mode_SP/unit_w10Mb/sliding.windows.10Mb.nochr.bed"
phased_snps="/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/all.vcf"
outdir="/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/mode_SP/unit_w10Mb"
cores=20

Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/get_xclone_input.R $single_cells_cn $unit $phased_snps $outdir $cores
