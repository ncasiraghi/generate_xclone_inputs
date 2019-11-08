#!/bin/sh 
#BSUB -J xclone_input
#BSUB -n 20
#BSUB -q long
#BSUB -e xclone_input.log 
#BSUB -o xclone_input.txt 

module load R/3.6.0 bedtools/2.24.0

single_cells_cn="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/single_cells_cn.txt"
genes="/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/data/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.bed"
phased_snps="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/eagle2_phasing/rawdata/EAGLE2_20190909_vcfs/all.vcf"
outdir="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_SP/unit_genes"
cores=20

Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/get_xclone_input_SP_mode.R $single_cells_cn $genes $phased_snps $outdir $cores
