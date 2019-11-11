#!/bin/sh 
#BSUB -J xclone_input
#BSUB -n 20
#BSUB -q long
#BSUB -e xclone_input.log 
#BSUB -o xclone_input.txt 

module load R/3.6.0 bedtools/2.24.0

single_cells_cn="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/single_cells_cn.txt"
haplotype_blocks="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.bed"
phased_snps="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf"
outdir="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_RP/unit_haplotype_blocks"
cores=20

Rscript /icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/get_xclone_input_RP_mode.R $single_cells_cn $haplotype_blocks $phased_snps $outdir $cores
