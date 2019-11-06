# Require: module load bedtools/2.24.0

# copy number profiles for each cell, or one single file for bulk data (chr, start, end, cn). No header. 
single_cells_cn <- readLines("/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/single_cells_cn.txt")

# haplotype blocks (chr, start, end). No header.
haplotype_blocks <- "/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/haploblocks.bed"

# List of (only) phased SNPs, VCF format.
phased_snps <- "/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K/all.vcf"

# output folder
outdir <- "/icgc/dkfzlsdf/analysis/B260/users/n790i/hipo_K08K"

# number of cells (or samples) in parallel
mc.cores <- 2
