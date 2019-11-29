# Generate inputs for the XClone pipeline

The R script `get_xclone_input.R` intersects (`intersectBed`) assign copy number status to a genomic unit defined by the user (i.e. Genes, Haplotype Blocks, ATAC peaks, custom bin). Each phased SNPs is then assigned to each genomic unit.<br/>
Each entry reported in the output corresponds to a SNP annotated with phasing information, coordiantes and copy number status of the genomic unit where it is located. 

## Requirements

If you are working on the DKFZ `odcf server` you can load the needed modules (`R`, `bedtools`, `bedops`) with:
```
module load R/3.6.0  
module load bedtools/2.24.0
module load bedops/2.4.14
```
> A `bsub` example script is saved in this repository.  

## Inputs

```R
# Usage example
Rscript get_xclone_input.R [single_cells_cn] [unit] [phased_snps] [outdir] [mc.cores]
```

### Copy number profiles of each cell (or sample) [`single_cells_cn`]
A **txt file** listing the full path to the copy number profile of each cell (or each sample, i.e. bulk data).

Each **copy number profile** is a **tab-delimited** file reporting coordinates of genomic segments and their copy number status. It has **4 columns** corresponding to chromosome, start position , end position and copy number status. **No header**.<br/>
Take a look at the standard [BED](http://genome.ucsc.edu/FAQ/FAQformat#format1) format for more details).

*Example of **copy number profile** file format:*
```
1	144010127	144019970	2
1	144920087	145739956	3
1	146540101	147829895	1
1	149815060	187128297	1
1	187128298	206257984	2
```
### The genomic unit [`unit`]

### SNP with phasing information [`phased_snps`]

### Output folder [`outdir`]

### Number jobs to run in parallel [`mc.cores`] 

## Outputs

# Get summary stats about XClone inputs
