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

Usage example:
```R

Rscript get_xclone_input.R [single_cells_cn] [unit] [phased_snps] [outdir] [mc.cores]

```

### Copy number profiles of each cell (or sample) [`single_cells_cn`]
A **txt file** listing the full path to the copy number profile of each cell (or each sample, i.e. bulk data).

Each **copy number profile** is a **tab-delimited** file reporting coordinates of genomic segments and their copy number status.<br/>
It has **4 columns** corresponding to chromosome, start position , end position and copy number status. **No header**.<br/>
> Take a look at the standard [BED](http://genome.ucsc.edu/FAQ/FAQformat#format1) format for more details.

*Example of **copy number profile** file:*
```
1	144010127	144019970	2
1	144920087	145739956	3
1	146540101	147829895	1
1	149815060	187128297	1
1	187128298	206257984	2
```
### The genomic unit [`unit`]
A **tab-delimited** file reporting coordinates of genomic segments (the genomic unit, i.e. genes).<br/>
It has **3 columns** corresponding to chromosome, start position , end position. **No header**.

### SNP with phasing information [`phased_snps`]
A **VCF** file listing SNPs of interest annotated with phasing information specified as `0|1` or `1|0`.
> Take a look at the standard [VCF]((https://samtools.github.io/hts-specs/VCFv4.2.pdf)) format for more details.

*Example of **phased_snps** file:*
```
5	180708328	.	T	G	1539.6	.	AN=2;AC=1	GT	0|1
5	180708937	.	T	A	1098.6	.	AN=2;AC=1	GT	0|1
5	180710261	.	T	C	1518.6	.	AN=2;AC=1	GT	0|1
5	180711595	.	G	A	1534.6	.	AN=2;AC=1	GT	0|1
5	180712169	.	C	T	1463.6	.	AN=2;AC=1	GT	0|1
```

### Output folder [`outdir`]
The full-path to the folder where all outputs will be saved.

### Number jobs to run in parallel [`mc.cores`] 
The script can be run on multiple samples in parallel.
Note:
```
mc.cores ≤ wc -l single_cells_cn
```
## Outputs

# Get summary stats about XClone inputs
