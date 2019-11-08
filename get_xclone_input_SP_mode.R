args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=5){
  message("\n[USAGE] Rscript get_xclone_input_SP_mode.R [single_cells_cn] [genes] [phased_snps] [outdir] [mc.cores]\n")
  message("[1] single_cells_cn \tA txt file listing the full path to each copy number profiles for each cell or each sample (bulk data); 4 cols: chr, start, end, cn. No header.")
  message("[2] genes           \tBED-like file with genomic coordinates of genes; 3 cols: chr, start, end, HUGO_name, ENSG_id. No header. [optional, write NA to ignore this param, output will be based on copy number blocks]")
  message("[3] phased_snps     \tVCF file with only phased SNPs; VCF format. This is the output from EAGLE2")
  message("[4] outdir          \tFull path to the output folder.")
  message("[5] mc.cores        \tNumber of cells (or samples) to run in parallel.\n")
  quit()
}

single_cells_cn <- readLines(args[1])
genes <- args[2]
phased_snps <- args[3]
outdir <- args[4]
mc.cores <- as.numeric(args[5])

library( tidyr )
library( data.table )
library( parallel )

setwd( outdir )

if(genes=="NA"){
  message("running in mode: Statistical Phasing & unit: cn_blocks")
  ## unit: cn_blocks 
  quit()
  
} else {
  ## unit: genes
  GetAnnotedSNPs <- function(i,single_cells_cn,genes,phased_snps){
    bed_file <- single_cells_cn[i]
    CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
    
    # step 4 : intersect haplotype blocks with step3 
    step4 <- paste0('step4_',CELL_ID,'.bed')
    cmd <- paste('intersectBed -a',genes,'-b',bed_file,'-wa -wb >',step4)
    system(cmd)
    
    m <- fread(file = step4,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:9))
    m[,9] <- round(m[,9])
    m <- unite(m,col = geneblock,seq(1,3),sep = ":",remove=FALSE)
    
    genes_to_exclude <- c(which(duplicated(m$geneblock,fromLast = FALSE)),which(duplicated(m$geneblock,fromLast = TRUE)))
    
    message(paste("removing",length(unique(m[genes_to_exclude,1])),"genes out of",length(unique(m[,1]))))
    
    if(length(genes_to_exclude)>0){
      tab <- m[-genes_to_exclude,]
    } else {
      tab <- m
    }
    
    # step 5 : write genes with associated copy number
    step5 <- paste0('step5_',CELL_ID,'.bed')
    write.table(tab[,c(2:6,10)],file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    # step 7 : add phased SNPs to each gene
    step7 <- paste0('step7_',CELL_ID,'.bed')
    cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
    system(cmd)
    file.remove(step4)
    
    # step 8 : reduce table columns and write the final output
    step8 <- paste0('step8_',CELL_ID,'.bed')
    
    m <- fread(file = step7,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:6,7,8,10,11,16))
    header <- c("GENE_CHROM","GENE_START","GENE_END","GENE_HUGO","GENE_ID","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
    colnames(m) <- header
    
    m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","GENE_CHROM","GENE_START","GENE_END","GENE_HUGO","GENE_ID","COPY_NUMBER")]
    write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    file.remove(step7)
    
  }
  
  message("running in mode: Statistical Phasing & unit: genes")
  mclapply(seq(single_cells_cn),GetAnnotedSNPs,single_cells_cn=single_cells_cn,genes=genes,phased_snps=phased_snps,mc.cores = mc.cores)

  # keep only gene that are present in all cells
  ks <- list.files(path = file.path(outdir),pattern = 'step5_',full.names = TRUE)
  
  a <- fread(file = ks[1],stringsAsFactors = FALSE,header = FALSE, data.table = FALSE)
  a <- unite(a,col = group,seq(1,3),sep = ":",remove=FALSE)
  a <- a[,c(1,7)]
  colnames(a)[2] <- basename(ks[1])
  
  for(k in seq(2,length(ks))){
    b <- fread(file = ks[k],stringsAsFactors = FALSE,header = FALSE, data.table = FALSE)
    b <- unite(b,col = group,seq(1,3),sep = ":",remove=FALSE)
    b <- b[,c(1,7)]
    colnames(b)[2] <- basename(ks[k])
    a <- merge(x = a,y = b,by='group',all = TRUE)
  }
  
  row.has.na <- apply(a, 1, function(x){any(is.na(x))})
  message(paste("excluding",sum(row.has.na),"genes out of",nrow(a)))
  
  genes_to_keep <- a[!row.has.na,]
  
  # step9 : filtering step8
  ks <- list.files(path = file.path(outdir),pattern = 'step8_',full.names = TRUE)
  
  for(k in seq_len(length(ks))){
    message(paste("filtering:",basename(ks[k])))
  
    # step 9
    step9 <- gsub(basename(ks[k]),pattern = "step8_",replacement = "step9_")
    
    m <- fread(file = ks[k],stringsAsFactors = FALSE,header = TRUE, data.table = FALSE)
    m <- unite(m,col = group,seq(6,8),sep = ":",remove=FALSE)
    m <- m[which(m$group %in% unique( genes_to_keep$group )),]
    
    write.table(m[,-6],file = step9,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  }
  message("done.")
}













# GetAnnotedSNPs <- function(i,single_cells_cn,haplotype_blocks,phased_snps){
#   bed_file <- single_cells_cn[i]
#   message(basename(bed_file))
#   CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
#   
#   # step 4 : intersect haplotype blocks with step3 
#   step4 <- paste0('step4_',CELL_ID,'.bed')
#   cmd <- paste('intersectBed -a',haplotype_blocks,'-b',bed_file,'-wa -wb >',step4)
#   system(cmd)
#   
#   m <- fread(file = step4,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:7))
#   m[,7] <- round(m[,7])
#   m <- unite(m,col = haploblock,seq(1,3),sep = ":",remove=FALSE)
#   
#   hb_to_exclude <- c(which(duplicated(m$haploblock,fromLast = FALSE)),which(duplicated(m$haploblock,fromLast = TRUE)))
#   
#   message(paste("removing",length(unique(m[hb_to_exclude,1])),"haplotype blocks out of",length(unique(m[,1]))))
#   
#   if(length(hb_to_exclude)>0){
#     tab <- m[-hb_to_exclude,]
#   } else {
#     tab <- m
#   }
#   
#   # step 5 : write haplo type blocks with associated copy number
#   step5 <- paste0('step5_',CELL_ID,'.bed')
#   write.table(tab[,c(2:4,8)],file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
#   
#   # step 7 : add phased SNPs to each block
#   step7 <- paste0('step7_',CELL_ID,'.bed')
#   cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
#   system(cmd)
#   file.remove(step4,step5)
#   
#   # step 8 : reduce table columns and write the final output
#   step8 <- paste0('step8_',CELL_ID,'.bed')
# 
#   m <- fread(file = step7,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:6,8,9,14))
#   header <- c("HB_CHROM","HB_START","HB_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
#   colnames(m) <- header
#   
#   reduce_info <- function(info){
#     return(unlist(strsplit(info,split = ":"))[1])
#   }
#   
#   m$SNP_PHASE_INFO <- unlist(lapply(X = m$SNP_PHASE_INFO,FUN = reduce_info))
#   
#   m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","HB_CHROM","HB_START","HB_END","COPY_NUMBER")]
#   write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
#   
#   file.remove(step7)
#   
# }
