args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  warning("\n[USAGE] Rscript get_xclone_input.R get_xclone_input_configure.R\n")
  quit()
}

configure_file <- args[1]

library( tidyr )
library( data.table )
library( parallel )

source( configure_file )

setwd( outdir )

GetAnnotedSNPs <- function(i,single_cells_cn,haplotype_blocks,phased_snps){
  bed_file <- single_cells_cn[i]
  message(basename(bed_file))
  CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
  
  # step 4 : intersect haplotype blocks with step3 
  step4 <- paste0('step4_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',haplotype_blocks,'-b',bed_file,'-wa -wb >',step4)
  system(cmd)
  
  m <- fread(file = step4,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:7))
  m[,7] <- round(m[,7])
  m <- unite(m,col = haploblock,seq(1,3),sep = ":",remove=FALSE)
  
  hb_to_exclude <- c(which(duplicated(m$haploblock,fromLast = FALSE)),which(duplicated(m$haploblock,fromLast = TRUE)))
  
  message(paste("removing",length(unique(m[hb_to_exclude,1])),"haplotype blocks out of",length(unique(m[,1]))))
  
  if(length(hb_to_exclude)>0){
    tab <- m[-hb_to_exclude,]
  } else {
    tab <- m
  }
  
  # step 5 : write haplo type blocks with associated copy number
  step5 <- paste0('step5_',CELL_ID,'.bed')
  write.table(tab[,c(2:4,8)],file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  # step 7 : add phased SNPs to each block
  step7 <- paste0('step7_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
  system(cmd)
  file.remove(step4,step5)
  
  # step 8 : reduce table columns and write the final output
  step8 <- paste0('step8_',CELL_ID,'.bed')

  m <- fread(file = step7,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:6,8,9,14))
  header <- c("HB_CHROM","HB_START","HB_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
  colnames(m) <- header
  
  reduce_info <- function(info){
    return(unlist(strsplit(info,split = ":"))[1])
  }
  
  m$SNP_PHASE_INFO <- unlist(lapply(X = m$SNP_PHASE_INFO,FUN = reduce_info))
  
  m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","HB_CHROM","HB_START","HB_END","COPY_NUMBER")]
  write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  file.remove(step7)
  
}

mclapply(seq(single_cells_cn),GetAnnotedSNPs,single_cells_cn=single_cells_cn,haplotype_blocks=haplotype_blocks,phased_snps=phased_snps,mc.cores = mc.cores)

message("alright, alright, alright.")
