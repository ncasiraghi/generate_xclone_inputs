args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript get_xclone_input.R get_xclone_input_configure.R\n")
  quit()
}

configure_file <- args[1]
configure_file <- "get_xclone_input_configure.R"

library( data.table )
library( parallel )

source(configure_file)

## RUN
setwd(outdir)

GetAnnotedSNPs <- function(i,single_cells_cn,haplotype_blocks,phased_snps){
  bed_file <- single_cells_cn[i]
  message(basename(bed_file))
  CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
  
  # step 4 : intersect haplotype blocks with step3 
  step4 <- paste0('step4_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',haplotype_blocks,'-b',bed_file,'-wa -wb >',step4)
  system(cmd)
  
  m <- read.delim(file = step4,header = FALSE,as.is = TRUE,stringsAsFactors = FALSE)
  m <- m[,1:7]
  m$haploblock <- paste(m[,1],m[,2],m[,3],sep = ":")
  
  hb <- unique(m$haploblock)
  
  hb_to_exclude <- c()
  tab <- c()
  for(x in hb){
    check <- length(unique(m[which(m$haploblock == x),7]))
    if( check > 1 ){
      hb_to_exclude <- c(hb_to_exclude, x)
    } else {
      tab <- rbind(tab, unique(m[which(m$haploblock == x),c(1:3,7)]))
    }
  }
  
  length(hb_to_exclude)/length(hb) # fraction of hb with multiple spanning CNs 
  
  # step 5 : write haplo type blocks with associated copy number
  step5 <- paste0('step5_',CELL_ID,'.bed')
  write.table(tab,file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  # step 7 : add phased SNPs to each block
  step7 <- paste0('step7_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
  system(cmd)
  file.remove(step4,step5)
  
  # step 8 : reduce table columns and wrinte final output
  step8 <- paste0('step8_',CELL_ID,'.bed')
  cmd <- paste('cut -f 1-6,8,9,14',step7,'>',step8)
  system(cmd)
  file.remove(step7)
  
  m <- fread(file = step8,data.table = FALSE)
  header <- c("HB_CHROM","HB_START","HB_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
  colnames(m) <- header
  
  reduce_info <- function(info){
    return(unlist(strsplit(info,split = ":"))[1])
  }
  
  m$SNP_PHASE_INFO <- unlist(lapply(X = m$SNP_PHASE_INFO,FUN = reduce_info))
  
  m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","HB_CHROM","HB_START","HB_END","COPY_NUMBER")]
  write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
}

mclapply(seq(single_cells_cn),GetAnnotedSNPs,single_cells_cn=single_cells_cn,haplotype_blocks=haplotype_blocks,phased_snps=phased_snps,mc.cores = mc.cores)