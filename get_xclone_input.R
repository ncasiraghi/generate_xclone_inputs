args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=5){
  message("\n[USAGE] Rscript get_xclone_input.R [single_cells_cn] [unit] [phased_snps] [outdir] [mc.cores]\n")
  message("[1] single_cells_cn \tA txt file listing the full path to each copy number profiles for each cell or each sample (bulk data); 4 cols: chr, start, end, cn. No header.")
  message("                    \tif unit = 'NA' (unit_cn_blocks) required 5 cols: chr, start, end, cn, cell/sample_id. No header.")
  message("[2] unit            \tBED-like file with genomic coordinates of genomic segment of interest; 3 cols: chr, start, end. No header. [optional: write 'NA' to ignore this param, output will be based on copy number blocks]")
  message("[3] phased_snps     \tVCF file with only phased SNPs; VCF format. This is the output from EAGLE2")
  message("[4] outdir          \tFull path to the output folder.")
  message("[5] mc.cores        \tNumber of cells (or samples) to run in parallel.\n")
  quit()
}

single_cells_cn <- readLines(args[1])
unit <- args[2]
phased_snps <- args[3]
outdir <- args[4]
mc.cores <- as.numeric(args[5])

library( tidyr )
library( data.table )
library( parallel )

setwd( outdir )

if(unit=="NA"){
  ## unit: cn_blocks 
  message("running in mode: Statistical Phasing & unit: cn_blocks")
  
  # step 1 : intersect all CN profiles from cells belonging to the same cluster
  step1 <- paste0('step1.bed')
  multiple_beds <- paste(single_cells_cn,collapse = ' ')
  cmd <- paste('bedops --partition',multiple_beds,'>',step1)
  system(cmd)
  
  # step 2 : intersect step1 with all BEDs to assign CN to each genomic segment 
  step2 <- paste0('step2.bed')
  cmd <- paste('intersectBed -a',step1,'-b',multiple_beds,'-wa -wb | cut -f1,2,3,8,9 >',step2)
  system(cmd)
  file.remove(step1)
  
  # step 3 : write consensus copy number blocks in each cell
  m <- fread(file = step2,header = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
  m[,4] <- round(m[,4])
  out <- split(m, f = m[,5] )
  
  ## make sure that all cells have same number of copy number blocks
  a <- out[[1]]
  a <- unite(a,col = group,seq(1,3),sep = ":",remove=FALSE)
  CELL_ID <- a$V5[1]
  a <- a[,c(1,5)]
  colnames(a)[2] <- CELL_ID
  
  for(k in seq(2,length(out))){
    b <- out[[k]]
    b <- unite(b,col = group,seq(1,3),sep = ":",remove=FALSE)
    CELL_ID <- b$V5[1]
    b <- b[,c(1,5)]
    colnames(b)[2] <- CELL_ID
    a <- merge(x = a,y = b,by='group',all = TRUE)
  }
  
  file.remove(step2)
  row.has.na <- apply(a, 1, function(x){any(is.na(x))})
  message(paste("excluding",sum(row.has.na),"copy number blocks out of",nrow(a)))
  
  cn_blocks_to_keep <- a[!row.has.na,]
  
  # step4 : filtering out
  for(k in seq(length(out))){
    m <- out[[k]]
    m <- unite(m,col = group,seq(1,3),sep = ":",remove=FALSE)
    CELL_ID <- m$V5[1]
    message(paste("filtering:",CELL_ID))
    m <- m[which(m$group %in% unique( cn_blocks_to_keep$group )),]
    
    step5 <- paste0("step5_",CELL_ID,'.bed')
    write.table(m[,2:5],file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    # step 8 : add phased SNPs to each copy number block
    step8 <- paste0('step8_',CELL_ID,'.bed')
    cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step8)
    system(cmd)
    
    # step 9 : reduce table columns and write the final output
    step9 <- paste0('xci_',CELL_ID,'.bed')
    
    m <- fread(file = step8,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:4,5,6,8,9,14))
    header <- c("CNB_CHROM","CNB_START","CNB_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
    colnames(m) <- header
    
    m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","CNB_CHROM","CNB_START","CNB_END","COPY_NUMBER")]
    write.table(m,file = step9,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    file.remove(step5)
    file.remove(step8)
  }
  message("done.")
  
} else {
  
  message(paste("running in mode: Statistical Phasing & unit based on:",unit))
  
  format.unit  <- fread(unit,data.table = FALSE,select = 1:3)
  checked.unit <- file.path(outdir,'unit.bed')
  write.table(format.unit,file = checked.unit,quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t')
  
  GetAnnotedSNPs <- function(i,single_cells_cn,unit,phased_snps){
    bed_file <- single_cells_cn[i]
    CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
    
    # step 4 : intersect haplotype blocks with step3 
    step4 <- paste0('step4_',CELL_ID,'.bed')
    cmd <- paste('intersectBed -a',unit,'-b',bed_file,'-wa -wb >',step4)
    system(cmd)
    
    m <- fread(file = step4,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:7))
    if(is.numeric(m[,7])){
      m[,7] <- round(m[,7])
    }
    m <- unite(m,col = UNIT,seq(1,3),sep = ":",remove=FALSE)
    
    unit_to_exclude <- c(which(duplicated(m$UNIT,fromLast = FALSE)),which(duplicated(m$UNIT,fromLast = TRUE)))
    
    message(paste("removing",length(unique(m[unit_to_exclude,1])),"units out of",length(unique(m[,1]))))
    
    if(length(unit_to_exclude)>0){
      tab <- m[-unit_to_exclude,]
    } else {
      tab <- m
    }
    
    # step 5 : write unit with associated copy number
    step5 <- paste0('step5_',CELL_ID,'.bed')
    write.table(tab[,c(2:4,8)],file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    # step 7 : add phased SNPs to each gene
    step7 <- paste0('step7_',CELL_ID,'.bed')
    cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
    system(cmd)
    file.remove(step4)
    
    # step 8 : reduce table columns and write the final output
    step8 <- paste0('step8_',CELL_ID,'.bed')
    
    m <- fread(file = step7,sep = "\t",header = FALSE, stringsAsFactors = FALSE, data.table = FALSE,select = c(1:6,8,9,14))
    header <- c("UNIT_CHROM","UNIT_START","UNIT_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
    colnames(m) <- header
    
    reduce_info <- function(info){
      return(unlist(strsplit(info,split = ":"))[1])
    }
    
    m$SNP_PHASE_INFO <- unlist(lapply(X = m$SNP_PHASE_INFO,FUN = reduce_info))
    
    m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","UNIT_CHROM","UNIT_START","UNIT_END","COPY_NUMBER")]
    write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
    
    file.remove(step7)
    
  }
  
  mclapply(seq(single_cells_cn),GetAnnotedSNPs,single_cells_cn=single_cells_cn,unit=checked.unit,phased_snps=phased_snps,mc.cores = mc.cores)
  
  # keep only units that are present in all cells
  ks <- list.files(path = file.path(outdir),pattern = 'step5_',full.names = TRUE)
  
  a <- fread(file = ks[1],stringsAsFactors = FALSE,header = FALSE, data.table = FALSE)
  a <- unite(a,col = group,seq(1,3),sep = ":",remove=FALSE)
  a <- a[,c(1,5)]
  colnames(a)[2] <- basename(ks[1])
  
  for(k in seq(2,length(ks))){
    b <- fread(file = ks[k],stringsAsFactors = FALSE,header = FALSE, data.table = FALSE)
    b <- unite(b,col = group,seq(1,3),sep = ":",remove=FALSE)
    b <- b[,c(1,5)]
    colnames(b)[2] <- basename(ks[k])
    a <- merge(x = a,y = b,by='group',all = TRUE)
  }
  
  row.has.na <- apply(a, 1, function(x){any(is.na(x))})
  message(paste("excluding",sum(row.has.na),"units out of",nrow(a)))
  
  unit_to_keep <- a[!row.has.na,]
  
  # step9 : filtering step8
  ks <- list.files(path = file.path(outdir),pattern = 'step8_',full.names = TRUE)
  
  for(k in seq_len(length(ks))){
    message(paste("filtering:",basename(ks[k])))
    
    # step 9
    step9 <- gsub(basename(ks[k]),pattern = "step8_",replacement = "xci_")
    
    m <- fread(file = ks[k],stringsAsFactors = FALSE,header = TRUE, data.table = FALSE)
    m <- unite(m,col = group,seq(6,8),sep = ":",remove=FALSE)
    m <- m[which(m$group %in% unique( unit_to_keep$group )),]
    
    write.table(m[,-6],file = step9,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
    
  }
  remove.files <- list.files(path = outdir,pattern = "step5_|step8_",full.names = TRUE)
  do.call(file.remove,list(remove.files))
  file.remove(checked.unit)
  message("done.")
}