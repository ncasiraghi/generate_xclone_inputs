#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

if(length(args)!=4){
  message("\n\tError!\n\tUsage: studio_xclone_input.R [xci_folder or unit] [cellsnp] [outdir] [cores]\n")
  quit()
}

input_data <- args[1]
cellsnp <- args[2]
outdir <- args[3]
cores <- args[4]

library( tidyverse )
library( parallel )
library( Matrix )
library( data.table )
library( vioplot )

## [ FUNCTIONS ]

DataFilter <- function(i,bed,minDP=2,maxDP=Inf,minAD=0,cn=NA,maxCN=5){
  message(i)
  x <- fread(bed[i],data.table = FALSE)
  x <- x[which(x$SNP_DP >= minDP & x$SNP_DP <= maxDP),]
  x <- x[which(x$SNP_AD >= minAD),]
  x <- x[which(x$COPY_NUMBER <= maxCN),]
  if(unique(!is.na(cn))){
    x <- x[which(x$COPY_NUMBER %in% cn),,drop=FALSE]
  }
  if(nrow(x)==0){
    return(NULL)
  } else {
    x$SNP_ASR <- x$SNP_AD/x$SNP_DP
    x$CELL_ID <- gsub(basename(bed[i]),pattern = "xci_dp_lane1DNA|_sequence.cbs.nochr.bed",replacement = "")
    return(x)
  }
}

plotres <- function(df,title){
  
  boxplot(df$SNP_DP,df$SNP_AD,outline = F,names = c("DP","AD"),ylab="n. reads",main=title)
  text(1:2,y = c(median(df$SNP_DP),median(df$SNP_AD)), labels = c(median(df$SNP_DP),median(df$SNP_AD)),pos = 3)
  
  hist(df$SNP_ASR,20,col="grey",border = "grey",xlab = "ASR",xlim=c(0,1),main="")
  abline(v = median(df$SNP_ASR),col="orangered")
  
  plot(density(df$SNP_ASR),col="black",lwd=4,xlab = "ASR",main=title)
  
}

# [ GENERATE DATA ] 
xci_info_folder <- file.path(outdir,'xci_info')
dir.create(xci_info_folder) 
setwd(xci_info_folder)

if(file_test("-d",input_data)){
  
  dp.mtx.file  <- file.path(cellsnp,"cellSNP.tag.DP.mtx")
  ad.mtx.file  <- file.path(cellsnp,"cellSNP.tag.AD.mtx")
  samples.file <- file.path(cellsnp,"cellSNP.samples.tsv")
  snps.file <- file.path(cellsnp,"cellSNP.base.vcf.gz")
  
  xci.files <- list.files(input_data,pattern = "xci_",full.names = TRUE)
  
  dp.mtx <- readMM(file = dp.mtx.file)
  ad.mtx <- readMM(file = ad.mtx.file)
  snps <- fread(input = snps.file,data.table = FALSE,stringsAsFactors = FALSE,skip = 1,select = c(1,2),verbose = FALSE)
  snps[,1] <- gsub(snps[,1],pattern = "chr",replacement = "")
  snps <- unite(snps,col = SNP,seq(1,2),sep = ":",remove=TRUE)
  
  samples <- readLines(samples.file)
  
  getDP <- function(i,xci.files,snps,dp.mtx,gsub.pattern=NA){
    xc <- xci.files[i]
    m <- fread(xc,data.table = F)
    m <- unite(m,col = SNP,seq(1,2),sep = ':',remove=TRUE)
    s <- gsub(basename(xc),pattern = gsub.pattern, replacement = '')
    s <- gsub(s,pattern = '\\.bed$',replacement = '')
    dpc <- data.frame(SNP=snps$SNP,
                      SNP_DP=as.numeric(dp.mtx[,which(samples == s)]),
                      SNP_AD=as.numeric(ad.mtx[,which(samples == s)]),
                      stringsAsFactors = F)
    m <- merge(m,dpc,by = "SNP",all.x = TRUE)
    write.table(m,file = gsub(basename(xc),pattern = "xci_",replacement = "xci_info_"),col.names = T,row.names = F,sep = "\t",quote = F)
  }
  
  mclapply(seq(length(xci.files)),getDP,xci.files=xci.files,snps=snps,dp.mtx=dp.mtx,gsub.pattern='xci_cellbarcode_',mc.cores = cores)
  
}

if(file_test("-f",input_data)){
  
  dp.mtx.file  <- file.path(cellsnp,"cellSNP.tag.DP.mtx")
  ad.mtx.file  <- file.path(cellsnp,"cellSNP.tag.AD.mtx")
  samples.file <- file.path(cellsnp,"cellSNP.samples.tsv")
  snps.file <- file.path(cellsnp,"cellSNP.base.vcf.gz")
  
  format.unit  <- fread(input_data,data.table = FALSE,select = 1:3)
  checked.unit <- file.path(xci_info_folder,'unit.bed')
  write.table(format.unit,file = checked.unit,quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t')
  
  cmd <- paste('intersectBed -a',checked.unit,'-b',snps.file,'-wa -wb > snps_unit.bed')
  system(cmd)

}
