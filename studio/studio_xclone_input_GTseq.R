library( tidyverse )
library( parallel )
library( Matrix )
library( data.table )
library( vioplot )

setwd('/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/')

cores <- 5
## [ input ]

GenomicUnit <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_SP/unit_genes/"

omic <- "rna"
#omic <- "dna"


## [ functions ]

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

# generate data 
if(TRUE){

dir.create(paste(omic,basename(GenomicUnit),sep = '_'))

## GTseq - DNA
if(omic=="dna"){
  setwd("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/dna_unit_genes")
  cellsnp <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scase/cellSNP_test/scDNA_data/sample_3/sparseVCF"
}

## GTseq - RNA
if(omic=="rna"){
  setwd("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/rna_unit_genes")
  cellsnp <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/scase/cellSNP_test/scRNA_data/sample_3/sparseVCF"
}

getwd()

dp.mtx.file  <- file.path(cellsnp,"cellSNP.tag.DP.mtx")
ad.mtx.file  <- file.path(cellsnp,"cellSNP.tag.AD.mtx")
samples.file <- file.path(cellsnp,"cellSNP.samples.tsv")
snps.file <- file.path(cellsnp,"cellSNP.base.vcf.gz")

xci.files <- list.files(GenomicUnit,pattern = "xci_",full.names = TRUE)

dp.mtx <- readMM(file = dp.mtx.file)
ad.mtx <- readMM(file = ad.mtx.file)
snps <- fread(input = snps.file,data.table = FALSE,stringsAsFactors = FALSE,skip = 1,select = c(1,2),verbose = FALSE)
snps[,1] <- gsub(snps[,1],pattern = "chr",replacement = "")
snps <- unite(snps,col = SNP,seq(1,2),sep = ":",remove=TRUE)

samples <- readLines(samples.file)

if(omic=="dna"){
  samples <- gsub(samples,pattern = "lane1DNA|_sequence.bam",replacement = "")
}
if(omic=="rna"){
  samples <- gsub(samples,pattern = ".Aligned.sortedByCoord.out.bam",replacement = "")
}

cat(samples)

getDP <- function(i,xci.files,snps,dp.mtx){
  xc <- xci.files[i]
  m <- fread(xc,data.table = F)
  m <- unite(m,col = SNP,seq(1,2),sep = ":",remove=TRUE)
  s <- gsub(basename(xc),pattern = "xci_lane1DNA|_sequence.cbs.nochr.bed",replacement = "")
  dpc <- data.frame(SNP=snps$SNP,
                    SNP_DP=as.numeric(dp.mtx[,which(samples == s)]),
                    SNP_AD=as.numeric(ad.mtx[,which(samples == s)]),
                    stringsAsFactors = F)
  m <- merge(m,dpc,by = "SNP",all.x = TRUE)
  write.table(m,file = gsub(basename(xc),pattern = "xci_",replacement = "xci_dp_"),col.names = T,row.names = F,sep = "\t",quote = F)
}

mclapply(seq(length(xci.files)),getDP,xci.files=xci.files,snps=snps,dp.mtx=dp.mtx,mc.cores = cores)

}

##################################################
# SNP Analysis
##################################################

dir.create('png')

bed.dna <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/dna_unit_genes",pattern = "\\.bed$",full.names = T)
bed.rna <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/rna_unit_genes",pattern = "\\.bed$",full.names = T)

# SNPs with DP > th in DNA or RNA - cell annotated
df.all.dna <- do.call(rbind,mclapply(seq(length(bed.dna)),DataFilter,minDP=1,bed=bed.dna,cn=NA,mc.cores = cores))
df.all.rna <- do.call(rbind,mclapply(seq(length(bed.rna)),DataFilter,minDP=1,bed=bed.rna,cn=NA,mc.cores = cores))

df.all.dna <- unite(df.all.dna,col = UNIT,seq(5,7),sep = ":",remove=FALSE)
df.all.rna <- unite(df.all.rna,col = UNIT,seq(5,7),sep = ":",remove=FALSE)

# Number of SNPs with DP in DNA and RNA - per cell
m <- merge(data.frame(table(df.all.dna$CELL_ID),stringsAsFactors = F),
           data.frame(table(df.all.rna$CELL_ID),stringsAsFactors = F),
           by = "Var1",suffixes = c('_dna','_rna'))

# Plot1 | n. of SNPs per cell in DNA and RNA
png("png/Plot1.png",width = 12,height = 4,units = 'in', res = 300)

layout(matrix(c(1,1,1,3,3,
                2,2,2,3,3),2,5,byrow = T))

barplot(m$Freq_dna,col = 'grey60',border = 'grey60',names.arg = m$Var1,cex.names = .6,las=2,main='n. of SNPs per cell in scDNA')
barplot(m$Freq_rna,col = 'grey60',border = 'grey60',names.arg = m$Var1,cex.names = .6,las=2,main='n. of SNPs per cell in scRNA')
vioplot(m$Freq_dna,m$Freq_rna,names = c("scDNA", "scRNA"),ylab="n. of SNPs per cell",pchMed = 20,col="white",rectCol = "white",colMed = "black",main=paste(median(m$Freq_dna),median(m$Freq_rna),sep = ' | '))

dev.off()

# Plot2 | DP of SNPs in DNA and RNA 
png("png/Plot2.png",width = 14,height = 7,units = 'in', res = 300)

par(pty='s',mfrow=c(1,2))

boxplot(df.all.dna$SNP_DP,df.all.rna$SNP_DP,names = c("scDNA", "scRNA"),ylab="DP",outline=F,varwidth = T,main=paste(median(df.all.dna$SNP_DP),median(df.all.rna$SNP_DP),sep = ' | '))
vioplot(df.all.dna$SNP_DP,df.all.rna$SNP_DP,names = c("scDNA", "scRNA"),pchMed = 20,ylab="DP",col="white",rectCol = "white",colMed = "black")

dev.off()

# SNPs with DP in both DNA and RNA - annotated per cell
dfm <- merge(x = df.all.dna[,c(1,4,5,11:15)],df.all.rna[,c(1,4,5,11:15)],by = c('SNP','SNP_PHASE_INFO','CELL_ID','COPY_NUMBER',"UNIT"),all = FALSE,suffixes = c('_dna','_rna'))

# Plot3 | SNPs with DP in both DNA and RNA per cell 
png("png/Plot3.png",width = 12,height = 6,units = 'in', res = 300)

layout(matrix(c(1,1,1,1,2,2,
                1,1,1,1,2,2,
                1,1,1,1,3,3,
                1,1,1,1,3,3),4,6,byrow = T))

barplot(sort(table(dfm$CELL_ID),decreasing = T),cex.names = 0.6,las=2,xlab='cells',ylab='n. of SNPs covered both in DNA and RNA',col = 'grey60',border = 'grey60')
abline(h = median(table(dfm$CELL_ID)),lwd=0.6)
boxplot(dfm$SNP_DP_dna,dfm$SNP_DP_rna,names = c("scDNA", "scRNA"),ylab="DP",outline=F,varwidth = T,main=paste(median(dfm$SNP_DP_dna),median(dfm$SNP_DP_rna),sep = ' | '))
vioplot(dfm$SNP_DP_dna,dfm$SNP_DP_rna,names = c("scDNA", "scRNA"),pchMed = 20,ylab="DP",col="white",rectCol = "white",colMed = "black")

dev.off()

# Plot4a | ASR absolute deviance from 0.5 for SNPs in copy number 2 segments
png("png/Plot4.png",width = 15,height = 7,units = 'in', res = 300)

par(pty='s')

k.dna <- df.all.dna[which(df.all.dna$COPY_NUMBER==2),]
k.rna <- df.all.rna[which(df.all.rna$COPY_NUMBER==2),]

vioplot(abs(0.5-k.dna$SNP_ASR),abs(0.5-k.rna$SNP_ASR),names = c("scDNA", "scRNA"),pchMed = 20,ylab="|0.5-(AD/DP)|",col="white",rectCol = "white",colMed = "black")

dev.off()

# Plot4b | ASR absolute deviance from 0.5 for SNPs - considering only SNPs with DP in both DNA and RNA - in copy number 2 segments
png("png/Plot4_selected_snps.png",width = 9,height = 6,units = 'in', res = 300)

layout(matrix(c(1,3,3,
                2,3,3),2,3,byrow = T))

k <- dfm[which(dfm$COPY_NUMBER==2),]

plot(k$SNP_DP_dna,k$SNP_DP_rna,pch=20,col='grey60',xlab='DP in DNA',ylab='DP in RNA',main=paste('n. SNPs = ',nrow(k),' in CN = 2'))
plot(abs(0.5-k$SNP_ASR_dna),abs(0.5-k$SNP_ASR_rna),pch=20,col='grey60',xlab='|0.5-(AD/DP)| in DNA',ylab='|0.5-(AD/DP)| in RNA')
vioplot(abs(0.5-k$SNP_ASR_dna),abs(0.5-k$SNP_ASR_rna),names = c("scDNA", "scRNA"),pchMed = 20,ylab="|0.5-(AD/DP)|",col="white",rectCol = "white",colMed = "black")

dev.off()

# Plot5_cn2 | ASR absolute deviance from 0.5 for GENEs in copy number 2 segments

getCountPerGene <- function(tab,cn=NA){
  
  if(!is.na(cn)){
    x <- tab %>% filter(COPY_NUMBER == cn)
    out <- split(x,f = x$UNIT)
  } else {
    out <- split(tab,f = tab$UNIT)
  }
  
  sumDP <- function(df){
    return(sum(df$SNP_DP))
  }

  sum.maternal <- function(df){
    mat.ad <- df$SNP_AD[which(df$SNP_PHASE_INFO == '1|0')]
    pat.ad <- df$SNP_AD[which(df$SNP_PHASE_INFO == '0|1')]
    pat.dp <- df$SNP_DP[which(df$SNP_PHASE_INFO == '0|1')]
    return(sum(mat.ad,(pat.dp-pat.ad)))
  }
  
  mat.reads <- unlist(lapply(out,FUN = sum.maternal))
  dp <- unlist(lapply(out,FUN = sumDP))
  asr <- mat.reads/dp
  return(list(dp=dp,ad=mat.reads,asr=asr))
  
}

png("png/Plot5_cn2.png",width = 8,height = 4,units = 'in', res = 300)

genes.dna <- getCountPerGene(tab = df.all.dna,cn = 2)
genes.rna <- getCountPerGene(tab = df.all.rna,cn = 2)

layout(matrix(c(1,2,3,3,
                1,2,3,3),2,4,byrow = T))

vioplot(genes.dna$dp,genes.rna$dp,names = c("scDNA", "scRNA"),pchMed = 20,ylab="DP per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F)
boxplot(genes.dna$dp,genes.rna$dp,outline = F,names = c("scDNA", "scRNA"),ylab="DP per gene",frame.plot=F,main='(no outliers)')
vioplot(abs(0.5-(genes.dna$asr)),abs(0.5-(genes.rna$asr)),names = c("scDNA", "scRNA"),ylim=c(0,0.5),pchMed = 20,ylab="|0.5-(AD/DP)| per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F,main='Genes in CN = 2')

dev.off()

png("png/Plot5_cn1.png",width = 8,height = 4,units = 'in', res = 300)

genes.dna <- getCountPerGene(tab = df.all.dna,cn = 1)
genes.rna <- getCountPerGene(tab = df.all.rna,cn = 1)

layout(matrix(c(1,2,3,3,
                1,2,3,3),2,4,byrow = T))

vioplot(genes.dna$dp,genes.rna$dp,names = c("scDNA", "scRNA"),pchMed = 20,ylab="DP per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F)
boxplot(genes.dna$dp,genes.rna$dp,outline = F,names = c("scDNA", "scRNA"),ylab="DP per gene",frame.plot=F,main='(no outliers)')
vioplot(abs(0.5-(genes.dna$asr)),abs(0.5-(genes.rna$asr)),names = c("scDNA", "scRNA"),ylim=c(0,0.5),pchMed = 20,ylab="|0.5-(AD/DP)| per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F,main='Genes in CN = 1')

dev.off()


# Plot6 | ASR absolute deviance from 0.5 per GENEs - considering only SNPs with DP in both DNA and RNA - in copy number 2 segments
d.dp <- c()
d.ad <- c()
r.dp <- c()
r.ad <- c()

sum.maternal <- function(df,omic){
  if(omic == 'dna'){
    mat.ad <- df$SNP_AD_dna[which(df$SNP_PHASE_INFO == '1|0')]
    pat.ad <- df$SNP_AD_dna[which(df$SNP_PHASE_INFO == '0|1')]
    pat.dp <- df$SNP_DP_dna[which(df$SNP_PHASE_INFO == '0|1')]
  }
  if(omic == 'rna'){
    mat.ad <- df$SNP_AD_rna[which(df$SNP_PHASE_INFO == '1|0')]
    pat.ad <- df$SNP_AD_rna[which(df$SNP_PHASE_INFO == '0|1')]
    pat.dp <- df$SNP_DP_rna[which(df$SNP_PHASE_INFO == '0|1')]
  }
  return(sum(mat.ad,(pat.dp-pat.ad)))
}

for( u in unique(dfm$UNIT[which(dfm$COPY_NUMBER == 2 )]) ){
  d.dp <- c(d.dp, sum(dfm$SNP_DP_dna[which(dfm$UNIT == u)]))
  d.ad <- c(d.ad, sum.maternal(df = dfm[which(dfm$UNIT == u),],omic = 'dna'))
  r.dp <- c(r.dp, sum(dfm$SNP_DP_rna[which(dfm$UNIT == u)]))
  r.ad <- c(r.ad, sum.maternal(df = dfm[which(dfm$UNIT == u),],omic = 'rna'))
}

png("png/Plot6.png",width = 9,height = 7,units = 'in', res = 300)

layout(matrix(c(1,3,3,
                2,3,3),2,3,byrow = T))

plot(d.dp,r.dp,pch=20,col='grey60',xlab='DP in DNA',ylab='DP in RNA',main=paste('n. genes = ',length(d.dp),' in CN = 2'))
plot(abs(0.5-(d.ad/d.dp)),abs(0.5-(r.ad/r.dp)),pch=20,col='grey60',xlab='|0.5-(AD/DP)| in DNA',ylab='|0.5-(AD/DP)| in RNA')
vioplot(abs(0.5-(d.ad/d.dp)),abs(0.5-(r.ad/r.dp)),names = c("scDNA", "scRNA"),pchMed = 20,ylab="|0.5-(AD/DP)|",col="white",rectCol = "white",colMed = "black")

dev.off()

# Plot7 | Select genes that are diploid in 50% of cells and loss in the other 50%

xci.files <- list.files(GenomicUnit,pattern = "xci_",full.names = TRUE)

tab <- c()
for(i in seq(length(xci.files))){
  message(i,"\t",basename(xci.files[i]))
  if(i==1){
    xc <- xci.files[i]
    m <- unique(fread(xc,data.table = F,select = c(6:8,11)))
    m <- unite(m,col = UNIT,seq(1,3),sep = ":",remove=TRUE)
    names(m)[2] <- gsub(basename(xc),pattern = "xci_lane1DNA|_sequence.cbs.nochr.bed",replacement = "")
    tab <- m
  } else {
    xc <- xci.files[i]
    m <- unique(fread(xc,data.table = F,select = c(6:8,11)))
    m <- unite(m,col = UNIT,seq(1,3),sep = ":",remove=TRUE)
    names(m)[2] <- gsub(basename(xc),pattern = "xci_lane1DNA|_sequence.cbs.nochr.bed",replacement = "")
    tab <- merge(tab,m,by = "UNIT",all = TRUE)
  }
}

cn1 <- rowSums(tab == 1)/ncol(tab)
cn2 <- rowSums(tab == 2)/ncol(tab)

selected_genes <- intersect(which(cn1 > 0.4 & cn1 < 0.6),which(cn2 > 0.4 & cn2 < 0.6))

png("png/Plot7.png",width = 7,height = 7,units = 'in', res = 300)

par(pty='s')

col <- rep("grey60",length(cn1))
col[selected_genes] <- "orangered"

plot(cn1,cn2,pch=20,
     xlab = "fraction of cells with this gene in CN = 1",
     ylab = "fraction of cells with this gene in CN = 2",ylim=c(0,1),xlim=c(0,1),
     col=col,main=paste("selected genes = ",length(selected_genes)))
abline(h = c(0.4,0.6),v = c(0.4,0.6),lwd=0.4)

dev.off()

# Plot8 | plot ASR in genes that have cn=1 and cn=2 in about 50% of cells

sum.maternal <- function(df){
  mat.ad <- df$SNP_AD[which(df$SNP_PHASE_INFO == '1|0')]
  pat.ad <- df$SNP_AD[which(df$SNP_PHASE_INFO == '0|1')]
  pat.dp <- df$SNP_DP[which(df$SNP_PHASE_INFO == '0|1')]
  return(sum(mat.ad,(pat.dp-pat.ad)))
}

cate <- function(bed, selected_genes,minDP,filter_snps=FALSE,selected_snps=NA){
  
  df <- do.call(rbind,mclapply(seq(length(bed)),DataFilter,bed=bed,minDP=minDP,cn=c(1,2),mc.cores = cores))
  df <- unite(df,col = UNIT,seq(5,7),sep = ":",remove=TRUE)
  
  if(filter_snps){
    df <- df[which(df$SNP %in% selected_snps),,drop=FALSE]
  }

  genes_cn1 <- c()
  genes_cn2 <- c()
  for(gene in tab$UNIT[selected_genes]){
    
    u <- df[which(df$UNIT == gene),,drop=F]
    if(nrow(u) > 0){
      summaryGene <- function(u,cn){
        g <- u[which(u$COPY_NUMBER == cn),,drop=FALSE]
        if(nrow(g) > 0){
          this <- data.frame(UNIT=unique(g$UNIT),
                             CN=unique(g$COPY_NUMBER),
                             SNPS=nrow(g),
                             CELLS=length(unique(g$CELL_ID)),
                             DP=sum(g$SNP_DP),
                             AD=sum.maternal(g),
                             ASR=sum(g$SNP_AD)/sum(g$SNP_DP),stringsAsFactors = F)
          return(this)
        } else {
          return(NULL)
        }
      }
      
      genes_cn1 <- rbind(genes_cn1,summaryGene(u=u,cn=1))
      genes_cn2 <- rbind(genes_cn2,summaryGene(u=u,cn=2))
    }
  }
  
  gcn <- merge(x = genes_cn1,genes_cn2,by = "UNIT",all = FALSE,suffixes = c("_cn1","_cn2"))
  
  return(gcn)
  
}

png("png/Plot8.png",width = 7,height = 7,units = 'in', res = 300)

gcn.dna <- cate(bed.dna,selected_genes,minDP = 1,filter_snps = FALSE,selected_snps = unique(dfm$SNP))
gcn.rna <- cate(bed.rna,selected_genes,minDP = 1,filter_snps = FALSE,selected_snps = unique(dfm$SNP))

par(pty="s",mfrow=c(2,2))

plot(abs(0.5-gcn.dna$ASR_cn1),abs(0.5-gcn.dna$ASR_cn2),xlab="abs(0.5-ASR) in CN=1",ylab="abs(0.5-ASR) in CN=2",xlim=c(0,0.5),ylim=c(0,0.5),main='DNA')
vioplot(abs(0.5-gcn.dna$ASR_cn1),abs(0.5-gcn.dna$ASR_cn2),names = c("CN = 1", "CN = 2"),ylab="abs(0.5-ASR)",col="white",rectCol = "white",colMed = "black",pchMed = 20,main='DNA')

plot(abs(0.5-gcn.rna$ASR_cn1),abs(0.5-gcn.rna$ASR_cn2),xlab="abs(0.5-ASR) in CN=1",ylab="abs(0.5-ASR) in CN=2",xlim=c(0,0.5),ylim=c(0,0.5),main='RNA')
vioplot(abs(0.5-gcn.rna$ASR_cn1),abs(0.5-gcn.rna$ASR_cn2),names = c("CN = 1", "CN = 2"),ylab="abs(0.5-ASR)",col="white",rectCol = "white",colMed = "black",pchMed = 20,main='RNA')

dev.off()

png("png/Plot8_selected_snps.png",width = 7,height = 7,units = 'in', res = 300)

gcn.dna <- cate(bed.dna,selected_genes,minDP = 1,filter_snps = TRUE,selected_snps = unique(dfm$SNP))
gcn.rna <- cate(bed.rna,selected_genes,minDP = 1,filter_snps = TRUE,selected_snps = unique(dfm$SNP))

par(pty="s",mfrow=c(2,2))

plot(abs(0.5-gcn.dna$ASR_cn1),abs(0.5-gcn.dna$ASR_cn2),xlab="abs(0.5-ASR) in CN=1",ylab="abs(0.5-ASR) in CN=2",xlim=c(0,0.5),ylim=c(0,0.5),main='DNA')
vioplot(abs(0.5-gcn.dna$ASR_cn1),abs(0.5-gcn.dna$ASR_cn2),names = c("CN = 1", "CN = 2"),ylab="abs(0.5-ASR)",col="white",rectCol = "white",colMed = "black",pchMed = 20,main='DNA')

plot(abs(0.5-gcn.rna$ASR_cn1),abs(0.5-gcn.rna$ASR_cn2),xlab="abs(0.5-ASR) in CN=1",ylab="abs(0.5-ASR) in CN=2",xlim=c(0,0.5),ylim=c(0,0.5),main='RNA')
vioplot(abs(0.5-gcn.rna$ASR_cn1),abs(0.5-gcn.rna$ASR_cn2),names = c("CN = 1", "CN = 2"),ylab="abs(0.5-ASR)",col="white",rectCol = "white",colMed = "black",pchMed = 20,main='RNA')

dev.off()

# Plot9 | check for allelic imbalance on chromosome 3p 

genes.dna <- getCountPerGene(tab = df.all.dna[which( df.all.dna$GENE_CHROM == 3 & df.all.dna$GENE_START < 88000000),],cn=1)
genes.rna <- getCountPerGene(tab = df.all.rna[which( df.all.rna$GENE_CHROM == 3 & df.all.rna$GENE_START < 88000000),],cn=1)

png("png/Plot9.png",width = 8,height = 4,units = 'in', res = 300)

layout(matrix(c(1,2,3,3,
                1,2,3,3),2,4,byrow = T))

vioplot(genes.dna$dp,genes.rna$dp,names = c("scDNA", "scRNA"),pchMed = 20,ylab="DP per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F)
boxplot(genes.dna$dp,genes.rna$dp,outline = F,names = c("scDNA", "scRNA"),ylab="DP per gene",frame.plot=F,main='(no outliers)')
vioplot(abs(0.5-(genes.dna$asr)),abs(0.5-(genes.rna$asr)),names = c("scDNA", "scRNA"),pchMed = 20,ylab="|0.5-(AD/DP)| per gene",col="white",rectCol = "white",colMed = "black",frame.plot = F,main='Genes in 3p')

dev.off()
