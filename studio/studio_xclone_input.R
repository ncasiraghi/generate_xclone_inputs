library( tidyr )
library( parallel )
library( Matrix )
library( data.table )

cores <- 10

omic <- "rna"

# GenomicUnit <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_SP/unit_cn_blocks/"
GenomicUnit <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_SP/unit_genes/"

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


## [ functions ]

DataFilter <- function(i,bed,minDP=1,maxDP=Inf,minAD=0,cn=NA,maxCN=5){
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


# SNP analysis
if( FALSE ){

#bed <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/dna_unit_genes",pattern = "\\.bed$",full.names = T)
bed <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/rna_unit_genes",pattern = "\\.bed$",full.names = T)
  
df.all <- do.call(rbind,mclapply(seq(length(bed)),DataFilter,bed=bed,cn=2,mc.cores = cores))

summary(df.all$COPY_NUMBER)
quantile(df.all$SNP_DP)
quantile(df.all$SNP_DP,probs = seq(0,1,0.01))


plotres <- function(df,title){
  
  boxplot(df$SNP_DP,df$SNP_AD,outline = F,names = c("DP","AD"),ylab="n. reads",main=title)
  text(1:2,y = c(median(df$SNP_DP),median(df$SNP_AD)), labels = c(median(df$SNP_DP),median(df$SNP_AD)),pos = 3)
  
  hist(df$SNP_ASR,20,col="grey",border = "grey",xlab = "ASR",xlim=c(0,1),main="")
  abline(v = median(df$SNP_ASR),col="orangered")
  
  plot(density(df$SNP_ASR),col="black",lwd=4,xlab = "ASR",main=title)
  
}

par(pty="s",mfrow=c(1,2))

# check in diploid segments
minDP <- quantile(df.all$SNP_DP,probs = 0.99)
minDP <- 20

df <- do.call(rbind,mclapply(seq(length(bed)),DataFilter,bed=bed,minDP=minDP,maxDP=Inf,minAD=0,cn=2,mc.cores = cores))

summary(df$COPY_NUMBER)
summary(df$SNP_DP)
summary(df$SNP_AD)

plotres(df,title = "SNPs in CN = 2")


# check in segments with CN = 1
df <- do.call(rbind,mclapply(seq(length(bed)),DataFilter,bed=bed,minDP=minDP,maxDP=Inf,minAD=0,cn=1,mc.cores = cores))

summary(df$COPY_NUMBER)
summary(df$SNP_DP)
summary(df$SNP_AD)

plotres(df,title = "SNPs in CN = 1")

}

# Gene analysis

GenomicUnit <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/xclone_inputs/GTseq/STP-PDX/mode_SP/unit_genes/"
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

par(pty='s')
col <- rep("grey60",length(cn1))

selected_genes <- intersect(which(cn1 > 0.4 & cn1 < 0.6),which(cn2 > 0.4 & cn2 < 0.6))

col[selected_genes] <- "orangered"

plot(cn1,cn2,pch=20,
     xlab = "fraction of cells with this gene in CN = 1",
     ylab = "fraction of cells with this gene in CN = 2",ylim=c(0,1),xlim=c(0,1),
     col=col,main=paste("selected genes = ",length(selected_genes)))


#bed <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/dna_unit_genes",pattern = "\\.bed$",full.names = T)
#bed <- list.files("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/rna_unit_genes",pattern = "\\.bed$",full.names = T)

df <- do.call(rbind,mclapply(seq(length(bed)),DataFilter,bed=bed,minDP=1,cn=c(1,2),mc.cores = cores))

summary(df$COPY_NUMBER)
summary(df$SNP_DP)
summary(df$SNP_AD)

df <- unite(df,col = UNIT,seq(5,7),sep = ":",remove=TRUE)

genes_cn1 <- c()
genes_cn2 <- c()
for(gene in tab$UNIT[selected_genes]){
  message(gene)
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
                           AD=sum(g$SNP_AD),
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

# Add boxplots to a scatterplot
library(vioplot)

par(pty="s",mfrow=c(1,2))
plot(abs(0.5-gcn$ASR_cn1),abs(0.5-gcn$ASR_cn2),xlab="abs(0.5-ASR) in CN=1",ylab="abs(0.5-ASR) in CN=2",xlim=c(0,0.5),ylim=c(0,0.5))
vioplot(abs(0.5-gcn$ASR_cn1),abs(0.5-gcn$ASR_cn2),names = c("CN = 1", "CN = 2"),ylab="abs(0.5-ASR)",col="white",rectCol = "white",colMed = "black")
#vioplot(gcn$DP_cn1,gcn$DP_cn2,names = c("CN = 1", "CN = 2"),ylab="DP",col="white",rectCol = "white",colMed = "black")
mtext("scRNA", side=3, outer=TRUE, line=-3)

# old
if( FALSE ){

# scDNA
if(F){
ad.dna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scDNA/for_nicola/ad.csv",header = TRUE,stringsAsFactors = FALSE)
dp.dna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scDNA/for_nicola/dp.csv",header = TRUE,stringsAsFactors = FALSE)
block.dna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scDNA/for_nicola/block_cnv_per_cell.csv",header = TRUE,stringsAsFactors = FALSE)
block.dna <- unite(block.dna,col = group,seq(1,3),sep = ":",remove=TRUE)
ad.dna$group <- block.dna$group
dp.dna$group <- block.dna$group

# scRNA
ad.rna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scRNA/for_nicola/ad.csv",header = TRUE,stringsAsFactors = FALSE)
dp.rna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scRNA/for_nicola/dp.csv",header = TRUE,stringsAsFactors = FALSE)
block.rna <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/xclone/data/processed/STP_G&T/scRNA/for_nicola/block_cnv_per_cell.csv",header = TRUE,stringsAsFactors = FALSE)
block.rna <- unite(block.rna,col = group,seq(1,3),sep = ":",remove=TRUE)
ad.rna$group <- block.rna$group
dp.rna$group <- block.rna$group

# keep only intersected blocks
keep <- intersect(block.dna$group,block.rna$group)

block.dna <- block.dna[which(block.dna$group %in% keep),]
block.rna <- block.rna[which(block.rna$group %in% keep),]

ad.dna <- melt(ad.dna[which(ad.dna$group %in% keep),],id=c("group"))
ad.rna <- melt(ad.rna[which(ad.rna$group %in% keep),],id=c("group"))
ad.dna$variable <- gsub(ad.dna$variable,pattern = "_ad",replacement = "")
ad.rna$variable <- gsub(ad.rna$variable,pattern = "_ad",replacement = "")

ad <- merge(ad.dna,ad.rna,by = c("group","variable"),suffixes = c("_dna","_rna"))
colnames(ad) <- gsub(colnames(ad),pattern = "value_",replacement = "AD_")
head(ad)

dp.dna <- melt(dp.dna[which(dp.dna$group %in% keep),],id=c("group"))
dp.rna <- melt(dp.rna[which(dp.rna$group %in% keep),],id=c("group"))
dp.dna$variable <- gsub(dp.dna$variable,pattern = "_dp",replacement = "")
dp.rna$variable <- gsub(dp.rna$variable,pattern = "_dp",replacement = "")

dp <- merge(dp.dna,dp.rna,by = c("group","variable"),suffixes = c("_dna","_rna"))
colnames(dp) <- gsub(colnames(dp),pattern = "value_",replacement = "DP_")
head(dp)

cn <- melt(block.dna, id=c("group"))
colnames(cn)[3] <- "COPY_NUMBER"

df <- merge(cn,ad,by = c("group","variable"))
df <- merge(df,dp,by = c("group","variable"))
df$ASR_dna <- df$AD_dna/df$DP_dna
df$ASR_rna <- df$AD_rna/df$DP_rna
}

m <- df[which(df$COPY_NUMBER %in% c(1:3)),]
m <- m[grep(m$group,pattern = "^3:"),]

barplot(table(m$COPY_NUMBER))

m <- m[which(m$DP_dna > 0 & m$DP_rna > 0),]

par(pty="s",mfrow=c(3,2))
boxplot(DP_dna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)
boxplot(DP_rna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)

boxplot(AD_dna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)
boxplot(AD_rna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)

boxplot(ASR_dna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)
abline(h = 0.5)
boxplot(ASR_rna ~ COPY_NUMBER,data = m,varwidth=TRUE,outline=FALSE)
abline(h = 0.5)

##

m <- m[which(m$COPY_NUMBER == 2),]

par(pty="s",mfrow=c(2,2))

rbPal <- colorRampPalette(c('#edf8b1','#225ea8'))
m$Col <- rbPal(10)[as.numeric(cut(log(m$DP_dna),breaks = 10))]

plot(m$ASR_dna,m$ASR_rna,xlim=c(0,1),ylim=c(0,1),pch = 20,col = m$Col)

plot(density(m$ASR_dna),col="orangered",lwd=4,xlab = "ASR",main="")
lines(density(m$ASR_rna),col="black",lwd=4)

par(pty="s",mfrow=c(3,3))
for(x in unique(m$variable)){
  a <- m[which(m$variable == x),]
  plot(a$ASR_dna,a$ASR_rna,main=x,xlim=c(0,1),ylim=c(0,1))
  abline(h = 0.5,v = 0.5,col="red",lwd=0.4)
}


}
