library( tidyr )
setwd("/icgc/dkfzlsdf/analysis/B260/users/n790i/generate_xclone_inputs/studio/")

## GTseq

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

