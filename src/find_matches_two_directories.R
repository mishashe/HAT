# For github
# cd /home/msheinman/Development/HAT/
# git add ./src/*.*
# git commit --all -m "added empirical"
# git push origin

# To avoid problems run
# export OPENBLAS_NUM_THREADS=1
# export GOTO_NUM_THREADS=1
# export OMP_NUM_THREADS=1

library(base)
library(foreach)
require(doParallel)
library(stringr)
library(seqinr)
library(tseries)
library(data.table)
library(seqinr)
library(dplyr)
library(stringi)
library(R.utils)
library(plotrix)
registerDoParallel(40)



TablesCombine <- function(tab1,tab2)
{
  r1 <- paste0(tab1[,1],"_str")
  r2 <- paste0(tab2[,1],"_str")
  rownames(tab1) <- r1
  rownames(tab2) <- r2
  r <- unique(c(r1,r2))
  tab <- data.frame(r=unique(c(tab1[,1],tab2[,1])),Freq=rep(0,length(r)))
  rownames(tab) <- r
  
  r1and2 <- intersect(r1,r2)
  r1not2 <- setdiff(r1,r2)
  r2not1 <- setdiff(r2,r1)
  tab[r1and2,2] <- tab1[r1and2,2] + tab2[r1and2,2]
  tab[r1not2,2] <- tab1[r1not2,2]
  tab[r2not1,2] <- tab2[r2not1,2]
  return(tab)
}

genus1 <- "Escherichia"
sp1 <- "coli"
genus2 <- "Buchnera"
sp2 <- "aphidicola"
dir1 <- paste0("/home/m.sheinman/Development/HGTnew/data/external/",genus1,"_",sp1,"/")
dir2 <- paste0("/home/m.sheinman/Development/HGTnew/data/external/",genus2,"_",sp2,"/")
dirOut <- paste0("/home/m.sheinman/Development/HGTnew/data/external/",genus1,"_",sp1,"_VS_",genus2,"_",sp2,"/")
system(paste0("mkdir -p ",dirOut))
system(paste0("mkdir -p ",dirOut,"temp/"))
setwd(paste0(dirOut))
write.table(data.frame("S1",	"E1",	"S2",	"E2",	"LEN1",	"LEN2",	"PER_IDY",	"LEN_R",	"LEN_Q",	"COV_R",	"COV_Q",	"FRM",	"TAGS","genome1","genome2","AccNum1","AccNum2"),
            file = paste0(dirOut,"algns.csv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)  

system(paste0("rm -fr ",dirOut,"temp/*/"))[1:1]
files1 <- Sys.glob(file.path(dir1, "*.fna"))[1:5]
files2 <- Sys.glob(file.path(dir2, "*.fna"))
setwd(paste0(dirOut,"temp/"))

# add lengths
if (1==0)
{
  for (dirt in c(dir1,dir2))
  {
    filest <- Sys.glob(file.path(dirt, "*.fna"))
    Length <- foreach (file = filest, .inorder=FALSE, .combine=c) %do%
    {
      sum(getLength(read.fasta(file = file)))
    }
    write.table(Length,file = paste0(dirt,"Length.csv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  }
}

rT <- foreach (i = 1:length(files1), .inorder=FALSE, .combine=TablesCombine) %do%
{
  file1 <- files1[i]
  AccNum1 <- strsplit(file1,"/")[[1]]
  AccNum1 <- AccNum1[length(AccNum1)]
  AccNum1 <- str_replace(AccNum1,"_genomic.fna","")
  rT <- foreach (j = 1:length(files2), .inorder=FALSE, .combine=TablesCombine) %dopar%
  {
    setwd(paste0(dirOut,"temp/"))
    file2 <- files2[j]
    AccNum2 <- strsplit(file2,"/")[[1]]
    AccNum2 <- AccNum2[length(AccNum2)]
    AccNum2 <- str_replace(AccNum2,"_genomic.fna","")
    system(paste0("/home/m.sheinman/Development/Software/mummer-4.0.0rc1/bin/nucmer --threads 1 --maxmatch -l 8 --prefix=ref_qry_",AccNum1,"_",AccNum2,"_unfiltered ",file1," ",file2))
    system(paste0("/home/m.sheinman/Development/Software/mummer-4.0.0rc1/bin/delta-filter -l 50 -i 30 -q -r ref_qry_",AccNum1,"_",AccNum2,"_unfiltered.delta > ref_qry_",AccNum1,"_",AccNum2,".delta"))
    system(paste0("/home/m.sheinman/Development/Software/mummer-4.0.0rc1/bin/show-coords -c -d -r -l -o -T  ref_qry_",AccNum1,"_",AccNum2,".delta > ref_qry_",AccNum1,"_",AccNum2,".coords"))
    system(paste0("/home/m.sheinman/Development/Software/mummer-4.0.0rc1/bin/show-snps -l -T ref_qry_",AccNum1,"_",AccNum2,".delta > ref_qry_",AccNum1,"_",AccNum2,".snps"))

    algns <- read.table(paste0(dirOut,"temp/ref_qry_",AccNum1,"_",AccNum2,".coords"), header = FALSE, sep = "\t", skip = 3, stringsAsFactors = FALSE,fill=TRUE)
    algns <- unique(algns)
    algns[1,c(14,15)] <- c("genome1","genome2")
    algns <- algns[,c(-16)]
    colnames(algns) <- algns[1,]
    algns <- algns[-1,]
    algns$'[S1]' <- as.numeric(algns$'[S1]')
    algns$'[E1]' <- as.numeric(algns$'[E1]')
    algns$'[S2]' <- as.numeric(algns$'[S2]')
    algns$'[E2]' <- as.numeric(algns$'[E2]')
    algns$'[LEN R]' <- as.numeric(algns$'[LEN R]')
    algns$'[LEN Q]' <- as.numeric(algns$'[LEN Q]')
    # algns <- algns[which(algns$'[LEN R]' > 1e5 & algns$'[LEN Q]' > 1e5),]
    if (nrow(algns)>0)
    {
      snps <- read.table(paste0(dirOut,"temp/ref_qry_",AccNum1,"_",AccNum2,".snps"), header = FALSE, sep = "\t", skip = 3, stringsAsFactors = FALSE,fill=TRUE)
      snps[1,c(13,14)] <- c("genome1","genome2")
      colnames(snps) <- snps[1,]
      snps <- snps[-1,]
      snps$'[P1]' <- as.numeric(snps$'[P1]')
      snps$'[P2]' <- as.numeric(snps$'[P2]')
      
      algns$AccNum1 <- AccNum1
      algns$AccNum2 <- AccNum2
      rT <- foreach (algnsI = 1:nrow(algns), .inorder=FALSE, .combine=TablesCombine) %do%
      {
        S1 <- algns$'[S1]'[algnsI]
        E1 <- algns$'[E1]'[algnsI]
        S2 <- algns$'[S2]'[algnsI]
        E2 <- algns$'[E2]'[algnsI]
        if (S1>E1) {foo <- E1; E1 <- S1; S1 <- foo;}
        if (S2>E2) {foo <- E2; E2 <- S2; S2 <- foo;}
        Ind <- which(snps$genome1==algns$genome1[algnsI] & snps$genome2==algns$genome2[algnsI] & snps$'[P1]'<=E1 & snps$'[P1]'>=S1)
        r1 <- diff(sort(c(S1-1,snps$'[P1]'[Ind],E1+1),decreasing=FALSE))-1
        Ind <- which(snps$genome1==algns$genome1[algnsI] & snps$genome2==algns$genome2[algnsI] & snps$'[P2]'<=E2 & snps$'[P2]'>=S2)
        r2 <- diff(sort(c(S2-1,snps$'[P2]'[Ind],E2+1),decreasing=FALSE))-1
        r <- c(r1,r2)
        r <- r[r>0]
        rT <- as.data.frame(table(r),stringsAsFactors = FALSE)
        colnames(rT) <- c("r","Freq")
        rT$Freq <- rT$Freq/2.0
        return(rT)
      }
      if (i %in% c(1,2,3) & j %in% c(1,2,3))
      {
        system(paste0("/home/m.sheinman/Development/Software/mummer-4.0.0rc1/bin/mummerplot -f -S -l --layout ref_qry_",AccNum1,"_",AccNum2,".delta"))
        system(paste0("mkdir -p ",dirOut,"/plots/dotplot/",AccNum1,"_",AccNum2,"/"))
        system(paste0("mv ",dirOut,"/temp/out.* ",dirOut,"/plots/dotplot/",AccNum1,"_",AccNum2,"/"))
        pdf(paste0(dirOut,"/plots/dotplot/",AccNum1,"_",AccNum2,"/snps.pdf"),width=10,height=20)
        plot(c(0,nrow(algns)+1),c(0,max(algns$'[E1]'-algns$'[S1]')),pch=20, cex=0.001)
        algns <- algns[order( algns$'[E1]'- algns$'[S1]'),]
        for (algnsI in 1:nrow(algns))
        {
          S1 <- algns$'[S1]'[algnsI]
          E1 <- algns$'[E1]'[algnsI]
          if (S1>E1) {foo <- E1; E1 <- S1; S1 <- foo;}
          Ind <- which(snps$genome1==algns$genome1[algnsI] & snps$genome2==algns$genome2[algnsI] & snps$'[P1]'<=E1 & snps$'[P1]'>=S1)
          pos <- snps$'[P1]'[Ind]
          lines(c(algnsI,algnsI),c(0,E1-S1),col="red",lwd=0.2)
          points(rep(algnsI,length(pos)),pos-S1,pch='-', cex=0.2,col="blue")
        }
        dev.off()
      }
    }
    return(rT)
  }
  algns <- foreach (j = 1:length(files2), .combine=rbind) %do%
  {
    file2 <- files2[j]
    AccNum2 <- strsplit(file2,"/")[[1]]
    AccNum2 <- AccNum2[length(AccNum2)]
    AccNum2 <- str_replace(AccNum2,"_genomic.fna","")
    algns <- read.table(paste0(dirOut,"temp/ref_qry_",AccNum1,"_",AccNum2,".coords"), header = FALSE, sep = "\t", skip = 3, stringsAsFactors = FALSE,fill=TRUE)
    algns$AccNum1 <- AccNum1
    algns$AccNum2 <- AccNum2
    system(paste0("rm -fr ",dirOut, "temp/*",AccNum1,"_",AccNum2,"*.*"))
    return(algns[-1,])
  }
  write.table(algns,file = paste0(dirOut,"algns.csv"),append=TRUE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  print(paste0(i," from ",length(files1)))
  return(rT)
}
write.table(rT,file = paste0(dirOut,"rT.csv"),append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)
system(paste0("rm -fr ",dirOut, "temp/*"))
print(paste0(" complete ", dirOut))


IDs <- read.table(paste0(dirOut,"algns.csv"), header = TRUE, sep = "\t",row.names=NULL)
pdf(paste0(dirOut,"plots/","ID.pdf"))
p <- weighted.hist(x=IDs$PER_IDY,w=IDs$LEN1,plot=TRUE)
p <- weighted.hist(x=IDs$PER_IDY,w=IDs$LEN1*0+1,plot=TRUE)
dev.off()

options(digits=2)
Length1 <- read.table(paste0(dir1,"Length.csv"), header = FALSE, sep = "\t",row.names=NULL)$V1
Length2 <- read.table(paste0(dir2,"Length.csv"), header = FALSE, sep = "\t",row.names=NULL)$V1
rT <- read.table(paste0(dirOut,"rT.csv"), header = TRUE, sep = "\t",row.names=NULL)
d <- 0.1
rV <- c(seq(0.5,10.5,1),10^seq(log10(11.5),1*d+log10(max(rT$r)),d))
p <- weighted.hist(x=rT$r,w=rT$Freq,breaks=rV,plot=FALSE)
while(sum(p$counts[2:(length(p$counts)-1)]==0)>0 & d<0.2)
{
  d <- d*1.1
  rV <- c(seq(0.5,10.5,1),10^seq(log10(11.5),1*d+log10(max(rT$r)),d))
  p <- weighted.hist(x=rT$r,w=rT$Freq,breaks=rV,plot=FALSE)
}
knee <- 10
p$counts <- p$counts/sum(Length1)/sum(Length2)
pdf(paste0(dirOut,"plots/","m.pdf"))
plot(log10(p$mids),log10(p$counts/diff(p$breaks)), pch = 19,
     main=paste0("1/<r>=",round(1/sum(rT$r*rT$Freq)*sum(rT$Freq),3)))
lines(log10(p$mids),log10(p$counts/diff(p$breaks))[10]+3*log10(p$mids[10]+knee)-3*log10(p$mids+knee))
dev.off()

