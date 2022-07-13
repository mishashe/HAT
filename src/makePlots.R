library(plotrix)

Comparison <- "Escherichia_albertii_vs_Escherichia_fergusonii"

rT <- read.table(paste0("/home/msheinman/Development/HAT/data/processed/",Comparison,".r"), 
                 header = FALSE, sep = "\t",row.names=NULL)

L <- read.table(paste0("/home/msheinman/Development/HAT/data/processed/",Comparison,".L"), 
                 header = FALSE, sep = "\t",row.names=NULL)
colnames(rT) <- c("r","Freq")

d <- 0.1
rV <- c(seq(0.5,10.5,1),10^seq(log10(11.5),1*d+log10(max(rT$r)),d))
p <- weighted.hist(x=rT$r,w=rT$Freq,breaks=rV,plot=FALSE)
while(sum(p$counts[2:(length(p$counts)-1)]==0)>0 & d<0.2)
{
  d <- d*1.1
  rV <- c(seq(0.5,10.5,1),10^seq(log10(11.5),1*d+log10(max(rT$r)),d))
  p <- weighted.hist(x=rT$r,w=rT$Freq,breaks=rV,plot=FALSE)
}

K <- 3e10
N <- 50
rho <- 10
rho12 <- 0.002
rV <- p$mids
mTh <- ifelse(rV<3/2*N*rho,2*K*rho/(rV+rho)^3,ifelse(rV<3/2*N*rho*rho/rho12,2*K*N*rho^2/rV^4,2*K*rho12/rV^3))
mTh0 <- 2*K*rho/(rV+rho)^3
knee <- 10
pdf(paste0("/home/msheinman/Development/HAT/plots/",Comparison,"_m.pdf"))
plot(log10(p$mids),log10(p$counts/diff(p$breaks)), pch = 19)
lines(log10(rV),log10(mTh))
lines(log10(rV),log10(mTh0),col="red")
lines(log10(rV),18-5*log10(rV),col="blue")

plot((p$mids),log10(p$counts/diff(p$breaks)), pch = 19)

dev.off()
