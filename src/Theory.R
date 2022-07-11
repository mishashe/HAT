getMatrixR <- function(rateMatrix)
{
  N <- dim(rateMatrix)[1]
  I <- 1:N
  J <- 1:N
  Q <- list()
  q <- 0
  for (i in I)
  {
    for (j in J)
    {
      if (i < j )
      {
        q <- q+1
        Q[[q]] <- c(i,j)
      }
    }
  }
  R <- matrix(0,nrow=N*(N-1)/2,ncol=N*(N-1)/2)
  
  for (q1 in 1:length(Q))
  {
    i1 <- Q[[q1]][1]
    j1 <- Q[[q1]][2]
    for (q2 in 1:length(Q))
    {
      i2 <- Q[[q2]][1]
      j2 <- Q[[q2]][2]
      if (q1==q2) 
      {
        R[q1,q2] <- -2*rateMatrix[i1,j1]
        for (k in 1:N)
        {
          if (k != i1 & k!=j1)
          {
            R[q1,q2] <- R[q1,q2] - rateMatrix[i1,k] - rateMatrix[j1,k]
          }
        }
      }
      else if (i1==i2) R[q1,q2] <- rateMatrix[j1,j2]
      else if (j1==j2) R[q1,q2] <- rateMatrix[i1,i2]
      else if (i1==j2) R[q1,q2] <- rateMatrix[j1,i2]
      else if (j1==i2) R[q1,q2] <- rateMatrix[i1,j2]
    }
  }
  return(R)
}


getRhoQ <- function(rateMatrix)
{
  N <- dim(rateMatrix)[1]
  I <- 1:N
  J <- 1:N
  Q <- list()
  q <- 0
  rhoQ <- data.frame()
  for (i in I)
  {
    for (j in J)
    {
      if (i < j )
      {
        q <- q+1
        Q[[q]] <- c(i,j)
        rhoQ <- rbind(rhoQ, data.frame(i=i,j=j,rho=rateMatrix[i,j]))
      }
    }
  }
  return(rhoQ)
}

getEMrate <- function(rateMatrix)
{
  N <- dim(rateMatrix)[1]
  ToOptim <- function(rho)
  {
    sumEM <- 0
    for (i in 1:(N-1))
    {
      for (j in (i+1):N)
      {
        sumEM <- sumEM + 2.0/((N-2)*rho+2*rateMatrix[i,j])/(N-1)
      }
    }
    return((1/sumEM-rho)^2)
  }

  rho <- mean(rateMatrix)*10
  return(optim(rho,fn=ToOptim,lower = 0.0001, upper = 100*rho,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                 "Brent")[4],control=c(reltol=1e-20,abstol=1e-20,maxit=1e5,ndeps=1e-6))$par)
}


library(Matrix)
N <- 4
meanlog <- log(0.1)
sdlog <- log(50)
rateMatrix <- matrix(0,nrow=N,ncol=N)

# Random matrix
# for (i in 1:(N-1))
# {
#   for (j in (i+1):N)
#   {
#     rateMatrix[i,j] <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
#     rateMatrix[j,i] <- rateMatrix[i,j]
#   }
# }

rateMatrix[1,2] <- 1
rateMatrix[2,3] <- 1
rateMatrix[3,4] <- 1
rateMatrix[2,1] <- 1
rateMatrix[3,2] <- 1
rateMatrix[4,3] <- 1

R <- getMatrixR(rateMatrix)

RhoQ <- getRhoQ(rateMatrix)
tauQ <- 4*solve(R %*% R) %*% as.matrix(RhoQ$rho)
# plot(log10(1/RhoQ$rho),log10(tauQ))

K <- 1e5
rV <- 10^(seq(0,5,0.1))
m <- matrix(0,nrow=nrow(RhoQ),ncol=length(rV))
library(Matrix)
for (ir in 1:length(rV))
{
  foo <- solve(diag(rV[ir],dim(R)[1])-R/2)
  m[,ir] <- 2*(diag(K,dim(R)[1])-R/2) %*% foo %*% foo %*% foo %*% as.matrix(RhoQ$rho)[,1]
}


rho <- getEMrate(rateMatrix)
library(latex2exp)
pdf("/home/misha/Documents/Development/GeneralTheoryHGT/plots/rateHist.pdf",width=3,height=3.5)
p <- hist(log10(RhoQ$rho),20,plot=FALSE)
plot(p$mids,p$density,type="o", xlab = "", ylab = "",cex.axis=0.5)
foo <- max(p$density)
p <- hist(log10(1/tauQ),20,plot=FALSE)
points(p$mids,p$density/max(p$density)*foo,col="red",type="l", lwd = 2, lty = 1)
axis(1, at = log10(rho),cex.axis=0.75,labels=TeX('$\\rho$'))
dev.off()


pdf("/home/misha/Documents/Development/GeneralTheoryHGT/plots/mEM.pdf",width=3,height=3.5)
plot(log10(rV),log10(2*K*rho/(rho+rV)^3),col="grey",type="l",lwd=5,ylim=c(-12,4),xlab=TeX('$log_{10} r$'), xaxt = "n", yaxt = "n", ylab=TeX('$log_{10} m(r)$'))
axis(1, at = 0:5,cex.axis=0.5)
axis(2, at = seq(-12,4,2),cex.axis=0.5, labels = seq(-12,4,2))
points(log10(rV),log10(m[which.min(RhoQ$rho),]),pch=".",col="blue",cex=3)
rho12 <- RhoQ$rho[which.min(RhoQ$rho)]
mth <- ifelse(rV<3/2*N*rho,2*K*rho/(rV+rho)^3,ifelse(rV<3/2*N*rho*rho/rho12,3*K*N*rho^2/rV^4,2*K*rho12/rV^3))
lines(log10(rV),log10(mth),col="blue")
points(log10(rV),log10(m[which.max(RhoQ$rho),]),pch=".",col="red",cex=3)
rho12 <- RhoQ$rho[which.max(RhoQ$rho)]
mth <- ifelse(rV<rho12*(N*rho^2/rho12^2/2)^(1/3),N*rho/rho12*K*rho/(rV+rho)^3,2*K*rho12/(rV+rho12)^3)
lines(log10(rV),log10(mth),col="red")
points(log10(rV),log10(m[which.min(abs(RhoQ$rho-rho)),]),pch=".",col="black",cex=3)
dev.off()


pdf("/home/misha/Documents/Development/GeneralTheoryHGT/plots/mEM14.pdf",width=3,height=3.5)
plot(log10(rV),log10(2*K*rho/(rho+rV)^3),col="white",type="l",lwd=5,ylim=c(-20,4),xlab=TeX('$log_{10} r$'), xaxt = "n", yaxt = "n", ylab=TeX('$log_{10} m(r)$'))
axis(1, at = 0:5,cex.axis=0.5)
axis(2, at = seq(-18,4,2),cex.axis=0.5, labels = seq(-18,4,2))
points(log10(rV),log10(m[3,]),pch=".",col="black",cex=3)
lines(log10(rV[20:length(rV)]),6.5-3*log10(rV[20:length(rV)]),col="red")
points(log10(rV),log10(m[1,]),pch=".",col="red",cex=3)
lines(log10(rV[20:length(rV)]),7.5-5*log10(rV[20:length(rV)]),col="black")
dev.off()




