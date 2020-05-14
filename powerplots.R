
setwd("C:\\Users\\ullahi\\OneDrive - Queensland University of Technology\\Desktop\\plos-latex-template")
setwd("C:\\VSPCA")
# Required packages for generating data and performing MTT
library(MASS)
library(limma)
# install.packages("devtools")
# devtools::install_github("andrewthomasjones/EMMIXcontrasts")
library(EMMIXcontrasts2)

# number of blocks 
nb <- 6
# number of variables/genes in per block
veb <- 500
# total number of variables/genes
p <- veb*nb

##################################################################
##################################################################

rho <- c(0,.5,.8)
rhon <- c(0,5,8)
cntrho <- 1:length(rho)
phi <- seq(0.2,0.8,0.2)
phin <- c(2,4,6,8)
cntphi <- 1:length(phi)
for (rhol in cntrho){
  # I use exchangeable covariance structure for both groups
  Sigma <- matrix(rho[rhol],veb,veb);diag(Sigma)=1;
  # Sigmax =Sigma + diag(c(rgamma(p,1,1)))
  # Sigmay=Sigma + diag(c(rgamma(p,1,1)))
  #Sigma <- matrix(rep(1,3))%*%t(matrix(rep(1,3)))*rho+diag(3)*(1-rho)
  
  Sigmax  <- Sigmay <- kronecker(diag(nb), Sigma)
  
  # level of significance 
  alpha <- 0.01
  #null mean
  NLV <- 0
  Y_m <- mvrnorm(n=10000, mu= rep(NLV,p), Sigma=Sigmax)
  
  # below
  # shifts sizes are randomly chosen for 600 of variables
  mu <- rep(0,p)
  ones <- sample(1:p, 600, replace=F)
  mu[ones]=1
  
  # different levels of shift
  
  c0=0 # No shft. Null is true
  c1=0.6 # smallest shift
  c2=0.9 # middle level shift
  c3=1.2 # largest shift
  
  mu <- NLV + mu*sample(c(-c3,-c2,-c1,c1,c2,c3), p, replace=T)
  
  #mu <- NLV + rbinom(p,1,prob=.2)*sample(c(-c3,-c2,-c1,c1,c2,c3), p, replace=T)
  Y_p <- mvrnorm(n=10000, mu=mu, Sigma=Sigmay)
  # We repeat 1000 independent experiment of similar nature
  
  for (n1 in c(10,20,30,50,100)){
    for (n2 in c(10,20,30)){
      for (phil in cntphi) {
        
        res <- replicate(1000,{
          
          # generate n1 X p matrix of data. This is the training data
          Y_minus <- Y_m[sample(1:10000, n1, replace=F), ] #mvrnorm(n=n1, mu= rep(0,p), Sigma=Sigmax)
          
          # generate n2 X p matrix of data. This is the cases data.
          Y_plus <- Y_p[sample(1:10000, n2, replace=F), ] #mvrnorm(n=n2, mu=mu, Sigma=Sigmay)
          
          # standardize columns of training data
          X_minus <- scale(Y_minus)
          # standardize columns of cases data using the mean 
          # and standard deviations of the columns of training data.
          # R is clever! It stores means and sds for you in X_minus 
          X_plus <- scale(Y_plus, center = attr(X_minus,"scaled:center"), 
                          scale=attr(X_minus,"scaled:scale"))
          
          # perform SVD using training data
          SVDec <- svd(X_minus)
          V <- SVDec$v
          # scree plot to decide how many components to use to recover cases data.
          # plot(SVDec$d)
          
          # number of components to use to recover cases data.
          # We specify phi=.80
          nc=1
          while (cumsum(sum(SVDec$d[1:nc])/sum(SVDec$d)) < phi[phil]){
            nc=nc+1
          }
          #nc
          
          X_plus_tilde <- X_plus%*%V[,1:nc]%*%t(V[,1:nc])
          # induced error
          err=(X_plus-X_plus_tilde)
          # average error 
          tjplus <- apply(err,2,mean) # assumed to follow a normal distribution
          pvalues <- p.adjust(2*pnorm(abs(tjplus/mad(tjplus)),
                                      lower.tail = F), method = "fdr")
          
          # perform MTT
          YY <- t(rbind(Y_minus, Y_plus))
          design <- cbind(Grp1=1,Grp2vs1=rep(c(0,1),c(n1,n2)))
          fit <- lmFit(YY,design)
          efit <- eBayes(fit)
          de.table <- topTable( efit, coef=2, number=Inf, sort.by="none", 
                                adjust.method="fdr" )
          
          # finaly results
          c(mean(pvalues[abs(mu) %in% (c0+NLV)] < alpha), # null with VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c0+NLV)] < alpha), # null with MTT
            mean(pvalues[abs(mu) %in% (c1+c(NLV,-NLV))] < alpha), # identify shift c1 with VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c1+c(NLV,-NLV))] < alpha), # identify shift c1 with MTT
            mean(pvalues[abs(mu) %in% (c2+c(NLV,-NLV))] < alpha), # identify shift c2 with VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c2+c(NLV,-NLV))] < alpha), # identify shift c2 with MTT
            mean(pvalues[abs(mu) %in% (c3+c(NLV,-NLV))] < alpha), # identify shift c3 with VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c3+c(NLV,-NLV))] < alpha)) # identify shift c3 with MTT
        })
        
        # compare results for different levels of shifts
        png(file=paste("n1",n1,"n2",n2,"rho",rhon[rhol],"phi",phin[phil],"a.png",sep = ""))
        par(mar=c(5, 5, 2, 1))
        bcol=rep(c("gray","red"), 4)
        boxplot(res~row(res), xlab="", ylab="", main="", cex.axis=1.5, las=1, ylim=c(0,1), 
                xaxt='n', col=bcol, at=c(1,2,4,5,7,8,10,11))
        axis(1, at=c(1.5,4.5,7.5,10.5),labels=c(c0,c1,c2,c3),line=0, tck=-.01, cex.axis=1.5)
        
        axis(1,at=c(1,2),labels=NA,line=-0.1,tck=0.01, cex.axis=1)
        axis(1,at=c(4,5),labels=NA,line=-0.1,tck=0.01, cex.axis=1)
        axis(1,at=c(7,8),labels=NA,line=-0.1,tck=0.01, cex.axis=1)
        axis(1,at=c(10,11),labels=NA,line=-0.1,tck=0.01, cex.axis=1)
        title(xlab=expression(paste("|", delta, "|")), cex.lab=1.5, line = 3.5)
        title(ylab="power", cex.lab=1.5, line=3.5)
        title(main=bquote(paste("(",n[1],", ",n[2], ", ",rho, ", ",phi, ")", " = ","(",.(n1),", ",.(n2), ", ",.(rho[rhol]),", ",.(phi[phil]),")" )), cex.main=1.7)
        legend("topleft", 
               legend=c("VSPCA", "MTT"), horiz = F , 
               cex = 1.2, bg = "white",box.lty = 1,
               fill = unique(bcol))
        
        dev.off()
        
      }
    }
  }
}