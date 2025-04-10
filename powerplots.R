# Code Description: Simulation Study for Variable Selection in High-Dimensional Data
# 
# This R script conducts a simulation study to compare the performance of a new 
# variable selection method, Variable Selection based on PCA (VSPCA), with 
# existing methods, including:
#   
#   Moderated t-test (MTT) from the limma package,
#   Lasso Regression using the glmnet package,
#   Significance Analysis of Microarrays (SAM) from the samr package.
# 
# The study evaluates these methods under different conditions of sample size, 
# correlation structure, and proportion of variance explained for VSPCA only. 
# The primary objective is to assess the statistical power and false discovery 
# rates of these methods for detecting differentially expressed genes.

# Set working directory (update this path as needed)
setwd("C:\\VSPCA\\NewEXP\\UnequalVar")

# Load required packages for data simulation and statistical testing
library(MASS)    # For generating multivariate normal data
library(limma)   # For Moderated t-test (MTT)
library(glmnet)  # For Lasso regression
library(samr)    # For Significance Analysis of Microarrays (SAM)

# Define simulation parameters
nb <- 6          # Number of independent blocks
veb <- 500       # Number of variables/genes per block
p <- veb * nb    # Total number of variables/genes across all blocks

# Define correlation structures for covariance matrices
rho <- c(0, 0.5, 0.8)  # Levels of correlation
rhon <- c(0, 5, 8)     # Labels for output files
cntrho <- 1:length(rho)

# Define proportion of variance explained in PCA-based adjustment
phi <- seq(0.2, 0.8, 0.2)  # Variance thresholds
phin <- c(2, 4, 6, 8)       # Labels for output files
cntphi <- 1:length(phi)

# Loop over different correlation structures
for (rhol in cntrho) {
  
  # Define an exchangeable covariance structure for both groups
  Sigma <- matrix(rho[rhol], veb, veb)  # Create correlation matrix
  diag(Sigma) <- 1  # Set diagonal elements to 1 (ensuring variance normalization)
  
  Sigmax =Sigma + diag(c(rgamma(veb,1,4))) # Covariance for control group (X)
  Sigmay =Sigma + diag(c(rgamma(veb,1,4))) # Covariance for case group (Y)
  
  # Expand the covariance matrices for all blocks
  Sigmax <- kronecker(diag(nb), Sigmax)
  Sigmay <- kronecker(diag(nb), Sigmay)
  
  # Define significance level for hypothesis testing
  alpha <- 0.01
  NLV <- 0  # Null hypothesis mean level
  
  # Generate multivariate normal data for control group (Y_m)
  Y_m <- mvrnorm(n = 10000, mu = rep(NLV, p), Sigma = Sigmax)
  
  # Generate differential expression effects
  mu <- rep(0, p)  # Initialize all mean effects to 0 (null)
  ones <- sample(1:p, 600, replace = FALSE)  # Select 600 differentially expressed variables
  mu[ones] <- 1  # Assign differential expression
  
  # Define effect size categories
  c0 <- 0
  c1 <- 0.6
  c2 <- 0.9
  c3 <- 1.2
  
  # Assign effect sizes randomly to differentially expressed variables
  mu <- NLV + mu * sample(c(-c3, -c2, -c1, c1, c2, c3), p, replace = TRUE)
  
  # Generate multivariate normal data for the case group (Y_p)
  Y_p <- mvrnorm(n = 10000, mu = mu, Sigma = Sigmay)
  
  # Loop over different sample sizes
  for (n1 in c(10, 20, 30, 50, 100, 200)) {  # Control group sample sizes
    for (n2 in c(10, 20, 30)) {             # Case group sample sizes
      for (phil in cntphi) {                # PCA variance thresholds
        
        # Repeat simulation 1000 times
        res <- replicate(1000, {
          
          # Randomly sample observations for both groups
          Y_minus <- Y_m[sample(1:10000, n1, replace = FALSE), ]
          Y_plus <- Y_p[sample(1:10000, n2, replace = FALSE), ]
          
          # Standardize datasets
          X_minus <- scale(Y_minus)
          X_plus <- scale(Y_plus, center = attr(X_minus, "scaled:center"), 
                          scale = attr(X_minus, "scaled:scale"))
          
          # Perform Singular Value Decomposition (SVD)
          SVDec <- svd(X_minus)
          V <- SVDec$v
          
          # Determine number of principal components explaining phi[phil] variance
          nc <- 1
          while (cumsum(sum(SVDec$d[1:nc]) / sum(SVDec$d)) < phi[phil]) {
            nc <- nc + 1
          }
          
          # Project case group onto the control group PCA space
          X_plus_tilde <- X_plus %*% V[, 1:nc] %*% t(V[, 1:nc])
          err <- (X_plus - X_plus_tilde)  # Compute reconstruction errors
          tjplus <- apply(err, 2, mean)   # Compute mean error per variable
          
          # Compute p-values using FDR correction
          pvalues <- p.adjust(2 * pnorm(abs(tjplus / mad(tjplus)), 
                                        lower.tail = FALSE), method = "fdr")
          
          # Perform Moderated t-test (MTT)
          YY <- t(rbind(Y_minus, Y_plus))  # Combine data for MTT
          design <- cbind(Grp1 = 1, Grp2vs1 = rep(c(0, 1), c(n1, n2)))  # Design matrix
          fit <- lmFit(YY, design)
          efit <- eBayes(fit)
          de.table <- topTable(efit, coef = 2, number = Inf, sort.by = "none", 
                               adjust.method = "fdr")
          
          # Perform Lasso regression for variable selection
          X_combined <- rbind(Y_minus, Y_plus)
          y_combined <- c(rep(0, n1), rep(1, n2))
          lasso_fit <- cv.glmnet(X_combined, y_combined, family = "binomial", 
                                 alpha = 1, standardize = TRUE)
          selected_genes <- which(coef(lasso_fit, s = "lambda.min")[-1] != 0)
          
          # Perform Significance Analysis of Microarrays (SAM)
          data <- list(x = YY, y = rep(c(1, 2), c(n1, n2)), geneid = as.character(1:nrow(YY)),
                       genenames = paste("g", as.character(1:nrow(YY)), sep = ""), logged2 = TRUE)
          samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 100)
          pvsam <- p.adjust(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar), method = "fdr")
          
          # Compute detection power for different methods across different effect sizes
          c(
            # Effect size category c0 (smallest effect)
            mean(pvalues[abs(mu) %in% (c0 + NLV)] < alpha),  # VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c0 + NLV)] < alpha),  # MTT
            ifelse(length(selected_genes) == 0, 0, 
                   mean(which(abs(mu) %in% (c0 + NLV)) %in% selected_genes)),  # Lasso
            mean(pvsam[abs(mu) %in% (c0 + NLV)] < alpha),  # SAM
            
            # Effect size category c1
            mean(pvalues[abs(mu) %in% (c1 + c(NLV, -NLV))] < alpha),  # VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c1 + c(NLV, -NLV))] < alpha),  # MTT
            ifelse(length(selected_genes) == 0, 0, 
                   mean(which(abs(mu) %in% (c1 + c(NLV, -NLV))) %in% selected_genes)),  # Lasso
            mean(pvsam[abs(mu) %in% (c1 + c(NLV, -NLV))] < alpha),  # SAM
            
            # Effect size category c2
            mean(pvalues[abs(mu) %in% (c2 + c(NLV, -NLV))] < alpha),  # VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c2 + c(NLV, -NLV))] < alpha),  # MTT
            ifelse(length(selected_genes) == 0, 0, 
                   mean(which(abs(mu) %in% (c2 + c(NLV, -NLV))) %in% selected_genes)),  # Lasso
            mean(pvsam[abs(mu) %in% (c2 + c(NLV, -NLV))] < alpha),  # SAM
            
            # Effect size category c3 (largest effect)
            mean(pvalues[abs(mu) %in% (c3 + c(NLV, -NLV))] < alpha),  # VSPCA
            mean(de.table$adj.P.Val[abs(mu) %in% (c3 + c(NLV, -NLV))] < alpha),  # MTT
            ifelse(length(selected_genes) == 0, 0, 
                   mean(which(abs(mu) %in% (c3 + c(NLV, -NLV))) %in% selected_genes)),  # Lasso
            mean(pvsam[abs(mu) %in% (c3 + c(NLV, -NLV))] < alpha)  # SAM
          )
        })
        
        # Save the output as a PNG file with relevant parameter labels
        png(file=paste("n1", n1, "n2", n2, "rho", rhon[rhol], "phi", 
                       phin[phil], "_UnEqVar.png", sep=""))
        
        # Adjust margins for the plot layout
        par(mar=c(5, 5, 2, 1))
        
        # Define colors for different methods
        bcol <- rep(c("gray", "red", "lightblue", "darkgoldenrod1"), 4)
        
        # Create a boxplot to compare detection power across methods
        boxplot(res ~ row(res), xlab="", ylab="", main="", cex.axis=1.5, las=1, ylim=c(0,1), 
                xaxt='n', col=bcol, at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19))
        
        # X-axis labels for different effect size categories (|delta| values)
        axis(1, at=c(2.5,7.5,12.5,17.5), labels=c(c0, c1, c2, c3), line=0, tck=-.01, cex.axis=1.5)
        
        # X-axis ticks without labels for better alignment
        axis(1, at=c(1,2,3,4), labels=NA, line=-0.1, tck=0.01, cex.axis=1)
        axis(1, at=c(6,7,8,9), labels=NA, line=-0.1, tck=0.01, cex.axis=1)
        axis(1, at=c(11,12,13,14), labels=NA, line=-0.1, tck=0.01, cex.axis=1)
        axis(1, at=c(16,17,18,19), labels=NA, line=-0.1, tck=0.01, cex.axis=1)
        
        # Add axis labels
        title(xlab=expression(paste("|", delta, "|")), cex.lab=1.5, line=3.5)
        title(ylab="Power", cex.lab=1.5, line=3.5)
        
        # Add a title showing parameter values dynamically
        title(main=bquote(paste("(", n[1], ", ", n[2], ", ", rho, ", ", phi, ")", " = ", 
                                "(", .(n1), ", ", .(n2), ", ", .(rho[rhol]), ", ", .(phi[phil]), ")")), 
              cex.main=1.7)
        
        # Add a legend explaining the color coding for different methods
        legend("topleft", 
               legend=c("VSPCA", "MTT", "Lasso", "SAM"), 
               horiz=FALSE, cex=1.2, bg="white", box.lty=1, fill=unique(bcol))
        
        # Close the PNG device to save the file
        dev.off()
      }
    }
  }
}


