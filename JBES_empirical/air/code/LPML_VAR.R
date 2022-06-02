###########################################################################
## Calculate LPML given posterior samples after burnin
###########################################################################

LPML <- function(Y_tilde, X_tilde, out) 
{
  #===========================================================#
  ## Input: 
  ##        Y_tilde, T*N matrix ##
  ##        X_tilde, T*N*(p+3) array ##
  ##        out = the result from GAGNAR      ##
  
  ## Output:
  ##         LPML = LPML value ##
  #===========================================================#
  
  LPML <- 0.0
  
  n.nodes <- dim(Y_tilde)[2]
  d <- dim(X_tilde)[3]  # dimension of theta
  t <- dim(X_tilde)[1]  # time length
  
  for (i in 1:n.nodes) {  # loop for each node
    Yi <- Y_tilde[2:t,i]  # [2:T, 1]
    Xi <- X_tilde[2:t,i,]  # [2:T, d]
    lik <- lapply(out$Iterates, function(y) { # loop for each iteration result
      sigma2out <- y$sigma2out
      thetaout <- y$thetaout
      clu <- y$zout  # this is a cluster assignment
      clu_factor <- as.factor(clu)
      zout <- model.matrix(~clu_factor - 1)  # no intercept
      
      # cat("dimension of zout, sigma2out and thetaout:", dim(zout), "\n", dim(sigma2out), "\n", dim(thetaout), "\n")
      
      sigma2zi <- zout %*% sigma2out
      thetazi <- zout %*% thetaout
      ratio <- zout %*% (table(clu)/n.nodes)
      ll <- exp((t(Yi - Xi %*% thetazi[i, ]) %*% (Yi - Xi %*% thetazi[i, ]))/(2*sigma2zi[i, ])) *
        (( ratio[i] / sqrt(2*pi* (sigma2zi[i, ])))^(t-1))
      # ll <- (- (t-1)*log(2*pi* (sigma2zi[i, ]))/2 -  
      #          (t(Yi - Xi %*% thetazi[i, ]) %*% (Yi - Xi %*% thetazi[i, ]))/(2*sigma2zi[i, ]))
      return(1/ll)
    })
    # cat("1/likelihood for node", i, ":", unlist(lik), "\n")
    log_cpo_i <- log(length(out$Iterates)) - log(sum(unlist(lik)))
    # cat(log_cpo_i,"\n")
    LPML <- LPML + log_cpo_i
  }
  return(LPML)
}

