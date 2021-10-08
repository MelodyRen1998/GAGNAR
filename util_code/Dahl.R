###########################################################################
## Dahl's method to summarize the samples from the MCMC
###########################################################################

getDahl <- function(out, burn)
{
  #===========================================================#
  ## Input: out = the result from GWCRP ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = clustering configuration, n by 1 vector ##
  ##         sigma2out = estimate sigma2, K by 1 vector ##
  ##         thetaout = estimate theta, K*d matrix ##
  #===========================================================#
  
  iters <- out$Iterates[-(1:burn)]
  # n <- length(iters[[1]][[1]])
  niters <- length(iters)
  # membership matrix
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x$zout
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  # square error
  SqError <- sapply(membershipMatrices, function(x, av) {sum((x - av)^2)},
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

###########################################################################
## Reorder the zout and betaout from getDahl function, to make zout increasing
###########################################################################

reorderDahl <- function(x)
{
  #===========================================================#
  ## Input: x = the output of getDahl function ##
  
  ## Output: 
  ##         zout = clustering configuration, n by 1 vector ##
  ##         sigma2out = estimate sigma2, K by 1 vector ##
  ##         thetaout = estimate theta, K*d matrix ##
  #===========================================================#
  
  zout <- x$zout
  sigma2out <- x$sigma2out
  thetaout <- x$thetaout
  
  for (i in 1:length(unique(x$zout))){  # loop for each cluster
    zout[x$zout == unique(x$zout)[i]] <- i
    sigma2out[i] <- x$sigma2out[unique(x$zout)[i]]
    thetaout[i, ] <- x$thetaout[unique(x$zout)[i], ]
  }
  Dahlout <- list(zout=zout, sigma2out=sigma2out, thetaout=thetaout)
  attr(Dahlout, "iterIndex") <- attr(x, "iterIndex")
  attr(Dahlout, "burnin") <- attr(x, "burnin")
  return(Dahlout)
}