source("./code/gagnar.R")
source("./code/Dahl.R")
source("./code/LPML_VAR.R")

# =========
# real data function 
# =========
runReal <- function(Y_tilde, X_tilde, dd, hyperpara, niterations, 
                    a0, b0, tau0, sigma0, initNClusters, verbose = F){
  simresult <- list()
  for (pp in 1:nrow(hyperpara)){
    cat("[ h:", hyperpara$h[pp], "alpha: ", hyperpara$alpha[pp], "] \n")
    set.seed(hyperpara$h[pp]*1000)
    out <- GAGNAR(Y_tilde, X_tilde, dd, hyperpara$alpha[pp], hyperpara$h[pp], 
                 niterations, a0, b0, tau0, sigma0, initNClusters, verbose)
    # out$Iterates[[100]]$thetaout
    Dahlout <- reorderDahl(getDahl(out, burn=burnin))
    LPML <- LPML(Y_tilde, X_tilde, list(Iterates = out$Iterates[(burnin+1):niterations]))
    simresult[[pp]] <- list(alpha = hyperpara$alpha[pp], h = hyperpara$h[pp], out = out, Dahlout = Dahlout, LPML = LPML)
  }
  return(simresult)
}


NAR.predict <-function(Ymat, W, theta, Z)
{
  # demean the responses
  Ymat1 = Ymat 
  N = nrow(Ymat); Time = ncol(Ymat)
  
  # covariates
  WYmat = W%*%Ymat1[,-Time]
  Ylag = as.vector(Ymat1[,-Time])
  zz = do.call(rbind, rep(list(Z), Time-1))
  X = cbind(rep(1, N*(Time-1)), as.vector(WYmat), Ylag, zz)
  
  # response
  Yvec = as.vector(Ymat1[, -1])
  
  # residuals
  resi1 = Yvec - X%*%theta
  
  return(list(rmse = sqrt(mean(resi1^2)),
              mese = median(resi1^2),
              resi_mat = matrix(resi1, nrow = N),
              pred = X%*%theta))
}

NAR.est<-function(Ymat, W, Z)
{
  # demean the responses
  Ymat1 = Ymat 
  N = nrow(Ymat); Time = ncol(Ymat)
  
  # covariates
  WYmat = W%*%Ymat1[,-Time]
  Ylag = as.vector(Ymat1[,-Time])
  zz = do.call(rbind, rep(list(Z), Time-1))
  X = cbind(rep(1, N*(Time - 1)), as.vector(WYmat), Ylag, zz)
  
  # response
  Yvec = as.vector(Ymat1[, -1])
  XX_inv = solve(crossprod(X))
  # NAR estimate
  theta = XX_inv%*%(t(X)%*%Yvec)
  
  # sigma estimate
  sig2 = as.numeric(t(Yvec) %*% (diag(1, N*(Time-1)) - X %*% XX_inv %*% t(X)) %*% Yvec) / N
  
  # cov of coef
  cov = sig2 * XX_inv
  
  return(list(theta = theta, mu = rowMeans(Ymat), sig2 = sig2, cov = cov))
}
