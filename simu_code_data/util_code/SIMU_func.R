# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("igraph", "mvtnorm", "MCMCpack", "parallel", "mclust", "plyr", "fcd")
ipak(packages)

source("./util_code/gwcrp.R")
source("./util_code/Dahl.R")
source("./util_code/LPML_VAR.R")
source("./util_code/estimator.R")

sim <- function(isim, Y_tilde_list, X_tilde_list, d_matrix, hyperpara,
                niterations, a0, b0, tau0, sigma0, initNClusters, burnin){
  Y_tilde <- Y_tilde_list[[isim]]
  X_tilde <- X_tilde_list[[isim]]
  dd <- d_matrix[[isim]]
  
  simresult <- list()
  for (pp in 1:nrow(hyperpara)){
    cat("[ h:", hyperpara$h[pp], "alpha: ", hyperpara$alpha[pp], "] \n")
    set.seed(hyperpara$h[pp]*1000 + isim)
    out <- GWCRP(Y_tilde, X_tilde, dd, hyperpara$alpha[pp], hyperpara$h[pp], 
                 niterations, a0, b0, tau0, sigma0, initNClusters, verbose = T)
    # out$Iterates[[100]]$thetaout
    Dahlout <- reorderDahl(getDahl(out, burn=burnin))
    LPML <- LPML(Y_tilde, X_tilde, list(Iterates = out$Iterates[(burnin+1):niterations]))
    simresult[[pp]] <- list(alpha = hyperpara$alpha[pp], h = hyperpara$h[pp], out = out, 
                            Dahlout = Dahlout, LPML = LPML)
  }
  return(simresult)
}


# ==========
# simulation with best h
# ==========

sim_best_h = function(isim, best_h, Y_tilde_list, X_tilde_list, d_matrix,
                      niterations, a0, b0, tau0, sigma0, initNClusters, burnin) {
  Y_tilde <- Y_tilde_list[[isim]]
  X_tilde <- X_tilde_list[[isim]]
  dd <- d_matrix[[isim]]
  set.seed(round(best_h))
  out <- GWCRP(Y_tilde, X_tilde, dd, alpha = 1, best_h, 
               niterations, a0, b0, tau0, sigma0, initNClusters, verbose = T)
  # out$Iterates[[100]]$thetaout
  Dahlout <- reorderDahl(getDahl(out, burn=burnin))
  LPML <- LPML(Y_tilde, X_tilde, list(Iterates = out$Iterates[(burnin+1):niterations]))
  simresult_besth <- list(out = out, Dahlout = Dahlout, LPML = LPML)
  return(simresult_besth)
}


# ==========
# EM algorithm
# ==========
simEM <- function(isim, kk, Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list, niterations) {
  Y_tilde <- Y_tilde_list[[isim]]
  X_tilde <- X_tilde_list[[isim]]
  dd <- d_matrix[[isim]]
  WW = W_list[[isim]]
  t = nrow(Y_tilde)
  n.nodes = ncol(Y_tilde)
  p = dim(X_tilde)[3] - 3
  
  Ymat <- Y_tilde %>% t() # N*T
  X <- aperm(X_tilde, c(2,1,3))  # by time
  dim(X) <- c(t*n.nodes, (p+3))  # (TN)*p
  Xmat_1 <- X[((t-1)*n.nodes+1):(t*n.nodes),4:(4+p-1)]  # covariates for all nodes
  true_label <- label_list[[isim]]
  ### Estimate the model
  EM_MSE <- c()
  aris <- c()
  SSE_para <- rep(0,5)
  
  em_theta = list()
  em_clu = list()
  em_sigma = list()
  EM_summary = list()
  
  for (iter in 1:niterations) {
    cat("replicates:", isim, "iter:",iter, "\r")
    em_res = EM.NAR(Ymat, WW, Xmat_1, K = kk, seed = iter) # EM estimation
    em_theta[[iter]] <- em_res$theta %>% t()
    em_clu[[iter]] <- em_res$ezK
    em_sigma[[iter]] <- em_res$sigma
    if (iter %% 300 == 0) {
      print(em_theta[[iter]])
    }
    EM_summary[[iter]] = list(theta = em_theta[[iter]], clu = em_clu[[iter]], sigma = em_sigma[[iter]])
  }
  return(EM_summary)
}

# ==========
# GNAR algorithm 
# ==========
simGNAR <- function(isim, kk, Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list) {
  Y_tilde <- Y_tilde_list[[isim]]
  X_tilde <- X_tilde_list[[isim]]
  dd <- d_matrix[[isim]]
  WW = W_list[[isim]]
  t = nrow(Y_tilde)
  n.nodes = ncol(Y_tilde)
  p = dim(X_tilde)[3] - 3
  
  Ymat <- Y_tilde %>% t() # N*T
  X <- aperm(X_tilde, c(2,1,3))  # by time
  dim(X) <- c(t*n.nodes, (p+3))  # (TN)*p
  Xmat_1 <- X[((t-1)*n.nodes+1):(t*n.nodes),4:(4+p-1)]  # covariates for all nodes, N*p, only covariates
  true_label <- label_list[[isim]]
  ### Estimate the model
  GNAR_MSE <- c()
  aris <- c()
  SSE_para <- rep(0,5)
  
  GNAR_theta = list()
  GNAR_clu = list()
  
  GNAR_summary = Cluster.NAR(Ymat, WW, Xmat_1, K = kk, method="kmeans", seed = isim) 
  GNAR_summary$theta = t(GNAR_summary$theta)
  return(GNAR_summary)
}
