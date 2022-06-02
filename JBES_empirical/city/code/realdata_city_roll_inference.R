# =========
# load packages
# =========
setwd("~/renyimeng/GAGNAR_140/city")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("scales", "igraph", "mvtnorm", "MCMCpack", "parallel", "mclust", "plyr")
ipak(packages)

source("./code/gagnar.R")
source("./code/Dahl.R")
source("./code/LPML_VAR.R")
# source("./code/GNAR_code-master/estimator.R")
source("./code/real_estimator.R")
# 自定义GNAR函数
Cluster.NAR<-function(Ymat, W, Z, K, method = "complete", memb = NULL,  
                      plot = F, group = NULL, seed = F)
{
  N = nrow(Ymat)
  Time = ncol(Ymat)
  Ymat1 = W%*%Ymat     
  if (seed)
    set.seed(1234)
  ### obtain the regression dataset
  if (is.null(Z))
    yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                        inter = 1,
                        net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                        lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                        Y = as.vector(t(Ymat[,-1]))) 
  else
    yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                        inter = 1,  # intercept
                        net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                        lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                        Z[rep(1:N, each = Time-1),],
                        Y = as.vector(t(Ymat[,-1])))
  
  ### for each node obtain the estimates for b_i
  paras = ddply(yy_dat, .(id), function(x){
    X = as.matrix(x[,2:4])
    invXX = ginv(crossprod(X))                                                                                   ### the response vector
    thetaEst = invXX%*%colSums(X*x$Y)   
    df = data.frame(matrix(c(thetaEst), nrow = 1))
  })
  
  colnames(paras)[-1] = c("intercept", "network", "momentum")
  
  ### scale the network and momentum parameter and calculate the nodal distances
  para_scale = apply(paras[,3:4], 2, function(x) x/max(abs(x)))#scale(paras[,3:4])
  para_dist = dist(as.matrix(para_scale))
  
  ### conduct the cluster algorithm
  if (method=="kmeans")
  {
    # nar_para = betaOLS(Ymat, W, Z)
    qts = seq(0.1, 0.9, length.out = K)
    for(ss in 1:nseeds){
      ini_theta = apply(para_scale, 2, 
                        function(x) {set.seed(ss); sample(quantile(x, qts))})
      k_res = kmeans(para_scale, centers = ini_theta)
      memb = k_res$cluster
    }
  }
  else if (method == "complete")
  {
    hcl = hclust(para_dist, method = method)
    memb = cutree(hcl, k = K)
  }
  
  alpha = table(memb)/length(memb)
  if (plot==T)
  {
    par(mfrow = c(1,2))
    plot(paras$network, paras$momentum, col = group)
    plot(paras$network, paras$momentum, col = memb)
  }
  yy_dat$group = rep(memb, each = Time - 1) # obtain the group for each node
  
  
  theta_est = matrix(0, nrow = ncol(Z)+3, ncol = K)
  sig2Est = numeric(K)
  covEst = list()
  for (g in 1:max(memb)) {
    x = yy_dat[which(yy_dat$group == g),]
    X = as.matrix(x[,2:(ncol(x)-2)])
    invXX = ginv(crossprod(X))                                                                                   ### the response vector
    theta_est[,g] = invXX%*%colSums(X*x$Y)
    sig2Est[g] = as.numeric(t(x$Y) %*% (diag(1, nrow(X)) - X%*%invXX%*%t(X)) %*% x$Y) / nrow(X)
    covEst[[g]] = sig2Est[g] * invXX
  }
  
  theta_scale = apply(theta_est, 2, function(x) x/max(abs(x)))#scale(paras[,3:4])
  return(list(theta = t(theta_scale), alpha = alpha, group = memb, sig2 = sig2Est, cov = covEst))
}

# =========
# load data
# =========
load("./city_data/Ymat.rda")
load("./city_data/Xmat.rda")
load("./city_data/Adjmat.rda")
cat("network density is", sum(Adjmat)/(nrow(Adjmat)^2 - nrow(Adjmat)), "\n")
load("./city_data/dmat.rda")  # d.matrix

cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]

# =========
# experiment setting
# =========
train_window = 7
test_window = 7
n_window = 1
train_start = 5
n_test = 1
flag = TRUE
test_days = NULL
pred_rmse = matrix(0, nrow = 7, ncol = 100)
p_values = list()
# basic setting & tuning param
niterations <- 1500 # number of iterations for MCMC
initNClusters <- 10
burnin <- 500 # number of burnin iterations for MCMC
# h <- c(0,2) # crp and optimal h
h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)

tryCatch({
  while(flag) {
    train_end = train_start + train_window - 1
    test_start = train_start + 1
    test_end = test_start + test_window - 1
    
    # 1. split train and test set
    Ymat_train <- Ymat[,(train_start:train_end)]
    Ymat_test <- Ymat[,((test_start-1):test_end)]  # use the last of trainset as start point
    
    Ymat_train <- t(Ymat_train)  # T by N
    Ymat_test <- t(Ymat_test)
    
    N = ncol(Ymat_train)
    TT = nrow(Ymat_train)
    
    # 2. build X matrix
    p <- ncol(Zmat)
    X_cov <- array(0, dim = c(train_window, N, (p+3)))  # T by N by p
    for (t in 1:train_window) {
      X_cov[t,,1] <- 1  # constant
      X_cov[t,,2] <- 0  # network effect
      X_cov[t,,3] <- 0  # momentum effect
      X_cov[t,,4:(p+3)] <- Zmat
    }
    
    W <- Adjmat/rowSums(Adjmat)
    W[rowSums(Adjmat)==0,] = 0
    for (t in 2:train_window) {
      X_cov[t,,2] <- W %*% Ymat_train[(t-1),]  # network effect
      X_cov[t,,3] <- Ymat_train[(t-1),]  # network effect
    }
    
    # 3. estimate
    # (1) initialize
    a0 <- 0.01
    b0 <- 0.01
    tau0 <- rep(0, p + 3)
    sigma0 <- diag(100, p + 3)
    # (2) est
    real_res <- runReal(Ymat_train, X_cov, d.matrix, hyperpara, niterations,
                        a0, b0, tau0, sigma0, initNClusters, verbose = T)
    # save(real_res, file = "./report/city_res.rda")
    # besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
    besth_ind <- 2  # index of the optimal h
    best_h <- real_res[[besth_ind]]$h
    # (3) predict
    ZZ <- real_res[[besth_ind]]$Dahlout$zout %>% as.factor()
    ZZ1 <- model.matrix( ~ ZZ - 1)  # N*K
    Theta <- real_res[[besth_ind]]$Dahlout$thetaout  # K*(p+3)
    node_est <- ZZ1 %*% Theta
    Sigma2 <- real_res[[1]]$Dahlout$sigma2out
    
    YY_pred <- matrix(0, nrow = test_window, ncol = N)
    for (tt in 2:(test_window+1)) {
      XX_const <- rep(1,N)
      XX_net <- W %*% Ymat_test[(tt-1),]
      XX_mom <- Ymat_test[(tt-1),]
      XX_covar <- Zmat
      XX <- cbind(XX_const, XX_net, XX_mom, XX_covar)
      YY_t1 <- diag(XX %*% t(ZZ1 %*% Theta))
      YY_pred[tt-1,] <- YY_t1
    }
    
    # ===============
    # inference by GNAR
    # ===============
    
    # 1. split train and test set
    Ymat_train <- Ymat[,(train_start:train_end)]
    Ymat_test <- Ymat[,((test_start-1):test_end)]  # use the last of trainset as start point
    
    TT = ncol(Ymat_train)
    N = nrow(Ymat_train)
    
    # 2. build X matrix
    p <- ncol(Zmat)
    W <- Adjmat/rowSums(Adjmat)
    W[rowSums(Adjmat)==0,] = 0
    
    # 3. estimate
    # (2) est
    memb = as.numeric(ZZ)
    K = max(memb)
    GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K, method="other", memb, seed = n_test)
    Theta = GNAR_summary$theta
    Sigma2 = GNAR_summary$sig2
    p_values_gnar = matrix(0, nrow = K, ncol = p+3)
    for (kk in 1:K) {
      p_values_gnar[kk, ] = (1 - pnorm(abs(Theta[kk, ])/sqrt(diag(GNAR_summary$cov[[kk]]))))*2
    }
    
    p_values[[n_test]] = p_values_gnar

    cat( p_values[[n_test]], "\n")
    n_test = n_test + 1
    train_start = train_start + n_window
    test_days = c(test_days, test_start)
    if ((train_start + train_window - 1) > ncol(Ymat)){flag = F}
  }
} , error = function(err) {
  print(paste0("rolling ERROR: ",err))
  save(p_values, file = "./report/city_p_value_train7_test7.rda")
})
