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
train_window = 12
test_window = 5
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
# h <- c(seq(0,2,0.2),3,4,5)
h <- c(2)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)

tryCatch({
  
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
  save(real_res, file = "./report/city_res_12_5.rda")
  
  besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
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
} , error = function(err) {
  print(paste0("ERROR: ",err))
})


# analysis for four cities
ZZ_best[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "北京市")])]
p_value_node[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "北京市")]), ]
ZZ_best[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "临沂市")])]
p_value_node[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "临沂市")]), ]
ZZ_best[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "抚顺市")])]
p_value_node[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "抚顺市")]), ]
ZZ_best[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "泰州市")])]
p_value_node[which(rownames(Adjmat) == dijishiID$市代码[which(dijishiID$市 == "泰州市")]), ]
