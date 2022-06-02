# =========
# load packages
# =========
setwd("~/renyimeng/GAGNAR_140/JBES_empirical/stock")
# setwd("../stock")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("scales", "igraph", "mvtnorm", "MCMCpack", "parallel", "mclust", "plyr", "ggplot2")
ipak(packages)

source("./code/gagnar.R")
source("./code/Dahl.R")
source("./code/LPML_VAR.R")
source("./code/GNAR_code-master/estimator.R")
source("./code/real_estimator.R")
# =========
# load data
# =========
load("./stock_data/rda_full_sparse/Ymat.rda")
load("./stock_data/rda_full_sparse/Zmat.rda")
load("./stock_data/rda_full_sparse/Adjmat.rda")
cat("network density is", sum(adj.mt)/(nrow(adj.mt)^2 - nrow(adj.mt)), "\n")
load("./stock_data/rda_full_sparse/dmat.rda")  # d.matrix
Zmat <- scale(Zmat)

# visualize the histogram of stock return
# yy_dat <- data.frame(ID = 1:nrow(Ymat), value = rowMeans(Ymat))
# pdf(width=4, height=4, file = "./report/stk_des_2020.pdf")
# ggplot(yy_dat, aes(x = value)) +
#   geom_histogram(binwidth = 0.004, color = "black", fill = "grey") +
#   theme_bw() +
#   labs(x = "Average Weekly Stock Return Rate", y = "Frequency") +
#   theme(axis.text = element_text(size = 10))
# dev.off()
# mean(yy_dat$value)

# ==========
# experiment setting
# ==========

train_window = 52
test_window = 6
n_window = 1  # rolling window size
train_start = 1
n_test = 1
flag = TRUE
test_days = NULL
pred_rmse = matrix(0, nrow = 7, ncol = 100)


while (flag) {
  train_end = train_start + train_window - 1
  # test_start = train_end
  # test_end = test_start + test_window
  
  cat("GAGNAR", train_start, "\n")
  # 1. split train and test set
  Ymat_train <- Ymat[,(train_start:train_end)]
  # Ymat_test <- Ymat[,(test_start:test_end)]
  
  Ymat_train <- t(Ymat_train)  # T by N
  # Ymat_test <- t(Ymat_test)
  
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
  adj.mt <- as.matrix(adj.mt)
  W <- adj.mt/rowSums(adj.mt)
  W[rowSums(adj.mt)==0,] = 0
  
  for (t in 2:train_window) {
    X_cov[t,,2] <- W %*% Ymat_train[(t-1),]  # network effect
    X_cov[t,,3] <- Ymat_train[(t-1),]  # network effect
  }
  
  # 3. estimate parameters
  # (1) initialize
  a0 <- 0.02
  b0 <- 0.02
  tau0 <- rep(0, p + 3)
  sigma0 <- diag(100, p + 3)
  
  # (2) estimate
  niterations <- 1500 # number of iterations for MCMC
  initNClusters <- 10
  burnin <- 500 # number of burnin iterations for MCMC
  h <- seq(0,2,0.4)
  # h <- c(0.6)  # crp and optimal h full & sparse > 3
  alpha <- c(0.1)
  hyperpara <- expand.grid(h = h, alpha = alpha)
  res_topstk <- runReal(Ymat_train, X_cov, d.matrix, hyperpara, niterations,
                        a0, b0, tau0, sigma0, initNClusters, verbose = T)
  if(train_start == 1) {
    save(res_topstk, file = "./report/stock_res_supp_full_train52.rda")
    break
    }
}