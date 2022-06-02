# =========
# load packages
# =========
# setwd("./case_study/stock/")
setwd("~/Desktop/renyimeng/复旦项目/GAGNAR/GAGNAR_code_JBES/case_study/stock")
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
source("./code/GNAR_code-master/estimator.R")
source("./code/real_estimator.R")
# =========
# load data
# =========
load("./stock_data/Ymat_top.rda")
load("./stock_data/Zmat_top.rda")
load("./stock_data/Adjmat_top.rda")
cat("network density is", sum(adj.mt)/(nrow(adj.mt)^2 - nrow(adj.mt)))
load("./stock_data/dmat_top.rda")  # d.matrix
Zmat <- scale(Zmat)

# visualize the histogram of stock return
yy_dat <- data.frame(ID = 1:nrow(Ymat), value = rowMeans(Ymat))
pdf(width=4, height=4, file = "./report/stk_top_histogram.pdf") 
ggplot(yy_dat, aes(x = value)) +
  geom_histogram(binwidth = 0.005, color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "Average Weekly Stock Return Rate", y = "Frequency") +
  theme(axis.text = element_text(size = 10))
dev.off()
# mean(yy_dat$value)

# ==========
# experiment setting
# ==========

train_window = 40
test_window = 12
n_window = 5
train_start = 1
n_test = 1
flag = TRUE
test_days = NULL
train_end = train_start + train_window - 1
test_start = train_end
test_end = test_start + test_window

# 1. split train and test set
Ymat_train <- Ymat[,(train_start:train_end)]
Ymat_test <- Ymat[,(test_start:test_end)]

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
h <- c(seq(0,2,0.2),3,4,5) 
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)

res_topstk <- runReal(Ymat_train, X_cov, d.matrix, hyperpara, niterations, 
                      a0, b0, tau0, sigma0, initNClusters, verbose = T)
# save(res_topstk, file = "./report/stock_res.rda")
besth_ind <- which.max(unlist(lapply(res_topstk, function(x){x$LPML})))
besth <- res_topstk[[besth_ind]]$h
ZZ_best <- res_topstk[[besth_ind]]$Dahlout$zout %>% as.factor()
ZZ1 <- model.matrix( ~ ZZ_best - 1)  # N*K
Theta <- res_topstk[[besth_ind]]$Dahlout$thetaout  # K*(p+3)
Sigma2 <- res_topstk[[besth_ind]]$Dahlout$sigma2out

# (3) predict
pred_rmse = matrix(0, nrow = 7, ncol = 100)
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

# (4) save predict relative RMSE
pred_rmse[1,1:test_window] <- rowMeans((Ymat_test[-1,] - YY_pred)^2)/rowMeans((Ymat_test[-1,] - colMeans(Ymat_train))^2)

# predict with CRP
# (3) predict
ZZ_crp <- res_topstk[[1]]$Dahlout$zout %>% as.factor()
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K
Theta_crp <- res_topstk[[1]]$Dahlout$thetaout  # K*(p+3)
Sigma2_crp <- res_topstk[[1]]$Dahlout$sigma2out
node_est <- ZZ1_crp %*% Theta_crp

YY_pred <- matrix(0, nrow = test_window, ncol = N)
for (tt in 2:(test_window+1)) {
  XX_const <- rep(1,N)
  XX_net <- W %*% Ymat_test[(tt-1),]
  XX_mom <- Ymat_test[(tt-1),]
  XX_covar <- Zmat
  XX <- cbind(XX_const, XX_net, XX_mom, XX_covar)
  YY_t1 <- diag(XX %*% t(ZZ1_crp %*% Theta_crp))
  YY_pred[tt-1,] <- YY_t1
}
pred_rmse[2,1:test_window] <- rowMeans((Ymat_test[-1,] - YY_pred)^2)/rowMeans((Ymat_test[-1,] - colMeans(Ymat_train))^2)

# ===============
# compare with GNAR
# ===============
train_window = 40
test_window = 12
n_window = 5
train_start = 1
n_test = 1

# ==========
# estimation
# ==========
train_end = train_start + train_window - 1
test_start = train_end + 1
test_end = test_start + test_window - 1

# 1. split train and test set
Ymat_train <- Ymat[,(train_start:train_end)]
Ymat_test <- Ymat[,((test_start-1):test_end)]  # use the last of trainset as start point

TT = ncol(Ymat_train)
N = nrow(Ymat_train)
# 2. build XX matrix
p <- ncol(Zmat)
W <- adj.mt/rowSums(adj.mt)
W[rowSums(adj.mt)==0,] = 0
# 3. estimate
# (2) est
GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 5, method="complete", seed = n_test)
Theta = t(GNAR_summary$theta)
Sigma2 <- (GNAR_summary$sigma)^2
# (3) predict
ind = order(GNAR_summary$alpha)
group_cl = order(ind)[GNAR_summary$group]
ZZ <- group_cl %>% as.factor()
ZZ1 <- model.matrix( ~ ZZ - 1)  # N*K

YY_pred <- matrix(0, nrow = N, ncol = test_window)
for (tt in 2:(test_window+1)) {
  cat(tt,"\n")
  XX_const <- rep(1,N)
  XX_net <- W %*% Ymat_test[,(tt-1)]
  XX_mom <- Ymat_test[,(tt-1)]
  XX_covar <- Zmat
  XX <- cbind(XX_const, XX_net, XX_mom, XX_covar)
  YY_t1 <- diag(XX %*% t(ZZ1 %*% Theta))
  YY_pred[,tt-1] <- YY_t1
}

# (4) save predict relative RMSE
pred_rmse[3,1:test_window] <- colMeans((Ymat_test[,-1] - YY_pred)^2)/colMeans((Ymat_test[,-1] - rowMeans(Ymat_train))^2)

write.csv(pred_rmse[1:3,1:test_window], "./report/pred_rmse_stock.csv", row.names = F)
