# =========
# load packages
# =========
setwd("~/renyimeng/GAGNAR_140/JBES_empirical/city/")
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
load("./city_data/Ymat_GDP_growth.rda")
load("./city_data/Xmat.rda")
load("./city_data/Adjmat.rda")
load("./city_data/dmat.rda")  # d.matrix
# visualize the histogram of stock return
yy_dat <- data.frame(ID = 1:nrow(Ymat), value = rowMeans(Ymat))
pdf(width=4, height=4, file = "./report/city_des.pdf") 
ggplot(yy_dat, aes(x = value)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "Average FR/GDP", y = "Frequency") +
  theme(axis.text = element_text(size = 10))
dev.off()

cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]

# 第37个城市在2015年是异常值，删除这个城市
# Ymat = Ymat[-37,]
# Adjmat = Adjmat[-37, -37]
# d.matrix = d.matrix[-37, -37]
# Zmat = Zmat[-37,]
# =========
# experiment setting
# =========
train_window = 12
test_window = 5
n_window = 1
train_start = 1
n_test = 1
# set data frame for result
pred_rmse = matrix(0, nrow = 7, ncol = 100)
# basic setting & tuning param
niterations <- 1500 # number of iterations for MCMC
initNClusters <- 10
burnin <- 500 # number of burnin iterations for MCMC
h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)
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
save(real_res, file = "./report/city_res_GDP_growth.rda")
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

# (4) save predict relative RMSE
pred_rmse[1,1:test_window] <- rowMeans((Ymat_test[-1,] - YY_pred)^2)/rowMeans((Ymat_test[-1,] - colMeans(Ymat_train))^2)

# predict with CRP
# (3) predict
ZZ_crp <- real_res[[1]]$Dahlout$zout %>% as.factor()
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K
Theta_crp <- real_res[[1]]$Dahlout$thetaout  # K*(p+3)
node_est <- ZZ1_crp %*% Theta_crp
Sigma2_crp <- real_res[[1]]$Dahlout$sigma2out

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
train_window = 12
test_window = 5
n_window = 1
train_start = 1
n_test = 1
train_end = train_start + train_window - 1
test_start = train_start + 1
test_end = test_start + test_window - 1

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
GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 4, method="kmeans", seed = n_test)
Theta = GNAR_summary$theta
Sigma2 = (GNAR_summary$sigma)^2
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

# write.csv(pred_rmse[1:3,1:test_window], "./report/pred_rmse_city.csv", row.names = F)