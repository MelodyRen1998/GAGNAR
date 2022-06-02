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
load("./stock_data/Ymat_top.rda")
load("./stock_data/Zmat_top.rda")
load("./stock_data/Adjmat_top.rda")
cat("network density is", sum(adj.mt)/(nrow(adj.mt)^2 - nrow(adj.mt)), "\n")
load("./stock_data/dmat_top.rda")  # d.matrix
Zmat <- scale(Zmat)

# visualize the histogram of stock return
# yy_dat <- data.frame(ID = 1:nrow(Ymat), value = rowMeans(Ymat))
# pdf(width=4, height=4, file = "./report/stk_top_histogram.pdf") 
# ggplot(yy_dat, aes(x = value)) +
#   geom_histogram(binwidth = 0.005, color = "black", fill = "grey") +
#   theme_bw() +
#   labs(x = "Average Weekly Stock Return Rate", y = "Frequency") +
#   theme(axis.text = element_text(size = 10))
# dev.off()
# mean(yy_dat$value)

# ==========
# experiment setting
# ==========

train_window = 35
test_window = 10
n_window = 1  # rolling window size
train_start = 1s
n_test = 1
flag = TRUE
test_days = NULL
pred_rmse = matrix(0, nrow = 7, ncol = 100)


while (flag) {
  train_end = train_start + train_window - 1
  test_start = train_end
  test_end = test_start + test_window
  
  if ((train_end<ncol(Ymat)) & (test_end>ncol(Ymat))){
    test_end = ncol(Ymat)
    flag = F}
  if (test_start>ncol(Ymat)){
    flag = F}
  
  cat("GAGNAR", train_start, "\n")
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
  # h <- seq(0,2,0.2)
  h <- c(0, 0.6)  # crp and optimal h
  alpha <- c(1)
  hyperpara <- expand.grid(h = h, alpha = alpha)
  
  res_topstk <- runReal(Ymat_train, X_cov, d.matrix, hyperpara, niterations, 
                        a0, b0, tau0, sigma0, initNClusters, verbose = T)
  # save(res_topstk, file = "./report/stock_res.rda")
  # besth_ind <- which.max(unlist(lapply(res_topstk, function(x){x$LPML})))
  besth_ind <- 2
  besth <- res_topstk[[besth_ind]]$h
  ZZ_best <- res_topstk[[besth_ind]]$Dahlout$zout %>% as.factor()
  ZZ1 <- model.matrix( ~ ZZ_best - 1)  # N*K
  Theta <- res_topstk[[besth_ind]]$Dahlout$thetaout  # K*(p+3)
  Sigma2 <- res_topstk[[besth_ind]]$Dahlout$sigma2out
  
  # (3) predict
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
  pred_rmse[1,n_test] <- mean((Ymat_test[-1,] - YY_pred)^2)/mean((Ymat_test[-1,] - colMeans(Ymat_train))^2)
  
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
  pred_rmse[2,n_test] <- mean((Ymat_test[-1,] - YY_pred)^2)/mean((Ymat_test[-1,] - colMeans(Ymat_train))^2)
  
  # ===============
  # compare with GNAR
  # ===============
  # ==========
  # estimation
  # ==========
  # 1. split train and test set
  Ymat_train <- Ymat[,(train_start:train_end)]
  Ymat_test <- Ymat[,(test_start:test_end)]  # use the last of trainset as start point
  
  TT = ncol(Ymat_train)
  N = nrow(Ymat_train)
  # 2. build XX matrix
  p <- ncol(Zmat)
  W <- adj.mt/rowSums(adj.mt)
  W[rowSums(adj.mt)==0,] = 0
  # 3. estimate
  # (2) est
  # GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 5, method="complete", seed = n_test)
  # Theta = t(GNAR_summary$theta)
  # Sigma2 <- (GNAR_summary$sigma)^2
  K = 5
  GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 5, method="complete", seed = n_test)
  Theta = GNAR_summary$theta
  Sigma2 = GNAR_summary$sig2
  # p_values_gnar = matrix(0, nrow = K, ncol = p+3)
  # for (kk in 1:K) {
  #   p_values_gnar[kk, ] = (1 - pnorm(abs(Theta[kk, ])/sqrt(diag(GNAR_summary$cov[[kk]]))))*2
  # }
  
  # (3) predict
  ind = order(GNAR_summary$alpha)
  group_cl = order(ind)[GNAR_summary$group]
  ZZ <- group_cl %>% as.factor()
  ZZ1 <- model.matrix( ~ ZZ - 1)  # N*K
  
  YY_pred <- matrix(0, nrow = N, ncol = test_window)
  for (tt in 2:(test_window+1)) {
    # cat(tt,"\r")
    XX_const <- rep(1,N)
    XX_net <- W %*% Ymat_test[,(tt-1)]
    XX_mom <- Ymat_test[,(tt-1)]
    XX_covar <- Zmat
    XX <- cbind(XX_const, XX_net, XX_mom, XX_covar)
    YY_t1 <- diag(XX %*% t(ZZ1 %*% Theta))
    YY_pred[,tt-1] <- YY_t1
  }
  
  # (4) save predict relative RMSE
  pred_rmse[3,n_test] <- mean((Ymat_test[,-1] - YY_pred)^2)/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
  
  # ===============
  # compare with NAR
  # ===============
  cat("NAR" ,train_start, "\n")
  nar_est = NAR.est(Ymat_train, W, Zmat)
  nar_pred = NAR.predict(Ymat_test, W, theta = nar_est$theta, Zmat)
  pred_rmse[4, n_test] = nar_pred$rmse^2/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
  
  # ===============
  # compare with ARMA(1,1)
  # ===============
  arma_pred = matrix(0, nrow = nrow(Ymat_test), ncol = ncol(Ymat_test)-1)
  for (r in 1:N) {
    cat("ARMA" ,train_start, r, "\r")
    arma_pred[r,] = tryCatch({
      if (sum(Ymat_train[r,]==0)>(ncol(Ymat_train)/2))
      {
        rep(0, ncol(Ymat_test)-1)
      }else{
        fit = arima(Ymat_train[r,], order = c(1,0,1), method = "ML")
        coef(fit)['intercept'] + coef(fit)['ar1'] * Ymat_test[r,-ncol(Ymat_test)]
      }
    }, error = function(err)
    {
      print(c(r, paste("ARMA ERROR:  ",err)))
      rep(0, ncol(Ymat_test)-1)
    })
  }
  pred_rmse[5, n_test] = mean((Ymat_test[,-1] - arma_pred)^2)/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
  
  # ===============
  # compare with AR(1)
  # ===============
  ar_pred = matrix(0, nrow = nrow(Ymat_test), ncol = ncol(Ymat_test)-1)
  for (r in 1:N) {
    cat("AR" ,train_start, r, "\r")
    ar_pred[r,] = tryCatch({
      if (sum(Ymat_train[r,]==0)>(ncol(Ymat_train)/2))
      {
        rep(0, ncol(Ymat_test)-1)
      }else{
        fit = arima(Ymat_train[r,], order = c(1,0,0), method = "ML")
        coef(fit)['intercept'] + coef(fit)['ar1'] * Ymat_test[r,-ncol(Ymat_test)]
      }
    }, error = function(err)
    {
      print(c(r, paste("AR ERROR:  ",err)))
      rep(0, ncol(Ymat_test)-1)
    })
  }
  pred_rmse[6, n_test] = mean((Ymat_test[,-1] - ar_pred)^2)/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
  
  # pred_rmse每一行分别为 GAGNAR, CRP, GNAR, NAR, ARMA, AR
  cat(pred_rmse[, n_test], "\n")
  n_test = n_test + 1
  train_start = train_start + n_window
  test_days = c(test_days, test_start)
  if ((train_start + train_window - 1) > ncol(Ymat)){flag = F}
}


write.csv(pred_rmse[,1:(n_test - 1)], "./report/pred_rmse_stock_train30_test10.csv", row.names = F)
# write.csv(p_values_gnar, "./report/gnar_p_values.csv", row.names = F)
# 
# # RMSE line plot
# start_time_point <- 41:47
# plot(x = start_time_point, y = pred_rmse[1,1:(n_test - 1)], 
#      xlab = "Starting Time Point of Test Data", ylab = "ReMSPE", 
#      col = "red", type = "l", lwd = 2, xaxt="n", lty = 1, ylim = c(0.8,1.5))
# axis(1, at = start_time_point)
# lines(x = start_time_point, y = pred_rmse[2,1:(n_test - 1)], col = "#96cac1", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[3,1:(n_test - 1)], col = "#6A93CC", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[4,1:(n_test - 1)], col = "#facb4c", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[5,1:(n_test - 1)], col = "#df9a96", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[6,1:(n_test - 1)], col = "#7f7f7f", type = "l", lwd = 2)
# legend("topleft", c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"), 
#        col = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f"),
#        lty = c(1,1,1,1,1,1), lwd = 2)
# 
# df_pred = data.frame(ID = rep(start_time_point,6), 
#                      value = c(t(pred_rmse[c(1,2,3,4,5,6),1:(n_test - 1)])),
#                      type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
#                               rep("NAR", n_test - 1), rep("ARMA", n_test - 1), rep("AR", n_test - 1)))
# df_pred$type <- factor(df_pred$type, levels = c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"))
# 
# 
# 
# pdf(width=5, height=4, file = "./report/stk_pred_rmse.pdf") 
# ggplot(df_pred, aes(x = as.integer(ID), y = value, group = type)) +
#   geom_line(aes(color = type), lwd = 0.8) +
#   scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
#   theme_bw() +
#   ylim(0.8,1.5) +
#   theme(legend.position = c(0.12,0.78),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         text = element_text(size = 12)) +
#   scale_x_continuous(breaks = c(start_time_point)) +
#   labs(x = "", y = "ReMSPE")
# dev.off() 
# 
# # 从本地读入结果数据
# pred_rmse = read.csv("./report/pred_rmse_stock.csv")
# df_pred = data.frame(ID = rep(start_time_point,6), 
#                      value = c(t(pred_rmse[c(1,2,3,4,5,6),1:(n_test - 1)])),
#                      type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
#                               rep("NAR", n_test - 1), rep("ARMA", n_test - 1), rep("AR", n_test - 1)))
# df_pred$type <- factor(df_pred$type, levels = c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"))
# 
# # 不考虑NAR的
# df_pred1 = df_pred[which(df_pred$type != "NAR"),]
# pdf(width=5, height=4, file = "./report/stk_pred_rmse_roll.pdf") 
# ggplot(df_pred1, aes(x = as.integer(ID), y = value, group = type)) +
#   geom_line(aes(color = type), lwd = 0.8) +
#   scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96")) +
#   theme_bw() +
#   ylim(0.8,1.5) +
#   theme(legend.position = c(0.12,0.83),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         text = element_text(size = 12)) +
#   scale_x_continuous(breaks = c(start_time_point)) +
#   labs(x = "Start Time Points", y = "ReMSPE")
# dev.off() 
