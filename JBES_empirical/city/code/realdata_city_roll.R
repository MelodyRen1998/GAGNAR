# =========
# load packages
# =========
setwd("~/renyimeng/GAGNAR_140/JBES_empirical/city")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("scales", "igraph", "mvtnorm", "MCMCpack", 
              "parallel", "mclust", "plyr", "ggplot2")
ipak(packages)

source("./code/gagnar.R")
source("./code/Dahl.R")
source("./code/LPML_VAR.R")
source("./code/GNAR_code-master/estimator.R")
source("./code/real_estimator.R")
# =========
# load data
# =========
load("./city_data/Ymat.rda")
load("./city_data/Xmat.rda")
load("./city_data/Adjmat.rda")
cat("network density is", sum(Adjmat)/(nrow(Adjmat)^2 - nrow(Adjmat)), "\n")
load("./city_data/dmat.rda")  # d.matrix

# visualize the histogram of stock return
Ymat1 = Ymat[,c(5:16)]
yy_dat <- data.frame(ID = 1:nrow(Ymat1), value = rowMeans(Ymat1))
pdf(width=4, height=4, file = "./report/city_des_2005_2016.pdf")
ggplot(yy_dat, aes(x = value)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "Average FR/GDP", y = "Frequency") +
  theme(axis.text = element_text(size = 10))
dev.off()

yy_dat2 <- data.frame(YY = colnames(Ymat1), value = colMeans(Ymat1))
pdf(width=4, height=4, file = "./report/city_line_2005_2016.pdf")
ggplot(yy_dat2, aes(x = YY, y = value, group = 1)) +
  geom_line() +
  theme_bw() +
  labs(x = "Year", y = "Average FR/GDP") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylim(0,0.1)
dev.off()



mean(Ymat1)

cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]

# Ymat = Ymat[-c(32,119),]
# Zmat = Zmat[-c(32,119),]
# Adjmat = Adjmat[-c(32,119),-c(32,119)]
# d.matrix = d.matrix[-c(32,119),-c(32,119)]
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
# p_value = matrix(0, nrow = 7, ncol = p+3)
# basic setting & tuning param
niterations <- 1500 # number of iterations for MCMC
initNClusters <- 10
burnin <- 500 # number of burnin iterations for MCMC
h <- c(0,2) # crp and optimal h
# h <- c(seq(0,2,0.2),3,4,5)
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
    
    # (4) save predict relative RMSE
    pred_rmse[1,n_test] <- mean((Ymat_test[-1,] - YY_pred)^2)/mean((Ymat_test[-1,] - colMeans(Ymat_train))^2)
    # p_value[1,] <- rep(NA, p+3)
    
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
    pred_rmse[2,n_test] <- mean((Ymat_test[-1,] - YY_pred)^2)/mean((Ymat_test[-1,] - colMeans(Ymat_train))^2)
    # p_value[2,] <- rep(NA, p+3)
    
    # ===============
    # compare with GNAR
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
    K = 4
    GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 4, method="kmeans", seed = n_test)
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
      # cat(tt,"\n")
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
    (1 - pnorm(abs(nar_est$theta)/sqrt(diag(nar_est$cov))))*2
    nar_pred = NAR.predict(Ymat_test, W, theta = nar_est$theta, Zmat)
    pred_rmse[4, n_test] = nar_pred$rmse^2/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
    
    # ===============
    # compare with AR(1) * do not report p-values for ARMA/AR model
    # ===============
    arma_res = list()
    arma_pred = matrix(0, nrow = nrow(Ymat_test), ncol = ncol(Ymat_test)-1)
    arma_p_value = matrix(0, nrow = nrow(Ymat_test), ncol = 2)
    for (r in 1:N) {
      cat("AR" ,train_start, r, "\r")
      tryCatch({
        if (sum(Ymat_train[r,]==0)>(ncol(Ymat_train)/2)) { 
          arma_pred[r, ] = rep(0, ncol(Ymat_test)-1) 
          arma_p_value[r, ] = rep(0, 2) } else {
            fit = arima(Ymat_train[r,], order = c(1,0,0), method = "ML")
            arma_pred[r, ] = coef(fit)['intercept'] + coef(fit)['ar1'] * Ymat_test[r,-ncol(Ymat_test)]  # vector
            arma_p_value[r, ] = (1 - pnorm(abs(fit$coef)/sqrt(diag(fit$var.coef))))*2  # vector
          }}, error = function(err) { 
            print(c(r, paste("ARMA ERROR:  ",err)))
            arma_pred[r, ] = rep(0, ncol(Ymat_test)-1)
            arma_p_value[r, ] = rep(0, 2)})
    }
    
    pred_rmse[5, n_test] = mean((Ymat_test[,-1] - arma_pred)^2)/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
    # p_values[5, 2] = colMeans(arma_p_value)
    # ===============
    # compare with ARMA(1,1) 注意这里的变量命名和AR(1)反了
    # ===============
    ar_pred = matrix(0, nrow = nrow(Ymat_test), ncol = ncol(Ymat_test)-1)
    ar_p_value = matrix(0, nrow = nrow(Ymat_test), ncol = 3)
    
    for (r in 1:N) {
      cat("AR" ,train_start, r, "\r")
      tryCatch({
        if (sum(Ymat_train[r,]==0)>(ncol(Ymat_train)/2)){ 
          ar_pred[r,] = rep(0, ncol(Ymat_test)-1) 
          ar_p_value[r, ] = rep(0, 3)} else {
            fit = arima(Ymat_train[r,], order = c(1,0,1), method = "ML")
            ar_pred[r,] = coef(fit)['intercept'] + coef(fit)['ar1'] * Ymat_test[r,-ncol(Ymat_test)]
            ar_p_value[r, ] = (1 - pnorm(abs(fit$coef)/sqrt(diag(fit$var.coef))))*2  # vector
          }}, error = function(err) { 
            print(c(r, paste("ARMA ERROR:  ",err)))
            ar_pred[r,] = rep(0, ncol(Ymat_test)-1)
            ar_p_value[r, ] = rep(0, 3)})
    }
    pred_rmse[6, n_test] = mean((Ymat_test[,-1] - ar_pred)^2)/mean((Ymat_test[,-1] - rowMeans(Ymat_train))^2)
    # p_value[6, 3] = colMeans(ar_p_value)
    
    # pred_rmse每一行分别为 GAGNAR, CRP, GNAR, NAR, AR, ARMA
    cat(pred_rmse[, n_test], "\n")
    n_test = n_test + 1
    train_start = train_start + n_window
    test_days = c(test_days, test_start)
    if ((train_start + train_window - 1) > ncol(Ymat)){flag = F}
  }
} , error = function(err) {
  print(paste0("rolling ERROR: ",err))})

write.csv(pred_rmse[,1:(n_test - 1)], "./report/pred_rmse_city_train7_test7_remove2.csv", row.names = F)

# write.csv(p_values_gnar, "./report/gnar_p_values_train12.csv", row.names = F)
# 
# 
# # RMSE line plot
# start_time_point <- 2002:2006
# plot(x = start_time_point, y = pred_rmse[1,1:(n_test - 1)], 
#      xlab = "Starting Time Point of Test Data", ylab = "ReMSPE", 
#      col = "red", type = "l", lwd = 2, xaxt="n", lty = 1, ylim = c(0.1,12))
# axis(1, at = start_time_point)
# lines(x = start_time_point, y = pred_rmse[2,1:(n_test - 1)], col = "#96cac1", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[3,1:(n_test - 1)], col = "#6A93CC", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[4,1:(n_test - 1)], col = "#facb4c", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[6,1:(n_test - 1)], col = "#df9a96", type = "l", lwd = 2)
# lines(x = start_time_point, y = pred_rmse[5,1:(n_test - 1)], col = "#7f7f7f", type = "l", lwd = 2)
# legend("topleft", c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"), 
#        col = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f"),
#        lty = c(1,1,1,1,1,1), lwd = 2)
# 
# df_pred = data.frame(ID = rep(start_time_point,6), 
#                      value = c(t(pred_rmse[c(1,2,3,4,6,5),1:(n_test - 1)])),
#                      type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
#                               rep("NAR", n_test - 1), rep("ARMA", n_test - 1), rep("AR", n_test - 1)))
# df_pred$type <- factor(df_pred$type, levels = c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"))
# 
# 
# pdf(width=5, height=4, file = "./report/stk_pred_rmse.pdf") 
# ggplot(df_pred, aes(x = as.integer(ID), y = value, group = type)) +
#   geom_line(aes(color = type), lwd = 0.8) +
#   scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
#   theme_bw() +
#   ylim(0.1,12) +
#   theme(legend.position = c(0.12,0.77),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         text = element_text(size = 12)) +
#   scale_x_continuous(breaks = c(start_time_point)) +
#   labs(x = "", y = "ReMSPE")
# dev.off() 
# 
# 
# # 从本地读入结果数据
# pred_rmse_city_train12 <- read.csv("~/Desktop/renyimeng/复旦项目/GAGNAR/GAGNAR_code_JBES/case_study/city/report/pred_rmse_city_train12.csv")
# df_pred = data.frame(ID = rep(start_time_point,6), 
#                      value = c(t(pred_rmse_city_train12[c(1,2,3,4,6,5),1:(n_test - 1)])),
#                      type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
#                               rep("NAR", n_test - 1), rep("ARMA", n_test - 1), rep("AR", n_test - 1)))
# df_pred$type <- factor(df_pred$type, levels = c("GAGNAR", "CRP", "GNAR", "NAR", "ARMA", "AR"))
# df_pred2 = df_pred[which(df_pred$type != "NAR"),]
# 
# 
# # 不考虑NAR
# pdf(width=5, height=4, file = "./report/city_pred_rmse_roll.pdf") 
# ggplot(df_pred2, aes(x = as.integer(ID), y = value, group = type)) +
#   geom_line(aes(color = type), lwd = 0.8) +
#   scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
#   theme_bw() +
#   ylim(0.1,11.5) +
#   theme(legend.position = c(0.12,0.83),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         text = element_text(size = 12)) +
#   scale_x_continuous(breaks = c(start_time_point)) +
#   labs(x = "Start Time Points", y = "ReMSPE")
# dev.off() 
