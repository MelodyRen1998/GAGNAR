# 1. city
city_theta = read.csv("./city_report/city_theta.csv")
load("./city_report/city_sigma2.rda")
load("./city_report/city_label.rda")

# given sigma_k^2, theta_k | sigma_k^2 ~ N(tau_0, sigma_k^2 Sigma_0)
K = 4
p = 4
p_values_gagnar = matrix(0, nrow = K, ncol = p+3)
for (kk in 1:K) {
  p_values_gagnar[kk, ] = (1 - pnorm(abs(as.numeric(city_theta[kk, ]))/sqrt(Sigma2[kk] * 100)))*2
}
# 使用NIG的方差公式计算
a0 <- 0.01
b0 <- 0.01
tau0 <- rep(0, p + 3)
sigma0 <- diag(100, p + 3)
sigma0_inv = solve(sigma0)

load("./city/city_data/Ymat.rda")
load("./city/city_data/Xmat.rda")
load("./city/city_data/Adjmat.rda")
cat("network density is", sum(Adjmat)/(nrow(Adjmat)^2 - nrow(Adjmat)), "\n")
load("./city/city_data/dmat.rda")  # d.matrix
cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]

train_window = 12
train_start = 5
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

X_tilde = X_cov
clusterAssign = ZZ_best
Y_tilde = Ymat_train
t = nrow(Ymat_train)
var_theta = list()
for (r in 1:K){
  cat("group:", r,"\n")
  ## sigma*
  X_group_r <- X_tilde[2:t, which(clusterAssign == r), ]
  
  if (sum(clusterAssign == r) == 1) {  # if there is single node, then X_group_r appears to be a matrix
    X_tilde_t_sum <- t(X_group_r) %*% X_group_r
  } else if (sum(clusterAssign == r) > 1) { 
    X_list_bygroup <- lapply(1:dim(X_group_r)[2], function(rr) {X_group_r[,rr,]})
    X_tilde_t <- lapply(X_list_bygroup, function(x){t(x) %*% x})
    X_tilde_t_sum <- Reduce('+', X_tilde_t)
  }
  
  post_sigma_star <- ginv(sigma0_inv + X_tilde_t_sum)
  
  ## tau*
  Y_group_r <- Y_tilde[2:t, clusterAssign == r]  # dim(Y_group_r) = [T, Nr]
  
  if (sum(clusterAssign == r) == 1) {  # if there is single node, then Y_group_r appears to be a vector
    X_tilde_t_Y_sum <- t(X_group_r) %*% Y_group_r  # dim = [d, 1]
  } else if (sum(clusterAssign == r) > 1) {
    X_tilde_t_Y <- lapply(1:dim(Y_group_r)[2], function(x) {t(X_list_bygroup[[x]]) %*% Y_group_r[,x]})
    X_tilde_t_Y_sum <- Reduce('+', X_tilde_t_Y)
  }
  
  post_tau_star <- post_sigma_star %*% (sigma0_inv %*% tau0 + X_tilde_t_Y_sum)
  
  ## a*
  post_a_star <- a0 + (t-1) * sum(clusterAssign == r)/2
  
  ## b*
  Y_tilde_sum <- sum(apply(as.matrix(Y_tilde[2:t, clusterAssign == r]), 2, function(x){t(x) %*% x}))
  post_b_star <- (b0 + (t(tau0) %*% (sigma0_inv) %*% tau0 + Y_tilde_sum - t(post_tau_star) %*% (sigma0_inv + X_tilde_t_sum) %*% post_tau_star)/2) %>% as.numeric()
  
  var_theta[[r]] = ginv(post_sigma_star) * post_b_star/(post_a_star - 1)
  (p_values_gagnar[r, ] = (1 - pnorm(abs(as.numeric(city_theta[r, ]))/sqrt(diag(var_theta[[r]]))))*2)
}

# 使用GAGNAR估计的Z_i代入GNAR计算p-value


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
ZZ_best = as.numeric(ZZ_best)
GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 4, method="other", memb = ZZ_best, seed = n_test)
Theta = GNAR_summary$theta
Sigma2 = GNAR_summary$sig2
p_values_gagnar_gnar = matrix(0, nrow = K, ncol = p+3)
for (kk in 1:K) {
  p_values_gagnar_gnar[kk, ] = (1 - pnorm(abs(Theta[kk, ])/sqrt(diag(GNAR_summary$cov[[kk]]))))*2
}
