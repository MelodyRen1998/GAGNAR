# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("reshape2", "ggplot2")
ipak(packages)
h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)
nsims = 100
exp_id = 3
sce_id = 2
idx_design = sce_id + (exp_id - 1)*2
K = ifelse(exp_id == 1, 3, 5)
# data & para
load(paste0("./simulation_pack/data/data_design_", idx_design, ".RData")) # data
load(paste0("./simulation_pack/data/true_para_design_", idx_design, ".RData")) # para
NN <- dim(all_simu_dat[[1]][[1]])[2]
# result
res_path <- paste0("./simulation_pack/result/example",exp_id,"/")
load(paste0(res_path, "/sim_res_", idx_design, ".rda"))
load(paste0(res_path, "/sim_res_h_", idx_design, ".rda"))
load(paste0(res_path, "/sim_res_EM_", idx_design, ".rda"))
load(paste0(res_path, "/sim_res_GNAR_", idx_design, ".rda"))

all_simu_dat[[5]][[1]]
all_simu_dat[[5]][[2]]

# ===========
# histogram
# ===========
df1 = matrix(0, nrow = nsims, ncol = nrow(hyperpara))
for (reppp in 1:nsims) {
  for (hpp in 1:nrow(hyperpara)) {
    df1[reppp, hpp] <- sim_res[[reppp]][[hpp]]$Dahlout$zout %>% unique() %>% length()
  }
}
df1 = data.frame(df1)
colnames(df1) = paste0('h=', h, sep = "")
best_hs = c()
for (reppp in 1:nsims) {
  lpml <- lapply(sim_res[[reppp]], function(each_rep) {each_rep$LPML}) %>% unlist()
  best_hs[reppp] <- hyperpara$h[which.max(lpml)]
}
best_h = mean(best_hs)
cat("best h is", best_h)
df2 = c()
for (reppp in 1:nsims) {
  df2[reppp] <- sim_res_h[[reppp]]$Dahlout$zout %>% unique() %>% length()
}
df1$`LPML` = df2
df_vis <- df1[,c(paste0("h=", c(0,0.4,1,2,5)), "LPML")]
df22 <- melt(df_vis)
df22$flag = df22$value == K
df22 <- df22[df22$value <= 7,]
pdf(width=5, height=4, file = paste0("./add_simu_result/report/exp_", exp_id, "sce_", sce_id, "_k.pdf"))
ggplot(df22, aes(x = value, fill = flag)) + 
  geom_histogram(bins=4, binwidth = 0.5) +
  facet_wrap(~variable) +
  scale_x_continuous(name ="Number of Clusters", breaks=c(3,4,5,6,7,8,9)) +
  ylab("Number of Replicates") +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c("#709fd4", "#CD6155"))
dev.off()

# ==========
# ari
# ==========
ari = matrix(0, nrow = nsims, ncol = nrow(hyperpara))
acc_df = matrix(0, nrow = nsims, ncol = nrow(hyperpara))
for (reppp in 1:nsims) {
  for (hypp in 1:nrow(hyperpara)) {
    model_label = sim_res[[reppp]][[hypp]]$Dahlout$zout
    true_label = all_simu_dat[[5]][[reppp]]
    ariii = adjustedRandIndex(model_label, true_label)
    acc_ = 1 - classError(model_label, true_label)$errorRate
    ari[reppp, hypp] = ariii
    acc_df[reppp, hypp] = acc_
    # cat(ariii, "\n")
  }
}

ari_h = c()
acc_h = c()
for (reppp in 1:nsims) {
  model_label = sim_res_h[[reppp]]$Dahlout$zout
  true_label = all_simu_dat[[5]][[reppp]]
  ariii = adjustedRandIndex(model_label, true_label)
  acc_ = 1 - classError(model_label, true_label)$errorRate
  ari_h[reppp] = ariii
  acc_h[reppp] = acc_
}
ari = data.frame(ari)
acc_df = data.frame(acc_df)
colnames(ari) = paste0('h=', h, sep = "")
colnames(acc_df) = paste0('h=', h, sep = "")
ari$LPML = ari_h
acc_df$LPML = acc_h


NNiter = 1500
ari_em = c()
acc_em = c()
for (reppp in 1:nsims) {
  if (length(sim_res_EM[[reppp]]) > 0) {
    em_clu = sim_res_EM[[reppp]][[NNiter]]$clu
    em_label = apply(em_clu, 1, function(x){which.max(x)}) 
    true_label = all_simu_dat[[5]][[reppp]]
    ariii = adjustedRandIndex(em_label, true_label)
    acc = 1 - classError(em_label, true_label)$errorRate
    ari_em[reppp] = ariii
    acc_em[reppp] = acc
  }
}
ari$EM = ari_em
acc_df$EM = acc_em

ari_vis = ari[,c(paste0("h=", 0), "LPML","EM")]
acc_df_vis = acc_df[,c(paste0("h=", 0), "LPML","EM")]
ari1 = melt(ari_vis)
acc_ratio = colMeans(acc_df_vis, na.rm = T) %>% round(4)
pdf(width=5, height=4, file = paste0("./add_simu_result/report_subset/exp_", exp_id, "sce_", sce_id, "_ari.pdf")) 
ggplot(ari1, aes(x = variable, y = value)) +
  geom_boxplot() +
  # ylim(0.5,1.12) +
  # stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               # width = 1, linetype = "dashed", col = "red")+
  geom_text(data=data.frame(), aes(x=names(acc_ratio), y=0.85, label=acc_ratio), col='red', size=3, fontface = "italic") +
  geom_text(data=data.frame(), aes(x="h=0", y=0.8, label="Accuracy"), col='grey', size=4, fontface = "italic") +
  xlab("") + ylab("ARI") +theme_bw()
dev.off()

# ==========
# RMSE
# ==========
true_theta = matrix(unlist(true_para[[2]]), nrow = K, byrow = T)
true_sigma2 = unlist(true_para[[1]])

coef_RMSE = function(est, true, needSort = F, scale = F) {
  # two vec SUM of MSE
  if (needSort) {
    est = sort(est)
    true = sort(true)
  }
  err2 = (est - true)^2
  if (scale) {
    for (kk in 1:K) {
      if (true[kk] != 0) {
        err2[kk] = err2[kk]/(true[kk]^2)
      }
    }
  }
  return(err2)  # by group
}
coef_SE <- function(est, true) {
  err2 = (est - true)^2
  return(err2)
}

RMSE_sim = matrix(0, nrow = nrow(hyperpara) + 3, ncol = length(true_theta[1,]))  # number of h +LPML+ EM+ GNAR
### BETA
for (hyppp in 1:(nrow(hyperpara))) {
  each_rep_beta0 = 0
  each_rep_beta1 = 0
  each_rep_beta2 = 0
  each_rep_gamma = 0
  each_rep_s2 = 0
  for (reppp in 1:nsims) {
    true_ZZ = model.matrix(~ as.factor(all_simu_dat[[5]][[reppp]]) - 1)
    node_true = true_ZZ %*% true_theta
    node_true_s2 = true_ZZ %*% true_sigma2
    Zhat <- sim_res[[reppp]][[hyppp]]$Dahlout$zout %>% as.factor()
    ZZhat <- model.matrix(~ Zhat -1)
    thetahat <- sim_res[[reppp]][[hyppp]]$Dahlout$thetaout
    s2hat <- sim_res[[reppp]][[hyppp]]$Dahlout$sigma2out
    node_est <- ZZhat %*% thetahat
    node_est_s2 <- ZZhat %*% s2hat
    
    each_rep_beta0 = sum(coef_SE(node_est[,1], node_true[,1])) + each_rep_beta0
    each_rep_beta1 = sum(coef_SE(node_est[,2], node_true[,2])) + each_rep_beta1
    each_rep_beta2 = sum(coef_SE(node_est[,3], node_true[,3])) + each_rep_beta2
    each_rep_gamma = sum(coef_SE(node_est[,4:6], node_true[,4:6])) + each_rep_gamma
    each_rep_s2 = sum(coef_SE(node_est_s2, node_true_s2)) + each_rep_s2
  }
  RMSE_sim[hyppp, 1] = sqrt(each_rep_beta0/(NN*nsims))
  RMSE_sim[hyppp, 2] = sqrt(each_rep_beta1/(NN*nsims))
  RMSE_sim[hyppp, 3] = sqrt(each_rep_beta2/(NN*nsims))
  RMSE_sim[hyppp, 4] = sqrt(each_rep_gamma/(NN*nsims*3))
  RMSE_sim[hyppp, 5] = sqrt(each_rep_s2/(NN*nsims))
  print(RMSE_sim[hyppp,])
}

# LPML
each_rep_beta0 = 0
each_rep_beta1 = 0
each_rep_beta2 = 0
each_rep_gamma = 0
each_rep_s2 = 0
for (reppp in 1:nsims) {
  true_ZZ = model.matrix(~ as.factor(all_simu_dat[[5]][[reppp]]) - 1)
  node_true = true_ZZ %*% true_theta
  node_true_s2 = true_ZZ %*% true_sigma2
  Zhat <- sim_res_h[[reppp]]$Dahlout$zout %>% as.factor()
  ZZhat <- model.matrix(~ Zhat -1)
  thetahat <- sim_res_h[[reppp]]$Dahlout$thetaout
  s2hat <- sim_res_h[[reppp]]$Dahlout$sigma2out
  node_est <- ZZhat %*% thetahat
  node_est_s2 <- ZZhat %*% s2hat
  
  each_rep_beta0 = sum(coef_SE(node_est[,1], node_true[,1])) + each_rep_beta0
  each_rep_beta1 = sum(coef_SE(node_est[,2], node_true[,2])) + each_rep_beta1
  each_rep_beta2 = sum(coef_SE(node_est[,3], node_true[,3])) + each_rep_beta2
  each_rep_gamma = sum(coef_SE(node_est[,4:6], node_true[,4:6])) + each_rep_gamma
  each_rep_s2 = sum(coef_SE(node_est_s2, node_true_s2)) + each_rep_s2
}
RMSE_sim[nrow(hyperpara)+1, 1] = sqrt(each_rep_beta0/(NN*nsims))
RMSE_sim[nrow(hyperpara)+1, 2] = sqrt(each_rep_beta1/(NN*nsims))
RMSE_sim[nrow(hyperpara)+1, 3] = sqrt(each_rep_beta2/(NN*nsims))
RMSE_sim[nrow(hyperpara)+1, 4] = sqrt(each_rep_gamma/(NN*nsims*3))
RMSE_sim[nrow(hyperpara)+1, 5] = sqrt(each_rep_s2/(NN*nsims))
print(RMSE_sim[nrow(hyperpara)+1,])

# EM
each_rep_beta0 = 0
each_rep_beta1 = 0
each_rep_beta2 = 0
each_rep_gamma = 0
each_rep_s2 = 0
for (reppp in 1:nsims) {
  if(length(sim_res_EM[[reppp]]) > 0) {
    em_clu = sim_res_EM[[reppp]][[NNiter]]$clu
    em_label = apply(em_clu, 1, function(x){which.max(x)}) %>% as.factor()
    ZZhat = model.matrix(~ em_label - 1)
    thetahat = sim_res_EM[[reppp]][[NNiter]]$theta
    s2hat = sim_res_EM[[reppp]][[NNiter]]$sigma
    node_est = ZZhat %*% thetahat
    node_est_s2 = ZZhat %*% s2hat
    ZZtrue = model.matrix(~ as.factor(all_simu_dat[[5]][[reppp]]) - 1)
    node_true = ZZtrue %*% true_theta
    node_true_s2 = true_ZZ %*% true_sigma2
    
    each_rep_beta0 = sum(coef_SE(node_est[,1], node_true[,1])) + each_rep_beta0
    each_rep_beta1 = sum(coef_SE(node_est[,2], node_true[,2])) + each_rep_beta1
    each_rep_beta2 = sum(coef_SE(node_est[,3], node_true[,3])) + each_rep_beta2
    each_rep_gamma = sum(coef_SE(node_est[,4:6], node_true[,4:6])) + each_rep_gamma
    each_rep_s2 = sum(coef_SE(node_est_s2, node_true_s2)) + each_rep_s2
  }
}  
RMSE_sim[nrow(hyperpara)+2, 1] = sqrt(sum(each_rep_beta0)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+2, 2] = sqrt(sum(each_rep_beta1)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+2, 3] = sqrt(sum(each_rep_beta2)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+2, 4] = sqrt(sum(each_rep_gamma)/(NN*nsims*3))
RMSE_sim[nrow(hyperpara)+2, 5] = sqrt(sum(each_rep_s2)/(NN*nsims))

## GNAR
each_rep_beta0 = 0
each_rep_beta1 = 0
each_rep_beta2 = 0
each_rep_gamma = 0
each_rep_s2 = 0
for (reppp in 1:nsims) {
  Zhat = sim_res_GNAR[[reppp]]$group %>% as.factor()
  ZZhat = model.matrix(~ Zhat - 1)
  thetahat = sim_res_GNAR[[reppp]]$theta
  s2hat = sim_res_GNAR[[reppp]]$sigma
  node_est = ZZhat %*% thetahat
  node_est_s2 = ZZhat %*% s2hat
  ZZtrue = model.matrix(~ as.factor(all_simu_dat[[5]][[reppp]]) - 1)
  node_true = ZZtrue %*% true_theta
  node_true_s2 = ZZtrue %*% true_sigma2
  
  each_rep_beta0 = sum(coef_SE(node_est[,1], node_true[,1])) + each_rep_beta0
  each_rep_beta1 = sum(coef_SE(node_est[,2], node_true[,2])) + each_rep_beta1
  each_rep_beta2 = sum(coef_SE(node_est[,3], node_true[,3])) + each_rep_beta2
  each_rep_gamma = sum(coef_SE(node_est[,4:6], node_true[,4:6])) + each_rep_gamma
  each_rep_s2 = sum(coef_SE(node_est_s2, node_true_s2)) + each_rep_s2
}  
RMSE_sim[nrow(hyperpara)+3, 1] = sqrt(sum(each_rep_beta0)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+3, 2] = sqrt(sum(each_rep_beta1)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+3, 3] = sqrt(sum(each_rep_beta2)/(NN*nsims))
RMSE_sim[nrow(hyperpara)+3, 4] = sqrt(sum(each_rep_gamma)/(NN*nsims*3))
RMSE_sim[nrow(hyperpara)+3, 5] = sqrt(sum(each_rep_s2)/(NN*nsims))

rownames(RMSE_sim) = c(paste0("h=", h),'LPML','EM','GNAR')
write.csv(round(RMSE_sim[,1:5],4), paste0("./add_simu_result/node_rmse/exp_", exp_id, "_sce_", sce_id, "_RMSE.csv"), row.names = T)

# credible interval
burnin <- 500
reppp = 1
true_labels <- all_simu_dat[[5]][[reppp]] %>% as.factor()
idx1 <- which(true_labels == 1)[1]
idx2 <- which(true_labels == 2)[1]
idx3 <- which(true_labels == 3)[1]
length(sim_res[[reppp]][[1]])
length(sim_res[[reppp]][[hyppp]]$out$Iterates[[1]])
for (hyppp in 1:length(sim_res[[reppp]])) {
  est_all_iter_b0 <- matrix(0, nrow = NNiter - burnin, ncol = K)
  est_all_iter_b1 <- matrix(0, nrow = NNiter - burnin, ncol = K)
  est_all_iter_b2 <- matrix(0, nrow = NNiter - burnin, ncol = K)
  est_all_para <- list()
  for (iterrr in (burnin + 1:(NNiter - burnin))) {
    zhat <- sim_res[[reppp]][[hyppp]]$out$Iterates[[iterrr]]$zout
    zhat <- as.factor(zhat)
    zz <- model.matrix(~ zhat - 1)
    thetahat <- sim_res[[reppp]][[hyppp]]$out$Iterates[[iterrr]]$thetaout
    node_meb <- zz %*% thetahat
    est_all_iter_b0[(iterrr - burnin),] <- node_meb[c(idx1, idx2, idx3), 1]
    est_all_iter_b1[(iterrr - burnin),] <- node_meb[c(idx1, idx2, idx3), 2]
    est_all_iter_b2[(iterrr - burnin),] <- node_meb[c(idx1, idx2, idx3), 3]
    # est_all_iter_b0[(iterrr - burnin),] <- node_meb[c(idx1, idx2, idx3), 1]
  }
  est_all_para[[1]] <- est_all_iter_b0
  est_all_para[[2]] <- est_all_iter_b1
  est_all_para[[3]] <- est_all_iter_b2
}


# estimate sigma2 in GNAR
gnar_s2hat <- numeric(6)
for (sce in 1:6) {
  cat("sce:", sce, "\n")
  if (sce %in% c(1,2)) {
    kk = 3
  } else {kk = 5}
  load(paste0("./add_simu_result_1/sim_res_GNAR_", sce, ".rda"))
  load(paste0("./data/data_design_", sce, ".RData"))
  load(paste0("./data/true_para_design_", sce, ".RData"))
  true_theta = matrix(unlist(true_para[[2]]), nrow = kk, byrow = T)
  true_sigma2 = unlist(true_para[[1]])
  
  each_rep_beta0 = 0
  each_rep_beta1 = 0
  each_rep_beta2 = 0
  each_rep_gamma = 0
  each_rep_s2 = 0
  for (reppp in 1:nsims) {
    Zhat = sim_res_GNAR[[reppp]]$group %>% as.factor()
    ZZhat = model.matrix(~ Zhat - 1)
    thetahat = sim_res_GNAR[[reppp]]$theta
    s2hat = sim_res_GNAR[[reppp]]$sigma
    node_est = ZZhat %*% thetahat
    node_est_s2 = ZZhat %*% s2hat
    ZZtrue = model.matrix(~ as.factor(all_simu_dat[[5]][[reppp]]) - 1)
    node_true = ZZtrue %*% true_theta
    node_true_s2 = ZZtrue %*% true_sigma2
    
    each_rep_beta0 = sum(coef_SE(node_est[,1], node_true[,1])) + each_rep_beta0
    each_rep_beta1 = sum(coef_SE(node_est[,2], node_true[,2])) + each_rep_beta1
    each_rep_beta2 = sum(coef_SE(node_est[,3], node_true[,3])) + each_rep_beta2
    each_rep_gamma = sum(coef_SE(node_est[,4:6], node_true[,4:6])) + each_rep_gamma
    each_rep_s2 = sum(coef_SE(node_est_s2, node_true_s2)) + each_rep_s2
  }
  cat("sce:", sqrt(sum(each_rep_s2)/(NN*nsims)), "\n")
  gnar_s2hat[sce] <- sqrt(sum(each_rep_s2)/(NN*nsims))
}
