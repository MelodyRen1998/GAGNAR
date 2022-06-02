# analyze the group result by GAGNAR, CRP, and GNAR

# =========
# load packages
# =========
setwd("../stock/")
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
load("./stock_data/rda_full_sparse/Ymat.rda")
load("./stock_data/rda_full_sparse/Zmat.rda")
load("./stock_data/rda_full_sparse/Adjmat.rda")
cat("network density is", sum(adj.mt)/(nrow(adj.mt)^2 - nrow(adj.mt)), "\n")
load("./stock_data/rda_full_sparse/dmat.rda")  # d.matrix
Zmat <- scale(Zmat)

# =========
# experiment setting
# =========

train_window = 40
test_window = 6
n_window = 1  # rolling window size
train_start = 1
n_test = 1
flag = TRUE
test_days = NULL

train_end = train_start + train_window - 1
test_start = train_end
test_end = test_start + test_window

# ===== GAGNAR & CRP =====
load("./report/stock_est_label.rda")  # ZZ_best
ZZ1_GAGNAR <- model.matrix( ~ ZZ_best- 1)  # N*K

load("./report/stock_est_label_crp.rda")  # ZZ_crp
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K

# ===== GNAR =====
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
K = 5
GNAR_summary = Cluster.NAR(Ymat_train, W, Zmat, K = 5, method="complete", seed = n_test)
ind = order(GNAR_summary$alpha)
group_cl = order(ind)[GNAR_summary$group]
ZZ_gnar <- group_cl %>% as.factor()
ZZ1_GNAR <- model.matrix( ~ ZZ_gnar - 1)  # N*K

save(ZZ_gnar, file = "./report/stock_est_label_gnar.rda")


# ===== SBM & GAGNAR boxplot =====
load("./report/stk_est_label_sbm.rda")  # sbm_label
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
ZZ1_sbm <- model.matrix( ~ sbm_label - 1)  # N*K


return_ave <- data.frame(ID = stk_ID_names, ret = rowMeans(Ymat), group = 1)

return_ave$group[which(sbm_label == 2)] <- 2
return_ave$group[which(sbm_label == 3)] <- 3
return_ave$group[which(sbm_label == 4)] <- 4
return_ave$group[which(sbm_label == 5)] <- 5
return_ave$group[which(sbm_label == 6)] <- 6
boxplot(return_ave$ret ~ return_ave$group)
return_ave$group = as.factor(return_ave$group)
pdf(width=5, height=4, file = "./report/stock_sbm_box.pdf")
ggplot(return_ave, aes(x = group, y = ret, group = group)) +
  geom_boxplot() +
  xlab("Group") + ylab("Average Stock Return") +theme_bw()
dev.off()


return_ave$group[which(ZZ_best == 2)] <- 2
return_ave$group[which(ZZ_best == 3)] <- 3
return_ave$group[which(ZZ_best == 4)] <- 4
return_ave$group[which(ZZ_best == 5)] <- 5
return_ave$group[which(ZZ_best == 6)] <- 6
boxplot(return_ave$ret ~ return_ave$group)
pdf(width=5, height=4, file = "./report/stk_group_box_full.pdf")
return_ave$group = as.factor(return_ave$group)
ggplot(return_ave, aes(x = group, y = ret, group = group)) +
  geom_boxplot() +
  xlab("Group") + ylab("Average Stock Return") +theme_bw()
dev.off()

# ===== 从建模结果保存参数 =====

load("./report/stock_res_supp_full.rda")
write.csv(round(res_topstk[[2]]$Dahlout$thetaout*10, 3), "./report/stk_theta.csv", row.names = F)
write.csv(round(res_topstk[[2]]$Dahlout$sigma2out*10, 3), "./report/stk_s2.csv", row.names = F)

load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/res/est_res_gnar_alpha_0_1_sparse_3.rda")
gnar_theta = GNAR_summary$theta
gnar_s2 = GNAR_summary$sig2
gnar_memb = table(GNAR_summary$group)
gnar_res = rbind(gnar_memb, t(gnar_theta), gnar_s2)
write.csv(round(gnar_res*10, 3), "./report/stk_est_gnar.csv", row.names = F)

load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/res/est_res_crp_alpha_0_1_sparse_3.rda")
crp_memb = table(crp_res[[1]])
crp_est_res = rbind(crp_memb, t(crp_res[[2]]), crp_res[[3]])
write.csv(round(crp_est_res*10, 3), "./report/stk_est_crp.csv", row.names = F)

# ===== GAGNAR, CRP, GNAR 对应股票 =====
library(readxl)
company_data <- read_excel("./stock_data/TRD_Co.xlsx")
company_data$Stkcd <- as.numeric(company_data$Stkcd)
load("./stock_data/rda_full_sparse/stockid.rda")
stock_type <- company_data$Nindnme[match(stk_ID_names, company_data$Stkcd)]

load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
load("./res/est_res_crp_alpha_0_1_sparse_3.rda")
ZZ_crp = crp_res[[1]]
load("./res/est_res_gnar_alpha_0_1_sparse_3.rda")
ZZ_gnar = GNAR_summary$group %>% as.factor()

for (type_id in 1:max(as.numeric(ZZ_best))) {
  print(sort(table(stock_type[ZZ_best == type_id]), decreasing = T)[1:5])
  print("=======")
}
for (type_id in 1:max(as.numeric(ZZ_crp))) {
  print(sort(table(stock_type[ZZ_crp == type_id]), decreasing = T)[1:5])
  print("=======")
}
for (type_id in 1:max(as.numeric(ZZ_gnar))) {
  print(sort(table(stock_type[ZZ_gnar == type_id]), decreasing = T)[1:5])
  print("=======")
}


# ===== graph vis =====
load("./stock_data/rda_full_sparse/g_stock.rda")
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
V(g.top2)$frame.color <- "white"
pdf(width=4, height=4, file = "./report/stk_group_res.pdf") 
plot(g.top2, vertex.size=6,
     vertex.color=c("#fed439", "steelblue","darkgray", "#fd7446", "skyblue", "#7ca729")[ZZ_best], vertex.label = NA)
dev.off()

# ===== interval =====
# ===== calculate interval =====
library(HDInterval)
load("./report/stock_res_supp_full.rda")
# load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stock_res_supp_full_train40_test6_alpha_1.rda")
# load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stock_res_supp_full_train52.rda")
besth_ind <- which.max(unlist(lapply(res_topstk, function(x){x$LPML})))
# besth <- res_topstk[[3]]$h
Theta_stk <- res_topstk[[2]]$Dahlout$thetaout
ZZ_best <- res_topstk[[2]]$Dahlout$zout
para_int <- list()
paramt_sub_list = list()

post_mean <- matrix(0, nrow = nrow(Ymat), ncol = ncol(Theta_stk))  # N by p
for (pp in 1:ncol(Theta_stk)) {
  cat(pp, "\n")# for each parameter
  para_mt <- matrix(0, nrow = nrow(Ymat), ncol = 1500)  # N by Niter
  for (ittt in ((burnin+1):niterations)) {
    thetahat <- res_topstk[[2]]$out$Iterates[[ittt]]$thetaout
    zhat <- res_topstk[[2]]$out$Iterates[[ittt]]$zout %>% as.factor()
    zzhat <- model.matrix(~ zhat - 1)
    node_est <- zzhat %*% thetahat
    para_mt[, ittt] <- node_est[,pp]
  }
  paramt_sub <- para_mt[, (burnin+1):niterations]
  paramt_sub_list[[pp]] = paramt_sub
  post_mean[,pp] <- rowMeans(paramt_sub)
  node_interval <- t(apply(paramt_sub, 1, hdi, credMass = 0.95))
  para_int[[pp]] <- node_interval
  # print(para_int[[pp]][c(1,50,100),])  # select three nodes
}

# 逐个样本查找0.95interval 不包含0的
## 结论1：1000次迭代算出来的 全部network都是包含0的
good_nodes = c()
for (ii in 1:nrow(Ymat)) {
  good_est_net = ((para_int[[2]][ii,]["upper"] < 0) | (para_int[[2]][ii,]["lower"] > 0))
  good_est_mom = ((para_int[[3]][ii,]["upper"] < 0) | (para_int[[3]][ii,]["lower"] > 0))
  if ((good_est_net > 0) | (good_est_mom > 0)) {
  # if ((good_est_net > 0)) {
    good_nodes = c(good_nodes, ii)
  }
}

# for (gg in 1:6) {
#   print(which(!is.na(match(good_nodes, which(ZZ_best == gg))))[1])
# }


idxs_inall <- c(which((stock_type == "房地产业") & (ZZ_best == 1))[1],
                which((stock_type == "信息技术业") & (ZZ_best == 2))[1],
                which((stock_type == "医药制造业") & (ZZ_best == 3))[1],
                which((stock_type == "交通运输设备制造业") & (ZZ_best == 4))[1],
                which((stock_type == "计算机应用服务业") & (ZZ_best == 5))[1],
                which((stock_type == "电力、蒸汽、热水的生产和供应业") & (ZZ_best == 6))[1])
idxs_inall = c()
for (gg in 1:6) {
  # set.seed(gg+1)
  if (sum(!is.na(match(good_nodes, which(ZZ_best == gg)))) == 1) {
    idxs_inall[gg] = good_nodes[which(!is.na(match(good_nodes, which(ZZ_best == gg))))]
  } else {
    idxs_inall[gg] = good_nodes[sample(c(which(!is.na(match(good_nodes, which(ZZ_best == gg))))), 1)]
  }
}
stock_type[idxs_inall]
load("./report/six_stock_idx.rda")
interval_mt <- matrix(0, nrow = ncol(Theta_stk), ncol = 6*3)  # six nodes
for (pp in 1:ncol(Theta_stk)) {
  for (ii in 1:length(idxs_inall)) {
    stock_id <- idxs_inall[ii]
    gg <- ZZ_best[stock_id] %>% as.numeric()
    interval_mt[pp,((gg-1)*3 + 1:3)] <- c(post_mean[idxs_inall[ii],pp], para_int[[pp]][idxs_inall[ii],])
  }
}
write.csv(round(interval_mt,4), file = "./report/stk_interval_hdi_95.csv", row.names = F)

interval_mt = read.csv("./report/stk_interval_hdi_95.csv")
interval_mt = round(interval_mt*10, 3)
for (gg in 1:6) {
  est_idx = (gg - 1)*3 + 1
  interval_mt[, est_idx] = paste0(interval_mt[, est_idx], 
                                  " (", interval_mt[, est_idx + 1], ", ", 
                                  interval_mt[, est_idx + 2], ")")
}
write.csv(interval_mt, file = "./report/stk_interval_hdi_95_process.csv", row.names = F)



# ===== HPD boxplot =====
load("./report/stock_res_supp_full.rda")
Theta_stk <- res_topstk[[2]]$Dahlout$thetaout
ZZ = res_topstk[[2]]$Dahlout$zout
plots_all_effects = list()
for (effect_id in 2:3) {
  # par(mfrow = c(2, 2))
  plots = list()
  for (MM_idx in 1:length(c(1000))){
    MM = c(1000)[MM_idx]
    node_net_long_allgroup = NULL
    
    for (gg in 1:6) {
      idxs = which(ZZ == gg)
      # cities_group1 = dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ == gg)])]
      # idxs <- which(cities %in% cities_group1)
      
      # network effect
      node_net = t(paramt_sub_list[[effect_id]][idxs,1:MM])
      # node_net_long = data.frame(node_net = c(node_net), param_idx = rep(1:5, each = MM))
      node_net_long = data.frame(node_net = c(node_net))
      node_net_long$group = gg
      node_net_long_allgroup = rbind(node_net_long_allgroup, node_net_long)
    }
    # node_net_long_allgroup$param_idx = as.factor(node_net_long_allgroup$param_idx)
    node_net_long_allgroup$group = as.factor(node_net_long_allgroup$group)
    plots[[MM_idx]] = ggplot(node_net_long_allgroup, aes(x = group, y = node_net)) +
      geom_boxplot(aes(fill = group)) +
      labs(x = "Group", y = ifelse(effect_id == 2, "Network Effect", "Momentum Effect"), 
           title = paste0("Iteration ", MM, " Times"))
    # boxplot(node_net_long$node_net ~ node_net_long$param_idx, xlab = "Node", 
    #         ylab = ifelse(effect_id == 2, "Network Effect", "Momentum Effect"), 
    #         main = paste0("Group ", gg, ", Iteration ", MM, " Times"))
  }
  plots_all_effects[[effect_id-1]] = ggplot(node_net_long_allgroup, aes(x = group, y = node_net)) +
    geom_boxplot() +
    # scale_fill_manual(values = c("#fed439", "steelblue","darkgray", "#fd7446", "skyblue", "#7ca729")) +
    labs(x = "Group", y = ifelse(effect_id == 2, "Network Effect", "Momentum Effect"), 
         title = paste0("Iteration ", MM, " Times after Burn-in")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text =element_text(size=12),
          axis.title = element_text(size = 12))
}
pdf(width=10, height=4, file = "./report/stk_HPD_boxplot.pdf")
ggarrange(plots_all_effects[[1]], plots_all_effects[[2]],
          # labels = c("Network Effect", "Momentum Effect"),
          ncol = 2, nrow = 1)
dev.off()

# ===== visualize interval =====
idxs_inall = c()
for (gg in 1:6) {
  # set.seed(gg+1)
  if (sum(!is.na(match(good_nodes, which(ZZ_best == gg)))) == 1) {
    idxs_inall[gg] = good_nodes[which(!is.na(match(good_nodes, which(ZZ_best == gg))))]
  } else {
    idxs_inall[gg] = good_nodes[sample(c(which(!is.na(match(good_nodes, which(ZZ_best == gg))))), 1)]
  }
}
stock_type[idxs_inall]
save(idxs_inall, file = "./report/six_stock_idx.rda")
load("./report/six_stock_idx.rda")
interval_mt <- matrix(0, nrow = ncol(Theta_stk), ncol = 6*3)  # six nodes
for (pp in 1:ncol(Theta_stk)) {
  for (ii in 1:length(idxs_inall)) {
    stock_id <- idxs_inall[ii]
    gg <- ZZ_best[stock_id] %>% as.numeric()
    interval_mt[pp,((gg-1)*3 + 1:3)] <- c(post_mean[idxs_inall[ii],pp], para_int[[pp]][idxs_inall[ii],])
  }
}

library(dplyr)
library(ggplot2)
stk_long <- rbind(interval_mt[,1:3],interval_mt[,4:6],
                   interval_mt[,7:9],interval_mt[,10:12],
                   interval_mt[,13:15], interval_mt[,16:18]) %>% as.data.frame()
colnames(stk_long) <- c("est", "lower", "upper")
stk_long$group <- c(rep(1, 9), rep(2, 9), rep(3, 9), rep(4, 9), rep(5, 9), rep(6, 9)) %>%
  factor(levels = c(6,5,4,3,2,1))
stk_long$industry <- c(rep("Social Services", 9), 
                       rep("IT (Group 2)", 9), 
                       rep("IT (Group 3)", 9), 
                       rep("Electrical", 9), 
                       rep("IT (Group 5)", 9), 
                       rep("Metal", 9)) %>% 
  factor(levels = c("Metal", "IT (Group 5)", "Electrical", 
                    "IT (Group 3)", "IT (Group 2)", "Social Services"))
stk_long$para <- rep(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4","gamma5","gamma6"), 6) %>%
  factor(levels = rev(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4","gamma5","gamma6")))

mycolor = c("#fed439", "steelblue","darkgray", "#fd7446", "skyblue", "#7ca729")
pdf(width=5, height=4, file = "./report/stk_interval_95_full.pdf")
ggplot(stk_long, aes(x = para, group = industry)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, color = industry), 
                position = position_dodge(0.6), size = 0.8) +
  scale_colour_manual(values = mycolor, limits = rev(c("Metal", "IT (Group 5)", "Electrical", 
                                                       "IT (Group 3)", "IT (Group 2)", "Social Services"))) +
  geom_point(aes(y = est, group = industry), color = "black", size = 0.6, 
             position = position_dodge(0.6)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.25),
        legend.background=element_blank(),
        # legend.title = element_text(face="italic"),
        legend.text = element_text(face = "italic", size =9),
        legend.title = element_blank(),
        axis.text = element_text(size=13)) +
  labs(x = "", y = "Estimates") + ylim(-0.2,0.12) +
  coord_flip() +
  scale_x_discrete(labels = rev(expression(beta[0], beta[net], beta[mom],
                                           gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], gamma[6])))
dev.off()


# ===== trace plot (two chains) =====
## trace plot for 4 cities
library(gridExtra)
library(ggplot2)
library(reshape2)
library(ggpubr)
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stock_res_supp_full.rda")
real_res = res_topstk
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stock_res_supp_full_chain2_a_0025_b_0015.rda")
real_res_chain2 = res_topstk

# dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_best == 4)])]
idxs_inall <- c(which((stock_type == "房地产业") & (ZZ_best == 1))[1],
                which((stock_type == "信息技术业") & (ZZ_best == 2))[1],
                which((stock_type == "医药制造业") & (ZZ_best == 3))[1],
                which((stock_type == "交通运输设备制造业") & (ZZ_best == 4))[1],
                which((stock_type == "计算机应用服务业") & (ZZ_best == 5))[1],
                which((stock_type == "电力、蒸汽、热水的生产和供应业") & (ZZ_best == 6))[1])
fig_title = c( "Real Estate",  "IT","Pharmaceutical","Transportation","Computer Service",  "Resources")
for (cc in 1:6) {print(stock_type[idxs_inall[cc]])}
plots <- list()
titles <- expression(beta[0], beta[net], beta[mom],
                     gamma[SIZE], gamma[BM], gamma[PR], 
                     gamma[AR], gamma[LEV], gamma[CASH], sigma^2)
for (cc in 1:length(idxs_inall)) {
  cat(stock_type[idxs_inall[cc]], "\n")
  
  # theta trace plot
  for (pp in 1:ncol(Theta_stk)) {
    myscale <- ifelse(pp %in% c(2,3,8,10), 0.5, 1)
    # myscale = 1
    cat(pp, "\n")# for each parameter
    para_mt_c1 <- matrix(0, nrow = nrow(Ymat), ncol = niterations)  # N by Niter
    para_mt_c2 <- matrix(0, nrow = nrow(Ymat), ncol = niterations)  # N by Niter
    for (ittt in ((burnin+1):niterations)) {
      # chain 1
      besth_ind = 2
      thetahat_c1 <- real_res[[besth_ind]]$out$Iterates[[ittt]]$thetaout
      zhat_c1 <- real_res[[besth_ind]]$out$Iterates[[ittt]]$zout %>% as.factor()
      zzhat_c1 <- model.matrix(~ zhat_c1 - 1)
      node_est_c1 <- zzhat_c1 %*% thetahat_c1
      para_mt_c1[, ittt] <- node_est_c1[,pp]*myscale
      # chain 2
      thetahat_c2 <- real_res_chain2[[1]]$out$Iterates[[ittt]]$thetaout
      zhat_c2 <- real_res_chain2[[1]]$out$Iterates[[ittt]]$zout %>% as.factor()
      zzhat_c2 <- model.matrix(~ zhat_c2 - 1)
      node_est_c2 <- zzhat_c2 %*% thetahat_c2
      para_mt_c2[, ittt] <- node_est_c2[,pp]*myscale
    }
    para_trace <- data.frame(ID = 1:1000+burnin, 
                             est_theta_c1 = para_mt_c1[idxs_inall[cc], (burnin+1):niterations],
                             est_theta_c2 = para_mt_c2[idxs_inall[cc], (burnin+1):niterations])
    para_trace <- melt(para_trace, id.vars = "ID", variable.name = "chain")
    # plot the trace for a parameter
    myplt <- ggplot(para_trace, aes(x = ID, y = value, color = chain)) +
      geom_line(aes(group = chain, color = chain), alpha = 0.7) + theme_bw() +
      labs(x = "", y = paste0("est*", myscale), title = titles[pp]) +
      scale_x_continuous(breaks = c(500,1000,1500)) +
      scale_color_discrete(labels = c("Chain 1", "Chain 2"), name = "") +
      ylim(-0.1, 0.1)+
      theme(legend.position = c(1, 1.05),
            legend.justification = c(1, 1),
            legend.background = element_blank())
    plots[[pp]] <- myplt
  }
  
  # sigma2 trace plot
  myscale <- 0.5
  s2_mt_c1 <- matrix(0, nrow = nrow(Ymat), ncol = niterations)  # N by Niter
  s2_mt_c2 <- matrix(0, nrow = nrow(Ymat), ncol = niterations)  # N by Niter
  for (ittt in ((burnin+1):niterations)) {
    s2hat_c1 <- real_res[[besth_ind]]$out$Iterates[[ittt]]$sigma2out
    zhat_c1 <- real_res[[besth_ind]]$out$Iterates[[ittt]]$zout %>% as.factor()
    zzhat_c1 <- model.matrix(~ zhat_c1 - 1)
    s2_mt_c1[, ittt] <- zzhat_c1 %*% s2hat_c1*myscale
    # chain 2
    s2hat_c2 <- real_res_chain2[[1]]$out$Iterates[[ittt]]$sigma2out
    zhat_c2 <- real_res_chain2[[1]]$out$Iterates[[ittt]]$zout %>% as.factor()
    zzhat_c2 <- model.matrix(~ zhat_c2 - 1)
    s2_mt_c2[, ittt] <- zzhat_c2 %*% s2hat_c2*myscale
  }
  s2_trace <- data.frame(ID = 1:1000+burnin, 
                         est_s2_c1 = s2_mt_c1[idxs_inall[cc], (burnin+1):niterations],
                         est_s2_c2 = s2_mt_c2[idxs_inall[cc], (burnin+1):niterations])
  s2_trace <- melt(s2_trace, id.vars = "ID", variable.name = "chain")
  
  myplt <- ggplot(s2_trace, aes(x = ID, y = value, color = chain)) +
    geom_line(aes(color = chain), alpha = 0.7) + theme_bw() +
    labs(x = "", y = paste0("est*", myscale), title = titles[10]) +
    ylim(0,0.01) +
    scale_x_continuous(breaks = c(500,1000,1500)) +
    theme(legend.position = "none")
  
  plots[[ncol(Theta_stk) + 1]] <- myplt

  g0 = ggarrange(plotlist = plots, nrow = 2 , ncol = 5,common.legend = TRUE, legend = "bottom")
  g_res = annotate_figure(g0, top = text_grob(fig_title[cc] ,size = 13))
  g_res
  ggsave(file=paste0("./report/trace_plt_", fig_title[cc], ".pdf"), g_res, width = 8, height = 4) #saves g
  # do.call("ggarrange", c(plots, ncol = 4,common.legend = TRUE, legend = "right", label.x = "Iterations", label.y = "Values", top = "Taizhou"))
  # dev.off()
}

# ===== p-values of GNAR =====
p = 6
K = 6
load("./res/est_res_gnar_alpha_0_1_sparse_3.rda")
theta_gnar = GNAR_summary$theta
p_values_gnar = matrix(0, nrow = K, ncol = p+3)
for (kk in 1:K) {
  p_values_gnar[kk, ] = (1 - pnorm(abs(theta_gnar[kk, ])/sqrt(diag(GNAR_summary$cov[[kk]]))))*2
}
write.csv(round(p_values_gnar, 4), "./report/gnar_p_values_train40.csv", row.names = F)



# ===== adj matrix with block =====
load("./stock_data/rda_full_sparse/Adjmat.rda")
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
load("./res/est_res_crp_alpha_0_1_sparse_3.rda")
ZZ_crp = crp_res[[1]]
load("./res/est_res_gnar_alpha_0_1_sparse_3.rda")
ZZ_gnar = GNAR_summary$group %>% as.factor()
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stk_est_label_sbm.rda")


# 组内平均密度
inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GAGNAR", inblock_edge_all/inblock_size_all, "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("CRP", inblock_edge_all/inblock_size_all, "\n")

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GNAR", inblock_edge_all/inblock_size_all, "\n")

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GAGNAR", inblock_edge_all/inblock_size_all, "\n")


# 组间平均密度
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(adj.mt)
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GAGNAR", (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("CRP", (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GNAR", (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("SBM", (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")

# ===== 计算每组和整体组内的平均密度 & 组间平均密度 =====

load("./stock_data/rda_full_sparse/Adjmat.rda")
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
load("./res/est_res_crp_alpha_0_1_sparse_3.rda")
ZZ_crp = crp_res[[1]]
load("./res/est_res_gnar_alpha_0_1_sparse_3.rda")
ZZ_gnar = GNAR_summary$group %>% as.factor()
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stk_est_label_sbm.rda")

group_dens = matrix(nrow = 12, ncol = 4)  # 4种方法（位于每列）/ 每行：6 groups(最多的CRP 分了10组) + in-block + between block
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(adj.mt)

for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  group_dens[gg, 1] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 1], decreasing = T))
group_dens[1:6, 1] <- sort(group_dens[1:6, 1], decreasing = T)
cat("GAGNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[11, 1] <- inblock_edge_all/inblock_size_all
group_dens[12, 1] <- (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:10) {
  inblock_edge = sum(adj.mt[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  group_dens[gg, 2] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:10, 2], decreasing = T))
group_dens[1:10, 2] <- sort(group_dens[1:10, 2], decreasing = T)
cat("CRP", inblock_edge_all/inblock_size_all, "\n")
group_dens[11, 2] <- inblock_edge_all/inblock_size_all
group_dens[12, 2] <- (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  group_dens[gg, 3] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 3], decreasing = T))
group_dens[1:6, 3] <- sort(group_dens[1:6, 3], decreasing = T)
cat("GNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[11, 3] <- inblock_edge_all/inblock_size_all
group_dens[12, 3] <- (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(adj.mt[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  group_dens[gg, 4] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 4], decreasing = T))
group_dens[1:6, 4] <- sort(group_dens[1:6, 4], decreasing = T)
cat("SBM", inblock_edge_all/inblock_size_all, "\n")
group_dens[11, 4] <- inblock_edge_all/inblock_size_all
group_dens[12, 4] <- (sum(adj.mt) - inblock_edge_all)/ (N^2 - inblock_size_all)

write.csv(round(group_dens, 4), "./report/block_density.csv")


# ===== adj matrix with block =====
library(plot.matrix)
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/stock_data/rda_full_sparse/stockid.rda")
rownames(adj.mt) = colnames(adj.mt) = stk_ID_names

load("./stock_data/rda_full_sparse/Adjmat.rda")
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
load("./res/est_res_crp_alpha_0_1_sparse_3.rda")
ZZ_crp = crp_res[[1]]
load("./res/est_res_gnar_alpha_0_1_sparse_3.rda")
ZZ_gnar = GNAR_summary$group %>% as.factor()
load("~/renyimeng/GAGNAR_140/JBES_empirical/stock/report/stk_est_label_sbm.rda")


# 注意分组的排序和由大到小的密度对应
node_order = c(which(ZZ_best == 6), which(ZZ_best == 1), 
               which(ZZ_best == 2), which(ZZ_best == 4),
               which(ZZ_best == 3), which(ZZ_best == 5))
Adjmat_gagnar = adj.mt[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_best == 6), 1:sum(ZZ_best == 6)] = 12*Adjmat_gagnar[1:sum(ZZ_best == 6), 1:sum(ZZ_best == 6)]
base = sum(ZZ_best == 6)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 1)), 
              (base + 1):(base + sum(ZZ_best == 1))] = 10*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 1)), 
                                                                       (base + 1):(base + sum(ZZ_best == 1))]
base = sum(ZZ_best == 6) + sum(ZZ_best == 1)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 2)), 
              (base + 1):(base + sum(ZZ_best == 2))] = 8*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 2)), 
                                                                       (base + 1):(base + sum(ZZ_best == 2))]
base = sum(ZZ_best == 6) + sum(ZZ_best == 1) + sum(ZZ_best == 2)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 4)), 
              (base + 1):(base + sum(ZZ_best == 4))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 4)), 
                                                                       (base + 1):(base + sum(ZZ_best == 4))]
                                                                       
base = sum(ZZ_best == 6) + sum(ZZ_best == 1) + sum(ZZ_best == 2) + sum(ZZ_best == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 3)), 
              (base + 1):(base + sum(ZZ_best == 3))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 3)), 
                                                                       (base + 1):(base + sum(ZZ_best == 3))]

base = sum(ZZ_best == 6) + sum(ZZ_best == 1) + sum(ZZ_best == 2) + sum(ZZ_best == 4) + sum(ZZ_best == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 5)), 
              (base + 1):(base + sum(ZZ_best == 5))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 5)), 
                                                                       (base + 1):(base + sum(ZZ_best == 5))]

pdf(width=10, height=10, file = "./report/block_heatmap_gagnar.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue","purple"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "GAGNAR",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

node_order = c(which(ZZ_crp == 10), which(ZZ_crp == 4), 
               which(ZZ_crp == 1), which(ZZ_crp == 7),
               which(ZZ_crp == 3), which(ZZ_crp == 6),
               which(ZZ_crp == 5), which(ZZ_crp == 2),
               which(ZZ_crp == 8), which(ZZ_crp == 9))
Adjmat_gagnar = adj.mt[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_crp == 10), 1:sum(ZZ_crp == 10)] = 20*Adjmat_gagnar[1:sum(ZZ_crp == 10), 1:sum(ZZ_crp == 10)]
base = sum(ZZ_crp == 10)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 4)), 
              (base + 1):(base + sum(ZZ_crp == 4))] = 18*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 4)), 
                                                                        (base + 1):(base + sum(ZZ_crp == 4))]
base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 1)), 
              (base + 1):(base + sum(ZZ_crp == 1))] = 16*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 1)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 1))]
base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 7)), 
              (base + 1):(base + sum(ZZ_crp == 7))] = 14*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 7)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 7))]

base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 3)), 
              (base + 1):(base + sum(ZZ_crp == 3))] = 12*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 3)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 3))]
base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7) + sum(ZZ_crp == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 6)), 
              (base + 1):(base + sum(ZZ_crp == 6))] = 10*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 6)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 6))]

base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7) + sum(ZZ_crp == 3) + sum(ZZ_crp == 6)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 5)), 
              (base + 1):(base + sum(ZZ_crp == 5))] = 8*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 5)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 5))]
base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7) + sum(ZZ_crp == 3) + sum(ZZ_crp == 6) + sum(ZZ_crp == 5)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 2)), 
              (base + 1):(base + sum(ZZ_crp == 2))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 2)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 2))]
base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7) + sum(ZZ_crp == 3) + sum(ZZ_crp == 6) + sum(ZZ_crp == 5) +
  sum(ZZ_crp == 2)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 8)), 
              (base + 1):(base + sum(ZZ_crp == 8))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 8)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 8))]

base = sum(ZZ_crp == 10) + sum(ZZ_crp == 4) + sum(ZZ_crp == 1) +
  sum(ZZ_crp == 7) + sum(ZZ_crp == 3) + sum(ZZ_crp == 6) + sum(ZZ_crp == 5) +
  sum(ZZ_crp == 2) + sum(ZZ_crp == 8)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 9)), 
              (base + 1):(base + sum(ZZ_crp == 9))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 9)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 9))]
pdf(width=10, height=10, file = "./report/block_heatmap_crp.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue","purple", 
                                     "orange", "green", "pink", "brown"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "CRP",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()


node_order = c(which(ZZ_gnar == 2), which(ZZ_gnar == 1), 
               which(ZZ_gnar == 3), which(ZZ_gnar == 4),
               which(ZZ_gnar == 6), which(ZZ_gnar == 5))
Adjmat_gagnar = adj.mt[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_gnar == 2), 1:sum(ZZ_gnar == 2)] = 12*Adjmat_gagnar[1:sum(ZZ_gnar == 2), 1:sum(ZZ_gnar == 2)]
base = sum(ZZ_gnar == 2)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 1)), 
              (base + 1):(base + sum(ZZ_gnar == 1))] = 10*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 1)), 
                                                                        (base + 1):(base + sum(ZZ_gnar == 1))]
base = sum(ZZ_gnar == 2) + sum(ZZ_gnar == 1)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 3)), 
              (base + 1):(base + sum(ZZ_gnar == 3))] = 8*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 3)), 
                                                                       (base + 1):(base + sum(ZZ_gnar == 3))]
base = sum(ZZ_gnar == 2) + sum(ZZ_gnar == 1) + sum(ZZ_gnar == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 4)), 
              (base + 1):(base + sum(ZZ_gnar == 4))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 4)), 
                                                                       (base + 1):(base + sum(ZZ_gnar == 4))]

base = sum(ZZ_gnar == 2) + sum(ZZ_gnar == 1) + sum(ZZ_gnar == 3) + sum(ZZ_gnar == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 6)), 
              (base + 1):(base + sum(ZZ_gnar == 6))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 6)), 
                                                                       (base + 1):(base + sum(ZZ_gnar == 6))]

base = sum(ZZ_gnar == 2) + sum(ZZ_gnar == 1) + sum(ZZ_gnar == 3) + sum(ZZ_gnar == 4) + sum(ZZ_gnar == 6)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 5)), 
              (base + 1):(base + sum(ZZ_gnar == 5))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 5)), 
                                                                       (base + 1):(base + sum(ZZ_gnar == 5))]

pdf(width=10, height=10, file = "./report/block_heatmap_gnar.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue","purple"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "GAGNAR",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

node_order = c(which(sbm_label == 6), which(sbm_label == 4), 
               which(sbm_label == 2), which(sbm_label == 3),
               which(sbm_label == 1), which(sbm_label == 5))
Adjmat_gagnar = adj.mt[node_order, node_order]
Adjmat_gagnar[1:sum(sbm_label == 6), 1:sum(sbm_label == 6)] = 12*Adjmat_gagnar[1:sum(sbm_label == 6), 1:sum(sbm_label == 6)]
base = sum(sbm_label == 6)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 4)), 
              (base + 1):(base + sum(sbm_label == 4))] = 10*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 4)), 
                                                                        (base + 1):(base + sum(sbm_label == 4))]
base = sum(sbm_label == 6) + sum(sbm_label == 4)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 2)), 
              (base + 1):(base + sum(sbm_label == 2))] = 8*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 2)), 
                                                                       (base + 1):(base + sum(sbm_label == 2))]
base = sum(sbm_label == 6) + sum(sbm_label == 4) + sum(sbm_label == 2)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 3)), 
              (base + 1):(base + sum(sbm_label == 3))] = 6*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 3)), 
                                                                       (base + 1):(base + sum(sbm_label == 3))]

base = sum(sbm_label == 6) + sum(sbm_label == 4) + sum(sbm_label == 2) + sum(sbm_label == 3)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 1)), 
              (base + 1):(base + sum(sbm_label == 1))] = 4*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 1)), 
                                                                       (base + 1):(base + sum(sbm_label == 1))]

base = sum(sbm_label == 6) + sum(sbm_label == 4) + sum(sbm_label == 2) + sum(sbm_label == 3) + sum(sbm_label == 1)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 5)), 
              (base + 1):(base + sum(sbm_label == 5))] = 2*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 5)), 
                                                                       (base + 1):(base + sum(sbm_label == 5))]


pdf(width=10, height=10, file = "./report/block_heatmap_sbm.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue","purple"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "SBM",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()
