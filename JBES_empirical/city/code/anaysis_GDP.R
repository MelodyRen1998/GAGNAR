# analyze the group result by GAGNAR, CRP, and GNAR

# =========
# load packages
# =========
setwd("../city/")
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

cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]

# ===== 从建模结果保存参数 =====
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
write.csv(round(real_res[[besth_ind]]$Dahlout$thetaout*10, 3), "./report/city_theta_GDP.csv", row.names = F)
write.csv(round(real_res[[besth_ind]]$Dahlout$sigma2out*10, 3), "./report/city_stk_s2_GDP.csv", row.names = F)
table(real_res[[besth_ind]]$Dahlout$zout)

# ==========
# graph vis
# ==========
load("./city_data/g_city.rda")  # g2
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
ZZ <- real_res[[besth_ind]]$Dahlout$zout
V(g2)$frame.color <- "white"
pdf(width=4, height=4, file = "./report/city_group_res_GDP.pdf") 
plot(g2, vertex.size=6,
     vertex.color=c("orange", "steelblue","gray", "salmon","gold","darkred")[ZZ], vertex.label = NA)
dev.off()


# ===== GAGNAR & CRP =====
load("./report/city_est_label.rda")  # ZZ_best
ZZ1_GAGNAR <- model.matrix( ~ ZZ_best- 1)  # N*K

load("./report/city_est_label_crp.rda")  # ZZ_crp
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K

# ===== SBM =====
load("./report/city_est_label_sbm.rda")  # sbm_label
ZZ1_sbm <- model.matrix( ~ sbm_label - 1)  # N*K




# ===== SBM & GAGNAR boxplot =====
load("./report/stk_est_label_sbm.rda")  # sbm_label
load("./report/stock_res_supp_full.rda")
ZZ_best = res_topstk[[2]]$Dahlout$zout %>% as.factor()
ZZ1_sbm <- model.matrix( ~ sbm_label - 1)  # N*K

# ===== GAGNAR boxplot =====
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
ZZ <- real_res[[besth_ind]]$Dahlout$zout

city_ID = rownames(Adjmat)
return_ave <- data.frame(ID = city_ID, ret = rowMeans(Ymat), group = 1)

return_ave$group[which(ZZ == 2)] <- 2
return_ave$group[which(ZZ == 3)] <- 3
return_ave$group[which(ZZ == 4)] <- 4
return_ave$group[which(ZZ == 5)] <- 5
return_ave$group[which(ZZ == 6)] <- 6
boxplot(return_ave$ret ~ return_ave$group)
return_ave$group = as.factor(return_ave$group)
pdf(width=5, height=4, file = "./report/city_box_GDP.pdf")
ggplot(return_ave, aes(x = group, y = ret, group = group)) +
  geom_boxplot() +
  xlab("Group") + ylab("Average GDP growth rate") +theme_bw()
dev.off()


load("./report/city_est_label_sbm_GDP.rda")
city_ID = rownames(Adjmat)
return_ave <- data.frame(ID = city_ID, ret = rowMeans(Ymat), group = 1)

return_ave$group[which(sbm_label == 2)] <- 2
return_ave$group[which(sbm_label == 3)] <- 3
return_ave$group[which(sbm_label == 4)] <- 4
return_ave$group[which(sbm_label == 5)] <- 5
return_ave$group[which(sbm_label == 6)] <- 6
boxplot(return_ave$ret ~ return_ave$group)
return_ave$group = as.factor(return_ave$group)
pdf(width=5, height=4, file = "./report/city_sbm_box_GDP.pdf")
ggplot(return_ave, aes(x = group, y = ret, group = group)) +
  geom_boxplot() +
  xlab("Group") + ylab("Average GDP growth rate") +theme_bw()
dev.off()


# ===== 对应城市名 =====
dijishiID <- read.csv("./city_data/dijishiID.csv")
dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ == 2)])]

# dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_best == 1)])]
# dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_crp == 1)])]
# dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_gnar == 1)])]
table(dijishiID$省[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_best == 1)])]) %>% sort(decreasing = T)
table(dijishiID$省[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_crp == 1)])]) %>% sort(decreasing = T)
table(dijishiID$省[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_gnar == 4)])]) %>% sort(decreasing = T)
table(dijishiID$省[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(sbm_label == 4)])]) %>% sort(decreasing = T)

dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(sbm_label == 1)])]

frgdp_avg <- data.frame(ID = rownames(Adjmat), ret = rowMeans(Ymat), group = 1)

frgdp_avg$group[which(sbm_label == 2)] <- 2
frgdp_avg$group[which(sbm_label == 3)] <- 3
frgdp_avg$group[which(sbm_label == 4)] <- 4
boxplot(frgdp_avg$ret ~ frgdp_avg$group)
pdf(width=5, height=4, file = "./report/city_sbm_box.pdf")
ggplot(frgdp_avg, aes(x = group, y = ret, group = group)) +
  geom_boxplot() +
  xlab("group") + ylab("Average FR/GDP") +theme_bw()
dev.off()



# ===== HPD interval =====
library(HDInterval)
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
Theta_city <- real_res[[besth_ind]]$Dahlout$thetaout
ZZ <- real_res[[besth_ind]]$Dahlout$zout

para_int <- list()
paramt_sub_list = list()
post_mean <- matrix(0, nrow = nrow(Ymat), ncol = ncol(Theta_city))  # N by p
for (pp in 1:ncol(Theta_city)) {
  cat(pp, "\n")# for each parameter
  para_mt <- matrix(0, nrow = nrow(Ymat), ncol = 1500)  # N by Niter
  for (ittt in ((burnin+1):niterations)) {
    thetahat <- real_res[[besth_ind]]$out$Iterates[[ittt]]$thetaout
    zhat <- real_res[[besth_ind]]$out$Iterates[[ittt]]$zout %>% as.factor()
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
cities <- dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat))]

idxs_inall <- c(which(cities == "无锡市"), which(cities == "杭州市"),
                which(cities == "廊坊市"), which(cities == "抚顺市"),
                which(cities == "沈阳市"), which(cities == "临沂市"))
interval_mt <- matrix(0, nrow = ncol(Theta_city), ncol = 6*3)  # six nodes
for (pp in 1:ncol(Theta_city)) {
  for (ii in 1:length(idxs_inall)) {
    stock_id <- idxs_inall[ii]
    gg <- ZZ[stock_id] %>% as.numeric()
    interval_mt[pp,((gg-1)*3 + 1:3)] <- c(post_mean[idxs_inall[ii],pp], para_int[[pp]][idxs_inall[ii],])
  }
}
write.csv(round(interval_mt,4), file = "./report/city_interval_hdi_95_GDP.csv", row.names = F)

interval_mt = read.csv("./report/city_interval_hdi_95_GDP.csv")
for (gg in 1:6) {
  est_idx = (gg - 1)*3 + 1
  interval_mt[, est_idx] = paste0(interval_mt[, est_idx], 
                                  " (", interval_mt[, est_idx + 1], ", ", 
                                  interval_mt[, est_idx + 2], ")")
}
write.csv(interval_mt, file = "./report/city_interval_hdi_95_GDP_process.csv", row.names = F)





# ===== HPD boxplot =====
library(ggpubr)
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
Theta_city <- real_res[[besth_ind]]$Dahlout$thetaout
ZZ <- real_res[[besth_ind]]$Dahlout$zout

plots_all_effects = list()
for (effect_id in 2:3) {
  # par(mfrow = c(2, 2))
  plots = list()
  for (MM_idx in 1:length(c(500))){
    MM = c(500)[MM_idx]
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
pdf(width=10, height=4, file = "./report/city_HPD_boxplot_GDP.pdf")
ggarrange(plots_all_effects[[1]], plots_all_effects[[2]],
          # labels = c("Network Effect", "Momentum Effect"),
          ncol = 2, nrow = 1)
dev.off()

# ===== visualize interval =====
library(dplyr)
library(ggplot2)
interval_mt = read.csv("./report/city_interval_hdi_95_GDP.csv")
interval_mt = as.matrix(interval_mt)
stk_long <- rbind(interval_mt[,1:3],interval_mt[,4:6],
                  interval_mt[,7:9],interval_mt[,10:12],
                  interval_mt[,13:15], interval_mt[,16:18]) %>% as.data.frame()
colnames(stk_long) <- c("est", "lower", "upper")
stk_long$group <- c(rep(1, 7), rep(2, 7), rep(3, 7), rep(4, 7), rep(5, 7), rep(6, 7)) %>%
  factor(levels = c(6,5,4,3,2,1))
stk_long$industry <- c(rep("Wuxi", 7), 
                       rep("Hangzhou", 7), 
                       rep("Langfang", 7), 
                       rep("Fushun", 7), 
                       rep("Shenyang", 7), 
                       rep("Linyi", 7)) %>% 
  factor(levels = c("Linyi", "Shenyang", "Fushun", 
                    "Langfang", "Hangzhou", "Wuxi"))
stk_long$para <- rep(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4"), 6) %>%
  factor(levels = rev(c("beta0","beta1","beta2","gamma1","gamma2","gamma3","gamma4")))

mycolor = c("#fed439", "steelblue","darkgray", "#fd7446", "skyblue", "#7ca729")
pdf(width=5, height=4, file = "./report/city_interval_95_GDP.pdf")
ggplot(stk_long, aes(x = para, group = industry)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, color = industry), 
                position = position_dodge(0.6), size = 0.8) +
  scale_colour_manual(values = mycolor, limits = rev(c("Linyi", "Shenyang", "Fushun", 
                                                       "Langfang", "Hangzhou", "Wuxi"))) +
  geom_point(aes(y = est, group = industry), color = "black", size = 0.6, 
             position = position_dodge(0.6)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.25),
        legend.background=element_blank(),
        # legend.title = element_text(face="italic"),
        legend.text = element_text(face = "italic", size =9),
        legend.title = element_blank(),
        axis.text = element_text(size=13)) +
  labs(x = "", y = "Estimates") + ylim(-1,1) +
  coord_flip() +
  scale_x_discrete(labels = rev(expression(beta[0], beta[net], beta[mom],
                                           gamma[1], gamma[2], gamma[3], gamma[4])))
dev.off()

# ===== RMSE =====
library(ggsci)

pred_rmse = read.csv("./report/pred_rmse_city_GDP.csv", header = T)
pred_rmse = pred_rmse[1:6,]
start_time_point <- 2006:2010
n_test = ncol(pred_rmse) + 1  # 和训练模型对应，n_test比实际测试集长度大1
df_pred = data.frame(ID = rep(start_time_point, nrow(pred_rmse)), 
                     value = c(t(pred_rmse[,1:(n_test - 1)])),
                     type = c(rep("GAGNAR", n_test - 1), rep("CRP", n_test - 1), rep("GNAR",n_test - 1),
                              rep("NAR", n_test - 1), rep("AR", n_test - 1), rep("ARMA", n_test - 1)))
df_pred$type <- factor(df_pred$type, levels = c("AR", "ARMA", "NAR", "GNAR","CRP","GAGNAR"))

# 由于 GNAR 在某一个test set下误差非常大（已讨论），因此这里不对 GNAR 做可视化
df_pred1 = df_pred[which(df_pred$type != "GNAR"),]
pdf(width=5, height=4, file = "./report/city_pred_rmse_train7_test7_roll_GDP.pdf")
ggplot(df_pred1, aes(x = as.integer(ID), y = value, group = type)) +
  geom_line(aes(color = type), lwd = 0.8) +
  # scale_color_manual(values = c("red", "#96cac1", "#6a93cc", "#facb4c", "#df9a96", "#7f7f7f")) +
  scale_color_simpsons() +
  # scale_color_manual(values = c("#9bc0ba", "")) +
  theme_bw() +
  ylim(0,5) +
  theme(legend.position = c(0.85,0.82),
        legend.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  scale_x_continuous(breaks = c(start_time_point)) +
  labs(x = "Start Time Point", y = "ReMSPE")
dev.off()









# ===== parameter estimate boxplot =====
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_res.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
Theta <- real_res[[besth_ind]]$Dahlout$thetaout
ZZ <- real_res[[besth_ind]]$Dahlout$zout %>% as.factor()
ZZ1 <- model.matrix( ~ ZZ - 1)  # N*K
Theta <- real_res[[besth_ind]]$Dahlout$thetaout  # K*(p+3)
node_est <- ZZ1 %*% Theta
node_net_mom = as.data.frame(node_est[,2:3])
node_net_mom$group = ZZ

boxplot(node_net_mom$V2 ~ node_net_mom$group)

# Sigma2 <- real_res[[1]]$Dahlout$sigma2out

# ===== adj matrix with block =====
library(plot.matrix)
load("./city_data/Adjmat.rda")
# load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label.rda")
# load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_label_crp.rda")
# load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label_gnar.rda")
load("./report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
Theta_city <- real_res[[besth_ind]]$Dahlout$thetaout
ZZ_best <- real_res[[besth_ind]]$Dahlout$zout
ZZ_crp <- real_res[[1]]$Dahlout$zout
load("./report/city_res_gnar_GDP.rda")
ZZ_gnar = GNAR_summary$group
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label_sbm_GDP.rda")


# 由各组密度从高到低的顺序
node_order = c(which(ZZ_best == 4), which(ZZ_best == 2), 
               which(ZZ_best == 3), which(ZZ_best == 1))
Adjmat_gagnar = Adjmat[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_best == 4), 1:sum(ZZ_best == 4)] = 8*Adjmat_gagnar[1:sum(ZZ_best == 4), 1:sum(ZZ_best == 4)]
base = sum(ZZ_best == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 2)), 
              (base + 1):(base + sum(ZZ_best == 2))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 2)), 
                                                                       (base + 1):(base + sum(ZZ_best == 2))]
base = sum(ZZ_best == 4) + sum(ZZ_best == 2)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 3)), 
              (base + 1):(base + sum(ZZ_best == 3))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 3)), 
                                                                       (base + 1):(base + sum(ZZ_best == 3))]
base = sum(ZZ_best == 4) + sum(ZZ_best == 2) + sum(ZZ_best == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 1)), 
              (base + 1):(base + sum(ZZ_best == 1))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_best == 1)), 
                                                                       (base + 1):(base + sum(ZZ_best == 1))]

pdf(width=10, height=10, file = "./report/block_heatmap_gagnar.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "GAGNAR",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

node_order = c(which(ZZ_crp == 4), which(ZZ_crp == 2), 
               which(ZZ_crp == 3), which(ZZ_crp == 1))
Adjmat_gagnar = Adjmat[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_crp == 4), 1:sum(ZZ_crp == 4)] = 8*Adjmat_gagnar[1:sum(ZZ_crp == 4), 1:sum(ZZ_crp == 4)]
base = sum(ZZ_crp == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 2)), 
              (base + 1):(base + sum(ZZ_crp == 2))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 2)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 2))]
base = sum(ZZ_crp == 4) + sum(ZZ_crp == 2)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 3)), 
              (base + 1):(base + sum(ZZ_crp == 3))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 3)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 3))]
base = sum(ZZ_crp == 4) + sum(ZZ_crp == 2) + sum(ZZ_crp == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 1)), 
              (base + 1):(base + sum(ZZ_crp == 1))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_crp == 1)), 
                                                                       (base + 1):(base + sum(ZZ_crp == 1))]
                                                                       
pdf(width=10, height=10, file = "./report/block_heatmap_crp.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "CRP",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

node_order = c(which(ZZ_gnar == 1), which(ZZ_gnar == 4), 
               which(ZZ_gnar == 3), which(ZZ_gnar == 2))
Adjmat_gagnar = Adjmat[node_order, node_order]
Adjmat_gagnar[1:sum(ZZ_gnar == 1), 1:sum(ZZ_gnar == 1)] = 8*Adjmat_gagnar[1:sum(ZZ_gnar == 1), 1:sum(ZZ_gnar == 1)]
base = sum(ZZ_gnar == 1)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 4)), 
              (base + 1):(base + sum(ZZ_gnar == 4))] = 6*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 4)), 
                                                                      (base + 1):(base + sum(ZZ_gnar == 4))]
base = sum(ZZ_gnar == 1) + sum(ZZ_gnar == 4)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 3)), 
              (base + 1):(base + sum(ZZ_gnar == 3))] = 4*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 3)), 
                                                                      (base + 1):(base + sum(ZZ_gnar == 3))]
base = sum(ZZ_gnar == 1) + sum(ZZ_gnar == 4) + sum(ZZ_gnar == 3)
Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 2)), 
              (base + 1):(base + sum(ZZ_gnar == 2))] = 2*Adjmat_gagnar[(base + 1):(base + sum(ZZ_gnar == 2)), 
                                                                       (base + 1):(base + sum(ZZ_gnar == 2))]
pdf(width=10, height=10, file = "./report/block_heatmap_gnar.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "GNAR",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

node_order = c(which(sbm_label == 1), which(sbm_label == 2), 
               which(sbm_label == 3), which(sbm_label == 4))
Adjmat_gagnar = Adjmat[node_order, node_order]
Adjmat_gagnar[1:sum(sbm_label == 1), 1:sum(sbm_label == 1)] = 8*Adjmat_gagnar[1:sum(sbm_label == 1), 1:sum(sbm_label == 1)]
base = sum(sbm_label == 1)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 2)), 
              (base + 1):(base + sum(sbm_label == 2))] = 6*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 2)), 
                                                                       (base + 1):(base + sum(sbm_label == 2))]
base = sum(sbm_label == 1) + sum(sbm_label == 2)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 3)), 
              (base + 1):(base + sum(sbm_label == 3))] = 4*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 3)), 
                                                                       (base + 1):(base + sum(sbm_label == 3))]
base = sum(sbm_label == 1) + sum(sbm_label == 2) + sum(sbm_label == 3)
Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 4)), 
              (base + 1):(base + sum(sbm_label == 4))] = 2*Adjmat_gagnar[(base + 1):(base + sum(sbm_label == 4)), 
                                                                         (base + 1):(base + sum(sbm_label == 4))]
                                                                       
pdf(width=10, height=10, file = "./report/block_heatmap_sbm.pdf")
plot(Adjmat_gagnar, border=NA, col=c("white", "grey","gold", "steelblue", "#fd7446", "darkred", "skyblue"),
     breaks = c(sort(unique(c(Adjmat_gagnar-0.5))),100),
     xlab="", ylab="", main = "SBM",
     key=NULL, axis.row = list(las = 2, tick=FALSE),
     axis.col = list(las = 2, tick=FALSE),
     cex.axis = 0.4, frame.plot=T)
dev.off()

# NetworkDensity = function(mat) {
#   dens = sum(mat)/(nrow(mat)^2 - nrow(mat))
#   cat(dens, "\n")
#   return(dens)
# }

# ===== 组内平均密度与各组分别的密度 =====
group_dens = matrix(nrow = 9, ncol = 9)  # 4种方法（位于每列）/ 每行：4 groups + in-block + between block
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(Adjmat)

for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  group_dens[gg, 1] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 1], decreasing = T))
group_dens[1:6, 1] <- sort(group_dens[1:6, 1], decreasing = T)
cat("GAGNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[7, 1] <- inblock_edge_all/inblock_size_all
group_dens[8, 1] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:7) {
  inblock_edge = sum(Adjmat[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  group_dens[gg, 2] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:7, 2], decreasing = T))
group_dens[1:7, 2] <- sort(group_dens[1:7, 2], decreasing = T)
cat("CRP", inblock_edge_all/inblock_size_all, "\n")
group_dens[8, 2] <- inblock_edge_all/inblock_size_all
group_dens[9, 2] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  group_dens[gg, 3] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 3], decreasing = T))
group_dens[1:6, 3] <- sort(group_dens[1:6, 3], decreasing = T)
cat("GNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[7, 3] <- inblock_edge_all/inblock_size_all
group_dens[8, 3] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  group_dens[gg, 4] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:6, 4], decreasing = T))
group_dens[1:6, 4] <- sort(group_dens[1:6, 4], decreasing = T)
cat("SBM", inblock_edge_all/inblock_size_all, "\n")
group_dens[7, 4] <- inblock_edge_all/inblock_size_all
group_dens[8, 4] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

write.csv(round(group_dens, 4), "./report/block_density_GDP.csv")

# 组间平均密度
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(Adjmat)
for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GAGNAR", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:7) {
  inblock_edge = sum(Adjmat[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("CRP", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GNAR", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:6) {
  inblock_edge = sum(Adjmat[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("SBM", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")






# # ===== inblock distance =====
# load("./city_data/dmat.rda")  # d.matrix
# ddmat = d.matrix
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(ZZ_best == gg), which(ZZ_best == gg)])
#   inblock_size = (sum(ZZ_best == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("GAGNAR", inblock_dis_all/inblock_size_all, "\n")
# 
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(ZZ_crp == gg), which(ZZ_crp == gg)])
#   inblock_size = (sum(ZZ_crp == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("CRP", inblock_dis_all/inblock_size_all, "\n")
# 
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(sbm_label == gg), which(sbm_label == gg)])
#   inblock_size = (sum(sbm_label == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("SBM", inblock_dis_all/inblock_size_all, "\n")
# 
# 
# # 组间
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(ZZ_best == gg), which(ZZ_best == gg)])
#   inblock_size = (sum(ZZ_best == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("GAGNAR", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")
# 
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(ZZ_crp == gg), which(ZZ_crp == gg)])
#   inblock_size = (sum(ZZ_crp == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("CRP", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")
# 
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
#   inblock_size = (sum(ZZ_gnar == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("GNAR", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")
# 
# 
# inblock_dis_all = 0
# inblock_size_all = 0
# for (gg in 1:4) {
#   inblock_dis = sum(d.matrix[which(sbm_label == gg), which(sbm_label == gg)])
#   inblock_size = (sum(sbm_label == gg))^2
#   inblock_dis_all = inblock_dis_all + inblock_dis
#   inblock_size_all = inblock_size_all + inblock_size
# }
# cat("SBM", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")
# 


# ===== GDP result GAGNAR =====
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_res_GDP.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
# besth_ind <- 2  # index of the optimal h
best_h <- real_res[[besth_ind]]$h
# (3) predict
ZZ <- real_res[[besth_ind]]$Dahlout$zout %>% as.factor()
ZZ1 <- model.matrix( ~ ZZ - 1)  # N*K
Theta <- real_res[[besth_ind]]$Dahlout$thetaout  # K*(p+3)
node_est <- ZZ1 %*% Theta
Sigma2 <- real_res[[1]]$Dahlout$sigma2out


# ===== crp =====
Theta_crp = real_res[[1]]$Dahlout$thetaout
ZZ_crp <- real_res[[1]]$Dahlout$zout %>% as.factor()
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K
node_est_crp <- ZZ1_crp %*% Theta_crp
Sigma2_crp <- real_res[[1]]$Dahlout$sigma2out

