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
load("./city_data/Ymat.rda")
load("./city_data/Xmat.rda")
load("./city_data/Adjmat.rda")
load("./city_data/dmat.rda")  # d.matrix

cov_scale <- apply(coviriates, 2, rescale)  # normalize to [0,1]
NN = nrow(Ymat)
Zmat = cov_scale[(nrow(cov_scale) - NN+1):nrow(cov_scale),]


# =========
# experiment setting
# =========
train_window = 12
test_window = 12
n_window = 1
train_start = 5
n_test = 1
flag = TRUE
test_days = NULL
pred_rmse = matrix(0, nrow = 7, ncol = 100)
p_value = matrix(0, nrow = 7, ncol = p+3)
# basic setting & tuning param
niterations <- 1500 # number of iterations for MCMC
initNClusters <- 10
burnin <- 500 # number of burnin iterations for MCMC
h <- c(0,2) # crp and optimal h
# h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)

train_end = train_start + train_window - 1
test_start = train_start + 1
test_end = test_start + test_window - 1

# ===== GAGNAR & CRP =====
load("./report/city_est_label.rda")  # ZZ_best
ZZ1_GAGNAR <- model.matrix( ~ ZZ_best- 1)  # N*K

load("./report/city_est_label_crp.rda")  # ZZ_crp
ZZ1_crp <- model.matrix( ~ ZZ_crp - 1)  # N*K

# ===== GNAR =====
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
ind = order(GNAR_summary$alpha)
group_cl = order(ind)[GNAR_summary$group]
ZZ_gnar <- group_cl %>% as.factor()
ZZ1_GNAR <- model.matrix( ~ ZZ_gnar - 1)  # N*K

save(ZZ_gnar, file = "./report/city_est_label_gnar.rda")

# ===== SBM =====
load("./report/city_est_label_sbm.rda")  # sbm_label
ZZ1_sbm <- model.matrix( ~ sbm_label - 1)  # N*K


# ===== 对应城市名 =====
dijishiID <- read.csv("./city_data/dijishiID.csv")
dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ_best == 1)])]

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


# ===== 计算每个节点的pvalue =====
setwd("~/renyimeng/GAGNAR_140/JBES_empirical/city")
load("./report/city_res_12_5.rda")  # h = 2 的估计结果
node_param_est = list()
for (iter in 501:1500) {
  ZZ = real_res[[1]]$out$Iterates[[iter]]$zout %>% as.factor()
  ZZ1 = model.matrix( ~ ZZ - 1)  # N*K
  Theta = real_res[[1]]$out$Iterates[[iter]]$thetaout
  node_param_est[[iter - 500]] = ZZ1 %*% Theta  # N*p
}
post_est_ZZ = real_res[[1]]$Dahlout$zout %>% as.factor()
post_est_ZZ1 = model.matrix( ~ post_est_ZZ - 1)  # N*K
post_est_Theta = real_res[[1]]$Dahlout$thetaout
post_node_param_est = post_est_ZZ1 %*% post_est_Theta

symbol_node_param = matrix(0, nrow = nrow(post_node_param_est), ncol = ncol(post_node_param_est))
symbol = matrix(0, nrow = 1000, ncol = ncol(post_node_param_est))
for (ii in 1:nrow(post_node_param_est)) {
  for (iter in 1:1000) {
    symbol[iter, ] = node_param_est[[iter]][ii,] * post_node_param_est[ii,]
  }
  symbol_node_param[ii,] = colSums(symbol < 0)
}
p_value_node = symbol_node_param / 1000


# ===== HPD interval =====
library(HDInterval)
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_res.rda")
besth_ind <- which.max(unlist(lapply(real_res, function(x){x$LPML})))
Theta_city <- real_res[[besth_ind]]$Dahlout$thetaout
ZZ <- real_res[[besth_ind]]$Dahlout$zout %>% as.factor()

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


plots_all_effects = list()
for (effect_id in 2:3) {
  # par(mfrow = c(2, 2))
  plots = list()
  for (MM_idx in 1:length(c(1000))){
    MM = c(1000)[MM_idx]
    node_net_long_allgroup = NULL
    
    for (gg in 1:4) {
      dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ == gg)])]
      cities_group1 = dijishiID$市[which(as.character(dijishiID$市代码) %in% rownames(Adjmat)[which(ZZ == gg)])]
      idxs <- which(cities %in% cities_group1)
      
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
      labs(x = "Group", y = ifelse(effect_id == 2, "Network Effect", "Momentum Effect"), title = paste0("Iteration ", MM, " Times"))
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
  # ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
  #           # labels = c("A", "B", "C"),
  #           ncol = 2, nrow = 2)
}
pdf(width=10, height=4, file = "./report/city_HPD_boxplot.pdf")
ggarrange(plots_all_effects[[1]], plots_all_effects[[2]],
          # labels = c("Network Effect", "Momentum Effect"),
          ncol = 2, nrow = 1)
dev.off()



interval_mt <- matrix(0, nrow = ncol(Theta_city), ncol = 10*1000)  # 10 nodes, 1000 iteration
for (pp in 1:ncol(Theta_city)) {
  for (ii in 1:length(idxs)) {
    stock_id <- idxs[ii]
    gg <- ZZ_best[stock_id] %>% as.numeric()
    interval_mt[pp,((gg-1)*3 + 1:3)] <- c(post_mean[idxs[ii],pp], para_int[[pp]][idxs[ii],])
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
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label.rda")
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_label_crp.rda")
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label_gnar.rda")
load("~/renyimeng/GAGNAR_140/JBES_empirical/city/report/city_est_label_sbm.rda")

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
group_dens = matrix(nrow = 6, ncol = 4)  # 4种方法（位于每列）/ 每行：4 groups + in-block + between block
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(Adjmat)

for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  group_dens[gg, 1] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:4, 1], decreasing = T))
group_dens[1:4, 1] <- sort(group_dens[1:4, 1], decreasing = T)
cat("GAGNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[5, 1] <- inblock_edge_all/inblock_size_all
group_dens[6, 1] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  group_dens[gg, 2] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:4, 2], decreasing = T))
group_dens[1:4, 2] <- sort(group_dens[1:4, 2], decreasing = T)
cat("CRP", inblock_edge_all/inblock_size_all, "\n")
group_dens[5, 2] <- inblock_edge_all/inblock_size_all
group_dens[6, 2] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  group_dens[gg, 3] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:4, 3], decreasing = T))
group_dens[1:4, 3] <- sort(group_dens[1:4, 3], decreasing = T)
cat("GNAR", inblock_edge_all/inblock_size_all, "\n")
group_dens[5, 3] <- inblock_edge_all/inblock_size_all
group_dens[6, 3] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  group_dens[gg, 4] <- inblock_edge/inblock_size
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat(order(group_dens[1:4, 4], decreasing = T))
group_dens[1:4, 4] <- sort(group_dens[1:4, 4], decreasing = T)
cat("SBM", inblock_edge_all/inblock_size_all, "\n")
group_dens[5, 4] <- inblock_edge_all/inblock_size_all
group_dens[6, 4] <- (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all)

write.csv(round(group_dens, 4), "./report/block_density.csv")

# 组间平均密度
inblock_edge_all = 0
inblock_size_all = 0
N = nrow(Adjmat)
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GAGNAR", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("CRP", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")

inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("GNAR", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_edge_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_edge = sum(Adjmat[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_edge_all = inblock_edge_all + inblock_edge
  inblock_size_all = inblock_size_all + inblock_size
  # NetworkDensity(Adjmat[which(ZZ_best == gg), which(ZZ_best == gg)])
}
cat("SBM", (sum(Adjmat) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")






# ===== inblock distance =====
load("./city_data/dmat.rda")  # d.matrix
ddmat = d.matrix

inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("GAGNAR", inblock_dis_all/inblock_size_all, "\n")


inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("CRP", inblock_dis_all/inblock_size_all, "\n")


inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("SBM", inblock_dis_all/inblock_size_all, "\n")


# 组间
inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(ZZ_best == gg), which(ZZ_best == gg)])
  inblock_size = (sum(ZZ_best == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("GAGNAR", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(ZZ_crp == gg), which(ZZ_crp == gg)])
  inblock_size = (sum(ZZ_crp == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("CRP", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(ZZ_gnar == gg), which(ZZ_gnar == gg)])
  inblock_size = (sum(ZZ_gnar == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("GNAR", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")


inblock_dis_all = 0
inblock_size_all = 0
for (gg in 1:4) {
  inblock_dis = sum(d.matrix[which(sbm_label == gg), which(sbm_label == gg)])
  inblock_size = (sum(sbm_label == gg))^2
  inblock_dis_all = inblock_dis_all + inblock_dis
  inblock_size_all = inblock_size_all + inblock_size
}
cat("SBM", (sum(d.matrix) - inblock_edge_all)/ (N^2 - inblock_size_all), "\n")



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

