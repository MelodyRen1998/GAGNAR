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
