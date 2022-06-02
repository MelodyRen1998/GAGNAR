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
load("./stock_data/rda_full_sparse/Ymat.rda")
load("./stock_data/rda_full_sparse/Zmat.rda")
load("./stock_data/rda_full_sparse/Adjmat.rda")
cat("network density is", sum(adj.mt)/(nrow(adj.mt)^2 - nrow(adj.mt)), "\n")
load("./stock_data/rda_full_sparse/dmat.rda")  # d.matrix
Zmat <- scale(Zmat)

K = 6
eig_res = eigen(adj.mt)
U = eig_res$vectors[, 1:K]
km = kmeans(U, 6, nstart = 10)
sbm_label = as.factor(km$cluster)
save(sbm_label, file = "./report/stk_est_label_sbm.rda")
