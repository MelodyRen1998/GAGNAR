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

K = 6
eig_res = eigen(Adjmat)
U = eig_res$vectors[, 1:K]
km = kmeans(U, 6, nstart = 10)
sbm_label = as.factor(km$cluster)
save(sbm_label, file = "./report/city_est_label_sbm_GDP.rda")

