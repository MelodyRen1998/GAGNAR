# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("igraph", "mvtnorm", "MCMCpack", "parallel", "mclust", "plyr", "fcd", "doParallel")
ipak(packages)

setwd("YOUR.WORK.PATH")

source("./util_code/gwcrp.R")
source("./util_code/Dahl.R")
source("./util_code/LPML_VAR.R")
# source("./util_code/estimator.R")
source("./util_code/SIMU_func.R")

h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)
# ==========
# run parallel
# ==========
# parallel setting
ncores <- detectCores() - 4
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores, type="FORK")

for (senar in 1:2) {
  idx_design <- senar+4
  cat("simulation design:", idx_design, "\n")
  # load data and true parameters
  load(paste0("./data_sbm/simu_data_design", idx_design, ".RData")) # data
  load(paste0("./data_sbm/true_para_design", idx_design, ".RData")) # para
  Y_tilde_list <- all_simu_dat[[1]]
  X_tilde_list <- all_simu_dat[[2]]
  d_matrix <- all_simu_dat[[3]]
  W_list <- all_simu_dat[[4]]
  label_list <- all_simu_dat[[5]]
  ## ==================basic setting================== ##
  n.nodes <- ncol(all_simu_dat[[1]][[1]])  # number of nodes
  p <- 3  # dimension of attributes
  kk <- 6  # number of clusters
  t <- 20  # time length
  nsims <- 100 # simulation replicates 
  group.sigma <- true_para[[1]]
  group.theta <- true_para[[2]]
  ## ==================hyper parameter of Design 6================ ##
  a0 <- 0.01
  b0 <- 0.01
  tau0 <- rep(0, p + 3)
  sigma0 <- diag(100, p + 3)
  niterations <- 1500 # number of iterations for MCMC
  initNClusters <- 10
  burnin <- 500 # number of burnin iterations for MCMC
  ## ======== 
  ## simulation 
  ## ========
  sim_res <- list()
  sim_res_EM = list()
  sim_res_GNAR = list()
  
  for(ss in 1:nsims) {
    sim_res <- sim(ss, Y_tilde_list, X_tilde_list, d_matrix, hyperpara,
                   niterations, a0, b0, tau0, sigma0, initNClusters, burnin)
  }
  
  sim_res <- clusterApply(cl, 1:nsims, sim, Y_tilde_list, X_tilde_list, d_matrix, hyperpara,
                          niterations, a0, b0, tau0, sigma0, initNClusters, burnin)
  save(sim_res, file = paste0("./result/sim_res_", idx_design, ".rda"))
}
