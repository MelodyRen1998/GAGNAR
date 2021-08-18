# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("mvtnorm", "mclust", "plyr")
ipak(packages)

# pls set the work directory as the path of `simulation_pack` folder
setwd("~/renyimeng/GAGNAR")
# setwd("YOUR_wd")

# source("./util_code/gwcrp.R")
# source("./util_code/Dahl.R")
# source("./util_code/LPML_VAR.R")
source("./util_code/estimator.R")
source("./util_code/SIMU_func.R")

h <- c(seq(0,2,0.2),3,4,5)
alpha <- c(1)
hyperpara <- expand.grid(h = h, alpha = alpha)
# ==========
# run parallel
# ==========
ncores <- detectCores() - 4
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores, type="FORK")

for (senar in 1:2) {
  idx_design <- senar
  # load data and true parameters
  load(paste0("./data/data_design_", idx_design, ".RData")) # data
  load(paste0("./data/true_para_design_", idx_design, ".RData")) # para
  Y_tilde_list <- all_simu_dat[[1]]
  X_tilde_list <- all_simu_dat[[2]]
  d_matrix <- all_simu_dat[[3]]
  W_list <- all_simu_dat[[4]]
  label_list <- all_simu_dat[[5]]
  ## ==================basic setting================== ##
  n.nodes <- ncol(all_simu_dat[[1]][[1]])  # number of nodes
  p <- 3  # dimension of attributes
  kk <- 3  # number of clusters
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
  sim_res_h <- list()
  sim_res_EM = list()
  sim_res_GNAR = list()
  
  ## 1.
  sim_res <- clusterApply(cl, 1:nsims, sim, Y_tilde_list, X_tilde_list, d_matrix, hyperpara,
                          niterations, a0, b0, tau0, sigma0, initNClusters, burnin)
  save(sim_res, file = paste0("./result/sim_res_", idx_design, ".rda"))
  ## find the best h
  best_hs = c()
  for (reppp in 1:nsims) {
    lpml <- lapply(sim_res[[reppp]], function(each_rep) {each_rep$LPML}) %>% unlist()
    best_hs[reppp] <- hyperpara$h[which.max(lpml)]
  }
  best_h = mean(best_hs)
  # cat("best h is", best_h)
  ## 2.
  sim_res_h <- clusterApply(cl, 1:nsims, sim_best_h, best_h, Y_tilde_list, X_tilde_list, d_matrix,
                            niterations, a0, b0, tau0, sigma0, initNClusters, burnin)
  save(sim_res_h, file = paste0("./result/sim_res_h_", idx_design, ".rda"))
  ## 3.
  sim_res_EM <- clusterApply(cl, 1:nsims, simEM, kk, Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list, niterations)
  save(sim_res_EM, file = paste0("./result/sim_res_EM_", idx_design, ".rda"))
  ## 4.
  sim_res_GNAR <- clusterApply(cl, 1:nsims, simGNAR, kk, Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list)
  save(sim_res_GNAR, file = paste0("./result/sim_res_GNAR_", idx_design, ".rda"))
}



nsims <- 100 # simulation replicates 
for (sce in 1:6) {
  if (sce %in% c(1,2)) {
    kk = 3
  } else {kk = 5}
  cat("sce:", sce, "\n")
  load(paste0("./data/data_design_", sce, ".RData"))
  load(paste0("./data/true_para_design_", sce, ".RData"))
  
  Y_tilde_list <- all_simu_dat[[1]]
  X_tilde_list <- all_simu_dat[[2]]
  d_matrix <- all_simu_dat[[3]]
  W_list <- all_simu_dat[[4]]
  label_list <- all_simu_dat[[5]]
  sim_res_GNAR = list()
  for(ss in 1:nsims) {
    cat("replicate:", ss, "\r")
    sim_res_GNAR[[ss]] <- simGNAR(ss, kk, Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list)
  }
  save(sim_res_GNAR, file = paste0("./add_simu_result_1/sim_res_GNAR_", sce, ".rda"))
}

