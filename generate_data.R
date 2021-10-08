# ======== generate data of scenario 1 in example 3 ====== #
# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("igraph", "mvtnorm", "MCMCpack", "parallel")
ipak(packages)
setwd("YOUR.WORK.PATH")
source("./util_code/gendata.R")

## ==================basic setting================== ##
p <- 3  # dimension of attributes
k <- 6  # number of clusters
t <- 20  # time length
nsims <- 100 # simulation replicates 

## ==================Parameter of Scenario 1 of Example 3================ ##
# set group parameters
group.sigma <- list()
group.theta <- list()
sigma0 <- diag(1, nrow = p + 3)  # identity matrix
## group 1
sigma1_true <- 2
theta1 <- rep(-3, p+3)
theta1[1] <- 5  # beta0k
theta1[2] <- 0.2  # beta1k
theta1[3] <- 0.1  # beta2k
theta1[4:6] <- c(0.5,0.7,1)
## group 2
sigma2_true <- 1
theta2 <- rep(1, p+3)
theta2[1] <- -5  # beta0
theta2[2] <- -0.4  # beta1
theta2[3] <- 0.2  # beta2
theta2[4:6] <- c(0.1,0.9,0.4)
## group 3
sigma3_true <- 3
theta3 <- rep(10, p+3)
theta3[1] <- 0  # beta0
theta3[2] <- 0.2  # beta1
theta3[3] <- 0.4  # beta2
theta3[4:6] <- c(0.2,-1,2)
## group 4
sigma4_true <- 4
theta4 <- rep(10, p+3)
theta4[1] <- 3  # beta0
theta4[2] <- 0.1  # beta1
theta4[3] <- 0.2  # beta2
theta4[4:6] <- c(1,-1,1.5)
## group 5
sigma5_true <- 2
theta5 <- rep(10, p+3)
theta5[1] <- -3 # beta0
theta5[2] <- 0.5  # beta1
theta5[3] <- 0.2  # beta2
theta5[4:6] <- c(0.8,0.5,-2)

## group 6
sigma6_true <- 3
theta6 <- rep(10, p+3)
theta6[1] <- 2  # beta0
theta6[2] <- -0.6  # beta1
theta6[3] <- -0.2  # beta2
theta6[4:6] <- c(-0.8,0.5,2)
## assemble group parameter list
group.sigma <- list(sigma1_true, sigma2_true, sigma3_true, 
                    sigma4_true, sigma5_true, sigma6_true)
group.theta <- list(theta1, theta2, theta3, theta4, theta5, theta6)

## ==================data of Example 3================ ##
# load data
load("./data/g3_simu.rda")
load("./data/label3_sbm_simu.rda")
load("./data/Adjmat_stk_simu.rda")

Y_tilde_list <- list()
X_tilde_list <- list()
d_matrix <- list()
W_list <- list()
label_list <- list()
n.nodes <- length(label_SBM)
sc_stk_km <- as.factor(label_SBM)
# visualize the graph
mycolor = c("#FFCE3E", "#203E5F","#81B214","#D56073","skyblue","#757270")
# pdf(width=4, height=4, file = "../simulation_realgraph/report/g_stk.pdf")
V(g.top2_simu)$frame.color <- "white"
plot(g.top2_simu, vertex.size=7,
     vertex.color=mycolor[1:length(unique(sc_stk_km))][sc_stk_km], 
     vertex.label = NA)
# dev.off()

# generate data
for (iisim in 1:nsims){
  cat("generate data of replicate", iisim, "\r")
  seed <- iisim + 123 # myseed same for different h, varies with sim replicates
  res <- gendata_K6(n.nodes, p, t, group.sigma, group.theta, myseed = seed,
                 adj.mt.simu, label_SBM, g.top2_simu)
  Y_tilde_list[[iisim]] <- res$Y_tilde
  X_tilde_list[[iisim]] <- res$X_tilde
  d_matrix[[iisim]] <- res$d_matrix
  W_list[[iisim]] <- res$W
  label_list[[iisim]] = label_SBM
}
all_simu_dat <- list(Y_tilde_list, X_tilde_list, d_matrix, W_list, label_list)
true_para <- list(group.sigma, group.theta)
# save(true_para, file = "./data_sbm/true_para_design5.RData")
# save(all_simu_dat, file = "./data_sbm/simu_data_design5.RData")
