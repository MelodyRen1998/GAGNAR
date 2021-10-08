# packages needed: "igraph", "mvtnorm", "MCMCpack"
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("igraph", "mvtnorm", "MCMCpack", "blockmodels")
ipak(packages)

###########################################################################
## Generate network data with 1 lag
## (no time-dependent covariates)
###########################################################################

gendata <- function(n.nodes, p, t, group.sigma, group.theta, myseed, AA, LL, g)
{
  #===========================================================#
  ## Input: n.nodes = number of nodes in network ##
  ##        n.edges = number of edges in network ##
  ##        p = attributes(covariates) dimension ##
  ##        t = observed time length ##
  ##        group.sigma = a list containing sigma^2 from all groups ##
  ##        group.theta = a list containing theta from all groups ##
  ##        z = matrix of cluster assignment, n.nodes by n.cluster, 1 indicates node belongs to the cluster, 0 else.
  ##        myseed = random seed
  ## Output:
  ##        a list containing:
  ##        [["X_tilde"]], [["Y_tilde"]], [["d_matrix"]]
  #===========================================================#
  
  ## generate graph with nodes and edges
  set.seed(myseed)
  ## SBM generate graph
  # set.seed(123)
  # node_label <- sample(1:3, n.nodes, replace = T)
  # inblock_prob <- 20/n.nodes
  # between_prob <- 2/n.nodes
  # adj.mt <- matrix(0, nrow = n.nodes, ncol = n.nodes)
  # for (k in 1:3) {  # for each block
  #   idx <- which(node_label == k)  # index for nodes in the group k
  #   n_within <- length(idx)
  #   n_between <- n.nodes - n_within
  #   
  #   adj_k <- sapply(idx, function(x) {
  #     adj_node <- rep(0, n.nodes)
  #     adj_node[idx] <- rbinom(n_within, 1, prob = inblock_prob)  # nodes within a block
  #     adj_node[setdiff(1:n.nodes, idx)] <- rbinom(n_between, 1, prob = between_prob)  # nodes not in a block
  #     # adj.mt[x, ] <- adj_node  # neighbors for each node in group k
  #     return(adj_node)
  #   })
  #   adj.mt[idx, ] <- t(adj_k)
  # }
  # diag(adj.mt) <- 0
  W <- AA/rowSums(AA)  # row-normalized adj matrix
  z <- model.matrix(~as.factor(LL) - 1)
  d.matrix <- distances(g, v = V(g), to = V(g))
  
  ## assume the nodes attributes follow MVN(0,1)
  sigma_attr <- diag(1, nrow = p)
  set.seed(myseed)
  Vattr <- rmvnorm(n.nodes, mean = rep(0, p), sigma = sigma_attr)
  
  ## generate data Y_tilde
  ### from t=2 to t=T, generate n.nodes at each time step. Follow the equation (2.3).
  Dk <- list(diag(z[, 1]), diag(z[, 2]), diag(z[, 3]), diag(z[, 4]), diag(z[, 5]))
  beta0k <- list(group.theta[[1]][1], group.theta[[2]][1], group.theta[[3]][1], group.theta[[4]][1], group.theta[[5]][1])
  beta1k <- list(group.theta[[1]][2], group.theta[[2]][2], group.theta[[3]][2], group.theta[[4]][2], group.theta[[5]][2])
  beta2k <- list(group.theta[[1]][3], group.theta[[2]][3], group.theta[[3]][3], group.theta[[4]][3], group.theta[[5]][3])
  gammak <- list(group.theta[[1]][4:(p+3)], group.theta[[2]][4:(p+3)], group.theta[[3]][4:(p+3)], group.theta[[4]][4:(p+3)], group.theta[[5]][4:(p+3)])
  node_degree <- 1/rowSums(AA)
  node_degree[which(node_degree == Inf)] <- 10e6  # process the Inf
  W <- diag(node_degree) %*% AA  # standerized adjacency matrix
  ### \cal{Beta}_0
  Beta0 <- numeric(n.nodes)
  for (i in 1:k) { Beta0 <- Beta0 + Dk[[i]] %*% (beta0k[[i]] * rep(1, n.nodes) + Vattr %*% gammak[[i]]) }
  ### \cal{Beta}_1
  Beta1 <- numeric(n.nodes)
  for (i in 1:k) { Beta1 <- Beta1 + Dk[[i]] %*% (beta1k[[i]] * diag(rep(1, n.nodes))) }
  ### \cal{Beta}_2
  Beta2 <- numeric(n.nodes)
  for (i in 1:k) { Beta2 <- Beta2 + Dk[[i]] %*% (beta2k[[i]] * diag(rep(1, n.nodes))) }
  G <- Beta1 %*% W + Beta2
  true_sigma_2 <- unlist(group.sigma)  # unlist to a vector
  delta <- z %*% sqrt(true_sigma_2)
  ## generate white noise (length of T) for each node
  wn <- sapply(1:n.nodes, FUN = function(x) {set.seed(x + myseed); ts(rnorm(t))})
  
  ## initialize Y
  Y <- matrix(0, nrow = t, ncol = n.nodes)
  ## generate Y from t=2 to t=T
  for (time in 2:t) {
    Y[time, ] <- Beta0 + G %*% Y[(time-1), ] + delta * wn[time,]
  }
  
  ## calculate X_tilde
  X.tilde <- array(0, dim = c(t, n.nodes, (p + 3)))
  for (time in 2:t) {
    sub_Xtilde <- matrix(c(rep(1, n.nodes), W %*% t(Y)[,(time-1)], t(Y)[,(time-1)]), byrow = F, nrow = n.nodes)
    X.tilde[time,,] <- cbind(sub_Xtilde, Vattr)  # Vattr: [n.nodes * p]
  }
  ## X_tilde at t=1
  X.tilde[1,,] <- cbind(matrix(c(rep(1, n.nodes), rep(0, n.nodes), rep(0, n.nodes)), byrow = F, nrow = n.nodes), Vattr)
  
  res <- list(X_tilde = X.tilde, Y_tilde = Y, d_matrix = d.matrix, node_label = LL, W = W)
  return(res)
}

gendata_K3 <- function(n.nodes, p, t, group.sigma, group.theta, myseed, AA, LL, g)
{
  #===========================================================#
  ## Input: n.nodes = number of nodes in network ##
  ##        n.edges = number of edges in network ##
  ##        p = attributes(covariates) dimension ##
  ##        t = observed time length ##
  ##        group.sigma = a list containing sigma^2 from all groups ##
  ##        group.theta = a list containing theta from all groups ##
  ##        z = matrix of cluster assignment, n.nodes by n.cluster, 1 indicates node belongs to the cluster, 0 else.
  ##        myseed = random seed
  ## Output:
  ##        a list containing:
  ##        [["X_tilde"]], [["Y_tilde"]], [["d_matrix"]]
  #===========================================================#
  
  ## generate graph with nodes and edges
  set.seed(myseed)
  W <- AA/rowSums(AA)  # row-normalized adj matrix
  z <- model.matrix(~as.factor(LL) - 1)
  d.matrix <- distances(g, v = V(g), to = V(g))
  
  ## assume the nodes attributes follow MVN(0,1)
  sigma_attr <- diag(1, nrow = p)
  set.seed(myseed)
  Vattr <- rmvnorm(n.nodes, mean = rep(0, p), sigma = sigma_attr)
  
  ## generate data Y_tilde
  ### from t=2 to t=T, generate n.nodes at each time step. Follow the equation (2.3).
  Dk <- list(diag(z[, 1]), diag(z[, 2]), diag(z[, 3]))
  beta0k <- list(group.theta[[1]][1], group.theta[[2]][1], group.theta[[3]][1])
  beta1k <- list(group.theta[[1]][2], group.theta[[2]][2], group.theta[[3]][2])
  beta2k <- list(group.theta[[1]][3], group.theta[[2]][3], group.theta[[3]][3])
  gammak <- list(group.theta[[1]][4:(p+3)], group.theta[[2]][4:(p+3)], group.theta[[3]][4:(p+3)])
  node_degree <- 1/rowSums(AA)
  node_degree[which(node_degree == Inf)] <- 10e6  # process the Inf
  W <- diag(node_degree) %*% AA  # standerized adjacency matrix
  ### \cal{Beta}_0
  Beta0 <- numeric(n.nodes)
  for (i in 1:k) { Beta0 <- Beta0 + Dk[[i]] %*% (beta0k[[i]] * rep(1, n.nodes) + Vattr %*% gammak[[i]]) }
  ### \cal{Beta}_1
  Beta1 <- numeric(n.nodes)
  for (i in 1:k) { Beta1 <- Beta1 + Dk[[i]] %*% (beta1k[[i]] * diag(rep(1, n.nodes))) }
  ### \cal{Beta}_2
  Beta2 <- numeric(n.nodes)
  for (i in 1:k) { Beta2 <- Beta2 + Dk[[i]] %*% (beta2k[[i]] * diag(rep(1, n.nodes))) }
  G <- Beta1 %*% W + Beta2
  true_sigma_2 <- unlist(group.sigma)  # unlist to a vector
  delta <- z %*% sqrt(true_sigma_2)
  ## generate white noise (length of T) for each node
  wn <- sapply(1:n.nodes, FUN = function(x) {set.seed(x + myseed); ts(rnorm(t))})
  
  ## initialize Y
  Y <- matrix(0, nrow = t, ncol = n.nodes)
  ## generate Y from t=2 to t=T
  for (time in 2:t) {
    Y[time, ] <- Beta0 + G %*% Y[(time-1), ] + delta * wn[time,]
  }
  
  ## calculate X_tilde
  X.tilde <- array(0, dim = c(t, n.nodes, (p + 3)))
  for (time in 2:t) {
    sub_Xtilde <- matrix(c(rep(1, n.nodes), W %*% t(Y)[,(time-1)], t(Y)[,(time-1)]), byrow = F, nrow = n.nodes)
    X.tilde[time,,] <- cbind(sub_Xtilde, Vattr)  # Vattr: [n.nodes * p]
  }
  ## X_tilde at t=1
  X.tilde[1,,] <- cbind(matrix(c(rep(1, n.nodes), rep(0, n.nodes), rep(0, n.nodes)), byrow = F, nrow = n.nodes), Vattr)
  
  res <- list(X_tilde = X.tilde, Y_tilde = Y, d_matrix = d.matrix, node_label = LL, W = W)
  return(res)
}


gendata_K4 <- function(n.nodes, p, t, group.sigma, group.theta, myseed, AA, LL, g)
{
  #===========================================================#
  ## Input: n.nodes = number of nodes in network ##
  ##        n.edges = number of edges in network ##
  ##        p = attributes(covariates) dimension ##
  ##        t = observed time length ##
  ##        group.sigma = a list containing sigma^2 from all groups ##
  ##        group.theta = a list containing theta from all groups ##
  ##        z = matrix of cluster assignment, n.nodes by n.cluster, 1 indicates node belongs to the cluster, 0 else.
  ##        myseed = random seed
  ## Output:
  ##        a list containing:
  ##        [["X_tilde"]], [["Y_tilde"]], [["d_matrix"]]
  #===========================================================#
  
  ## generate graph with nodes and edges
  set.seed(myseed)
  W <- AA/rowSums(AA)  # row-normalized adj matrix
  z <- model.matrix(~as.factor(LL) - 1)
  d.matrix <- distances(g, v = V(g), to = V(g))
  
  ## assume the nodes attributes follow MVN(0,1)
  sigma_attr <- diag(1, nrow = p)
  set.seed(myseed)
  Vattr <- rmvnorm(n.nodes, mean = rep(0, p), sigma = sigma_attr)
  
  ## generate data Y_tilde
  ### from t=2 to t=T, generate n.nodes at each time step. Follow the equation (2.3).
  Dk <- list(diag(z[, 1]), diag(z[, 2]), diag(z[, 3]), diag(z[, 4]))
  beta0k <- list(group.theta[[1]][1], group.theta[[2]][1], group.theta[[3]][1], group.theta[[4]][1])
  beta1k <- list(group.theta[[1]][2], group.theta[[2]][2], group.theta[[3]][2], group.theta[[4]][2])
  beta2k <- list(group.theta[[1]][3], group.theta[[2]][3], group.theta[[3]][3], group.theta[[4]][3])
  gammak <- list(group.theta[[1]][4:(p+3)], group.theta[[2]][4:(p+3)], group.theta[[3]][4:(p+3)],  group.theta[[4]][4:(p+3)])
  node_degree <- 1/rowSums(AA)
  node_degree[which(node_degree == Inf)] <- 10e6  # process the Inf
  W <- diag(node_degree) %*% AA  # standerized adjacency matrix
  ### \cal{Beta}_0
  Beta0 <- numeric(n.nodes)
  for (i in 1:k) { Beta0 <- Beta0 + Dk[[i]] %*% (beta0k[[i]] * rep(1, n.nodes) + Vattr %*% gammak[[i]]) }
  ### \cal{Beta}_1
  Beta1 <- numeric(n.nodes)
  for (i in 1:k) { Beta1 <- Beta1 + Dk[[i]] %*% (beta1k[[i]] * diag(rep(1, n.nodes))) }
  ### \cal{Beta}_2
  Beta2 <- numeric(n.nodes)
  for (i in 1:k) { Beta2 <- Beta2 + Dk[[i]] %*% (beta2k[[i]] * diag(rep(1, n.nodes))) }
  G <- Beta1 %*% W + Beta2
  true_sigma_2 <- unlist(group.sigma)  # unlist to a vector
  delta <- z %*% sqrt(true_sigma_2)
  ## generate white noise (length of T) for each node
  wn <- sapply(1:n.nodes, FUN = function(x) {set.seed(x + myseed); ts(rnorm(t))})
  
  ## initialize Y
  Y <- matrix(0, nrow = t, ncol = n.nodes)
  ## generate Y from t=2 to t=T
  for (time in 2:t) {
    Y[time, ] <- Beta0 + G %*% Y[(time-1), ] + delta * wn[time,]
  }
  
  ## calculate X_tilde
  X.tilde <- array(0, dim = c(t, n.nodes, (p + 3)))
  for (time in 2:t) {
    sub_Xtilde <- matrix(c(rep(1, n.nodes), W %*% t(Y)[,(time-1)], t(Y)[,(time-1)]), byrow = F, nrow = n.nodes)
    X.tilde[time,,] <- cbind(sub_Xtilde, Vattr)  # Vattr: [n.nodes * p]
  }
  ## X_tilde at t=1
  X.tilde[1,,] <- cbind(matrix(c(rep(1, n.nodes), rep(0, n.nodes), rep(0, n.nodes)), byrow = F, nrow = n.nodes), Vattr)
  
  res <- list(X_tilde = X.tilde, Y_tilde = Y, d_matrix = d.matrix, node_label = LL, W = W)
  return(res)
}

gendata_K6 <- function(n.nodes, p, t, group.sigma, group.theta, myseed, AA, LL, g)
{
  #===========================================================#
  ## Input: n.nodes = number of nodes in network ##
  ##        n.edges = number of edges in network ##
  ##        p = attributes(covariates) dimension ##
  ##        t = observed time length ##
  ##        group.sigma = a list containing sigma^2 from all groups ##
  ##        group.theta = a list containing theta from all groups ##
  ##        z = matrix of cluster assignment, n.nodes by n.cluster, 1 indicates node belongs to the cluster, 0 else.
  ##        myseed = random seed
  ## Output:
  ##        a list containing:
  ##        [["X_tilde"]], [["Y_tilde"]], [["d_matrix"]]
  #===========================================================#
  
  ## generate graph with nodes and edges
  set.seed(myseed)
  W <- AA/rowSums(AA)  # row-normalized adj matrix
  z <- model.matrix(~as.factor(LL) - 1)
  d.matrix <- distances(g, v = V(g), to = V(g))
  
  ## assume the nodes attributes follow MVN(0,1)
  sigma_attr <- diag(1, nrow = p)
  set.seed(myseed)
  Vattr <- rmvnorm(n.nodes, mean = rep(0, p), sigma = sigma_attr)
  
  ## generate data Y_tilde
  ### from t=2 to t=T, generate n.nodes at each time step. Follow the equation (2.3).
  Dk <- list(diag(z[, 1]), diag(z[, 2]), diag(z[, 3]), diag(z[, 4]), diag(z[, 5]), diag(z[, 6]))
  beta0k <- list(group.theta[[1]][1], group.theta[[2]][1], group.theta[[3]][1], 
                 group.theta[[4]][1], group.theta[[5]][1], group.theta[[6]][1])
  beta1k <- list(group.theta[[1]][2], group.theta[[2]][2], group.theta[[3]][2], 
                 group.theta[[4]][2], group.theta[[5]][2], group.theta[[6]][2])
  beta2k <- list(group.theta[[1]][3], group.theta[[2]][3], group.theta[[3]][3], 
                 group.theta[[4]][3], group.theta[[5]][3], group.theta[[6]][3])
  gammak <- list(group.theta[[1]][4:(p+3)], group.theta[[2]][4:(p+3)], group.theta[[3]][4:(p+3)],  
                 group.theta[[4]][4:(p+3)], group.theta[[5]][4:(p+3)], group.theta[[6]][4:(p+3)])
  node_degree <- 1/rowSums(AA)
  node_degree[which(node_degree == Inf)] <- 10e6  # process the Inf
  W <- diag(node_degree) %*% AA  # standerized adjacency matrix
  ### \cal{Beta}_0
  Beta0 <- numeric(n.nodes)
  for (i in 1:k) { Beta0 <- Beta0 + Dk[[i]] %*% (beta0k[[i]] * rep(1, n.nodes) + Vattr %*% gammak[[i]]) }
  ### \cal{Beta}_1
  Beta1 <- numeric(n.nodes)
  for (i in 1:k) { Beta1 <- Beta1 + Dk[[i]] %*% (beta1k[[i]] * diag(rep(1, n.nodes))) }
  ### \cal{Beta}_2
  Beta2 <- numeric(n.nodes)
  for (i in 1:k) { Beta2 <- Beta2 + Dk[[i]] %*% (beta2k[[i]] * diag(rep(1, n.nodes))) }
  G <- Beta1 %*% W + Beta2
  true_sigma_2 <- unlist(group.sigma)  # unlist to a vector
  delta <- z %*% sqrt(true_sigma_2)
  ## generate white noise (length of T) for each node
  wn <- sapply(1:n.nodes, FUN = function(x) {set.seed(x + myseed); ts(rnorm(t))})
  
  ## initialize Y
  Y <- matrix(0, nrow = t, ncol = n.nodes)
  ## generate Y from t=2 to t=T
  for (time in 2:t) {
    Y[time, ] <- Beta0 + G %*% Y[(time-1), ] + delta * wn[time,]
  }
  
  ## calculate X_tilde
  X.tilde <- array(0, dim = c(t, n.nodes, (p + 3)))
  for (time in 2:t) {
    sub_Xtilde <- matrix(c(rep(1, n.nodes), W %*% t(Y)[,(time-1)], t(Y)[,(time-1)]), byrow = F, nrow = n.nodes)
    X.tilde[time,,] <- cbind(sub_Xtilde, Vattr)  # Vattr: [n.nodes * p]
  }
  ## X_tilde at t=1
  X.tilde[1,,] <- cbind(matrix(c(rep(1, n.nodes), rep(0, n.nodes), rep(0, n.nodes)), byrow = F, nrow = n.nodes), Vattr)
  
  res <- list(X_tilde = X.tilde, Y_tilde = Y, d_matrix = d.matrix, node_label = LL, W = W)
  return(res)
}
