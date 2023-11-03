# Wrappers ----
# |||||||||||||

# Multivariate implementation

wrapper_PLS <- function(){
  
  # Results
  Y_hat_PLS <- list()
  X_hat_PLS <- list()
  B_hat_PLS <- list()
  for(h in 1:nComp) {
    model_PLS <- PLS(X_batch, Y_batch, A = h, deflation_Y = FALSE)
    Y_hat_PLS[[h]] <- model_PLS$Y_hat
    X_hat_PLS[[h]] <- model_PLS$X_hat
    B_hat_PLS[[h]] <- model_PLS$Beta
  }
  
  X_mean_PLS <- model_PLS$X.mean
  Y_mean_PLS <- model_PLS$Y.mean
  W_hat_PLS <- model_PLS$W
  C_hat_PLS <- model_PLS$C
  
  return(list(X_mean = X_mean_PLS,
              Y_mean = Y_mean_PLS,
              X_hat = X_hat_PLS,
              Y_hat = Y_hat_PLS,
              B_hat = B_hat_PLS,
              W_hat = W_hat_PLS,
              C_hat = C_hat_PLS))
}

# - fPLS 
# - C++ implementation
# - Lambda selection strategy: GCV

wrapper_fPLS_cpp <- function(){
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model_fPLS <- new(FPLSR_Laplacian_2D_GeoStatNodes, pde)
  model_fPLS$set_verbose(FALSE)
  model_fPLS$set_H(nComp)
  model_fPLS$set_full_functional(FALSE)
  model_fPLS$set_lambdas(lambdas_in)
  model_fPLS$set_smoothing_initialization_l(TRUE, lambdas_in)
  model_fPLS$set_smoothing_regression_l(TRUE, lambdas_in)
  
  # Set data
  model_fPLS$set_data(Y_batch, X_batch)
  
  # Solve
  model_fPLS$solve()
  
  # Results
  X_mean_fPLS <- model_fPLS$X_mean()
  Y_mean_fPLS <- model_fPLS$Y_mean()
  W_hat_fPLS <- model_fPLS$W()
  C_hat_fPLS <- model_fPLS$C()
  Y_hat_fPLS <- list()
  X_hat_fPLS <- list()
  B_hat_fPLS <- list()
  model_fPLS$set_verbose(FALSE)
  for(h in 1:nComp) {
    Y_hat_fPLS[[h]] <- model_fPLS$fitted(h)
    X_hat_fPLS[[h]] <- model_fPLS$reconstructed_field(h)
    B_hat_fPLS[[h]] <- model_fPLS$B(h)
  }
  
  return(list(X_mean = X_mean_fPLS,
              Y_mean = Y_mean_fPLS,
              X_hat = X_hat_fPLS,
              Y_hat = Y_hat_fPLS,
              B_hat = B_hat_fPLS,
              W_hat = W_hat_fPLS,
              C_hat = C_hat_fPLS))
  
}

# - fPLS 
# - R implementation
# - Lambda selection strategy: unique

wrapper_fPLS_R_unique <- function(){
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  registerDoParallel(cl)
  
  cv_r1fpls <- cv_unique_par(X = X_batch,
                             Y = Y_batch,
                             penalty_vec = lambdas_in,
                             ncomp = nComp,
                             folds = k_folds,
                             basisobj = FEM_basis,
                             method = "r1fpls_fem",
                             verbose = FALSE,
                             stripped = FALSE )
  
  stopCluster(cl)
  
  model_fPLS_R <- cv_r1fpls$final_model
  
  # Results
  X_mean_fPLS_R <- model_fPLS_R$X_mean
  Y_mean_fPLS_R <- model_fPLS_R$Y_mean
  W_hat_fPLS_R <- as.matrix(model_fPLS_R$W)
  C_hat_fPLS_R <- as.matrix(model_fPLS_R$C)
  TT <- as.matrix(model_fPLS_R$TT)
  Y_hat_fPLS_R <- list()
  X_hat_fPLS_R <- list()
  B_hat_fPLS_R <- list()
  for(h in 1:nComp) {
    Y_hat_fPLS_R[[h]] <- model_fPLS_R$fitted.values[,,h]
    X_hat_fPLS_R[[h]] <- TT[,1:h] %*% t(C_hat_fPLS_R[,1:h]) + rep(1, N) %*% t(X_mean_fPLS_R)
    B_hat_fPLS_R[[h]] <- model_fPLS_R$coefficient_function[,,h]
  }
  
  return(list(X_mean = X_mean_fPLS_R,
              Y_mean = Y_mean_fPLS_R,
              X_hat = X_hat_fPLS_R,
              Y_hat = Y_hat_fPLS_R,
              B_hat = B_hat_fPLS_R,
              W_hat = W_hat_fPLS_R,
              C_hat = C_hat_fPLS_R))
}

# - fPLS 
# - R implementation
# - Lambda selection strategy: sequential

wrapper_fPLS_R_seq <- function(){
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  registerDoParallel(cl)
  
  cv_r1fpls <- cv_seq_par(X = X_batch,
                          Y = Y_batch,
                          penalty_vec = lambdas_in,
                          ncomp = nComp,
                          folds = k_folds,
                          basisobj = FEM_basis,
                          method = "r1fpls_fem",
                          verbose = FALSE,
                          stripped = FALSE )
  
  stopCluster(cl)
  
  model_fPLS_R <- cv_r1fpls$final_model
  
  # Results
  X_mean_fPLS_R <- model_fPLS_R$X_mean
  Y_mean_fPLS_R <- model_fPLS_R$Y_mean
  W_hat_fPLS_R <- as.matrix(model_fPLS_R$W)
  C_hat_fPLS_R <- as.matrix(model_fPLS_R$C)
  TT <- as.matrix(model_fPLS_R$TT)
  Y_hat_fPLS_R <- list()
  X_hat_fPLS_R <- list()
  B_hat_fPLS_R <- list()
  for(h in 1:nComp) {
    Y_hat_fPLS_R[[h]] <- model_fPLS_R$fitted.values[,,h]
    X_hat_fPLS_R[[h]] <- TT[,1:h] %*% t(C_hat_fPLS_R[,1:h]) + rep(1, N) %*% t(X_mean_fPLS_R)
    B_hat_fPLS_R[[h]] <- model_fPLS_R$coefficient_function[,,h]
  }
  
  return(list(X_mean = X_mean_fPLS_R,
              Y_mean = Y_mean_fPLS_R,
              X_hat = X_hat_fPLS_R,
              Y_hat = Y_hat_fPLS_R,
              B_hat = B_hat_fPLS_R,
              W_hat = W_hat_fPLS_R,
              C_hat = C_hat_fPLS_R))
}

# Functional PLS based on the (2D) TPS representation of the data

wrapper_B_fPLS <- function(){
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  registerDoParallel(cl)
  
  cv_pTPS <- cv_unique_par(X = X_batch,
                           Y = Y_batch,
                           center = TRUE,
                           nodes = nodes,
                           nbasis = n_basis_tps,
                           penalty_vec = 0,
                           basisobj = NA,
                           ncomp = nComp,
                           folds = k_folds,
                           method = "fpls_tps",
                           verbose = FALSE,
                           stripped = FALSE,
                           R0 = R0_tps )
  
  stopCluster(cl)
  
  model_B_fPLS <- cv_pTPS$final_model
  
  # Results
  X_mean_B_fPLS  <- model_B_fPLS$X_mean
  Y_mean_B_fPLS  <- model_B_fPLS$Y_mean
  W_hat_B_fPLS  <- Psi_tps %*% as.matrix(model_B_fPLS$W)
  C_hat_B_fPLS  <- Psi_tps %*% as.matrix(model_B_fPLS$C)
  TT <- as.matrix(model_B_fPLS$TT)
  Y_hat_B_fPLS  <- list()
  X_hat_B_fPLS  <- list()
  B_hat_B_fPLS  <- list()
  for(h in 1:nComp) {
    Y_hat_B_fPLS[[h]] <- model_B_fPLS$fitted.values[,,h]
    X_hat_B_fPLS[[h]] <- TT[,1:h] %*% t(C_hat_B_fPLS[,1:h]) + rep(1, N) %*% t(X_mean_B_fPLS)
    B_hat_B_fPLS[[h]] <- model_B_fPLS$coefficient_function[,,h]
  }
  
  return(list(X_mean = X_mean_B_fPLS,
              Y_mean = Y_mean_B_fPLS,
              X_hat = X_hat_B_fPLS,
              Y_hat = Y_hat_B_fPLS,
              B_hat = B_hat_B_fPLS,
              W_hat = W_hat_B_fPLS,
              C_hat = C_hat_B_fPLS))
  
}

# Penalized functional PLS based on the (2D) TPS representation of the data

wrapper_PB_fPLS <- function(){
  
  nodes_CL = detectCores()   # Detect number of cores to use
  cl = makeCluster(nodes_CL) # Specify number of threads here
  registerDoParallel(cl)
  
  cv_pTPS <- cv_unique_par(X = X_batch,
                           Y = Y_batch,
                           center = TRUE,
                           nodes = nodes,
                           nbasis = n_basis_tps,
                           penalty_vec = lambdas_in,
                           basisobj = NA,
                           ncomp = nComp,
                           folds = k_folds,
                           method = "fpls_tps",
                           verbose = FALSE,
                           stripped = FALSE,
                           R0 = R0_tps )
  
  stopCluster(cl)
  
  model_PB_fPLS <- cv_pTPS$final_model
  
  # Results
  X_mean_PB_fPLS  <- model_PB_fPLS$X_mean
  Y_mean_PB_fPLS  <- model_PB_fPLS$Y_mean
  W_hat_PB_fPLS  <- Psi_tps %*% as.matrix(model_PB_fPLS$W)
  C_hat_PB_fPLS  <- Psi_tps %*% as.matrix(model_PB_fPLS$C)
  TT <- as.matrix(model_PB_fPLS$TT)
  Y_hat_PB_fPLS  <- list()
  X_hat_PB_fPLS  <- list()
  B_hat_PB_fPLS  <- list()
  for(h in 1:nComp) {
    Y_hat_PB_fPLS[[h]] <- model_PB_fPLS$fitted.values[,,h]
    X_hat_PB_fPLS[[h]] <- TT[,1:h] %*% t(C_hat_PB_fPLS[,1:h]) + rep(1, N) %*% t(X_mean_PB_fPLS)
    B_hat_PB_fPLS[[h]] <- model_PB_fPLS$coefficient_function[,,h]
  }
  
  return(list(X_mean = X_mean_PB_fPLS,
              Y_mean = Y_mean_PB_fPLS,
              X_hat = X_hat_PB_fPLS,
              Y_hat = Y_hat_PB_fPLS,
              B_hat = B_hat_PB_fPLS,
              W_hat = W_hat_PB_fPLS,
              C_hat = C_hat_PB_fPLS))

  
}

# - fPCA Regression
# - Lambda selection strategy: GCV

wrapper_fPCA <- function(){
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model_fPCA <- new(FPCA_Laplacian_2D_GeoStatNodes, pde)
  model_fPCA$set_lambdas(lambdas_in)
  model_fPCA$set_npc(nComp)
  
  # Set data
  model_fPCA$set_observations(X_batch - rep(1, N) %*% t(X_mean_fPLS))
  
  # Solve
  model_fPCA$solve()
  
  # Results
  F_hat_fPCA <- model_fPCA$loadings()
  S_hat_fPCA <- model_fPCA$scores()
  Y_mean_fPCA <- Y_mean_fPLS
  X_mean_fPCA <- X_mean_fPLS
  Y_hat_fPCA <- list()
  X_hat_fPCA <- list()
  B_hat_fPCA <- list()
  for(h in 1:nComp) {
    fit <- lm(Y_batch - rep(1, N) %*% t(Y_mean_fPLS) ~ 0 + S_hat_fPCA[,1:h])
    beta_hat <- matrix(fit$coefficients, ncol = 1)
    Y_hat_fPCA[[h]] <- fit$fitted.values + rep(1, N) %*% t(Y_mean_fPLS)
    X_hat_fPCA[[h]] <- S_hat_fPCA[,1:h] %*% t(F_hat_fPCA[,1:h]) + rep(1, N) %*% t(X_mean_fPLS)
    B_hat_fPCA[[h]] <- as.numeric(F_hat_fPCA[,1:h] %*% solve(t(F_hat_fPCA[,1:h]) %*% F_hat_fPCA[,1:h]) %*% beta_hat)
  }
  
  return(list(X_mean = X_mean_fPCA,
              Y_mean = Y_mean_fPCA,
              X_hat = X_hat_fPCA,
              Y_hat = Y_hat_fPCA,
              B_hat = B_hat_fPCA,
              F_hat = F_hat_fPCA))
  
}