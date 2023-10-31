# Wrappers ----
# |||||||||||||

## fPLS_R ----
## |||||||||||

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
                             folds = 5,
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
              C_hat = W_hat_fPLS_R))
}