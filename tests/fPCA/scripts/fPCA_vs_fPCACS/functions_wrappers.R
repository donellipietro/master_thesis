# Wrappers ----
# |||||||||||||

# The following functions wrap the statistical methods used during the analysis
# in order to proide an interface the is the same for all of them.

wrapper_PCA <- function() {
  
  tic()
  
  # PCA
  pca <- prcomp(X_batch, center = FALSE, rank. = nComp)
  
  # Results
  F_hat <- pca$rotation
  S_hat <- pca$x
  for(h in 1:nComp){
    f_norm <- as.numeric(sqrt(t(F_hat[,h]) %*% R0 %*% F_hat[,h]))
    F_hat[,h] <- F_hat[,h]/f_norm
    S_hat[,h] <- S_hat[,h]*f_norm
  }
  
  elapsed <- toc(quiet = TRUE)
  
  adjusted <- adjust_results(F_hat, S_hat)
  
  return(list(F_hat = adjusted$F_hat,
              S_hat = adjusted$S_hat,
              time = as.numeric(elapsed$toc - elapsed$tic)))
}

wrapper_fPCA <- function() {
  
  tic()
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model <- new(FPCA_Laplacian_2D_GeoStatNodes, pde)
  model$set_lambdas(lambdas)
  model$set_npc(nComp)
  model$init_regularization()
  
  # Set observations
  model$set_observations(X_batch)
  
  # Solve
  model$solve()
  
  # Results
  F_hat <- model$loadings()
  S_hat <- model$scores()
  
  elapsed <- toc(quiet = TRUE)
  
  adjusted <- adjust_results(F_hat, S_hat)
  
  return(list(F_hat = adjusted$F_hat,
              S_hat = adjusted$S_hat,
              time = as.numeric(elapsed$toc - elapsed$tic)))
}

wrapper_fPCACS <- function(ML= FALSE, IT = FALSE) {
  
  tic()
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model <- new(FPCA_CS_Laplacian_2D_GeoStatNodes, pde)
  model$set_lambda_s(lambdas_in)
  model$set_npc(nComp)
  model$set_mass_lumping(ML)
  model$set_iterative(IT)
  model$init_regularization()
  
  # Set observations
  model$set_observations(X_batch)
  
  # Solve
  model$solve()
  
  # Results
  F_hat <- model$loadings()
  S_hat <- model$scores()
  
  elapsed <- toc(quiet = TRUE)
  
  adjusted <- adjust_results(F_hat, S_hat)
  
  return(list(F_hat = adjusted$F_hat,
              S_hat = adjusted$S_hat,
              time = as.numeric(elapsed$toc - elapsed$tic)))
}

wrapper_SpatialPCA <- function(){

  # "Fake" Model initialization
  counts <- X_batch
  colnames(counts) <- paste("P", 1:S, sep = "")
  rownames(counts) <- paste("F", 1:N, sep = "")
  locations <- nodes
  rownames(locations) <- paste("P", 1:S, sep = "")
  model_SpatialPCA <- CreateSpatialPCAObject(counts = counts,
                                             location = locations,
                                             project='SpatialPCA',
                                             gene.type="none",
                                             sparkversion="sparkx",
                                             numCores_spark=5,
                                             gene.number=NULL,
                                             customGenelist=NULL,
                                             min.loctions = 20,
                                             min.features=20)
  # Model initialization
  model_SpatialPCA@params <- list()
  model_SpatialPCA@location <- locations
  model_SpatialPCA@normalized_expr <- counts
  
  tic()

  # Solve
  sink("log.txt")
  model_SpatialPCA <- SpatialPCA_buildKernel(model_SpatialPCA,
                                             kerneltype='gaussian',
                                             bandwidthtype="SJ",
                                             bandwidth.set.by.user=NULL)
  model_SpatialPCA <- SpatialPCA_EstimateLoading(model_SpatialPCA,
                                                 fast=FALSE,
                                                 SpatialPCnum=nComp)
  model_SpatialPCA <- SpatialPCA_SpatialPCs(model_SpatialPCA,
                                            fast=FALSE)
  sink()
  
  # Results
  F_hat <- t(model_SpatialPCA@SpatialPCs)
  S_hat <- model_SpatialPCA@W
  
  # Normalization
  for(h in 1:nComp){
    f_norm <- as.numeric(sqrt(t(F_hat[,h]) %*% R0 %*% F_hat[,h]))
    F_hat[,h] <- F_hat[,h]/f_norm
    S_hat[,h] <- S_hat[,h]*f_norm
  }

  elapsed <- toc(quiet = TRUE)
  
  adjusted <- adjust_results(F_hat, S_hat)
  
  return(list(F_hat = adjusted$F_hat,
              S_hat = adjusted$S_hat,
              time = as.numeric(elapsed$toc - elapsed$tic)))
  
}