# Wrappers ----
# |||||||||||||

wrapper_colMeans <- function(){
  return(list(x_hat = colMeans(X_batch)))
}


wrapper_FRPDE <- function(){
  
  # Set model
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  # Define and init model
  model_FRPDE <- new(FRPDE_Laplacian_2D_GeoStatNodes, pde)
  model_FRPDE$set_verbose(FALSE)
  
  # Set observations
  model_FRPDE$set_observations(X_batch)
  
  # Tune
  lambda_opt <- model_FRPDE$tune(lambdas)
  
  # Set optimal lambda
  model_FRPDE$set_lambda_s(lambda_opt)
  
  # Solve
  model_FRPDE$solve()
  
  # Results
  x_hat <- as.numeric(model_FRPDE$fitted())
  
  return(list(lambda_opt = lambda_opt,
              x_hat = x_hat))
}
