# fdaPDE ----
# |||||||||||

## FEM utilities ----
## ||||||||||||||||||

### PDE ----
### ||||||||

mesh_data_3D <- function(mesh) {
  
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$faces,
    "elements" = mesh$tetrahedrons,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  return(mesh_data)
}

pde.init <- function(mesh) {
  
  # Set model
  pde <- new(Laplacian_3D_Order1, mesh_data_3D(mesh))
  
  # Set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  return(pde)
}


### Visualization ----
### ||||||||||||||||||

field.eval <- function(grid, f_nodes, mesh) {
  
  FEMbasis = create.FEM.basis(mesh)
  
  FEMfunction = FEM(f_nodes, FEMbasis)
  f_grid <- eval.FEM(FEMfunction, grid)
  
  return(f_grid)
}


### Conversion ----
### |||||||||||||||

print_grid <- function(file, dim, p, nnodes, t, nelems){
  eltype = 10
  
  ## VTK-Points (mesh nodes)
  cat('<Points>\n', file = file, append = TRUE)
  cat('<DataArray type="Float64" Name="Array" NumberOfComponents="3" format="ascii">\n', file = file, append = TRUE)
  write.table(t(p),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('</Points>\n', file = file, append = TRUE)
  
  ## VTK-Cells (mesh elements)
  cat('<Cells>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="connectivity" format="ascii">\n', file = file, append = TRUE)
  write.table(format(t(t), scientific = FALSE),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="offsets" format="ascii">\n', file = file, append = TRUE)
  write.table(format(matrix(seq(dim+1,(dim+1)*nelems,by=(dim+1)),nrow = 1, byrow = TRUE),  scientific = FALSE),  
              file = file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('<DataArray type="Int32" Name="types" format="ascii">\n', file = file, append = TRUE)
  write.table(format(matrix(rep(eltype,nelems),nrow = 1, byrow = TRUE), scientific = FALSE),  file = file, append = TRUE, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('</DataArray>\n', file = file, append = TRUE)
  cat('</Cells>\n', file = file, append = TRUE)
}

print_data_points <- function(file, nodedata, nnodes){
  
  ## # of data to print in 
  ## <PointData> field
  nvdata = ncol(nodedata)  
  
  cat('<PointData>\n', file = file, append = TRUE)
  for (ind in 1:nvdata) {
    cat(paste0('<DataArray type="Float64" Name="Data',ind,'" NumberOfComponents="1" format="ascii">\n'), file = file, append = TRUE)
    write.table(matrix(nodedata[,ind],nrow = 1),  file = file, append = TRUE, row.names = FALSE, col.names = FALSE)
    cat('</DataArray>\n', file = file, append = TRUE)
  }
  cat('</PointData>\n', file = file, append = TRUE)
  
}

write.vtu <- function(FEMobj, file = ""){
  
  if(class(FEMobj) != "FEM"){
    stop("FEMobj is not of class FEM.")
  }
  
  if(!is.character(file)){
    stop("file is not of type character.")
  }
  
  ## Header
  cat('<?xml version="1.0"?>\n', file = file)
  cat('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n', file = file, append = TRUE)
  cat('<UnstructuredGrid>\n', file = file, append = TRUE)
  
  dim = 3
  p = t(FEMobj$FEMbasis$mesh$nodes)
  t = t(FEMobj$FEMbasis$mesh$tetrahedrons)
  t = t - 1
  
  nnodes = nrow(FEMobj$FEMbasis$mesh$nodes)
  nelems = nrow(FEMobj$FEMbasis$mesh$tetrahedrons)
  
  ## Header for <Piece>
  cat(paste0('<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',nelems,'">\n'), file = file, append = TRUE)
  
  ## Print grid
  print_grid(file, dim, p, nnodes, t, nelems);
  
  ## Print Data
  print_data_points(file, FEMobj$coeff, nnodes)
  
  ## Footer for <Piece>
  cat('</Piece>\n', file = file, append = TRUE)
  
  ## Footer
  cat('</UnstructuredGrid>\n', file = file, append = TRUE)
  cat('</VTKFile>', file = file, append = TRUE)
}


### Errors ----
### |||||||||||

RMSE <- function(x1, x2){
  sqrt(sum((x1 - x2)^2)/length(x1))
}

compute_residuals <- function(Y, X, Y_hat, X_hat) {
  
  residuals_X <- c()
  for(h in 1:length(X_hat)){
    residuals_X[h] <- RMSE(X, X_hat[[h]])
  }
  residuals_X <- c(1, residuals_X/RMSE(X, 0))
  
  if(!is.null(Y_hat)){
    residuals_Y <- c()
    for(h in 1:length(Y_hat)){
      residuals_Y[h] <- RMSE(Y, Y_hat[[h]])
    }
    residuals_Y <- c(1, residuals_Y/RMSE(Y, 0))
  } else {
    residuals_Y <- NULL
  }
  
  return(list(nComp = 0:length(X_hat),
              X = residuals_X,
              Y = residuals_Y))
}


## Wrappers ----
## |||||||||||||

# Column-wise mean

ColMean <- function(X, mesh = NULL, lambdas = NULL, verbose = FALSE){
  
  tic()
  
  fitted <- fitted_nodes <- colMeans(X)
  
  elapsed <- toc(quiet = TRUE)
  
  return(list(fitted = fitted,
              fitted_nodes = fitted_nodes, 
              execution_time = as.numeric(elapsed$toc-elapsed$tic)))
  
}

# Functional Regression with PDE regularization model

FRPDE <- function(X, mesh, lambdas, verbose = FALSE){
  
  tic()
  
  # Define and init model
  model_FRPDE <- new(FRPDE_Laplacian_3D_GeoStatNodes, pde.init(mesh))
  model_FRPDE$set_verbose(verbose)
  
  # Set observations
  model_FRPDE$set_observations(as.matrix(X))
  
  # Tune
  lambda_opt <- model_FRPDE$tune(lambdas)
  
  # Solve
  model_FRPDE$set_lambda_s(lambda_opt)
  model_FRPDE$solve()
  
  # Results
  fitted <- t(model_FRPDE$fitted())
  fitted_nodes <- t(model_FRPDE$f())
  
  elapsed <- toc(quiet = TRUE)
  
  return(list(lambda_opt = lambda_opt,
              fitted = fitted,
              fitted_nodes = fitted_nodes, 
              execution_time = as.numeric(elapsed$toc-elapsed$tic)))
}

# Multivariate Partial Least Squares Regression

MV_PLS <- function(Y, X, nComp, mesh = NULL, lambdas = NULL, verbose = FALSE){
  
  # Results
  Y_hat <- list()
  X_hat <- list()
  B_hat <- list()
  for(h in 1:nComp) {
    model_PLS <- PLS(X, Y, A = h, deflation_Y = FALSE)
    Y_hat[[h]] <- model_PLS$Y_hat
    X_hat[[h]] <- model_PLS$X_hat
    B_hat[[h]] <- model_PLS$Beta
  }
  
  X_mean <- model_PLS$X.mean
  Y_mean <- model_PLS$Y.mean
  W_hat <- model_PLS$W
  C_hat <- model_PLS$C
  T_hat <- model_PLS$TT
  
  return(list(X_mean = X_mean,
              Y_mean = Y_mean,
              X_hat = X_hat,
              Y_hat = Y_hat,
              B_hat = B_hat,
              W_hat = W_hat,
              C_hat = C_hat,
              T_hat = T_hat))
}


# Functional Partial Least Sqaures Regression

fPLS <- function(Y, X, nComp, mesh, lambdas, verbose = FALSE){
  
  tic()
  
  # Define and init model
  model_fPLS <- new(FPLSR_Laplacian_3D_GeoStatNodes, pde.init(mesh))
  model_fPLS$set_verbose(verbose)
  
  model_fPLS$set_H(nComp)
  model_fPLS$set_full_functional(FALSE)
  model_fPLS$set_lambdas(lambdas)
  model_fPLS$set_smoothing_initialization_l(TRUE, lambdas)
  model_fPLS$set_smoothing_regression_l(TRUE, lambdas)
  
  # Set data
  model_fPLS$set_data(as.matrix(Y), X)

  # Solve
  model_fPLS$solve()
  
  elapsed <- toc()
  
  # Results
  X_mean <- model_fPLS$X_mean()
  Y_mean <- model_fPLS$Y_mean()
  W_hat <- model_fPLS$W()
  C_hat <- model_fPLS$C()
  T_hat <- model_fPLS$T()
  Y_hat <- list()
  X_hat <- list()
  B_hat <- list()
  model_fPLS$set_verbose(FALSE)
  for(h in 1:nComp) {
    Y_hat[[h]] <- model_fPLS$fitted(h)
    X_hat[[h]] <- model_fPLS$reconstructed_field(h)
    B_hat[[h]] <- model_fPLS$B(h)
  }
  
  optimal_lambda <- list()
  optimal_lambda[["initialization"]] <- model_fPLS$get_lambda_initialization()
  optimal_lambda[["directions"]] <- model_fPLS$get_lambda_directions()
  optimal_lambda[["regression"]] <- model_fPLS$get_lambda_regression()

  return(list(X_mean = X_mean,
              Y_mean = Y_mean,
              X_hat = X_hat,
              Y_hat = Y_hat,
              B_hat = B_hat,
              W_hat = W_hat,
              C_hat = C_hat,
              T_hat = T_hat,
              optimal_lambda = optimal_lambda,
              execution_time = as.numeric(elapsed)))
  
}

# Functional Principal Components Analysis

fPCA <- function(mesh, X, lambdas, nComp, verbose = FALSE){
  
  tic()
  
  # Define and init model
  model_fPCA <- new(FPCA_Laplacian_3D_GeoStatNodes_GCV, pde.init(mesh))
  model_fPCA$set_npc(nComp)
  model_fPCA$set_lambdas(lambdas)
  
  # Set observations
  model_fPCA$set_observations(as.matrix(X))
  
  # Solve
  model_fPCA$solve()
  
  elapsed <- toc()
  
  # Optimal lambda
  lambdas_opt <- model_fPCA$lambda_opt()
  
  # Results
  scores <- model_fPCA$scores()
  loadings <-  model_fPCA$loadings()
  X_hat <- list()
  for(h in 1:nComp) {
    X_hat[[h]] <- scores[,1:h] %*% t(loadings[,1:h])
  }
  
  return(list(T_hat = scores,
              C_hat = loadings,
              X_hat = X_hat,
              lambdas_opt = lambdas_opt,
              execution_time = as.numeric(elapsed$toc-elapsed$tic)))
}
