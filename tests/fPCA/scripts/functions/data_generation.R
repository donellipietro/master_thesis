# Data generation ----
# ||||||||||||||||||||

eigen_functions_laplacian <- function(nodes, k1, k2){
  sin(k1*pi*nodes[,1])*sin(k2*pi*nodes[,2])
}

generate_2d_data <- function(num_grid_x_axes = 20,
                             num_grid_y_axes = 20,
                             N = 60,
                             H = 3,
                             NSR.X = 0.2^2) {
  
  # Nodes and 2D mesh:
  x <- seq(0, 1, length.out = num_grid_x_axes)
  y <- seq(0, 1, length.out = num_grid_y_axes)
  nodes <- expand.grid(x = x, y = y)
  mesh <- fdaPDE::create.mesh.2D(nodes = nodes)
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis)
  S = nrow(nodes)
  
  # Loading functions
  k1 <- c(1, 1, 2, 2, 1, 3, 2, 3, 1, 4, 2, 4, 3, 4, 4)
  k2 <- c(1, 2, 1, 2, 3, 2, 3, 3, 4, 1, 4, 2, 4, 3, 4)
  load_ex <- list()
  FF <- matrix(0, nrow = nrow(nodes), ncol = H)
  for(h in 1:H){
    FF[,h] <- eigen_functions_laplacian(nodes, k1[h], k2[h])
    L2norm <- sqrt(as.numeric(t(FF[,h]) %*% R0 %*% FF[,h]))
    FF[,h] <- FF[,h] / L2norm
  }
  
  # Latent variables sd.
  sd_s <- c(1, 0.8, 0.5, 0.3, 0.1)[1:H]
  sigma_s <- 2 * sd_s / (apply(FF, 2, max)-apply(FF, 2, min))
  
  # Noises sd.
  sigma_noise_x <- 0.5 * sqrt(NSR.X) * sum(sigma_s * 0.5 * (apply(FF, 2, max)-apply(FF, 2, min)))
  
  # Sampling
  Sampled <- mvrnorm(N, mu = rep(0, H + S), diag(c(sigma_s^2, rep(sigma_noise_x^2, S))) )
  SS <- as.matrix(Sampled[,1:H], ncol = H)
  EE <- Sampled[,(H+1):(H+S)]
  
  # Generating X:
  X_clean = SS %*% t(FF) # + rep(1,N) %*% t(W_star)
  
  # Centering:
  X_clean_center <- scale(X_clean, scale = F)
  X_mean <- attr(X_clean_center, "scaled:center")
  
  # Noisy observations:
  X = X_clean + EE
  X_center <- X - rep(1,N) %*% t(X_mean)
  
  return(list(X = X,
              X_center = X_center,
              X_clean = X_clean,
              X_clean_center = X_clean_center,
              X_mean = X_mean,
              basisobj = FEM_basis,
              mesh = mesh,
              F_true = FF,
              S_true = SS))
  
}