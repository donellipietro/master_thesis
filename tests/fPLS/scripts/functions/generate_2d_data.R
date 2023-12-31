# Generate 2d data ----
# |||||||||||||||||||||

# New implementation

comp_sin <- function(nodes, k1, k2, phi1, phi2){
  sin(k1*pi*(nodes[,1]+phi1))*sin(k2*pi*(nodes[,2]+phi2))
}

generate_2d_data <- function(num_grid_x_axes = 20,
                             num_grid_y_axes = 20,
                             N = 60,
                             H = 3,
                             L = 1,
                             B.index = 3,
                             NSR.X = 0.2^2,
                             NSR.Y = 0.2^2,
                             ORTH = FALSE) {
  
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
  if(ORTH){
    phi1 <- c(0, 0, 0, 0)
    phi2 <- c(0, 0, 0, 0)
  } else {
    phi1 <- c(0, 0.6, -0.3, 0.2)
    phi2 <- c(0.35, 0, 0.2, 0.45)
  }
  load_ex <- list()
  W <- matrix(0, nrow = nrow(nodes), ncol = H)
  for(h in 1:H){
    W[,h] <- comp_sin(nodes, k1[h], k2[h], phi1[h], phi2[h])
    L2norm <- sqrt(as.numeric(t(W[,h]) %*% R0 %*% W[,h]))
    W[,h] <- W[,h] / L2norm
  }
  W_star <- 0.5*comp_sin(nodes, k1[8], k2[8], 0, 0)
  
  # Latent variables sd.
  sd_s <- c(1, 0.8, 0.3, 0.2)
  sigma_s <- 2 * sd_s / (apply(W, 2, max)-apply(W, 2, min))
  
  # Loading values
  if(B.index == 1)
    D <- matrix(c(1, 0, 0, 0,
                  0, 1, 0, 0),
                ncol = 4, byrow = TRUE)
  if(B.index == 2)
    D <- matrix(c(0, 1, 0, 0,
                  0, 0, 1, 0),
                ncol = 4, byrow = TRUE)
  if(B.index == 3)
    D <- matrix(c(0, 0, 1, 0,
                  0, 0, 0, 1),
                ncol = 4, byrow = TRUE)
  if(B.index == 4)
    D <- matrix(c(0, 0, 0, 1,
                  1, 0, 0, 0),
                ncol = 4, byrow = TRUE)
  if(B.index == 5)
    D <- matrix(c(1, -3, 0, 8,
                  0, 0, 8, 0),
                ncol = 4, byrow = TRUE)
  if(B.index == 6)
    D <- matrix(c(0, 0, 8, 0,
                  1, -3, 0, 5),
                ncol = 4, byrow = TRUE)
  
  D <- as.matrix(D[1:L,], ncol = 4, byrow = TRUE) 
  if(L == 1) D <- t(D) # t because R is stupid
  
  # Noises sd.
  sigma_noise_x <- 0.5 * sqrt(NSR.X) * sum(sigma_s * 0.5 * (apply(W, 2, max)-apply(W, 2, min)))
  sigma_noise_y <- c()
  for(l in 1:L){
    sigma_noise_y[l] <- sqrt(NSR.Y * sum(sigma_s^2 * D[l,]^2))
    # cat("\n Variance explained by each component: ")
    # cat(sigma_s^2 * D[l,]^2)
  }
  
  # Sampling
  Sampled <- mvrnorm(N, mu = rep(0, H + L + S), diag(c(sigma_s^2, sigma_noise_y^2, rep(sigma_noise_x^2, S))) )
  SS <- Sampled[,1:H]
  FF <- matrix(Sampled[,(H+1):(H+L)], nrow = N, ncol = L)
  EE <- Sampled[,(H+L+1):(H+L+S)]
  
  # B
  B <- W %*% solve(t(W) %*% W) %*% t(D)
  
  # Generating X:
  X_clean = SS %*% t(W) + rep(1,N) %*% t(W_star)
  
  # Clean Y
  Y_clean <- SS %*% t(D) + rep(1,N) %*% t(c(2, 5)[1:L])
  
  # Centering:
  X_clean_center <- scale(X_clean, scale = F)
  X_mean <- attr(X_clean_center, "scaled:center")
  Y_clean_center <- scale(Y_clean, scale = F)
  Y_mean <- attr(Y_clean_center, "scaled:center")
  
  # Noisy observations:
  X = X_clean + EE
  X_center <- X - rep(1,N) %*% t(X_mean)
  Y <- Y_clean + FF
  Y_center <- Y - rep(1,N) %*% t(Y_mean)
  
  return(list(X = X,
              Y = Y,
              X_center = X_center,
              Y_center = Y_center,
              X_clean = X_clean,
              Y_clean = Y_clean,
              X_mean = X_mean,
              Y_mean = Y_mean,
              basisobj = FEM_basis,
              mesh = mesh,
              B_true = B,
              F_true = W,
              G_true = D))
  
}