# Data generation ----
# ||||||||||||||||||||

fun <- function(points){
  f <- function(x, y, z = 1){
    coe <- function(x,y){
      1/2*sin(5*pi*x)*exp(-x^2)+1
    }
    sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
  }
  z = f(points[,1], points[,2])
}

generate_2d_data <- function(num_grid_x_axes = 20,
                             num_grid_y_axes = 20,
                             N = 60,
                             NSR.X = 0.2^2) {
  
  # Nodes and 2D mesh:
  x <- seq(0, 1, length.out = num_grid_x_axes)
  y <- seq(0, 1, length.out = num_grid_y_axes)
  nodes <- expand.grid(x = x, y = y)
  mesh <- fdaPDE::create.mesh.2D(nodes = nodes)
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis)
  S = nrow(nodes)
  
  # Generating X:
  x_clean <- fun(nodes)
  X_clean <- matrix(x_clean, nrow = N, ncol = S, byrow = TRUE)
  
  # Noise
  # NSR.X = sigma_noise^2/Var(X)
  # => sigma_noise^2 = NSR.X*Var(X)
  sigma_noise <- sqrt(NSR.X*var(x_clean))
  # cat(sigma_noise)
  
  # Sampling
  EE <- mvrnorm(N, mu = rep(0, S), Sigma = diag(rep(sigma_noise^2, S)) )
  
  # Noisy observations:
  X = X_clean + EE
  
  return(list(X = X,
              x_clean = x_clean,
              basisobj = FEM_basis,
              mesh = mesh))
  
}