# Generate 2d data ----
# |||||||||||||||||||||

# penR1FPLS implementation
# Source: https://github.com/hhroig/penR1FPLS/blob/master/R/generate_2d_data.R
# Modifications: sd.X is now a parameter

generate_2d_data <- function(x,
                             y,
                             num_samples = 100,
                             beta_num = 3,
                             sd.X = 0.2,
                             Rsq = 0.95) {
  
  # Nodes and 2D mesh:
  nodes <- expand.grid(x = x, y = y)
  
  mesh <- fdaPDE::create.mesh.2D(nodes = nodes)
  
  # FEM basis:
  FEM_basis  <-  fdaPDE::create.FEM.basis(mesh)
  
  # Mass matrix:
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = FEM_basis)
  
  # Generate X:
  X = NULL
  
  for(ii in 1:num_samples){
    a1 = stats::rnorm(1, mean = 1, sd = 0.2)
    a2 = stats::rnorm(1, mean = 1, sd = 0.2)
    
    func_evaluation = numeric(nrow(mesh$nodes))
    
    for (i in 1:nrow(mesh$nodes)){
      
      func_evaluation[i] = a1* cos(2*pi*mesh$nodes[i,1]) +
        a2* cos(2*pi*mesh$nodes[i,2]) + 1
      
    }
    data = func_evaluation + stats::rnorm(nrow(mesh$nodes), mean = 0, sd = sd.X)
    X = rbind(X, data)
  }
  
  
  # Generate beta(x, y):
  if (beta_num == 1) {
    # centered at (0.5, 0.5):
    r  <-  0.4 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] - 0.5 )^2)/( 2*r^2 ))
    
  }else if (beta_num == 2) {
    # top right corner
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.2^2 ))
  }else if (beta_num == 3) {
    
    # bottom left + top right corner:
    z <- 5*exp(-((nodes[, 1] - 0.75)^2 + (nodes[, 2] - 0.75)^2)/( 2*0.25^2 )) +
      5*exp(-((nodes[, 1] - 0.1)^2 + (nodes[, 2] - 0.1)^2)/( 2*0.25^2 ))
    
  }else if (beta_num == 4) {
    
    # semi-circumference:
    r  <-  0.2 # set the r parameter
    z  <-  5*exp(-((nodes[, 1] - 0.5 )^2 + (nodes[, 2] )^2)/( 2*r^2 ))
    
  }else if (beta_num == 5) {
    # monkey saddle
    z = ((nodes[, 1]*4 - 2)^3 - 3*(nodes[, 1]*4 - 2)*((nodes[, 2]*4 - 2)^2))
    
  }else if (beta_num == 6) {
    # Test 1 - fdaPDE
    
    f = function(x, y, z = 1)
    {
      coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
      sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
    }
    
    # Exact solution (pointwise at nodes)
    z = f(nodes[,1], nodes[,2])
    
  }
  
  beta_true <- as.matrix(z)
  
  # Center X:
  Xc <- scale(X, scale = F)
  
  # No-noise Y:
  Y_clean <- as.matrix(Xc %*% R0 %*% beta_true)
  
  # Variance of errors:
  var_e <- (1/Rsq - 1)*stats::var(Y_clean)
  
  # Noisy Y:
  Y <- Y_clean +
    as.matrix(stats::rnorm(length(Y_clean), mean = 0, sd = sqrt(var_e)))
  
  return(list(X = X,
              Xc = Xc,
              Y = Y,
              Y_clean = Y_clean,
              basisobj = FEM_basis,
              mesh = mesh,
              coefficient_function = as.matrix(z)))
  
}
