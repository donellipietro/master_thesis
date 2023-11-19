# PLS, NIPALS algorithm

PLS <- function(X, Y, A, deflation_Y = FALSE) {
  
  # Centering
  Xc = scale(X, scale = FALSE)
  X.mean = attr(Xc, "scaled:center")
  Yc = scale(Y, scale = FALSE)
  Y.mean = attr(Yc, "scaled:center")
  
  # Dimensions
  N <- dim(Xc)[1]
  S <- dim(Xc)[2]
  L <- dim(Yc)[2]
  
  # Mean values
  Y.MEAN = matrix(Y.mean, ncol = L, nrow = N, byrow = TRUE)
  X.MEAN = matrix(X.mean, ncol = S, nrow = N, byrow = TRUE)
  
  # Room for solutions
  W <- matrix(0, nrow = S, ncol = A)
  V <- matrix(0, nrow = L, ncol = A)
  TT <- matrix(0, nrow = N, ncol = A)
  UU <- matrix(0, nrow = N, ncol = A)
  B <- matrix(0, nrow = A, ncol = A)
  
  # Extra
  C <- matrix(0, nrow = S, ncol = A)
  D <- matrix(0, nrow = L, ncol = A)
  
  # Initilization
  EE <- Xc
  FF <- Yc
  
  # Intermediate steps
  XX <- list()
  YY <- list()
  XX[[1]] <- EE
  YY[[1]] <- FF
  
  for(i in 1:A){
    
    # Decomposition
    C.YX <- t(FF) %*% EE
    SVD <- svd(C.YX, nu = 1, nv = 1)
    W[,i] <- SVD$v # x-direction
    V[,i] <- SVD$u # y-direction
    TT[,i] <- EE %*% W[,i] # x-component
    UU[,i] <- FF %*% V[,i] # y-component
    
    tt = sum(TT[,i]^2)
    
    # Regression
    C[,i] <- t(EE) %*% TT[,i] / tt
    D[,i] <- t(FF) %*% TT[,i] / tt
    B[i,i] <- t(UU[,i]) %*% TT[,i] / tt
    
    # Deflation
    EE <- EE - TT[,i] %*% t(C[,i])
    if(deflation_Y == TRUE)
      FF <- FF - B[i,i] * TT[,i] %*% t(V[,i])
    
    # Intermediate steps
    XX[[i+1]] = EE
    YY[[i+1]] = FF
    
  }
  
  Beta <- W %*% solve(t(C) %*% W, B %*% t(V))
  
  Y_hat <- Xc %*% Beta + Y.MEAN
  X_hat <- TT %*% t(C) + X.MEAN 
  
  return(list(W = W,
              V = V,
              TT = TT,
              UU = UU,
              B = B,
              C = C,
              D = D,
              EE = EE,
              FF = FF,
              XX = XX,
              YY = YY,
              X.mean = X.mean,
              Y.mean = Y.mean,
              Beta = Beta,
              Y_hat = Y_hat,
              X_hat = X_hat))
  
}

save(PLS, file = "scripts/functions/PLS.RData")