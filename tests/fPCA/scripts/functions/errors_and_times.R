# Errors and times ----
# |||||||||||||||||||||

## Utilities ----
## ||||||||||||||

# Computes:
# - The Euclidean norm foa a vector
# - The Frobenius norm for a matrix
# and normalized the results by the number of elements.

RMSE <- function(x){
  sqrt(sum(x^2)/length(x))
}

# Adjustes the sign of scores and loadings to be coherent with the one of
# the true ones.

adjust_results <- function(F_hat, S_hat) {
  
  for(h in 1:nComp){
    if (RMSE(F_hat[,h] + F_true[,h]) < RMSE(F_hat[,h] - F_true[,h])) {
      F_hat[,h] = -F_hat[,h]
      S_hat[,h] = -S_hat[,h]
    }
  }
  
  return(list(F_hat = F_hat, S_hat = S_hat))
}


## Errors calculation ----
## |||||||||||||||||||||||

# Computes the error made by a method with respect to the theoretical solutions.

compute_errors <- function(data, results) {
  
  # Reconstruction error
  error_X <- c()
  for(h in 1:nComp){
    error_X[h] <- RMSE(data$X_clean -
                         results$S_hat[,1:h] %*% t(results$F_hat[,1:h]))
  }
  
  # Loadings
  error_F <- c()
  for(h in 1:nComp){
    error_F[h] <- RMSE(data$F_true[,h] - results$F_hat[,h])
  }
  
  # Scores
  error_S <- c()
  for(h in 1:nComp){
    error_S[h] <- RMSE(data$S_true[,h] - results$S_hat[,h])
  }
  
  return(list(error_X = error_X,
              error_F = error_F,
              error_S = error_S))
  
}


## Errors reorganization ----
## ||||||||||||||||||||||||||

# Reorganizes errors by components given a specified K and N.
# Used to generate the plot "results_by_components"

errors_by_components <- function(name_K, name_N) {
  
  # Room for results
  errors_X <- NULL
  errors_F <- NULL
  errors_S <- NULL
  
  if(length(errors[[name_K]][[name_N]]) == 0 ){
    
    errors_X <- matrix(NaN, nrow = N.batches, ncol = 1+length(models))
    errors_F <- matrix(NaN, nrow = N.batches, ncol = 1+length(models))
    errors_S <- matrix(NaN, nrow = N.batches, ncol = 1+length(models))
    
  } else{
    
    for(i in 1:N.batches){
      
      labels <- paste("Comp.",1:nComp)
      temp_X <- temp_F <- temp_S <- NULL
      
      for(m in models){
        if(m %in% names(errors[[name_K]][[name_N]][[i]])){
          temp_X <- cbind(temp_X, errors[[name_K]][[name_N]][[i]][[m]]$error_X)
          temp_F <- cbind(temp_F, errors[[name_K]][[name_N]][[i]][[m]]$error_F)
          temp_S <- cbind(temp_S, errors[[name_K]][[name_N]][[i]][[m]]$error_S)
        } else {
          temp_X <- cbind(temp_X, rep(NaN, nComp))
          temp_F <- cbind(temp_F, rep(NaN, nComp))
          temp_S <- cbind(temp_S, rep(NaN, nComp))
        }
      }
      
      colnames(temp_X) <- models
      colnames(temp_F) <- models
      colnames(temp_S) <- models
      
      errors_X <- rbind(errors_X, data.frame(Group = labels, temp_X))
      errors_F <- rbind(errors_F, data.frame(Group = labels, temp_F))
      errors_S <- rbind(errors_S, data.frame(Group = labels, temp_S))
      
    }
    
  }
  
  return(list(errors_X = errors_X,
              errors_F = errors_F,
              errors_S = errors_S))
}

# Reorganizes times given a specified K and N.
# Used to generate the plot "results_by_components"

reorganize_times <- function(name_K, name_N){
  
  # Room for results
  times <- NULL
  
  if(length(errors[[name_K]][[name_N]]) == 0 ){
    
    times <- matrix(NaN, nrow = N.batches, ncol = 1+length(models))
    
  } else{
    
    for(i in 1:N.batches){
      
      temp_times <- NULL
      for(m in models){
        if(m %in% names(errors[[name_K]][[name_N]][[i]])){
          temp_times <- cbind(temp_times, time[[name_K]][[name_N]][[i]][[m]])
        } else {
          temp_times <- cbind(temp_times, NaN)
        }
      }
      
      colnames(temp_times) <- models
      times <- rbind(times, data.frame(Group = "", temp_times))
      
    }
    
  }
  
  return(times)
}

reorganize_times_nComp <- function(nComp){
  
  # Room for results
  times <- NULL
  
  if(length(errors[[nComp]]) == 0 ){
    
    times <- matrix(NaN, nrow = N.batches, ncol = 1+length(models))
    
  } else{
    
    for(i in 1:N.batches){
      
      temp_times <- NULL
      for(m in models){
        if(m %in% names(errors[[nComp]][[i]])){
          temp_times <- cbind(temp_times, time[[nComp]][[i]][[m]])
        } else {
          temp_times <- cbind(temp_times, NaN)
        }
      }
      
      colnames(temp_times) <- models
      times <- rbind(times, data.frame(Group = "", temp_times))
      
    }
    
  }
  
  return(times)
  
}

# Reorganizes errors by K.
# Used to generate the plot "accuracy"

errors_by_K <- function() {
  
  # Room for results
  errors_X <- list()
  errors_F <- list()
  errors_S <- list()
  
  for(K in K_vect){
    
    name_K <- paste("K", K, sep = "")
    
    errors_X[[name_K]] <- NULL
    errors_F[[name_K]] <- NULL
    errors_S[[name_K]] <- NULL
    
    labels <- paste("Comp.",1:nComp)
    
    for(N in N_vect){
      
      name_N <- paste("N", N, sep = "")
      
      if(length(errors[[name_K]][[name_N]]) == 0 ){
        
        temp_X <- matrix(NaN, nrow = N.batches*nComp, ncol = length(models))
        temp_F <- matrix(NaN, nrow = N.batches*nComp, ncol = length(models))
        temp_S <- matrix(NaN, nrow = N.batches*nComp, ncol = length(models))
        
        colnames(temp_X) <- c(models)
        colnames(temp_F) <- c(models)
        colnames(temp_S) <- c(models)
        
        errors_X[[name_K]] <- rbind(errors_X[[name_K]],
                                    data.frame(Group = N,
                                               Component = rep(labels, N.batches),
                                               temp_X))
        errors_F[[name_K]] <- rbind(errors_F[[name_K]],
                                    data.frame(Group = N,
                                               Component = rep(labels, N.batches),
                                               temp_F))
        errors_S[[name_K]] <- rbind(errors_S[[name_K]],
                                    data.frame(Group = N,
                                               Component = rep(labels, N.batches),
                                               temp_S))
        
      } else {
        
        for(i in 1:N.batches){
          
          temp_X <- temp_F <- temp_S <- NULL
          
          for(m in models){
            if(m %in% names(errors[[name_K]][[name_N]][[i]])){
              temp_X <- cbind(temp_X, errors[[name_K]][[name_N]][[i]][[m]]$error_X)
              temp_F <- cbind(temp_F, errors[[name_K]][[name_N]][[i]][[m]]$error_F)
              temp_S <- cbind(temp_S, errors[[name_K]][[name_N]][[i]][[m]]$error_S)
            } else {
              temp_X <- cbind(temp_X, rep(NaN, nComp))
              temp_F <- cbind(temp_F, rep(NaN, nComp))
              temp_S <- cbind(temp_S, rep(NaN, nComp))
            }
          }
          
          colnames(temp_X) <- models
          colnames(temp_F) <- models
          colnames(temp_S) <- models
          
          errors_X[[name_K]] <- rbind(errors_X[[name_K]],
                                      data.frame(Group = N,
                                                 Component = labels,
                                                 temp_X))
          errors_F[[name_K]] <- rbind(errors_F[[name_K]],
                                      data.frame(Group = N,
                                                 Component = labels,
                                                 temp_F))
          errors_S[[name_K]] <- rbind(errors_S[[name_K]],
                                      data.frame(Group = N,
                                                 Component = labels,
                                                 temp_S))
          
        }
        
      }
      
    }
    
  }
  
  return(list(errors_X = errors_X,
              errors_F = errors_F,
              errors_S = errors_S))
}
