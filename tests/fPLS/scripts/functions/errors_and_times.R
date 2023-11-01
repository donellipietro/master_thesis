# Errors and times ----
# |||||||||||||||||||||

## Utilites ----
## |||||||||||||

similarity <- function(f1, f2){
  as.numeric(2*(1 -                   t(f1) %*% R1 %*% f2 / 
                  sqrt( (t(f1) %*% R1 %*% f1)*(t(f2) %*% R1 %*% f2) )))
}

# Computes:
# - The Euclidean norm foa a vector
# - The Frobenius norm for a matrix
# and normalized the results by the number of elements.

RMSE <- function(x1, x2){
  sqrt(sum((x1-x2)^2)/length(x1))
}

## Errors ----
##||||||||||||

compute_errors <- function(results){
  
  error_X_mean <- RMSE(X_mean, results$X_mean)
  error_Y_mean <- RMSE(Y_mean, results$Y_mean)
  
  errors <- list()
  errors[["Y"]] <- c() 
  errors[["X"]] <- c() 
  errors[["B"]] <- c() 
  
  for(i in 1:nComp){
    errors[["Y"]][i] <- RMSE(Y_clean_batch, results$Y_hat[[i]])
    errors[["X"]][i] <- RMSE(X_clean_batch, results$X_hat[[i]])
    errors[["B"]][i] <- similarity(B_true, results$B_hat[[i]])
  }
  
  return(list(error_Y_mean = error_Y_mean,
              error_X_mean = error_X_mean,
              errors = errors))
  
}

errors_by_components <- function() {
  
  labels <- c("Init", paste("Comp.",1:nComp))
  
  errors_Y <- NULL
  errors_X <- NULL
  errors_B <- NULL
  
  for(i in 1:N.batches){
    
    temp_Y <- temp_X <- temp_B <- NULL
    
    for(m in models){
      
      temp_Y <- cbind(temp_Y, as.numeric(c(errors[[i]][[m]]$error_Y_mean, errors[[i]][[m]]$errors$Y)))
      temp_X <- cbind(temp_X, as.numeric(c(errors[[i]][[m]]$error_X_mean, errors[[i]][[m]]$errors$X)))
      temp_B <- cbind(temp_B, as.numeric(c(NaN, errors[[i]][[m]]$errors$B)))
      
    }
    
    colnames(temp_Y) <- models
    colnames(temp_X) <- models
    colnames(temp_B) <- models
    
    errors_Y <- rbind(errors_Y, data.frame(Group = labels, temp_Y))
    errors_X <- rbind(errors_X, data.frame(Group = labels, temp_X))
    errors_B <- rbind(errors_B, data.frame(Group = labels, temp_B))
    
  }
  
  return(list(errors_Y = errors_Y,
              errors_X = errors_X,
              errors_B = errors_B))
  
}


## Times ----
## ||||||||||

reorganize_times <- function(){
  
  # Room for results
  times <- NULL
  
  for(i in 1:N.batches){
    
    temp_times <- NULL
    for(m in models){
      temp_times <- cbind(temp_times, time[[i]][[m]])
    }
    
    colnames(temp_times) <- models
    times <- rbind(times, data.frame(Group = "", temp_times))
    
  }
  
  return(times)
}
