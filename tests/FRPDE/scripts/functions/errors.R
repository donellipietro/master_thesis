# Errors ----
# |||||||||||

## Utilities ----
## ||||||||||||||

# Computes:
# - The Euclidean norm foa a vector
# - The Frobenius norm for a matrix
# and normalized the results by the number of elements.

RMSE <- function(x1, x2){
  sqrt(sum((x1-x2)^2)/length(x1))
}


# Compute errors ----
# |||||||||||||||||||

compute_errors <- function(data, results){
  
  RMSE(data$x_clean, results$x_hat)
  
}


# Reorganize errors ----
# ||||||||||||||||||||||


reorganize_errors <- function(){
  
  # Room for results
  errors_ <- NULL
  
  for(NSR.X in NSR.X_vect){
    
    NSR.X_name <-  paste("NSR.X", sprintf("%2.4f", NSR.X), sep = "")
    NSR.X_final <-  paste("SNR = ", sprintf("%2.2f", 1/NSR.X), sep = "")
    
    for(i in 1:N.batches){
      
      temp_errors <- NULL
      for(m in models){
        temp_errors <- cbind(temp_errors, errors[[NSR.X_name]][[i]][[m]])
      }
      
      temp_errors <-data.frame(temp_errors)
      colnames(temp_errors) <- models
      errors_ <- rbind(errors_, cbind(data.frame(Group = NSR.X_final), temp_errors))
      
    }
    
  }
  
  return(errors_)
}
