# Performances assessment ----
# ||||||||||||||||||||||||||||

Y <- Yp

## ||||||||||||||||
## Directories ----
## ||||||||||||||||

# Results
directory.results_path <- paste(directory.results, name, "/", sep = "")
if (!file.exists(directory.results_path)){
  dir.create(directory.results_path)
}
directory.results_path <- paste(directory.results, name, "/CV/", sep = "")
if (!file.exists(directory.results_path)){
  dir.create(directory.results_path)
}


# Log
log_file <- paste(directory.results, "log", "_", name, format(Sys.time(), "_%Y%m%d_%H%M%S"), ".txt", sep = "")

# Create an empty text file
file.create(log_file)


## |||||||||||||||||||
## Generate folds ----
## |||||||||||||||||||

set.seed(123)
folds <- caret::createFolds(Y, k = nFolds)

# for(i in 1:nFolds){
#   cat(paste("\nFold ", i))
#   print(table(Ycl[folds[[i]]]))
# }


## |||||||||||||||
## Fit models ----
## |||||||||||||||

cat("\nModels fit\n\n")

# cl <- parallel::makeCluster(5)
# doParallel::registerDoParallel(cl)

# AUC_comp_fold <- foreach(i = 1:nFolds,
#                          .packages = c("pROC", "penR1FPLS", "tictoc", "fdaPDE2"),
#                          .combine = "cbind")  %dopar% {

for(i in 1:nFolds){
  
  # Update log
  # Sys.sleep(Sys.getpid()/10000) # To avoid overlaps in the log
  # cat(paste("Fold", i, "started on core", Sys.getpid()),
  #     file = log_file, append = TRUE, sep = "\n")
  
  # # Redirect console log in a file
  # sink(paste(directory.results_path, "log_fit_Prova", i, ".txt" , sep = ""), append = TRUE)    
  
  cat(paste("* Fit", i, "started")) # on core", Sys.getpid(), " "))
  
  # Build train
  Y_fold_train <- Y[-folds[[i]]]
  X_fold_train <- X[-folds[[i]],]
  
  tic()
  
  # PLS model:
  fit <- method(Y_fold_train, X_fold_train, nComp,
                mesh = mesh, lambdas = lambdas, VERBOSE)
  
  elapsed <- toc(quiet = TRUE)
  execution_time <- as.numeric(elapsed$toc - elapsed$tic)
  cat(paste(" - compleated after", execution_time, "seconds\n"))
  # sink()
  
  # Save fit
  save(fit, file = paste(directory.results_path, "AUC_CV_fit_", i, "_", name, ".RData" , sep = ""))
  rm(fit)
  
  # # Update log
  # cat(paste("Fold", i, "compleated by core", Sys.getpid(), "after", execution_time, "seconds"),
  #     file = log_file, append = TRUE, sep = "\n")
} 

# stopCluster(cl)
# rm(cl)

## Compute AUC ----
## ||||||||||||||||

cat("\nAUC indexes computation\n\n")

AUC_comp_fold <- matrix(NaN, nrow = nComp, ncol = nFolds)
colnames(AUC_comp_fold) <- paste0("fold_", 1:nFolds)
rownames(AUC_comp_fold) <- paste0("comp_", 1:nComp)

for(i in 1:nFolds){
  for(h in 1:nComp) {
    
    cat("* Fold", i,"comp.", h, "\n")
    
    # Loas fit on train
    load(paste(directory.results_path, "AUC_CV_fit_", i, "_", name,".RData", sep = ""))
    
    # Build test:
    Y_fold_test <- Y[folds[[i]]]
    Ycl_fold_test <- Ycl[folds[[i]]]
    X_fold_test <- X[folds[[i]], ]
    
    # Predict
    Y_hat_fold_test <- as.numeric(X_fold_test %*% fit$B_hat[[h]])
    
    # Compue AUC
    roc_object <- roc(Ycl_fold_test, Y_hat_fold_test, quiet = TRUE)
    AUC_comp_fold[h,i] <- as.numeric(auc(roc_object))
    
  }
}


## Save results ----
## |||||||||||||||||

save(folds, AUC_comp_fold,
     file = paste(directory.results, "AUC_CV_", name, ".RData" , sep = ""))

## Clean ----
## ||||||||||

rm(folds, AUC_comp_fold)
