# %%%%%%%%%%%%%%%%%%%%%%%%%
# %% fdaPDE installation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%


load("utils/functions/cat_utilities.RData")
cat.script_title("fdaPDE installation")


# |||||||||||||||||||||||||
# Install dependencies ----
# |||||||||||||||||||||||||

cat.section_title("Install dependencies")

# Rcpp
if (!require(Rcpp)) {
  install.packages("Rcpp")
} else{
  cat("Rcpp is already installed\n")
}

# RcppEigen
if (!require(RcppEigen)) {
  install.packages("RcppEigen")
} else{
  cat("RcppEigen is already installed\n")
}

# fdaPDE (old version)
if (!require(fdaPDE)) {
  install.packages("fdaPDE")
} else{
  cat("fdaPDE is already installed\n")
}


# |||||||||||||||||||||||||||
# Install fdaPDE library ----
# |||||||||||||||||||||||||||

cat.section_title("Install fdaPDE library")

library(Rcpp)

# Remove previous installations

if (require(fdaPDE2)) {
  remove.packages("fdaPDE2")
}

# New installation
path <- "libraries/fdaPDE/wrappers/R"
compileAttributes(path)
install.packages(path, type = "source", repos = NULL)
