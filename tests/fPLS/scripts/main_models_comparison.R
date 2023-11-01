# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPLS: Models comparison %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

setwd("~/master_thesis/tests/fPLS")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPLS: Models comparison")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Statistical utilities
library(MASS)

# Analysis
library(penR1FPLS)
library(doParallel)
library(fdaPDE2)

# Plots
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)
library(RColorBrewer)

# Time
library(tictoc)


# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

load("../../utils/functions/plot_utilities.RData")

load("scripts/functions/PLS.RData")

source("scripts/functions/generate_2d_data.R")
source("scripts/functions/wrappers.R")
source("scripts/functions/errors_and_times.R")
source("scripts/functions/plots.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Mesh options
num_grid_y_axes <- num_grid_x_axes <- 30

# Number of samples
N <- 200

# Data options
NSR.X <- 0.5^2
NSR.Y <- 0.5^2
H <- 4

# Tests
B.index <- 5

# Test options
lambdas_in <- 10^c(seq(-9, 1, by = 1), 4)
nComp <- 5

# Batches
N.batches <- 10


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

# Images directory
directory.images <- "images/models_comparison/"
if (!file.exists(directory.images)){
  dir.create(directory.images)
}
directory.images <- paste(directory.images, "B.index_", B.index, "/", sep = "")
if (!file.exists(directory.images)){
  dir.create(directory.images)
}

# Results directory
directory.results <- "results/models_comparison/"
if (!file.exists(directory.results)){
  dir.create(directory.results)
}

# Room for results
data <- list()
results <- list()
errors <- list()
time <- list()


# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")

# Dimensions
K <- S <- num_grid_y_axes * num_grid_x_axes


# ||||||||||
# Tests ----
# ||||||||||

cat.section_title("Tests")

for(i in 1:N.batches){
  
  cat(paste("\n\nBatch ", i, ":", sep = ""))
  
  # Room for results
  data[[i]] <- list()
  results[[i]] <- list()
  errors[[i]] <- list()
  time[[i]] <- list()
  
  # Data
  # Set seed for tests reproducibility 
  set.seed(i*100)
  
  # Data
  L <- generate_2d_data(num_grid_y_axes, num_grid_x_axes,
                        N = N,
                        H = H,
                        B.index = B.index,
                        NSR.X = NSR.X,
                        NSR.Y = NSR.Y)
  X_batch <- L$X
  Y_batch <- L$Y
  X_mean <- L$X_mean
  Y_mean <- L$Y_mean
  X_clean_batch <- L$X_clean
  Y_clean_batch <- L$Y_clean
  B_true <- L$B_true
  
  # FEM data
  FEM_basis <- L$basisobj
  mesh <- L$mesh
  nodes <- mesh$nodes
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEM_basis)
  R1 <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEM_basis)
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  data[[i]] <- L
  
  # PLS
  cat("\n - PLS: ")
  tic(quiet = TRUE)
  results[[i]][["PLS"]] <- wrapper_PLS()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["PLS"]] <- compute_errors(results[[i]][["PLS"]])
  time[[i]][["PLS"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["PLS"]], "sec elapsed"))
  
  # fPLS_cpp
  cat("\n - fPLS_cpp: ")
  tic(quiet = TRUE)
  results[[i]][["fPLS_cpp"]] <- wrapper_fPLS_cpp()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["fPLS_cpp"]] <- compute_errors(results[[i]][["fPLS_cpp"]])
  time[[i]][["fPLS_cpp"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPLS_cpp"]], "sec elapsed"))
  
  # Save mean results for the fPCA regression
  # (Data centering is not implemented yet in fPCA)
  Y_mean_fPLS <- results[[i]][["fPLS_cpp"]]$Y_mean
  X_mean_fPLS <- results[[i]][["fPLS_cpp"]]$X_mean
  
  # fPCA regression
  cat("\n - fPCA regression: ")
  tic(quiet = TRUE)
  results[[i]][["fPCA_regression"]] <- wrapper_fPCA()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["fPCA_regression"]] <- compute_errors(results[[i]][["fPCA_regression"]])
  time[[i]][["fPCA_regression"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPCA_regression"]], "sec elapsed"))
  
  # fPLS_R_unique
  cat("\n - fPLS_R_unique: ")
  tic(quiet = TRUE)
  results[[i]][["fPLS_R_unique"]] <- wrapper_fPLS_R_unique()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["fPLS_R_unique"]] <- compute_errors(results[[i]][["fPLS_R_unique"]])
  time[[i]][["fPLS_R_unique"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPLS_R_unique"]], "sec elapsed"))
  
  # fPLS_R_seq
  cat("\n - fPLS_R_seq: ")
  tic(quiet = TRUE)
  results[[i]][["fPLS_R_seq"]] <- wrapper_fPLS_R_seq()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["fPLS_R_seq"]] <- compute_errors(results[[i]][["fPLS_R_seq"]])
  time[[i]][["fPLS_R_seq"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPLS_R_seq"]], "sec elapsed"))
  
}

# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

save(num_grid_y_axes, num_grid_x_axes,
     N, N.batches,
     B.index,
     NSR.X, NSR.Y, H,
     nComp, lambdas_in,
     data, results, errors, time,
     file = paste(directory.results, "results_B.index_", B.index, format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))


# |||||||||||||||||
# Plot Results ----
# |||||||||||||||||

cat.section_title("Plot Results")

# Load data
# list.files(directory.results)
# load(paste(directory.results, tail(list.files(directory.results), n = 1), sep = ""))

# Models
models <- c("fPCA_regression",
            "fPLS_cpp",
            "fPLS_R_seq", "fPLS_R_unique",
            "PLS")

# Colors
m_colors <- c(brewer.pal(3, "Oranges")[3],
              brewer.pal(3, "Blues")[3],
              brewer.pal(3, "Greens")[3:2],
              brewer.pal(3, "Purples")[3])

## RMSE and Time ----
## ||||||||||||||||||

cat.subsection_title("RMSE and Time")

plot <- plot.results_comparison(errors_by_components(), reorganize_times())
ggsave(paste(directory.images, "comparison.jpg", sep = ""),
       plot = plot, width = 10, height = 13, dpi = 200)


## Qualitative results ----
## ||||||||||||||||||||||||

# Here the results for the first batch are considered

nodes <- data[[1]]$mesh$nodes

### X reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("X reconstruction")

# Images directory
directory.images_X <- paste(directory.images, "X/", sep = "")
if (!file.exists(directory.images_X)){
  dir.create(directory.images_X)
}

# Mean

X_mean <- data[[1]]$X_mean

fields <- list(X_mean)
for(m in models){
  fields <- c(fields, list(results[[1]][[m]]$X_mean))
}
plot <- plot.fields(nodes, fields, 3,
                    titles =  c("X_mean", paste("X_mean_", models, sep = "")))
ggsave(paste(directory.images_X, "mean.jpg", sep = ""),
       plot = plot, width = 3*3, height = 3*2, dpi = 200)

# X

X_clean <- data[[1]]$X_clean
X <- data[[1]]$X

n <- 4
set.seed(0)
indexes <- sort(sample(1:N, n, replace = FALSE))
fields <- NULL
titles <- NULL
for(id in indexes){
  fields <- rbind(fields, X_clean[id,], X[id,])
  titles <- c(titles, paste("X", id, c("_clean", ""), sep = ""))
  for(m in models){
    fields <- rbind(fields, results[[1]][[m]]$X_hat[[H]][id,])
    titles <- c(titles, paste("X", id, "_", m, sep = ""))
  }
}
fields <- split(fields, seq_len(length(indexes)*(2+length(models))))
plot <- plot.fields(nodes, fields, 7, titles = titles)
ggsave(paste(directory.images_X, "X.jpg", sep = ""),
       plot = plot, width = 3*(2+length(models)), height = 3*n, dpi = 200,
       device = "jpeg", path = NULL)

# X_c

X_c_clean <- data[[1]]$X_clean - rep(1,N) %*% t(data[[1]]$X_mean)
X_c <- data[[1]]$X_center

n <- 4
set.seed(0)
indexes <- sort(sample(1:N, n, replace = FALSE))
fields <- NULL
titles <- NULL
for(id in indexes){
  fields <- cbind(fields, X_c_clean[id,], X_c[id,])
  titles <- c(titles, paste("X", id, c("_clean", ""), sep = ""))
  for(m in models){
    fields <- cbind(fields, results[[1]][[m]]$X_hat[[H]][id,] - results[[1]][[m]]$X_mean)
    titles <- c(titles, paste("X", id, "_", m, sep = ""))
  }
}
fields <- split(t(fields), seq_len(length(indexes)*(2+length(models))))
plot <- plot.fields(nodes, fields, 7, titles = titles)
ggsave(paste(directory.images_X, "X_c.jpg", sep = ""),
       plot = plot, width = 3*(2+length(models)), height = 3*n, dpi = 200,
       device = "jpeg", path = NULL)


### Y reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("Y reconstruction")

# Images directory
directory.images_Y <- paste(directory.images, "Y/", sep = "")
if (!file.exists(directory.images_Y)){
  dir.create(directory.images_Y)
}

Y_clean <- data[[1]]$Y_clean
Y <- data[[1]]$Y

for(m in models){
  plot <- plot.Y_scatterplots(Y_clean, Y, results[[1]][[m]]$Y_hat)
  ggsave(paste(directory.images_Y, m, ".jpg", sep = ""),
         plot = plot, width = 3*2, height = 3*3, dpi = 200)
}


### B reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("Y reconstruction")

# Images directory
directory.images_B <- paste(directory.images, "B/", sep = "")
if (!file.exists(directory.images_B)){
  dir.create(directory.images_B)
}

B_true <- data[[1]]$B_true

for(m in models){
  n <- 3
  fields <- c(list(B_true), results[[1]][[m]]$B_hat)
  plot <- plot.fields(nodes, fields, n,
                             titles = c("B_true", paste("B_hat (", 1:nComp, " comp)", sep = "")))
  ggsave(paste(directory.images_B, m, ".jpg", sep = ""),
         plot = plot, width = 3*n, height = 3*n, dpi = 200)
}


### Directions ----
### |||||||||||||||

cat.subsection_title("Directions")

n <- nComp
fields <- NULL
titles <- c()
for(m in models[-1]){
  fields <- cbind(fields, results[[1]][[m]]$W_hat)
  titles <- c(titles, paste("W_hat", 1:nComp, " (", m, ")", sep = ""))
}
fields <- split(t(fields), seq_len((length(models)-1)*n))
plot <- plot.fields(nodes, fields, n, titles = titles)
ggsave(paste(directory.images, "W.jpg", sep = ""),
       plot = plot, width = 3*n, height = 3*(length(models)-1), dpi = 200)


### Loadings ----
### |||||||||||||

cat.subsection_title("Loadings")

n <- nComp
fields <- NULL
titles <- c()
for(m in models){
  if(m != models[1]){
    fields <- cbind(fields, results[[1]][[m]]$C_hat)
    titles <- c(titles, paste("C_hat", 1:nComp, " (", m, ")", sep = ""))
  } else{
    fields <- cbind(fields, results[[1]][[m]]$F_hat)
    titles <- c(titles, paste("F_hat", 1:nComp, " (", m, ")", sep = ""))
  }
}
fields <- split(t(fields), seq_len(length(models)*n))
plot <- plot.fields(nodes, fields, n, titles = titles)
ggsave(paste(directory.images, "C.jpg", sep = ""),
       plot = plot, width = 3*n, height = 3*length(models), dpi = 200)



