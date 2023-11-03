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
library(grid)
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
lambdas_in <- 10^c(seq(-9, 1, by = 0.2), 4)
nComp <- 5

# Batches
N.batches <- 10

# Number of TPS basis (for B_fPLS and PB_fPLS)
n_basis_tps <- 10

# Number K of folds to do cross-validation
k_folds <- 5


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
  
  
  # TPS data
  gam_fit <- mgcv::gam(X_batch[1, ] ~ s(nodes[ , 1],
                                        nodes[ , 2],
                                        bs = "tp",
                                        k = n_basis_tps))
  # Evaluate basis functions:
  # (rows corresponding to argument values and columns to basis functions)
  # X is approximated by A %*% t(Psi):
  Psi_tps <- stats::model.matrix(gam_fit)
  # Matrix of inner products (mass):
  R0_tps <- matrix(NA, nrow = ncol(Psi_tps), ncol = ncol(Psi_tps))
  # Numerical approx. of the inner products:
  for (ii in 1:nrow(R0_tps)) {
    for (jj in ii:ncol(R0_tps)) {
      df <- as.data.frame(nodes)
      df$z = as.numeric(Psi_tps[, ii]*Psi_tps[, jj])
      R0_tps[ii,jj] <-  penR1FPLS:::getVolume(df)
    } # jj
  } # ii
  R0_tps[lower.tri(R0_tps)] <- R0_tps[upper.tri(R0_tps)]
  
  # Save data
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
  
  # B-fPLS
  cat("\n - B-fPLS: ")
  tic(quiet = TRUE)
  results[[i]][["B_fPLS"]] <- wrapper_B_fPLS()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["B_fPLS"]] <- compute_errors(results[[i]][["B_fPLS"]])
  time[[i]][["B_fPLS"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["B_fPLS"]], "sec elapsed"))
  
  # PB-fPLS
  cat("\n - PB-fPLS: ")
  tic(quiet = TRUE)
  results[[i]][["PB_fPLS"]] <- wrapper_PB_fPLS()
  elapsed <- toc(quiet = TRUE)
  errors[[i]][["PB_fPLS"]] <- compute_errors(results[[i]][["PB_fPLS"]])
  time[[i]][["PB_fPLS"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["PB_fPLS"]], "sec elapsed"))
  
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
            "B_fPLS", "PB_fPLS",
            "PLS")

# Models names
models_names <- c(expression(fPCA[]),
                  expression(fPLS[cpp]),
                  expression(fPLS[R-seq]),  expression(fPLS[R-seq]),
                  expression(B-fPLS[]), expression(PB-fPLS[]),
                  expression(MV-PLS[]))
names(models_names) <- models

# Names
True <- expression(True[])
Noisy <- expression(Noisy[])

# Colors
m_colors <- c(brewer.pal(3, "Oranges")[3],
              brewer.pal(3, "Blues")[3],
              brewer.pal(3, "Greens")[3:2],
              brewer.pal(3, "Purples")[3],
              brewer.pal(3, "Reds")[3])

## RMSE and Time ----
## ||||||||||||||||||

cat.subsection_title("RMSE and Time")

# Figures
figures <- plot.results_comparison(errors_by_components(), reorganize_times())

# Title
title <- textGrob(expression(Comparison),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_comparison <- arrangeGrob(title, figures, heights = c(1, 19))
ggsave(paste(directory.images, "comparison.jpg", sep = ""),
       plot = plot_comparison, width = 21/2, height = 29.7/2, dpi = 200)


## Qualitative results ----
## ||||||||||||||||||||||||

# Here the results for the first batch are considered

nodes <- data[[1]]$mesh$nodes

### X reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("X reconstruction")

# # Images directory
# directory.images_X <- paste(directory.images, "X/", sep = "")
# if (!file.exists(directory.images_X)){
#   dir.create(directory.images_X)
# }

#### X_mean ----
#### |||||||||||

# True
X_mean <- data[[1]]$X_mean

# Title
title <- textGrob(expression(X[mean] ~ reconstruction),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Figures
nrow <- 1
ncol <- length(models)+1
fields <- list(X_mean)
for(m in models){
  fields <- c(fields, list(results[[1]][[m]]$X_mean))
}
figures <- plot.fields(nodes, fields, ncol,
                       titles = c(True, models_names, sep = ""))

# Plot
plot_X_mean <- arrangeGrob(title, figures, heights = c(1, 4))
# ggsave(paste(directory.images_X, "mean.jpg", sep = ""),
#        plot = plot_X_mean, width = 3*ncol, height = 1 + 4*nrow, dpi = 200)
# 

#### X ----
#### ||||||

# True
X_clean <- data[[1]]$X_clean

# Noisy
X <- data[[1]]$X

# Figures
nrow <- 4
ncol <- length(models) + 2
set.seed(0)
indexes <- sort(sample(1:N, nrow, replace = FALSE))
fields <- NULL
titles <- NULL
for(id in indexes){
  fields <- rbind(fields, X_clean[id,], X[id,])
  titles <- c(titles, True, Noisy)
  for(m in models){
    fields <- rbind(fields, results[[1]][[m]]$X_hat[[H]][id,])
    titles <- c(titles, models_names[m])
  }
}
fields <- split(fields, seq_len(nrow*ncol))
figures <- plot.fields(nodes, fields, ncol, titles = titles)

# Title
title <- textGrob(expression(X ~ reconstruction),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_X <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
# ggsave(paste(directory.images_X, "X.jpg", sep = ""),
#        plot = plot_X, width = 3*ncol, height = 1 + 4*nrow, dpi = 200,
#        device = "jpeg", path = NULL)


#### X_centered ----
#### |||||||||||||||

# True
X_c_clean <- data[[1]]$X_clean - rep(1,N) %*% t(data[[1]]$X_mean)

# Noisy
X_c <- data[[1]]$X_center

# Figures
nrow <- 4
ncol <- length(models) + 2
set.seed(0)
indexes <- sort(sample(1:N, nrow, replace = FALSE))
fields <- NULL
titles <- NULL
for(id in indexes){
  fields <- cbind(fields, X_c_clean[id,], X_c[id,])
  titles <- c(titles, True, Noisy)
  for(m in models){
    fields <- cbind(fields, results[[1]][[m]]$X_hat[[H]][id,] - results[[1]][[m]]$X_mean)
    titles <- c(titles, models_names[m])
  }
}
fields <- split(t(fields), seq_len(nrow*ncol))
figures <- plot.fields(nodes, fields, ncol, titles = titles)

# Title
title <- textGrob(expression(X[centered] ~ reconstruction),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_X_c <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
# ggsave(paste(directory.images_X, "X_centered.jpg", sep = ""),
#        plot = plot_X_c, width = 3*ncol, height = 1 + 4*nrow, dpi = 200,
#        device = "jpeg", path = NULL)


#### Global plot ----
#### ||||||||||||||||

plot_X_global <- arrangeGrob(plot_X_mean,
                             plot_X,
                             plot_X_c,
                             heights = c(0.25+1, 0.25+4, 0.25+4))

ggsave(paste(directory.images, "X.jpg", sep = ""),
       plot = plot_X_global, width = 21, height = 29.7, dpi = 200,
       device = "jpeg", path = NULL)


### Y reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("Y reconstruction")

# # Images directory
# directory.images_Y <- paste(directory.images, "Y/", sep = "")
# if (!file.exists(directory.images_Y)){
#   dir.create(directory.images_Y)
# }

# True
Y_clean <- data[[1]]$Y_clean

# Noisy
Y <- data[[1]]$Y

# Figures
nrow <- 4
ncol <- length(models)
plots <- list()
top  <- list()
for(m in models){
  plots[[m]] <- plot.Y_scatterplots(Y_clean, Y, results[[1]][[m]]$Y_hat)
  top[[m]] <-  textGrob(models_names[[m]],
                        gp = gpar(fontsize = 18, fontface = 'bold'))
}
figures <- arrangeGrob(grobs = plots, ncol = ncol)
top <-  arrangeGrob(grobs = top, ncol = ncol)
figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))

# Title
title <- textGrob(expression(Y ~ reconstruction),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_Y <- arrangeGrob(title, figures, heights = c(1, 4*nrow, 5))
ggsave(paste(directory.images, "Y.jpg", sep = ""),
       plot = plot_Y, width = 21, height = 29.7, dpi = 200)


### B reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("B reconstruction")

# # Images directory
# directory.images_B <- paste(directory.images, "B/", sep = "")
# if (!file.exists(directory.images_B)){
#   dir.create(directory.images_B)
# }

# True
B_true <- data[[1]]$B_true

# Figures
nrow <- nComp + 1
ncol <- length(models)
fields <- NULL
titles <- NULL
top <- list()
for(m in models){
  fields <- rbind(fields, as.numeric(B_true))
  titles <- c(titles, "True")
  top[[m]] <-  textGrob(models_names[[m]],
                        gp = gpar(fontsize = 18, fontface = 'bold'))
}
for(h in 1:nComp){
  for(m in models){
    fields <- rbind(fields, as.numeric(results[[1]][[m]]$B_hat[[h]]))
    titles <- c(titles, paste(h, "comp."))
  }
}
fields <- split(fields, seq_len(nrow*ncol))
figures <- plot.fields(nodes, fields, ncol, titles = titles)
top <-  arrangeGrob(grobs = top, ncol = ncol)
figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))

# Title
title <- textGrob(expression(B ~ reconstruction),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_B <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
ggsave(paste(directory.images, "B.jpg", sep = ""),
       plot = plot_B, width = 21, height = 29.7, dpi = 200)


### Directions ----
### |||||||||||||||

cat.subsection_title("Directions")

# Figures
nrow <- nComp
ncol <- length(models) - 1
fields <- NULL
titles <- NULL
top <- list()
for(m in models[-1]){
  top[[m]] <-  textGrob(models_names[[m]],
                        gp = gpar(fontsize = 18, fontface = 'bold'))
}
for(h in 1:nComp){
  for(m in models[-1]){
    fields <- rbind(fields, as.numeric(results[[1]][[m]]$W_hat[,h]))
    titles <- c(titles, paste(h, "comp."))
  }
}
fields <- split(fields, seq_len(nrow*ncol))
figures <- plot.fields(nodes, fields, ncol, titles = titles)
top <-  arrangeGrob(grobs = top, ncol = ncol)
figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))

# Title
title <- textGrob(expression(Directions),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_W <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
ggsave(paste(directory.images, "directions.jpg", sep = ""),
       plot = plot_W, width = 21, height = 29.7, dpi = 200)


### Loadings ----
### |||||||||||||

cat.subsection_title("Loadings")

# Figures
nrow <- nComp
ncol <- length(models)
fields <- NULL
titles <- NULL
top <- list()
for(m in models){
  top[[m]] <-  textGrob(models_names[[m]],
                        gp = gpar(fontsize = 18, fontface = 'bold'))
}
for(h in 1:nComp){
  for(m in models){
    if(m == "fPCA_regression") fields <- rbind(fields, as.numeric(results[[1]][[m]]$F_hat[,h]))
    else fields <- rbind(fields, as.numeric(results[[1]][[m]]$C_hat[,h]))
    titles <- c(titles, paste(h, "comp."))
  }
}
fields <- split(fields, seq_len(nrow*ncol))
figures <- plot.fields(nodes, fields, ncol, titles = titles)
top <-  arrangeGrob(grobs = top, ncol = ncol)
figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))

# Title
title <- textGrob(expression(Directions),
                  gp = gpar(fontsize = 20, fontface = 'bold'))

# Plot
plot_W <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
ggsave(paste(directory.images, "loadings.jpg", sep = ""),
       plot = plot_W, width = 21, height = 29.7, dpi = 200)




