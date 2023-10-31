# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPLS: Data generation problems %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

setwd("~/master_thesis/tests/fPLS")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPLS: Data generation problems")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Analysis
library(penR1FPLS)
library(doParallel)

# Plots
library(ggplot2)
library(viridis)
library(gridExtra)


# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

load("../../utils/functions/plot_utilities.RData")

source("scripts/functions/generate_2d_data_old.R")
source("scripts/functions/wrappers.R")
source("scripts/functions/plots.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Mesh grid
x <- seq(0, 1, length.out = 20)
y <- seq(0, 1, length.out = 20)

# Number of samples
N <- 60

# Tests
B.index_vect <- c(3, 5)
sd.X_vect <- c(0.2, 0)

# Test options
lambdas_in <- 10^c(seq(-9, 1, by = 1), 4)
nComp <- 5


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

# Image folder
images.directory <- "images/data_generation_problems/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}


# ||||||||||
# Tests ----
# ||||||||||

cat.section_title("Tests")

for(B.index in B.index_vect){
  for(sd.X in sd.X_vect){
    
    cat("- Test: B.index = ", B.index, ", sd.X = ", sd.X, "\n", sep = "")
    
    # Results directory
    images.directory_path <- paste(images.directory, "B.index_", B.index, "/", sep = "")
    if (!file.exists(images.directory_path)){
      dir.create(images.directory_path)
    }
    
    # Set seed for tests reproducibility 
    set.seed(0)
    
    # Data
    L <- generate_2d_data(x, y, N, B.index, sd.X, 0.95)
    Y_batch <- L$Y
    X_batch <- L$X
    B_true <- L$coefficient_function
    X_c <- L$Xc
    FEM_basis <- L$basisobj
    mesh <- L$mesh
    
    # Results
    results <- wrapper_fPLS_R_unique()
    B_hat_fPLS <- results$B_hat
    W_hat_fPLS <- results$W_hat
    C_hat_fPLS <- results$W_hat
    
    # Plot Y distribution
    df = data.frame(x = Y_batch)
    plot <- ggplot(df, aes(x)) + 
      geom_histogram(color = "white", fill = "grey") +
      standard_plot_settings() + 
      xlab("") + ylab("") +
      ggtitle("Y distribution")
    ggsave(paste(images.directory_path, "Y_dist_sd.X", sd.X, ".jpg", sep = ""),
           plot = plot, width = 5, height = 5, dpi = 200)
    
    # Plot B
    plot <- plot.fields(mesh$nodes, c(list(B_true), B_hat_fPLS), 3,
                        titles = c("B_true", paste("B_hat_", 1:nComp, sep = "")))
    ggsave(paste(images.directory_path, "B_sd.X", sd.X, ".jpg", sep = ""),
           plot = plot, width = 3*3, height = 3*2, dpi = 200)
    
    # Plot W
    fields <- split(t(W_hat_fPLS), 1:ncol(W_hat_fPLS))
    plot <- plot.fields(mesh$nodes, fields, 3,
                        titles = paste("W", 1:ncol(W_hat_fPLS), "_hat", sep = ""))
    ggsave(paste(images.directory_path, "W_sd.X", sd.X, ".jpg", sep = ""),
           plot = plot, width = 3*3, height = 3*2, dpi = 300)

    # Plot C
    fields <- split(t(C_hat_fPLS), 1:ncol(C_hat_fPLS))
    plot <- plot.fields(mesh$nodes, fields, 3,
                        titles = paste("C", 1:ncol(C_hat_fPLS), "_hat", sep = ""))
    ggsave(paste(images.directory_path, "C_sd.X", sd.X, ".jpg", sep = ""),
           plot = plot, width = 3*3, height = 3*2, dpi = 300)
    
    # Plot X
    data_plot <- list()
    for(i in 1:6){
      data_plot[[i]] <- X_batch[i,]
    }
    plot <- plot.fields(mesh$nodes, data_plot, 3, 
                        titles = paste("X_", 1:6, sep = ""),
                        type = "points")
    ggsave(paste(images.directory_path, "X_sd.X",sd.X,".jpg", sep = ""),
           plot = plot, width = 3*3, height = 3*2, dpi = 300)
    
    # Plot Xc
    data_plot <- list()
    for(i in 1:6){
      data_plot[[i]] <- X_c[i,]
    }
    plot <- plot.fields(mesh$nodes, data_plot, 3,
                        titles = paste("Xc_", 1:6, sep = ""),
                        type = "points")
    ggsave(paste(images.directory_path, "Xc_sd.X", sd.X, ".jpg", sep = ""),
           plot = plot, width = 3*3, height = 3*2, dpi = 300)

    
  }
}