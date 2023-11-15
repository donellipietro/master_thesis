# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPCA: fPCA vs fPCACS %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

setwd("~/master_thesis/tests/fPCA")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPCA: fPCA vs fPCACS")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Statistical utilities
library(MASS)

# Statistical methods
library(fdaPDE)
library(fdaPDE2)
library(SpatialPCA)

# Times traking
library(tictoc)

# Visualization
library(dplyr)
library(tidyr)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
library(viridis)


# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

source("scripts/functions/data_generation.R")
source("scripts/functions/wrappers.R")
source("scripts/functions/errors_and_times.R")
source("scripts/functions/plots.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Number of statistical units
N <- c(400)

# Mesh dimension
num_grid_axes <- c(20)
K <- num_grid_axes^2

# Number of batches
N.batches <- 10

# Data options
NSR.X <- 0.5^2

# Model options
nComp_vect <- 1:5
lambdas <- 10^c(seq(-6, 2, by = 0.5), 4)


# ||||||||||
# Tests ----
# ||||||||||

cat.section_title("Tests")

# Room for results
results <- list()
errors <- list()
time <- list()
data <- list()


for(nComp in nComp_vect){
  
  cat(paste("\nTests nComp = ", nComp, sep = ""))
  
  # Room for results
  name_N <- paste("N", N, sep = "") 
  results[[nComp]] <- list()
  errors[[nComp]] <- list()
  time[[nComp]] <- list()
  data[[nComp]] <- list()
  
  # Data generation
  set.seed(0)
  data_batch <- generate_2d_data(sqrt(K), sqrt(K),
                                 N = N,
                                 H = nComp,
                                 NSR.X = NSR.X)
  X_batch <- data_batch$X_center
  X_clean_batch <- data_batch$X_clean_center
  S_true_batch <-  data_batch$S_true
  F_true <-  data_batch$F_true
  
  # FEM data
  FEM_basis <- data_batch$basisobj
  mesh <- data_batch$mesh
  nodes <- mesh$nodes
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  lambda_s_opt <- wrapper_fPCA_GCV(lambdas)$lambda_s_opt
  
  for(i in 1:N.batches){
    
    cat(paste("\n\nBatch ", i, ":", sep = ""))
    
    set.seed(i)
    
    # Data generation
    set.seed(0)
    data_batch <- generate_2d_data(sqrt(K), sqrt(K),
                                   N = N,
                                   H = nComp,
                                   NSR.X = NSR.X)
    X_batch <- data_batch$X_center
    X_clean_batch <- data_batch$X_clean_center
    S_true_batch <-  data_batch$S_true
    F_true <-  data_batch$F_true
    
    # FEM data
    FEM_basis <- data_batch$basisobj
    mesh <- data_batch$mesh
    nodes <- mesh$nodes
    mesh_data <- list(
      "nodes"    = mesh$nodes,
      "edges"    = mesh$edges,
      "elements" = mesh$triangles,
      "neigh"    = mesh$neighbors,
      "boundary" = mesh$nodesmarkers
    )
    R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = data_batch$basisobj)
    S <- nrow(nodes)
    
    # Room for results
    results[[nComp]][[i]] <- list()
    errors[[nComp]][[i]] <- list()
    time[[nComp]][[i]] <- list()
    
    # Save data
    data[[nComp]][[i]] <- list(nodes = nodes,
                               X = X_batch,
                               X_clean = X_clean_batch,
                               S_true = S_true_batch,
                               F_true = F_true)
    
    # PCA
    cat("\n - PCA: ")
    results[[nComp]][[i]][["PCA"]] <- wrapper_PCA()
    errors[[nComp]][[i]][["PCA"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["PCA"]])
    time[[nComp]][[i]][["PCA"]] <- results[[nComp]][[i]][["PCA"]]$time
    cat(paste(time[[nComp]][[i]][["PCA"]], "sec elapsed"))
    
    # fPCA
    cat("\n - fPCA: ")
    results[[nComp]][[i]][["fPCA"]] <- wrapper_fPCA_fixed(lambda_s_opt)
    errors[[nComp]][[i]][["fPCA"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["fPCA"]])
    time[[nComp]][[i]][["fPCA"]] <- results[[nComp]][[i]][["fPCA"]]$time
    cat(paste(time[[nComp]][[i]][["fPCA"]], "sec elapsed"))
    
    # fPCACS
    cat("\n - fPCACS: ")
    results[[nComp]][[i]][["fPCACS"]] <- wrapper_fPCACS(ML = FALSE, IT = FALSE, lambda_s_opt)
    errors[[nComp]][[i]][["fPCACS"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["fPCACS"]])
    time[[nComp]][[i]][["fPCACS"]] <- results[[nComp]][[i]][["fPCACS"]]$time
    cat(paste(time[[nComp]][[i]][["fPCACS"]], "sec elapsed"))
    
    # fPCACS_IT
    cat("\n - fPCACS_IT: ")
    results[[nComp]][[i]][["fPCACS_IT"]] <- wrapper_fPCACS(ML = FALSE, IT = TRUE, lambda_s_opt)
    errors[[nComp]][[i]][["fPCACS_IT"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["fPCACS_IT"]])
    time[[nComp]][[i]][["fPCACS_IT"]] <- results[[nComp]][[i]][["fPCACS_IT"]]$time
    cat(paste(time[[nComp]][[i]][["fPCACS_IT"]], "sec elapsed"))
    
    
    # fPCACS_ML
    cat("\n - fPCACS_ML: ")
    results[[nComp]][[i]][["fPCACS_ML"]] <- wrapper_fPCACS(ML = TRUE, IT = FALSE, lambda_s_opt)
    errors[[nComp]][[i]][["fPCACS_ML"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["fPCACS_ML"]])
    time[[nComp]][[i]][["fPCACS_ML"]] <- results[[nComp]][[i]][["fPCACS_ML"]]$time
    cat(paste(time[[nComp]][[i]][["fPCACS_ML"]], "sec elapsed"))
    
    # fPCACS_ML_IT
    cat("\n - fPCACS_ML_IT: ")
    results[[nComp]][[i]][["fPCACS_ML_IT"]] <- wrapper_fPCACS(ML = TRUE, IT = TRUE, lambda_s_opt)
    errors[[nComp]][[i]][["fPCACS_ML_IT"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["fPCACS_ML_IT"]])
    time[[nComp]][[i]][["fPCACS_ML_IT"]] <- results[[nComp]][[i]][["fPCACS_ML_IT"]]$time
    cat(paste(time[[nComp]][[i]][["fPCACS_ML_IT"]], "sec elapsed"))
    
    # SpatialPCA
    if(K <= 900){
      cat("\n - SpatialPCA: ")
      results[[nComp]][[i]][["SpatialPCA"]] <- wrapper_SpatialPCA()
      errors[[nComp]][[i]][["SpatialPCA"]] <- compute_errors(data[[nComp]][[i]], results[[nComp]][[i]][["SpatialPCA"]])
      time[[nComp]][[i]][["SpatialPCA"]] <- results[[nComp]][[i]][["SpatialPCA"]]$time
      cat(paste(time[[nComp]][[i]][["SpatialPCA"]], "sec elapsed"))
    }
    
  }
  cat("\n")
}

cat("\n")

# |||||||||||||||||
# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

results.directory <- "results/speed_up/"
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

save(# Parameters
  nComp_vect, N.batches, NSR.X,
  # Results
  data, results, errors, time,
  # Locations
  file = paste(results.directory,
               "results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData",
               sep = ""))


# |||||||||||||||||
# Plot Results ----
# |||||||||||||||||

cat.section_title("Plot Results")

# Load data
# results.directory <- "results/speed_up/"
# list.files(results.directory)
# load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))


# Models
models_all <- c("fPCA",
                "fPCACS", "fPCACS_IT",
                "fPCACS_ML", "fPCACS_ML_IT", 
                "SpatialPCA", 
                "PCA")

# Names
m_names_all <- c(expression("fPCA"[" "]),
                 expression("fPCA"[" RSVD"]), expression("fPCA"[" RSVD, Seq."]),
                 expression("fPCA"[" RSVD, ML"]), expression(" fPCA"[" RSVD, ML, Seq."]),
                 expression("SpatialPCA"[" "]),
                 expression("MV-PCA"[" "]))


# Colors
m_colors_all <- c(brewer.pal(3, "Oranges")[3],
                  brewer.pal(3, "Greens")[3:2],
                  brewer.pal(3, "Blues")[3:2],
                  brewer.pal(3, "Purples")[3],
                  brewer.pal(3, "Greys")[3])

# Directory
images.directory <- "images/speed_up/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

models <- models_all
m_names <- m_names_all
m_colors <- m_colors_all

times <- NULL
for(nComp in nComp_vect){
  temp <- colMeans(reorganize_times_nComp(nComp)[,-1])
  temp_fPCA <- temp["fPCA"]
  temp_times <- data.frame(Group = paste(nComp, "comp."),
                           fPCACS = temp_fPCA/temp["fPCACS"],
                           fPCACS_IT = temp_fPCA/temp["fPCACS_IT"],
                           fPCACS_ML = temp_fPCA/temp["fPCACS_ML"],
                           fPCACS_ML_IT = temp_fPCA/temp["fPCACS_ML_IT"])
  times <- rbind(times, temp_times)
}


models <- models_all[2:5]
m_names <- m_names_all[2:5]
m_colors <- m_colors_all[2:5]

plot <- plot.speed_up(times)
ggsave(paste(images.directory, "speed-up.pdf", sep = ""),
       plot = plot, width = 16, height = 7, dpi = "print", unit = "cm")

