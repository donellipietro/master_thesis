# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPCA: Models comparison %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

setwd("~/master_thesis/tests/fPCA")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPCA: Models comparison")


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
library(ggtext)
require(grid)
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
N_vect <- c(50, 100, 200, 400, 800, 1600)

# Mesh dimension
num_grid_axes_vect <- c(20, 30, 40, 50)
K_vect <- num_grid_axes_vect^2

# Number of batches
N.batches <- 10

# Data options
NSR.X <- 0.5^2

# Model options
nComp <- 3
lambdas <- 10^c(seq(-6, 2, by = 0.5), 4)

parameters <- list(N_vect = N_vect,
                   num_grid_axes_vect = num_grid_axes_vect,
                   K_vect = K_vect,
                   N.batches = N.batches,
                   NSR.X = NSR.X,
                   nComp = nComp,
                   lambdas = lambdas)
print(parameters)


# ||||||||||
# Tests ----
# ||||||||||

cat.section_title("Tests")

# Room for results
results <- list()
errors <- list()
time <- list()
data <- list()

for(K in K_vect){
  
  # Room for results
  name_K <- paste("K", K, sep = "") 
  results[[name_K]] <- list()
  errors[[name_K]] <- list()
  time[[name_K]] <- list()
  data[[name_K]] <- list()
  
  for(N in N_vect){
    
    cat(paste("\nTests K = ", K, " N = ", N, sep = ""))
    
    # Room for results
    name_N <- paste("N", N, sep = "") 
    results[[name_K]][[name_N]] <- list()
    errors[[name_K]][[name_N]] <- list()
    time[[name_K]][[name_N]] <- list()
    data[[name_K]][[name_N]] <- list()
    
    # Data generation
    set.seed(0)
    data_batch <- generate_2d_data(sqrt(K), sqrt(K),
                                   N = N,
                                   H = 3,
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
      
      # Data generation
      set.seed(N*K + i)
      data_batch <- generate_2d_data(sqrt(K), sqrt(K),
                                     N = N,
                                     H = 3,
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
      results[[name_K]][[name_N]][[i]] <- list()
      errors[[name_K]][[name_N]][[i]] <- list()
      time[[name_K]][[name_N]][[i]] <- list()
      
      # Save data
      data[[name_K]][[name_N]][[i]] <- list(nodes = nodes,
                                            X = X_batch,
                                            X_clean = X_clean_batch,
                                            S_true = S_true_batch,
                                            F_true = F_true)
      
      # PCA
      cat("\n - PCA: ")
      results[[name_K]][[name_N]][[i]][["PCA"]] <- wrapper_PCA()
      errors[[name_K]][[name_N]][[i]][["PCA"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["PCA"]])
      time[[name_K]][[name_N]][[i]][["PCA"]] <- results[[name_K]][[name_N]][[i]][["PCA"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["PCA"]], "sec elapsed"))
      
      # fPCA
      cat("\n - fPCA: ")
      results[[name_K]][[name_N]][[i]][["fPCA"]] <- wrapper_fPCA_fixed(lambda_s_opt)
      errors[[name_K]][[name_N]][[i]][["fPCA"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["fPCA"]])
      time[[name_K]][[name_N]][[i]][["fPCA"]] <- results[[name_K]][[name_N]][[i]][["fPCA"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCA"]], "sec elapsed"))
      
      # fPCACS
      cat("\n - fPCACS: ")
      results[[name_K]][[name_N]][[i]][["fPCACS"]] <- wrapper_fPCACS(ML = FALSE, IT = FALSE, lambda_s_opt)
      errors[[name_K]][[name_N]][[i]][["fPCACS"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["fPCACS"]])
      time[[name_K]][[name_N]][[i]][["fPCACS"]] <- results[[name_K]][[name_N]][[i]][["fPCACS"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS"]], "sec elapsed"))
      
      # fPCACS_IT
      cat("\n - fPCACS_IT: ")
      results[[name_K]][[name_N]][[i]][["fPCACS_IT"]] <- wrapper_fPCACS(ML = FALSE, IT = TRUE, lambda_s_opt)
      errors[[name_K]][[name_N]][[i]][["fPCACS_IT"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["fPCACS_IT"]])
      time[[name_K]][[name_N]][[i]][["fPCACS_IT"]] <- results[[name_K]][[name_N]][[i]][["fPCACS_IT"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS_IT"]], "sec elapsed"))
      
      
      # fPCACS_ML
      cat("\n - fPCACS_ML: ")
      results[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- wrapper_fPCACS(ML = TRUE, IT = FALSE, lambda_s_opt)
      errors[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["fPCACS_ML"]])
      time[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- results[[name_K]][[name_N]][[i]][["fPCACS_ML"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS_ML"]], "sec elapsed"))
      
      # fPCACS_ML_IT
      cat("\n - fPCACS_ML_IT: ")
      results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- wrapper_fPCACS(ML = TRUE, IT = TRUE, lambda_s_opt)
      errors[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]])
      time[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]], "sec elapsed"))
      
      # SpatialPCA
      if(K <= 900){
        cat("\n - SpatialPCA: ")
        results[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- wrapper_SpatialPCA()
        errors[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- compute_errors(data[[name_K]][[name_N]][[i]], results[[name_K]][[name_N]][[i]][["SpatialPCA"]])
        time[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- results[[name_K]][[name_N]][[i]][["SpatialPCA"]]$time
        cat(paste(time[[name_K]][[name_N]][[i]][["SpatialPCA"]], "sec elapsed"))
      }
      
    }
    cat("\n")
  }
}

cat("\n")

# |||||||||||||||||
# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

results.directory <- "results/models_comparison/"
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

save(# Parameters
  N_vect, K_vect, N.batches, NSR.X, nComp,
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
# results.directory <- "results/models_comparison/"
# list.files(results.directory)
# load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))

# Models
models_all <- c("fPCA",
                "fPCACS", "fPCACS_IT",
                "fPCACS_ML", "fPCACS_ML_IT", 
                "SpatialPCA", 
                "PCA")

m_names_all <- c(expression("  fPCA"[" "]),
                 expression("  fPCA"[" RSVD"]), expression(" fPCA"[" RSVD, Seq."]),
                 expression(" fPCA"[" RSVD, ML"]), expression("fPCA"[" RSVD, ML, Seq."]),
                 expression("  SpatialPCA"[" "]),
                 expression("   MV-PCA"[" "]))

# Colors
m_colors_all <- c(brewer.pal(3, "Oranges")[3],
                  brewer.pal(3, "Greens")[3:2],
                  brewer.pal(3, "Blues")[3:2],
                  brewer.pal(3, "Purples")[3],
                  brewer.pal(3, "Greys")[3])

K_vect_all <- K_vect

# Directory
images.directory <- "images/models_comparison/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}

## Qualitative Results ----
## ||||||||||||||||||||||||

# cat.subsection_title("Qualitative Results")
# 
# images.directory_path <- paste(images.directory, "qualitative_results/", sep = "")
# if (!file.exists(images.directory_path)){
#   dir.create(images.directory_path)
# }
# 
# for(K in K_vect){
#   name_K <- paste("K", K, sep = "")
#   for(N in N_vect){
#     name_N <- paste("N", N, sep = "")
#     for(m in models){
# 
#       plot <- plot.results(data[[name_K]][[name_N]][[1]],
#                            results[[name_K]][[name_N]][[1]][[m]])
# 
#       ggsave(paste(images.directory_path, name_K,"_", name_N, "_", m, ".jpg", sep = ""),
#              plot = plot, width = 10, height = 11, dpi = 200)
# 
#     }
#   }
# }


## Errors by components ----
## |||||||||||||||||||||||||

# cat.subsection_title("Errors by components")
# 
# images.directory_path <- paste(images.directory, "by_components/", sep = "")
# if (!file.exists(images.directory_path)){
#   dir.create(images.directory_path)
# }
# 
# for(K in K_vect){
#   name_K <- paste("K", K, sep = "") 
#   for(N in N_vect){
#     name_N <- paste("N", N, sep = "") 
#       
#       plot <- plot.results_by_components(errors_by_components(name_K, name_N),
#                                          reorganize_times(name_K, name_N))
#       ggsave(paste(images.directory_path, name_K,"_", name_N, ".jpg", sep = ""),
#              plot = plot, width = 10, height = 11, dpi = 200)
# 
#   }
# }


## Time complexity analysis ----
## |||||||||||||||||||||||||||||

cat.subsection_title("Time complexity analysis")

images.directory_path <- paste(images.directory, "time_analysis/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}


# Time analysis in N, all

models <- models_all
m_names <- m_names_all
m_colors <- m_colors_all

times <- list()
for(K in K_vect){
  name_K <- paste("K", K, sep = "") 
  times[[name_K]] <- NULL
  for(N in N_vect){
    name_N <- paste("N", N, sep = "")
    temp_times <- data.frame(N = N,
                             value = as.numeric(colMeans(reorganize_times(name_K, name_N)[,-1])))
    times[[name_K]] <- rbind(times[[name_K]], temp_times)
  }
  times[[name_K]]$Model <- rep(models, length(N_vect))
}

models <- models_all
m_names <- m_names_all
m_colors <- m_colors_all
K_vect <- K_vect_all[1:2]

plot <- plot.N_time_analysis(times, title = "Execution times", LOG = FALSE)
ggsave(paste(images.directory_path, "N_all.pdf", sep = ""),
       plot = plot, width = 16, height = 12, dpi = "print", units = "cm")


# Time analysis in N, only fPCA

models <- models_all[-c(6)]
m_names <- m_names_all[-c(6)]
m_colors <- m_colors_all[-c(6)]
K_vect <- K_vect_all

times <- list()
for(K in K_vect){
  name_K <- paste("K", K, sep = "") 
  times[[name_K]] <- NULL
  for(N in N_vect){
    name_N <- paste("N", N, sep = "")
    temp_times <- data.frame(N = N,
                             value = as.numeric(colMeans(reorganize_times(name_K, name_N)[,-1])))
    times[[name_K]] <- rbind(times[[name_K]], temp_times)
  }
  times[[name_K]]$Model <- rep(models, length(N_vect))
}

plot <- plot.N_time_analysis(times, title = "Time complexity in N", LOG = TRUE)
ggsave(paste(images.directory_path, "N_log.pdf", sep = ""),
       plot = plot, width = 16, height = 21, dpi = "print", units = "cm")

plot <- plot.N_time_analysis(times, title = "Execution times", LOG = FALSE)
ggsave(paste(images.directory_path, "N.pdf", sep = ""),
       plot = plot, width = 16, height = 21, dpi = "print", units = "cm")


# Time analysis in K

times <- NULL
for(K in K_vect){
  name_K <- paste("K", K, sep = "")
  for(N in N_vect){
    name_N <- paste("N", N, sep = "")
    temp_times <- data.frame(K = K,
                             value = as.numeric(colMeans(reorganize_times(name_K, name_N)[,-1])),
                             Model = models,
                             N = paste(N))
    times <- rbind(times, temp_times)
  }
}

plot <- plot.K_time_analysis(times)
ggsave(paste(images.directory_path, "K.jpg", sep = ""),
       plot = plot, width = 10, height = 10, dpi = 200)


## Accuracy ----
## |||||||||||||

cat.subsection_title("Accuracy")

images.directory_path <- paste(images.directory, "accuracy/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

# All

models <- models_all
m_names <- m_names_all
m_colors <- m_colors_all
K_vect <- K_vect_all[1:2]

images.directory_path <- paste(images.directory, "accuracy/", "all/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

plot.accuracy(errors_by_K(), images.directory_path, "_all")


# Only fPCA

models <- models_all[-c(6,7)]
m_names <- m_names_all[-c(6,7)]
m_colors <- m_colors_all[-c(6,7)]
K_vect <- K_vect_all

images.directory_path <- paste(images.directory, "accuracy/", "only_fPCA/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

plot.accuracy(errors_by_K(), images.directory_path, "_only_fPCA")




