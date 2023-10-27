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

source("scripts/fPCA_vs_fPCACS/functions_data_generation.R")
source("scripts/fPCA_vs_fPCACS/functions_wrappers.R")
source("scripts/fPCA_vs_fPCACS/functions_errors_and_times.R")
source("scripts/fPCA_vs_fPCACS/functions_plots.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Number of statistical units
N_vect <- c(50, 100, 200, 400, 800, 1600)

# Mesh dimension
num_grid_axes_vect <- c(10, 20, 20 , 30, 40)
K_vect <- num_grid_axes_vect^2

# Number of batches
N.batches <- 10
N.tot <- max(N_vect)*N.batches

# Data options
NSR.X <- 0.2^2

# Model options
nComp <- 3
lambdas_in <- 10^(-4)
lambdas <- 10^(-4)

parameters <- list(N_vect = N_vect,
                   num_grid_axes_vect = num_grid_axes_vect,
                   K_vect = K_vect,
                   N.batches = N.batches,
                   N.tot = N.tot,
                   NSR.X = NSR.X,
                   nComp = nComp,
                   lambdas_in = lambdas_in)
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
  
  # Data generation
  set.seed(0)
  L <- generate_2d_data(sqrt(K), sqrt(K),
                        N = N.tot,
                        H = 3,
                        NSR.X = NSR.X)
  X <- L$X_center
  X_clean <- L$X_clean_center
  S_true <-  L$S_true
  F_true <-  L$F_true
  
  # FEM data
  FEM_basis <- L$basisobj
  mesh <- L$mesh
  nodes <- mesh$nodes
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  R0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis = L$basisobj)
  S <- nrow(nodes)
  
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
    
    for(i in 1:N.batches){
      
      cat(paste("\n\nBatch ", i, ":", sep = ""))
      
      # Room for results
      results[[name_K]][[name_N]][[i]] <- list()
      errors[[name_K]][[name_N]][[i]] <- list()
      time[[name_K]][[name_N]][[i]] <- list()
      
      # Data
      X_batch <- X[1:N + (i-1)*N,]
      X_clean_batch <- X_clean[1:N + (i-1)*N,]
      S_batch_true <- S_true[1:N + (i-1)*N,]
      
      data[[name_K]][[name_N]][[i]] <- list(nodes = nodes,
                                            X_batch = X_batch,
                                            X_clean = X_clean_batch,
                                            F_true = F_true,
                                            S_true = S_batch_true)
      
      # PCA
      cat("\n - PCA: ")
      results[[name_K]][[name_N]][[i]][["PCA"]] <- wrapper_PCA()
      errors[[name_K]][[name_N]][[i]][["PCA"]] <- compute_errors(results[[name_K]][[name_N]][[i]][["PCA"]])
      time[[name_K]][[name_N]][[i]][["PCA"]] <- results[[name_K]][[name_N]][[i]][["PCA"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["PCA"]], "sec elapsed"))
      
      # fPCA
      cat("\n - fPCA: ")
      results[[name_K]][[name_N]][[i]][["fPCA"]] <- wrapper_fPCA()
      errors[[name_K]][[name_N]][[i]][["fPCA"]] <- compute_errors(results[[name_K]][[name_N]][[i]][["fPCA"]])
      time[[name_K]][[name_N]][[i]][["fPCA"]] <- results[[name_K]][[name_N]][[i]][["fPCA"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCA"]], "sec elapsed"))
      
      # fPCACS_ML
      cat("\n - fPCACS_ML: ")
      results[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- wrapper_fPCACS(ML = TRUE, IT = FALSE)
      errors[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- compute_errors(results[[name_K]][[name_N]][[i]][["fPCACS_ML"]])
      time[[name_K]][[name_N]][[i]][["fPCACS_ML"]] <- results[[name_K]][[name_N]][[i]][["fPCACS_ML"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS_ML"]], "sec elapsed"))
      
      # fPCACS_ML_IT
      cat("\n - fPCACS_ML_IT: ")
      results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- wrapper_fPCACS(ML = TRUE, IT = TRUE)
      errors[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- compute_errors(results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]])
      time[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]] <- results[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]]$time
      cat(paste(time[[name_K]][[name_N]][[i]][["fPCACS_ML_IT"]], "sec elapsed"))

      # # SpatialPCA
      # cat("\n - SpatialPCA: ")
      # results[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- wrapper_SpatialPCA()
      # errors[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- compute_errors(results[[name_K]][[name_N]][[i]][["SpatialPCA"]])
      # time[[name_K]][[name_N]][[i]][["SpatialPCA"]] <- results[[name_K]][[name_N]][[i]][["SpatialPCA"]]$time
      # cat(paste(time[[name_K]][[name_N]][[i]][["SpatialPCA"]], "sec elapsed"))
      
    }
    cat("\n")
  }
}

cat("\n")

# |||||||||||||||||
# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

results.directory <- "results/fPCA_vs_fPCACS/"
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

save(# Parameters
  N_vect, K_vect, N.batches, N.tot, NSR.X, nComp, lambdas_in,
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
# results.directory <- "results/fPCA_vs_fPCACS/"
# list.files(results.directory)
# load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))

# Models
models <- c("fPCA",
            "fPCACS_ML", "fPCACS_ML_IT", 
            "PCA") # ,
            #"SpatialPCA")

# Colors
m_colors <- c(brewer.pal(3, "Reds")[3],
              brewer.pal(3, "Greens")[3:2],
              brewer.pal(3, "Purples")[3]) # ,
              # brewer.pal(3, "Blues")[3])

# Directory
images.directory <- "images/fPCA_vs_fPCACS/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}


## Qualitative Results ----
## ||||||||||||||||||||||||

cat.subsection_title("Qualitative Results")

images.directory_path <- paste(images.directory, "qualitative_results/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

for(K in K_vect){
  name_K <- paste("K", K, sep = "")
  for(N in N_vect){
    name_N <- paste("N", N, sep = "")
    for(m in models){

      plot <- plot.results(data[[name_K]][[name_N]][[1]],
                           results[[name_K]][[name_N]][[1]][[m]])

      ggsave(paste(images.directory_path, name_K,"_", name_N, "_", m, ".jpg", sep = ""),
             plot = plot, width = 10, height = 11, dpi = 200)

    }
  }
}


## Errors by components ----
## |||||||||||||||||||||||||

cat.subsection_title("Errors by components")

images.directory_path <- paste(images.directory, "by_components/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

for(K in K_vect){
  name_K <- paste("K", K, sep = "") 
  for(N in N_vect){
    name_N <- paste("N", N, sep = "") 
      
      plot <- plot.results_by_components(errors_by_components(name_K, name_N),
                                         reorganize_times(name_K, name_N))
      ggsave(paste(images.directory_path, name_K,"_", name_N, ".jpg", sep = ""),
             plot = plot, width = 10, height = 11, dpi = 200)

  }
}


## Time complexity analysis ----
## |||||||||||||||||||||||||||||

cat.subsection_title("Time complexity analysis")

images.directory_path <- paste(images.directory, "time_analysis/", sep = "")
if (!file.exists(images.directory_path)){
  dir.create(images.directory_path)
}

# Time analysis in N

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

plot <- plot.N_time_analysis(times, LOG = TRUE)
ggsave(paste(images.directory_path, "N_log.jpg", sep = ""),
       plot = plot, width = 10, height = 8, dpi = 200)

plot <- plot.N_time_analysis(times, LOG = FALSE)
ggsave(paste(images.directory_path, "N.jpg", sep = ""),
       plot = plot, width = 10, height = 8, dpi = 200)


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

plot.accuracy(temp <- errors_by_K())


