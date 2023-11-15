# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test FRPDE: Noise suppression %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

setwd("~/master_thesis/tests/FRPDE")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test FRPDE: Noise suppression")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Statistical utilities
library(MASS)

# Statistical methods
library(fdaPDE)
library(fdaPDE2)

# Times traking
library(tictoc)

# Visualization
library(dplyr)
library(tidyr)
library(ggplot2)
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
source("scripts/functions/errors.R")
source("scripts/functions/plots.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Number of statistical units
N <- 10

# Mesh dimension
num_grid_axes <- num_grid_x_axes <- num_grid_y_axes <- 40

# Number of batches
N.batches <- 10

# Data options
NSR.X_vect <- c(1/3, 1/2, 1, 2)

# Model options
lambdas <- 10^c(seq(-6, -3, by = 0.2), 4)


# ||||||||||
# Tests ----
# ||||||||||

cat.section_title("Tests")

# Room for results
results <- list()
errors <- list()
data <- list()


for(NSR.X in NSR.X_vect){
  
  NSR.X_name <- paste("NSR.X", sprintf("%2.4f", NSR.X), sep = "")
  cat(paste("\nTests NSR.X = ", NSR.X, sep = ""))
  
  # Room for results
  name_N <- paste("N", N, sep = "") 
  data[[NSR.X_name]] <- list()
  results[[NSR.X_name]] <- list()
  errors[[NSR.X_name]] <- list()
  
  for(i in 1:N.batches){
    
    cat(paste("\n\nBatch ", i, ":", sep = ""))
    
    set.seed(i)
    
    # Data generation
    set.seed(i)
    data_batch <- generate_2d_data(num_grid_axes, num_grid_axes,
                                   N = N,
                                   NSR.X = NSR.X)
    X_batch <- data_batch$X
    
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
    results[[NSR.X_name]][[i]] <- list()
    errors[[NSR.X_name]][[i]] <- list()
    
    # Save data
    data[[NSR.X_name]][[i]] <- data_batch
    
    # Column-wise mean
    cat("\n- ColMeans")
    results[[NSR.X_name]][[i]][["colMeans"]] <- wrapper_colMeans()
    errors[[NSR.X_name]][[i]][["colMeans"]] <- compute_errors(data[[NSR.X_name]][[i]], results[[NSR.X_name]][[i]][["colMeans"]])
    
    # FRPDE
    cat("\n- FR-PDE: ")
    results[[NSR.X_name]][[i]][["FR-PDE"]] <- wrapper_FRPDE()
    errors[[NSR.X_name]][[i]][["FR-PDE"]] <- compute_errors(data[[NSR.X_name]][[i]], results[[NSR.X_name]][[i]][["FR-PDE"]])
    cat(results[[NSR.X_name]][[i]][["FR-PDE"]]$lambda_opt)
  }
  cat("\n")
}

cat("\n")


# |||||||||||||||||
# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

results.directory <- "results/noise_suppression/"
if (!file.exists(results.directory)){
  dir.create(results.directory)
}

save(# Parameters
  NSR.X_vect, N.batches, N,
  # Results
  data, results, errors,
  # Locations
  file = paste(results.directory,
               "results", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData",
               sep = ""))


# |||||||||||||||||
# Plot Results ----
# |||||||||||||||||

cat.section_title("Plot Results")

# Load data
# results.directory <- "results/noise_suppression/"
# list.files(results.directory)
# load(paste(results.directory, tail(list.files(results.directory), n = 1), sep = ""))


# Models
models <- c("colMeans", "FR-PDE")

# Names
m_names <- c(expression("ColMeans"), 
             expression("FR-PDE"))

# Colors
m_colors <- c(brewer.pal(3, "Oranges")[3],
              brewer.pal(3, "Greens")[3])


# Directory
images.directory <- "images/noise_supression/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}


## Qualitative results ----
## ||||||||||||||||||||||||

nodes <- data[[1]][[1]]$mesh$nodes

lim <- NaN
for(i in 1:length(NSR.X_vect)){
  lim <- range(c(lim, range(data[[i]][[1]]$X[1,])), na.rm = TRUE)
}

plots <- list()
for(NSR.X in NSR.X_vect){
  
  NSR.X_name <- paste("NSR.X", sprintf("%2.4f", NSR.X), sep = "")
  
  # Fields
  fields <- list()
  titles <- list()
  # True mean
  fields[[1]] <- data[[NSR.X_name]][[1]]$x_clean
  titles[[1]] <- expression(True~field)
  # Noisy sample 
  titles[[2]] <- expression(Noisy~sample)
  fields[[2]] <- data[[NSR.X_name]][[1]]$X[1,]
  # FR_PDE estimate
  titles[[3]] <- expression("FR-PDE"~estimate)
  fields[[3]] <- results[[NSR.X_name]][[1]][["FR-PDE"]]$x_hat
  # ColMeans estimate
  titles[[4]] <- expression(ColMeans~estimate)
  fields[[4]] <- results[[NSR.X_name]][[1]][["colMeans"]]$x_hat
  grid_plot <- plot.fields(fields, nodes, titles, ncol = 4, lim = lim, legend = FALSE)
  
  # Title
  title <- textGrob(paste("SNR = ", 1/NSR.X), gp = gpar(fontsize = 14, fontface = 'bold'))
  
  # Final plot 
  plots[[NSR.X_name]] <- arrangeGrob(title, grid_plot, heights = c(0.5, 4))
}
plot <- arrangeGrob(grobs = plots, ncol = 1)
ggsave(paste(images.directory, "qualitative_results.pdf", sep = ""),
       plot = plot, width = 16, height = 21, dpi = "print", unit = "cm")


## Quantitative results ----
## |||||||||||||||||||||||||

final_errors <- reorganize_errors()

plot <- plot.groups_boxplots(final_errors, "RMSE")
ggsave(paste(images.directory, "RMSE.pdf", sep = ""),
       plot = plot, width = 16, height = 10, dpi = "print", unit = "cm")



