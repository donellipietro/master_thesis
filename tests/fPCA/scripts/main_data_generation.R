# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPCA: Data generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

setwd("~/master_thesis/tests/fPCA")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPCA: Data generation")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Statistical utilities
library(MASS)

# Statistical methods
library(fdaPDE)

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
source("scripts/functions/plots.R")

# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Number of statistical units
N <- c(9)

# Mesh dimension
num_grid_axes <- c(40)
K <- num_grid_axes^2

# Data options
NSR.X <- 0.5^2

# Model options
nComp <- 3


# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")

set.seed(0)
data <- generate_2d_data(sqrt(K), sqrt(K),
                               N = N,
                               H = nComp,
                               NSR.X = NSR.X)
X <- data$X_center
F_true <-  data$F_true

# FEM data
FEM_basis <- data$basisobj
mesh <- data$mesh
nodes <- mesh$nodes
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)


# |||||||||
# Plot ----
# |||||||||

cat.section_title("Plot")

# Directory
images.directory <- "images/data_generation/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}


## Components ----
## |||||||||||||||


fields <- list()
titles <- list()
for(h in 1:nComp){
  fields[[h]] <- F_true[,h]
  titles[[h]] <- bquote(italic(f)[.(h)])
}
plot <- plot.fields(fields, nodes, titles, ncol = 3, FALSE)
ggsave(paste(images.directory, "true_components.pdf", sep = ""),
       plot = plot, width = 16, height = 5, dpi = "print", unit = "cm")


## Samples ----
## ||||||||||||


fields <- list()
titles <- list()
for(i in 1:9){
  fields[[i]] <- X[i,]
  titles[[i]] <- bquote(italic(x)[.(i)])
}
plot <- plot.fields(fields, nodes, titles, ncol = 3, FALSE)
ggsave(paste(images.directory, "sample_examples.pdf", sep = ""),
       plot = plot, width = 16, height = 15, dpi = "print", unit = "cm")







