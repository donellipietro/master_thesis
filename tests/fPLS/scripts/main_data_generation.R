# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Test fPLS: Data generation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

setwd("~/master_thesis/tests/fPLS")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Test fPLS: Data generation")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# Statistical utilities
library(MASS)

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

source("scripts/functions/wrappers.R")
source("scripts/functions/plots.R")
source("scripts/functions/generate_2d_data.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Mesh options
num_grid_y_axes <- num_grid_x_axes <- 30

# Number of samples
N <- 60

# Data options
NSR.X <- 0.5^2  # 1/Signal to Noise Ratio X
NSR.Y <- 0.5^2  # 1/Signal to Noise Ratio Y
H <- 4          # Number of components
ORTH <- FALSE   # Type of components
L <- 2          # Number of responses
B.index <- 5


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

# Image folder
images.directory <- "images/data_generation/"
if (!file.exists(images.directory)){
  dir.create(images.directory)
}


# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")


# Set seed for tests reproducibility 
set.seed(0)

# Data generation
data <- generate_2d_data(num_grid_y_axes, num_grid_x_axes,
                               N = N,
                               H = H,
                               L = L,
                               B.index = B.index,
                               NSR.X = NSR.X,
                               NSR.Y = NSR.Y,
                               ORTH = ORTH)
X_batch <- data$X
X_mean <- data$X_mean
B_true <- data$B_true
F_true <- data$F_true

# FEM data
FEM_basis <- data$basisobj
mesh <- data$mesh
nodes <- mesh$nodes


# Plot B
titles <- list(expression(beta[1]), expression(beta[2]))
plot <- plot.fields(c(list(B_true[,1], B_true[,2])), mesh$nodes,
                    titles = titles, 2, legend = FALSE)
ggsave(paste(images.directory, "B_true.pdf", sep = ""),
       plot = plot, width = 8, height = 5.5, dpi = "print", unit = "cm")

grid.arrange(plot)

# Plot F
for(h in 1:H) titles[[h]] <- bquote(italic(f)[.(h)])
fields <- split(t(F_true), 1:H)
plot <- plot.fields(fields, mesh$nodes, titles, 4, legend = FALSE)
ggsave(paste(images.directory, "F_true.pdf", sep = ""),
       plot = plot, width = 16, height = 5.5, dpi = "print", unit = "cm")


# Plot X_mean
titles <- list(expression(X[mean]))
plot <- plot.fields(c(list(X_mean)), mesh$nodes,
                    titles = titles, 1, legend = FALSE)
ggsave(paste(images.directory, "X_mean.pdf", sep = ""),
       plot = plot, width = 4, height = 5.5, dpi = "print", unit = "cm")
