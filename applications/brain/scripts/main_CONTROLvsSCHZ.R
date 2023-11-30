# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Application Brain: Analysis %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()
def.par = par()

setwd("~/master_thesis/applications/brain")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Application Brain: Control vs Schz")

VERBOSE <- TRUE

# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

# fdaPDE
library(fdaPDE)
library(fdaPDE2)

# Clustering
library(caret)
library(e1071)

# Visualization
library(grid)
library(gridExtra)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Linear algebra
library(pracma)

# Statistics utilities
library(pROC)

# Parallelization
library(doParallel)
library(foreach)

# Time
library(tictoc)

# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

source("scripts/functions/fdaPDE.R")
source("scripts/functions/plots.R")
source("scripts/functions/linear_algebra.R")

load("scripts/functions/PLS.RData")


# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

# Directories
directory.images <- "images/"
directory.results <- "results/"

# Code flow control
RUN <- list()
RUN[["Mean Estimation - ColMean"]] <- FALSE
RUN[["Mean Estimation - FR-PDE"]] <- FALSE
RUN[["MV-PLS"]] <- FALSE
RUN[["fPLS"]] <- FALSE
RUN[["fPCA"]] <- FALSE
RUN[["Cl - MV-PLS"]] <- FALSE
RUN[["Cl - fPLS"]] <- FALSE
RUN[["Cl - MV-PCA"]] <- FALSE
RUN[["Cl - fPCA"]] <- FALSE
RUN[["PA - MV-PLS"]] <- TRUE
RUN[["PA - fPLS"]] <- TRUE
RUN[["PA - MV-PCA"]] <- TRUE
RUN[["PA - fPCA"]] <- TRUE

# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")


## |||||||||||||||||||||||||||||
## Importing available data ----
## |||||||||||||||||||||||||||||

# Mesh
FEMbasis <- readRDS("data/FEMbasis.Rds")

# Brain connectivity maps
X <- readRDS("data/X_imputed.rds")

# Responses
Ycl <- factor(readRDS("data/Ycl.rds"))
levels(Ycl)

# Names
names.groups <- c("Control", "Schz")

# Colors
colors <- c('#0C4B8E', '#BF382A')
names(colors) <- levels(Ycl)

# Y01 <- readRDS("data/Y01.rds")
# Yp <- readRDS("data/Yp.rds")


## ||||||||||||||||||||||
## Support variables ----
## ||||||||||||||||||||||

# Indexes groups
index.CONTROL <- Ycl == "CONTROL"
index.SCHZ <- Ycl == "SCHZ"

# Check
any(!(index.CONTROL | index.SCHZ))
any((index.CONTROL & index.SCHZ))

# Response re-scaling
p1 <- sum(index.CONTROL)/length(Ycl)
p2 <- sum(index.SCHZ)/length(Ycl)
Yp <- rep(0, length(Ycl))
Yp[index.CONTROL] <- sqrt(p2/p1)
Yp[index.SCHZ] <- -sqrt(p1/p2)

# Ycl: CONTROL                      vs    SCHZ
# Yp:  -sqrt(P{CONTROL}/P{SCHZ})    vs    sqrt(P{SCHZ}/P{CONTROL})


## ||||||||||||||
## Mesh data ----
## ||||||||||||||

mesh <- FEMbasis$mesh


# ||||||||||||||||||||
# Mean estimation ----
# ||||||||||||||||||||

cat.section_title("Mean Estimation")


## |||||||||||||||||||||
## Column-wise mean ----
## |||||||||||||||||||||

cat.subsection_title("Column-wise mean")

name <- "ColMean"
method <- ColMean

if(RUN[["Mean Estimation - ColMean"]]){
  source("scripts/sections/mean_estimation.R")
}


## ||||||||||||||||||||
## Functional mean ----
## ||||||||||||||||||||

cat.subsection_title("Functional mean")

name <- "FR-PDE"
method <- FRPDE

lambdas <- 10^c(-2, 1, by = 0.2)
  
if(RUN[["Mean Estimation - FR-PDE"]]){
  source("scripts/sections/mean_estimation.R")
}


# |||||||||||||
# Analysis ----
# |||||||||||||

cat.section_title("Analysis")

nComp <- 6


## |||||||||||
## MV-PLS ----
## |||||||||||

cat.subsection_title("MV-PLS")

if(RUN[["MV-PLS"]]){
  results_PLS <- MV_PLS(Yp, X, nComp, VERBOSE)
  save(results_PLS, file = paste(directory.results, "results_PLS.RData", sep = ""))
}

load(paste(directory.results, "results_PLS.RData", sep = ""))


## |||||||||
## fPLS ----
## |||||||||

cat.subsection_title("fPLS")

lambdas <- 10^seq(-2, -1.5, by = 0.1)[3:5]

if(RUN[["fPLS"]]){
  results_fPLS <- fPLS(Yp, X, nComp, mesh, lambdas, VERBOSE)
  save(results_fPLS, file = paste(directory.results, "results_fPLS.RData", sep = ""))
}

load(paste(directory.results, "results_fPLS.RData", sep = ""))


## |||||||||
## fPCA ----
## |||||||||

cat.subsection_title("fPCA")

lambdas <- 10^-1 # seq(-2, -1.5, by = 0.1)[3:5]

if(RUN[["fPCA"]]){
  results_fPCA <- fPCA(Yp, X, nComp, mesh, lambdas, VERBOSE)
  save(results_fPCA, file = paste(directory.results, "results_fPCA.RData", sep = ""))
}

load(paste(directory.results, "results_fPCA.RData", sep = ""))


# ||||||||||||||||||||||||||
# Latent space analysis ----
# ||||||||||||||||||||||||||

cat.section_title("Latent space analysis")

## |||||||||||
## MV-PLS ----
## |||||||||||

cat.subsection_title("MV-PLS")

name <- "MV-PLS"
results <- results_PLS

# Exlorative analysis
# pairs(results$T_hat,
#       main = paste(name, "latent space"),
#       col = colors[as.numeric(index.SCHZ)+1])
# plot3d(results$T_hat[,1:3], col = colors[as.numeric(index.SCHZ)+1], size = 10)

# Latent space analysis
if(RUN[["Cl - MV-PLS"]]){
  source("scripts/sections/latent_space_analysis.R")
}

## |||||||||
## fPLS ----
## |||||||||

cat.subsection_title("fPLS")

name <- "fPLS"
results <- results_fPLS

# Exlorative analysis
# pairs(results$T_hat,
#       main = paste(name, "latent space"),
#       col = colors[as.numeric(index.SCHZ)+1])
# plot3d(results$T_hat[,1:3], col = colors[as.numeric(index.SCHZ)+1], size = 10)

# Latent space analysis
if(RUN[["Cl - fPLS"]]){
  source("scripts/sections/latent_space_analysis.R")
}

## |||||||||
## fPCA ----
## |||||||||

cat.subsection_title("fPCA")

name <- "fPCA"
results <- results_fPCA

# Exlorative analysis
# pairs(results_fPCA$T_hat,
#       main = paste(name, "latent space"),
#       col = colors[as.numeric(index.SCHZ)+1])
# plot3d(results_fPCA$T_hat[,1:3], col = colors[as.numeric(index.SCHZ)+1], size = 10)

# Latent space analysis
if(RUN[["Cl - fPCA"]]){
  source("scripts/sections/latent_space_analysis.R")
}


# ||||||||||||||||||||||||||||
# Performances assessment ----
# ||||||||||||||||||||||||||||

cat.section_title("Performances assessment")

nFolds <- 10 #nrow(X)
nComp <- 6

EVALUATE_ONLY <- TRUE

## |||||||||||
## MV-PLS ----
## |||||||||||

cat.subsection_title("MV-PLS")

name <- "MV-PLS"
method <- MV_PLS

if(RUN[["PA - MV-PLS"]]){
  source("scripts/sections/performances_assessment.R")
}


## |||||||||
## fPLS ----
## |||||||||

cat.subsection_title("fPLS")

name <- "fPLS"
method <- fPLS

lambdas <- 10^seq(-2, -1.5, by = 0.1)[3:5]

if(RUN[["PA - fPLS"]]){
  source("scripts/sections/performances_assessment.R")
}

## |||||||||||
## MV-PCA ----
## |||||||||||

cat.subsection_title("MV-PCA")

name <- "MV-PCA"
method <- MV_PCA

if(RUN[["PA - MV-PCA"]]){
  source("scripts/sections/performances_assessment.R")
}


## |||||||||
## fPCA ----
## |||||||||

cat.subsection_title("fPCA")

name <- "fPCA"
method <- fPCA

lambdas <- 10^(-1)

if(RUN[["PA - fPCA"]]){
  source("scripts/sections/performances_assessment.R")
}


## ||||||||||||||||||||||||||||
## Post-processing results ----
## ||||||||||||||||||||||||||||

models <- c("MV-PCA", "fPCA", "MV-PLS", "fPLS")
m_colors <- c(brewer.pal(3, "Greens")[2:3],
              brewer.pal(3, "Blues")[2:3])
m_names <- c("MV-PCA reg.", "fPCA reg.", "MV-PLS", "fPLS")


plot <- plot.performance_indexes(directory.results)
ggsave(paste(directory.images, "CV_comparison.pdf", sep = ""),
       plot = plot, width = 16, height = 21, dpi = "print", unit = "cm")


# for(name in models) {
#   
#   plots <- list()
#   for(i in 1:8){
#     
#     load(paste(directory.results, name, "/CV/", "AUC_CV_fit_", i, "_", name, ".RData", sep = ""))
#     
#     plots[[(i-1)*2+1]] <- plot.separability(fit$Y_hat, Ycl[-folds[[i]]])
#     
#     # # Prediction on training-set
#     # Y_hat_train <- list()
#     # for(h in 1:nComp){
#     #   Y_hat_train[[h]] <- X[-folds[[i]],] %*% fit$B_hat[[h]]
#     # }
#     # plots[[(i-1)*3+2]] <- plot.separability(Y_hat_train, Ycl[-folds[[i]]])
#     
#     # Prediction on test-set
#     Y_hat_test <- list()
#     for(h in 1:nComp){
#       Y_hat_test[[h]] <- X[folds[[i]],] %*% fit$B_hat[[h]]
#     }
#     plots[[(i-1)*2+2]] <- plot.separability(Y_hat_test, Ycl[folds[[i]]])
#     
#   }
#   
#   plot <- arrangeGrob(grobs = plots, ncol = 2)
#   ggsave(paste(directory.images, "separability_CV_", name, ".pdf", sep = ""),
#          plot = plot, width = 12*2, height = 16*10, dpi = "print", unit = "cm",
#          limitsize = FALSE)
# }
