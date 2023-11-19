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

# Time
library(tictoc)

# Parallelization
library(doParallel)
library(foreach)


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
RUN[["Mean Estimation - ColMean"]] <- TRUE
RUN[["Mean Estimation - FR-PDE"]] <- FALSE
RUN[["MV-PLS"]] <- FALSE
RUN[["fPLS"]] <- FALSE
RUN[["fPCA"]] <- FALSE
RUN[["Cl - MV-PLS"]] <- FALSE
RUN[["Cl - fPLS"]] <- FALSE
RUN[["Cl - fPCA"]] <- FALSE
RUN[["PA - MV-PLS"]] <- FALSE
RUN[["PA - fPLS"]] <- FALSE
RUN[["PA - fPCA"]] <- FALSE

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

lambdas <- 10^seq(-2, -1.5, by = 0.1)[3:5]

if(RUN[["fPCA"]]){
  results_fPCA <- fPLS(X, nComp, mesh, lambdas, VERBOSE)
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
  # source("scripts/sections/latent_space_analysis.R")
}


# ||||||||||||||||||||||||||||
# Performances assessment ----
# ||||||||||||||||||||||||||||

cat.section_title("Performances assessment")

nFolds <- 10
nComp <- 6


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


## ||||||||||||||||||||||||||||
## Post-processing results ----
## ||||||||||||||||||||||||||||

models <- c("MV-PLS", "fPLS")
m_colors <- c(brewer.pal(3, "Oranges")[3],
              brewer.pal(3, "Greens")[3])
m_names <- c("MV-PLS", "fPLS")

AUC_results <- data.frame(matrix(nrow = nFolds*nComp, ncol = 0))
for(name in models){
  # Load AUC data
  load(paste(directory.results, "AUC_CV_",name,".RData", sep = ""))
  # Reorganize data
  temp <- data.frame(AUC_comp_fold)
  temp$Group <- paste("Comp.", 1:nComp)
  temp <- data.frame(temp %>% pivot_longer(cols = -Group))
  temp$name <- name
  AUC_results <- cbind(temp, AUC_results)
  rm(AUC_comp_fold)
}
AUC_results <- data.frame(AUC_results$Group, AUC_results[, c(3,6)])
colnames(AUC_results) <- c("Group", models)

# Plot
plot.groups_boxplots(AUC_results, "AUC MV-PLS vs fPLS")



plots <- list()
for(i in 1:8){

  load(paste(directory.results, name, "/CV/", "AUC_CV_fit_", i, "_fPLS.RData", sep = ""))

  plots[[(i-1)*3+1]] <- plot.separability(fit$Y_hat, Ycl[-folds[[i]]])

  # Prediction on training-set
  Y_hat_train <- list()
  for(h in 1:nComp){
    Y_hat_train[[h]] <- X[-folds[[i]],] %*% fit$B_hat[[h]]
  }
  plots[[(i-1)*3+2]] <- plot.separability(Y_hat_train, Ycl[-folds[[i]]])

  # Prediction on test-set
  Y_hat_test <- list()
  for(h in 1:nComp){
    Y_hat_test[[h]] <- X[folds[[i]],] %*% fit$B_hat[[h]]
  }
  plots[[(i-1)*3+3]] <- plot.separability(Y_hat_test, Ycl[folds[[i]]])

}

plot <- arrangeGrob(grobs = plots, ncol = 3)
ggsave(paste(directory.images, "separability_CV.pdf", sep = ""),
       plot = plot, width = 12*3, height = 16*10, dpi = "print", unit = "cm",
       limitsize = FALSE)

