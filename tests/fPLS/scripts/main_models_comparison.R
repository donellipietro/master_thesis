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

source("scripts/functions/wrappers.R")
source("scripts/functions/errors_and_times.R")
source("scripts/functions/plots.R")
source("scripts/functions/generate_2d_data.R")


# |||||||||||||||
# Parameters ----
# |||||||||||||||

cat.section_title("Parameters")

# Mesh options
num_grid_y_axes <- num_grid_x_axes <- 20

# Number of samples
N <- 60

# Data options
NSR.X <- 0.5^2  # 1/Signal to Noise Ratio X
NSR.Y <- 0.5^2  # 1/Signal to Noise Ratio Y
H <- 4          # Number of components
ORTH <- FALSE   # Type of components
L <- 2          # Number of responses
B.index <- 5

# Model hyper-parameters
lambdas_in <- 10^c(seq(-9, 1, by = 0.2), 4)
nComp <- 5

# Batches
N.batches <- 100

# Number of TPS basis (for B_fPLS and PB_fPLS)
n_basis_tps <- 10

# Number K of folds to do cross-validation
k_folds <- 5


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

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
  data_batch <- generate_2d_data(num_grid_y_axes, num_grid_x_axes,
                                 N = N,
                                 H = H,
                                 L = L,
                                 B.index = B.index,
                                 NSR.X = NSR.X,
                                 NSR.Y = NSR.Y,
                                 ORTH = ORTH)
  X_batch <- data_batch$X
  Y_batch <- data_batch$Y
  X_mean <- data_batch$X_mean
  Y_mean <- data_batch$Y_mean
  X_clean_batch <- data_batch$X_clean
  Y_clean_batch <- data_batch$Y_clean
  B_true <- data_batch$B_true
  
  # FEM data
  FEM_basis <- data_batch$basisobj
  mesh <- data_batch$mesh
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
  data[[i]] <- data_batch
  
  # PLS
  cat("\n - PLS: ")
  tic(quiet = TRUE)
  results[[i]][["PLS"]] <- wrapper_PLS()
  elapsed <- toc(quiet = TRUE)
  time[[i]][["PLS"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["PLS"]], "sec elapsed"))
  
  # fPLS_cpp
  cat("\n - fPLS_cpp: ")
  tic(quiet = TRUE)
  results[[i]][["fPLS_cpp"]] <- wrapper_fPLS_cpp()
  elapsed <- toc(quiet = TRUE)
  time[[i]][["fPLS_cpp"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPLS_cpp"]], "sec elapsed"))
  
  # Save mean results for the fPCA regression
  # (Data centering is not implemented yet in fPCA)
  Y_mean_fPLS <- results[[i]][["fPLS_cpp"]]$Y_mean
  X_mean_fPLS <- results[[i]][["fPLS_cpp"]]$X_mean
  
  # fPCA regression
  cat("\n - fPCA regression: ")
  tic(quiet = TRUE)
  results[[i]][["fPCA_regression"]] <- wrapper_fPCA_regression()
  elapsed <- toc(quiet = TRUE)
  time[[i]][["fPCA_regression"]] <- elapsed$toc - elapsed$tic
  cat(paste(time[[i]][["fPCA_regression"]], "sec elapsed"))
  
  # # fPLS_R_unique
  # cat("\n - fPLS_R_unique: ")
  # tic(quiet = TRUE)
  # results[[i]][["fPLS_R_unique"]] <- wrapper_fPLS_R_unique()
  # elapsed <- toc(quiet = TRUE)
  # time[[i]][["fPLS_R_unique"]] <- elapsed$toc - elapsed$tic
  # cat(paste(time[[i]][["fPLS_R_unique"]], "sec elapsed"))
  
  # # fPLS_R_seq
  # cat("\n - fPLS_R_seq: ")
  # tic(quiet = TRUE)
  # results[[i]][["fPLS_R_seq"]] <- wrapper_fPLS_R_seq()
  # elapsed <- toc(quiet = TRUE)
  # time[[i]][["fPLS_R_seq"]] <- elapsed$toc - elapsed$tic
  # cat(paste(time[[i]][["fPLS_R_seq"]], "sec elapsed"))
  
  # # B-fPLS
  # cat("\n - B-fPLS: ")
  # tic(quiet = TRUE)
  # results[[i]][["B_fPLS"]] <- wrapper_B_fPLS()
  # elapsed <- toc(quiet = TRUE)
  # time[[i]][["B_fPLS"]] <- elapsed$toc - elapsed$tic
  # cat(paste(time[[i]][["B_fPLS"]], "sec elapsed"))
  
  # # PB-fPLS
  # cat("\n - PB-fPLS: ")
  # tic(quiet = TRUE)
  # results[[i]][["PB_fPLS"]] <- wrapper_PB_fPLS()
  # elapsed <- toc(quiet = TRUE)
  # time[[i]][["PB_fPLS"]] <- elapsed$toc - elapsed$tic
  # cat(paste(time[[i]][["PB_fPLS"]], "sec elapsed"))
  
}


## |||||||||||||||||||
## Compute errors ----
## |||||||||||||||||||

for(i in 1:N.batches){
  
  cat(paste("\n\nBatch ", i, ":", sep = ""))
  
  # PLS
  cat("\n - PLS: ")
  errors[[i]][["PLS"]] <- compute_errors(data[[i]], results[[i]][["PLS"]])
  
  # fPLS_cpp
  cat("\n - fPLS_cpp: ")
  errors[[i]][["fPLS_cpp"]] <- compute_errors(data[[i]], results[[i]][["fPLS_cpp"]])
  
  # fPCA regression
  cat("\n - fPCA regression: ")
  errors[[i]][["fPCA_regression"]] <- compute_errors(data[[i]], results[[i]][["fPCA_regression"]])
  
  # Add other methods
}


# Save results ----
# |||||||||||||||||

cat.section_title("Save results")

save(num_grid_y_axes, num_grid_x_axes,
     N, N.batches,
     B.index,
     NSR.X, NSR.Y, H, L,
     nComp, lambdas_in,
     data, results, errors, time,
     file = paste(directory.results, "results",
                  "_B.index", B.index,
                  "_L", L,
                  format(Sys.time(), "_%Y%m%d_%H%M%S"), ".RData", sep = ""))


# |||||||||||||||||
# Plot Results ----
# |||||||||||||||||

cat.section_title("Plot Results")

# Images directory
directory.images <- "images/models_comparison/"
if (!file.exists(directory.images)){
  dir.create(directory.images)
}
directory.images <- paste(directory.images, "B.index_", B.index, "_L_", L, "/", sep = "")
if (!file.exists(directory.images)){
  dir.create(directory.images)
}

# Load data
# list.files(directory.results)
# load(paste(directory.results, tail(list.files(directory.results), n = 1), sep = ""))

# Models
models <- c("fPCA_regression",
            "fPLS_cpp",
            # "fPLS_R_seq", "fPLS_R_unique",
            # "B_fPLS", # "PB_fPLS",
            "PLS"
)

# Models names
m_names <- c(expression(fPCA-reg.), # []
             expression(fPLS), # [cpp]
             # expression(fPLS[R-seq]), expression(fPLS[R-unique]),
             # expression(B-fPLS[]), # expression(PB-fPLS[]),
             expression(MV-PLS) # []
)
names(m_names) <- models

# Names
True <- expression(Truth[])
Noisy <- expression(Noisy[])

# Colors
m_colors <- c(brewer.pal(3, "Greens")[3],
              brewer.pal(3, "Blues")[3],
              # brewer.pal(3, "Greys")[3:2],
              # brewer.pal(3, "Purples")[3], #:2],
              brewer.pal(3, "Reds")[3]
)


## RMSE and Time ----
## ||||||||||||||||||

cat.subsection_title("RMSE and Time")

# Figures
figures <- plot.results_comparison(errors_by_components(), reorganize_times())

# Figures
figures_no_outliers <- plot.results_comparison(errors_by_components(), reorganize_times(),
                                               OUTLIERS = FALSE)


# # Title
# title <- textGrob(expression("Test univariate - 1"),
#                   gp = gpar(fontsize = 14, fontface = 'bold'))
# title <- textGrob(expression("Test univariate - 1"),
#                   gp = gpar(fontsize = 14, fontface = 'bold'))
title <- textGrob(expression("Test bivariate"),
                  gp = gpar(fontsize = 14, fontface = 'bold'))


# Plot
plot_comparison <- arrangeGrob(title, figures, heights = c(1, 19))
ggsave(paste(directory.images, "comparison.pdf", sep = ""),
       plot = plot_comparison, width = 16, height = 21, dpi = "print", units = "cm")

# Plot
plot_comparison <- arrangeGrob(title, figures_no_outliers, heights = c(1, 19))
ggsave(paste(directory.images, "comparison_no_outliers.pdf", sep = ""),
       plot = plot_comparison, width = 16, height = 21, dpi = "print", units = "cm")


## Data ----
## |||||||||
# 
# cat.subsection_title("Original components")
# 
# # Figures
# nrow <- 2
# ncol <- H/2
# fields <- NULL
# titles <- NULL
# for(h in 1:nComp){
#   titles <- c(titles, paste(h, "comp."))
# }
# fields <- split(t(data[[1]]$F_true), seq_len(nrow*ncol))
# figures <- plot.fields(nodes, fields, ncol, titles = titles)
# 
# # Scalar products matrix
# scalar_products <- abs(as.numeric(t(data[[1]]$F_true) %*% R0 %*%  data[[1]]$F_true))
# names <- expand.grid(1:H, 1:H)
# data_plot <- data.frame(names, scalar_products)
# colnames(data_plot) <- c("x", "y", "value")
# sp_matrix <- ggplot(data =  data_plot, mapping = aes(x = x, y = y)) +
#   geom_tile(aes(fill = value), colour = "white") +
#   geom_text(aes(label = sprintf("%.2f",value)), vjust = 1) +
#   ggtitle("Scalar products") +
#   standard_plot_settings() +
#   scale_x_continuous(position = "top") + 
#   scale_y_continuous(position = "left") +
#   scale_y_reverse() +
#   theme(axis.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_blank()) + 
#   coord_fixed() + guides(fill = "none") +
#   scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,1.1))
# 
# # Title
# title <- textGrob(expression(Original ~ components),
#                   gp = gpar(fontsize = 20, fontface = 'bold'))
# 
# # Plot
# plot_OC <- arrangeGrob(figures, sp_matrix, heights = c(2, 1))
# plot_OC <- arrangeGrob(title, plot_OC, heights = c(0.2, 3))
# ggsave(paste(directory.images, "original_components.jpg", sep = ""),
#        plot = plot_OC, width = 21/2, height = 29.7/2, dpi = 200)


## Qualitative results ----
## ||||||||||||||||||||||||

# Here the results for the first batch are considered

nodes <- data[[1]]$mesh$nodes

### X reconstruction ----
### |||||||||||||||||||||

cat.subsection_title("X reconstruction")

indexes <- c(1, 3, 44, 60)

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
title <- textGrob(expression(X[mean]),
                  gp = gpar(fontsize = 14, fontface = 'bold'))

# Figures
nrow <- 1
ncol <- length(models)+1
fields <- list(X_mean)
for(m in models){
  fields <- c(fields, list(results[[1]][[m]]$X_mean))
}
figures <- plot.fields(fields, nodes,
                       titles = c(list(True), as.list(m_names)),
                       ncol, legend = FALSE, size = 10)

# Plot
plot_X_mean <- arrangeGrob(title, figures, heights = c(1, 4))
# ggsave(paste(directory.images_X, "mean.jpg", sep = ""),
#        plot = plot_X_mean, width = 3*ncol, height = 1 + 4*nrow, dpi = 200)
# 


#### X ----
#### ||||||

# # True
# X_clean <- data[[1]]$X_clean
# 
# # Noisy
# X <- data[[1]]$X
# 
# # Figures
# nrow <- 4
# ncol <- length(models) + 2
# fields <- NULL
# titles <- list()
# for(id in indexes){
#   fields <- rbind(fields, X_clean[id,], X[id,])
#   titles <- c(titles, list(True), list(Noisy))
#   for(m in models){
#     fields <- rbind(fields, results[[1]][[m]]$X_hat[[H]][id,])
#     titles <- c(titles, as.list(m_names[m]))
#   }
# }
# fields <- split(fields, seq_len(nrow*ncol))
# figures <- plot.fields(fields, nodes, titles = titles, ncol,
#                        legend = FALSE, size = 10)
# 
# # Title
# title <- textGrob(expression(X ~ reconstruction),
#                   gp = gpar(fontsize = 14, fontface = 'bold'))
# 
# # Plot
# plot_X <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
# # ggsave(paste(directory.images_X, "X.jpg", sep = ""),
# #        plot = plot_X, width = 3*ncol, height = 1 + 4*nrow, dpi = 200,
# #        device = "jpeg", path = NULL)


#### X_centered ----
#### |||||||||||||||

# True
X_c_clean <- data[[1]]$X_clean - rep(1,N) %*% t(data[[1]]$X_mean)

# Noisy
X_c <- data[[1]]$X_center

# Figures
nrow <- 4
ncol <- length(models) + 2
fields <- NULL
titles <- list()
side <- list()
top <- list()
top[[1]] <- NULL
top[[2]] <- textGrob(True,
                       gp = gpar(fontsize = 10, fontface = 'bold'))
top[[3]] <- textGrob(Noisy,
                     gp = gpar(fontsize = 10, fontface = 'bold'))
for(m in models){
  top <- c(top, list(textGrob(m_names[[m]],
                         gp = gpar(fontsize = 10, fontface = 'bold'))))
}
for(id in indexes){
  fields <- cbind(fields, X_c_clean[id,], X_c[id,])
  titles <- c(titles, FALSE, FALSE) #list(True), list(Noisy))
  side[[as.character(id)]] <- textGrob(id, gp = gpar(fontsize = 10))
  for(m in models){
    fields <- cbind(fields, results[[1]][[m]]$X_hat[[H]][id,] - results[[1]][[m]]$X_mean)
    titles <- c(titles, FALSE) #as.list(m_names[m]))
  }
}
fields <- split(t(fields), seq_len(nrow*ncol))
figures <- plot.fields(fields, nodes, titles = titles, ncol,
                       legend = FALSE, size = 10)
side <- arrangeGrob(grobs = side, ncol = 1)
figures <- arrangeGrob(side, figures, ncol = 2, widths = c(0.3, 5))
top <- arrangeGrob(grobs = top, ncol = ncol+1, widths = c(0.3, 1,1,1,1,1))
figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 4*5))

# Title
title <- textGrob(expression(X[centered]),
                  gp = gpar(fontsize = 14, fontface = 'bold'))

# Plot
plot_X_c <- arrangeGrob(title, figures, heights = c(2, 4*nrow))
# ggsave(paste(directory.images_X, "X_centered.jpg", sep = ""),
#        plot = plot_X_c, width = 3*ncol, height = 1 + 4*nrow, dpi = 200,
#        device = "jpeg", path = NULL)

#### Global plot ----
#### ||||||||||||||||

plot_X_global <- arrangeGrob(plot_X_mean,
                             #plot_X,
                             plot_X_c,
                             heights = c(0.25+1, 0.25+4))#, 0.25+4))

ggsave(paste(directory.images, "X.pdf", sep = ""),
       plot = plot_X_global, width = 16, height = 21, dpi = "print", units = "cm")


### Y reconstruction ----
### |||||||||||||||||||||

# cat.subsection_title("Y reconstruction")
# 
# # # Images directory
# # directory.images_Y <- paste(directory.images, "Y/", sep = "")
# # if (!file.exists(directory.images_Y)){
# #   dir.create(directory.images_Y)
# # }
# 
# # True
# Y_clean <- data[[1]]$Y_clean
# 
# # Noisy
# Y <- data[[1]]$Y
# 
# for(l in 1:L){
#   
#   # Figures
#   nrow <- 4
#   ncol <- length(models)
#   plots <- list()
#   top  <- list()
#   Y_hat <- list()
#   for(m in models){
#     for(h in 1:H) Y_hat[[h]] <- as.matrix(results[[1]][[m]]$Y_hat[[h]], ncol = L)[,l]
#     plots[[m]] <- plot.Y_scatterplots(Y_clean[,l], Y[,l], Y_hat)
#     top[[m]] <-  textGrob(models_names[[m]],
#                           gp = gpar(fontsize = 18, fontface = 'bold'))
#   }
#   figures <- arrangeGrob(grobs = plots, ncol = ncol)
#   top <-  arrangeGrob(grobs = top, ncol = ncol)
#   figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))
#   
#   # Title
#   title <- textGrob(expression(Y[l] ~ reconstruction),
#                     gp = gpar(fontsize = 20, fontface = 'bold'))
#   
#   # Plot
#   plot_Y <- arrangeGrob(title, figures, heights = c(1, 4*nrow, 5))
#   ggsave(paste(directory.images, "Y",l,".jpg", sep = ""),
#          plot = plot_Y, width = 21, height = 29.7, dpi = 200)
#   
# }


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
G_true <- data[[1]]$G_true

for(l in 1:L){
  
  # Figures
  nrow <- nComp
  ncol <- length(models)
  fields <- NULL
  titles <- list()
  top <- list(NULL)
  side <- list()
  for(m in models){
    # fields <- rbind(fields, as.numeric(B_true[,l]))
    # titles <- c(titles, list(True))
    top[[m]] <- textGrob(m_names[[m]],
                         gp = gpar(fontsize = 10, fontface = 'bold'))
  }
  top <- c(top, list(textGrob(True,
                         gp = gpar(fontsize = 10, fontface = 'bold'))))
  for(h in 1:nComp){
    side[[h]] <- textGrob(bquote(.(h)~comp.),
                          gp = gpar(fontsize = 10, fontface = 'bold'), just = "center")
    for(m in models){
      fields <- rbind(fields, as.numeric(as.matrix(results[[1]][[m]]$B_hat[[h]], ncol = L)[,l]))
      titles <- c(titles, FALSE) #as.list(paste(h, "comp.")))
    }
  }
  fields <- split(fields, seq_len(nrow*ncol))
  figures <- plot.fields(fields, nodes, titles, ncol, legend = FALSE, size = 10)
  side <- arrangeGrob(grobs = side, ncol = 1)
  figures <- arrangeGrob(side, figures, ncol = 2, widths = c(0.5, 3))
  
  true_fied <- plot.fields(list(B_true[,l]), nodes, list(FALSE), 1, legend = FALSE, size = 10)
  true_fied <- c(as.list(rep(list(NULL), tail(which(G_true[l,]!=0), 1)-1)), list(true_fied))
  true_fied <- arrangeGrob(grobs = true_fied, nrow = nrow)
  figures <- arrangeGrob(figures, true_fied, ncol = 2, widths = c(3.5, 1))
  
  top <- arrangeGrob(grobs = top, ncol = 2 + ncol, widths = c(0.5, 1,1,1, 1))
  figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))
  
  # Title
  title <- textGrob(bquote(beta[.(l)]),
                    gp = gpar(fontsize = 14, fontface = 'bold'))
  
  # Plot
  plot_B <- arrangeGrob(title, figures, heights = c(1.5, 4*nrow))
  ggsave(paste(directory.images, "B", l, ".pdf", sep = ""),
         plot = plot_B, width = 16, height = 21, dpi = "print", units = "cm")
  
}


# ### Directions ----
# ### |||||||||||||||
# 
# cat.subsection_title("Directions")
# 
# # Figures
# nrow <- nComp
# ncol <- length(models) - sum(models == "fPCA_regression")
# fields <- NULL
# titles <- NULL
# top <- list()
# for(m in models){
#   if(m != "fPCA_regression"){
#     top[[m]] <-  textGrob(models_names[[m]],
#                           gp = gpar(fontsize = 18, fontface = 'bold'))
#   }
# }
# for(h in 1:nComp){
#   for(m in models){
#     if(m != "fPCA_regression"){
#       fields <- rbind(fields, as.numeric(results[[1]][[m]]$W_hat[,h]))
#       titles <- c(titles, paste(h, "comp."))
#     }
#   }
# }
# fields <- split(fields, seq_len(nrow*ncol))
# figures <- plot.fields(nodes, fields, ncol, titles = titles)
# top <-  arrangeGrob(grobs = top, ncol = ncol)
# figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))
# 
# # Title
# title <- textGrob(expression(Directions),
#                   gp = gpar(fontsize = 20, fontface = 'bold'))
# 
# # Plot
# plot_W <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
# ggsave(paste(directory.images, "directions.jpg", sep = ""),
#        plot = plot_W, width = 21, height = 29.7, dpi = 200)
# 
# 
# ### Loadings ----
# ### |||||||||||||
# 
# cat.subsection_title("Loadings")
# 
# # Figures
# nrow <- nComp
# ncol <- length(models)
# fields <- NULL
# titles <- NULL
# top <- list()
# for(m in models){
#   top[[m]] <-  textGrob(models_names[[m]],
#                         gp = gpar(fontsize = 18, fontface = 'bold'))
# }
# for(h in 1:nComp){
#   for(m in models){
#     if(m == "fPCA_regression") fields <- rbind(fields, as.numeric(results[[1]][[m]]$F_hat[,h]))
#     else fields <- rbind(fields, as.numeric(results[[1]][[m]]$C_hat[,h]))
#     titles <- c(titles, paste(h, "comp."))
#   }
# }
# fields <- split(fields, seq_len(nrow*ncol))
# figures <- plot.fields(nodes, fields, ncol, titles = titles)
# top <-  arrangeGrob(grobs = top, ncol = ncol)
# figures <- arrangeGrob(top, figures, nrow = 2, heights = c(1, 6*5))
# 
# # Title
# title <- textGrob(expression(Directions),
#                   gp = gpar(fontsize = 20, fontface = 'bold'))
# 
# # Plot
# plot_W <- arrangeGrob(title, figures, heights = c(1, 4*nrow))
# ggsave(paste(directory.images, "loadings.jpg", sep = ""),
#        plot = plot_W, width = 21, height = 29.7, dpi = 200)




