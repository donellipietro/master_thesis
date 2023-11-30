# Latent space analysis ----
# ||||||||||||||||||||||||||

## ||||||||||||||||
## Directories ----
## ||||||||||||||||

# Images
directory.images_path <- paste(directory.images, name, "/", sep = "")
if (!file.exists(directory.images_path)){
  dir.create(directory.images_path)
}

# Results
directory.results_path <- paste(directory.results, name, "/", sep = "")
if (!file.exists(directory.results_path)){
  dir.create(directory.results_path)
}


## |||||||||||||||||||||||||||||||||
## Optimal number of components ----
## |||||||||||||||||||||||||||||||||

residuals <- compute_residuals(Yp, X, results$Y_hat, results$X_hat)

### Elbow plots ----
### ||||||||||||||||

# X variable
plot <- plot.nComp_selection(residuals$X, nComp, 3, 
                             title = "X reconstruction",
                             subtitle1 = "Percentage data reconstructed",
                             subtitle2 = "Percentage gain")
ggsave(paste(directory.images_path, "elbox_X.pdf", sep = ""),
       plot = plot, width = 8, height = 8, dpi = "print", unit = "cm")

# Y variable
if(!is.null(results$Y_hat)){
  plot <- plot.nComp_selection(residuals$Y, nComp, 3, 
                               title = "Y reconstruction",
                               subtitle1 = "Accuracy",
                               subtitle2 = "Accuracy gain")
  ggsave(paste(directory.images_path, "elbox_Y.pdf", sep = ""),
         plot = plot, width = 8, height = 8, dpi = "print", unit = "cm")
}


### Groups separability ----
### ||||||||||||||||||||||||

if(!is.null(results$Y_hat)){
  plot <- plot.separability(results$Y_hat, Ycl)
  ggsave(paste(directory.images_path, "separability.pdf", sep = ""),
         plot = plot, width = 16, height = 10, dpi = "print", unit = "cm")
}


## |||||||||||||||||
## Latent space ----
## |||||||||||||||||

# Latent scores data
T_hat <- data.frame(results$T_hat[,1:3])
colnames(T_hat) <- c("t1", "t2", "t3")
T_hat$group[index.CONTROL] <- 'Control'
T_hat$group[index.SCHZ] <- 'Schz'
T_hat$group <- as.factor(T_hat$group)


# SVM fit
svm_model <- svm(group ~ t1 + t2 + t3, data = T_hat, kernel = "linear", cost = 1e15, scale = FALSE)

# Decision boundary plane
intercept <- as.numeric(coef(svm_model)[1])
a <- as.numeric(coef(svm_model)[2:4])

# Plot
fig <- plot.latent_space(T_hat, c(intercept, a))
saveWidget(fig, file = paste(directory.images_path, "latent_space.html"))


# Compute transformation
R <- rotation_matrix(a)
translation <- - intercept / (R%*%a)[1]

# Apply transformation
T_hat_transf <- T_hat
T_hat_transf[,1:3] <- as.matrix(T_hat[,1:3]) %*% t(R)
T_hat_transf[,1] <- T_hat_transf[,1] - translation


fig <- plot.latent_space(T_hat_transf, NULL)
saveWidget(fig, file = paste(directory.images_path, "latent_space_rotated.html"))



## |||||||||||||||||||||||
## Effect on the mean ----
## |||||||||||||||||||||||

# X <- T_hat %*% t(C_hat) + x_mean
# X <- (T_hat %*% t(R)) %*% (R %*% t(C_hat)) + x_mean
# X <- T_hat_rotated %*% (R %*% t(C_hat)) + x_mean
# X1 <- t1_hat_rotated %*% t(c1_hat_rotated) + x_mean
# X1 <- (t1_hat_transf + translation*1) %*% t(c1_hat_rotated) + x_mean
# X1 <- t1_hat_transf %*% t(c1_hat_rotated) + (x_mean + translation*1 %*% t(c1_hat_rotated))
# X1 <- s %*% t(f) + x_critical

x_mean <- results$X_mean
s <- T_hat_transf[,1]
f <- as.matrix((results$C_hat[,1:3] %*% t(R))[,1])
x_critical <- x_mean + translation * f

data_plot <- data.frame(Y_hat = s, group = Ycl)

ggplot(data = data_plot) +
  geom_histogram(aes(x = Y_hat, fill = group)) +
  scale_fill_manual(values = colors, labels = names.groups) +
  ggtitle("Reduced latent space") +
  standard_plot_settings() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10, face = "plain"),
        text = element_text(size = 8))

s_control <- mean(s[index.CONTROL])
s_schz <-  mean(s[index.SCHZ])

x_control <- s_control * f + x_critical
x_schz <- s_schz * f + x_critical

range(x_critical)
range(x_control)
range(x_schz)


# Save f_critical
FEMobj <- FEM(f, FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "f.vtu", sep = ""))

# Save x_critical
FEMobj <- FEM(x_critical, FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "x_critical.vtu", sep = ""))

# Save x_control
FEMobj <- FEM(x_control, FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "x_control.vtu", sep = ""))

# Save x_schz
FEMobj <- FEM(x_schz, FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "x_schz.vtu", sep = ""))


