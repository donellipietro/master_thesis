# Plots ----
# ||||||||||

## Settings ----
## |||||||||||||

standard_plot_settings <- function(){
  standard_plot_settings <- theme_bw() +
    theme(text = element_text(size = 12),
          plot.title = element_text(
            color = "black",
            face = "bold",
            size = 14,
            hjust = 0.5,
            vjust = 1),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid = element_line(linewidth = 0.1),
          panel.border = element_rect(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.2))
}


## Utilities ----
## ||||||||||||||

# Extracts the legend (guide-box) from a given ggplot object.
# It operates on a ggplot object and returns the legend component.

g_legend<-function(a.gplot){
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


## Plots ----
## ||||||||||

plot.nComp_selection <- function(residuals_norm, nComp, nComp_opt, add = "") {
  
  residuals_norm <- as.list(residuals_norm)
  names(residuals_norm) <- c("N", as.character(1:nComp))
  
  # Percentage of data explained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          residuals_norm = unlist(residuals_norm))
  plot1 <- ggplot(data = data_plot) +
    geom_hline(yintercept = min(data_plot$residuals_norm), linetype = "dashed", color = "red", size = 0.5) +
    geom_line(aes(x = id, y = residuals_norm, group = 1)) +
    geom_point(aes(x = id, y = residuals_norm), size = 1.5) +
    geom_point(aes(x = id[1+nComp_opt], y = residuals_norm[1+nComp_opt]), size = 1.5, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Percentage explained") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings() +
    theme(plot.title = element_text(size = 10, face = "plain"),
          text = element_text(size = 8))
  
  # Information gained
  data_plot <- data.frame(id = 0:nComp,
                          name = factor(names(residuals_norm), levels = names(residuals_norm)),
                          gain = - unlist(residuals_norm) + c(NaN, unlist(residuals_norm)[1:(nComp)]))
  plot2 <- ggplot(data = data_plot) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    geom_line(aes(x = id, y = gain, group = 1)) +
    geom_point(aes(x = id, y = gain), size = 1.5) +
    geom_point(aes(x = id[1+nComp_opt], y = gain[1+nComp_opt]), size = 1.5, col = "blue") +
    xlab("Number of components") + ylab("") + ggtitle("Percentage gain") +
    scale_x_continuous(breaks = data_plot$id, labels = data_plot$name) +
    standard_plot_settings() +
    theme(plot.title = element_text(size = 10, face = "plain"),
          text = element_text(size = 8))
  
  # Title
  title <- textGrob(paste("Optimal nComp selection", add),
                    gp = gpar(fontsize = 14, fontface = 'bold'))
  
  # Plot
  plot <- arrangeGrob(title, plot1, plot2, nrow = 3, heights = c(1, 4, 4))
  
  return(plot)
}


plot.separability <- function(Y_hat, Ycl) {
  # Figures
  plots <- list()
  for(h in 1:nComp){
    
    # Data
    data_plot <- data.frame(Y_hat = Y_hat[[h]], group = Ycl)
    
    # Histogram
    plots[[h]] <- ggplot(data = data_plot) +
      geom_histogram(aes(x = Y_hat, fill = group)) +
      scale_fill_manual(values = colors, labels = names.groups) +
      ggtitle(paste("nComp =", h)) +
      standard_plot_settings() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(size = 10, face = "plain"),
            text = element_text(size = 8))
    
    # Legend
    legend <- g_legend(plots[[h]])
    plots[[h]] <-  plots[[h]] + guides(fill = "none")
    
  }
  figures <- arrangeGrob(grobs = plots, nclol = 2)
  
  # Title
  title <- textGrob("Groups identification varying nComp",
                    gp = gpar(fontsize = 14, fontface = 'bold'))
  
  # Plot
  plot <- arrangeGrob(title, figures, legend, heights = c(1, 12, 1))
}


plot.latent_space <- function(svm_model, latent_scores) {
  
  # Support vectors
  # index.support_vectors <- svm_model$index
  # support_vectors <- latent_scores[index.support_vectors, ]
  # support_vectors$group <- "Support vectors"
  # data <- rbind(latent_scores, support_vectors)
  
  # Decision boundary coefficients
  intercept <- coef(svm_model)[1]
  a <- coef(svm_model)[2:4]
  
  # Mesh-grid for the decision boundary plane
  t1_range <- seq(min(latent_scores$t1), max(latent_scores$t1), length.out = 20)
  t2_range <- seq(min(latent_scores$t2), max(latent_scores$t2), length.out = 20)
  meshgrid <- expand.grid(t1 = t1_range, t2 = t2_range)
  meshgrid$t3 <-  -(a[1] * meshgrid$t1 + a[2] * meshgrid$t2 + intercept) / a[3]
  
  
  # Scatter-plot latent scores
  fig <- plot_ly(latent_scores,
                 x = ~t1, y = ~t2, z = ~t3,
                 color = ~group, 
                 colors = c('#0C4B8E', '#BF382A'), #, "black"),
                 marker = list(size = 5))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 't1'),
                                     yaxis = list(title = 't2'),
                                     zaxis = list(title = 't3')),
                        title = "Latent space")
  
  # Plot the decision boundary plane
  fig <- fig %>% add_trace(data = meshgrid,
                           x = ~t1, y = ~t2, z = ~t3,
                           type = "mesh3d",
                           opacity = 0.5,
                           name = 'Boundary plane',
                           inherit = FALSE)
  return(fig)
}

plot.ROC <- function(response, predictors, nComp) {
  
  roc_data <- NULL
  optimal_points <- NULL
  for(h in 1:nComp){
    roc_object <- roc(response, predictors)
    
    # Compute FPR and TPR
    roc_data <- rbind(roc_data,
                      data.frame(
                        group = paste(h, "comp."),
                        FPR = rev(1-roc_object$specificities),
                        TPR = rev(roc_object$sensitivities)))
    
    # Compute optimal treshold
    optimal_threshold <- coords(roc_object, "best")
    optimal_points <-  rbind(optimal_points,
                             data.frame(
                               group = paste(h, "comp."),
                               FPR = 1-optimal_threshold$specificity,
                               TPR = optimal_threshold$sensitivity, 
                               treshold = sprintf("%1.2f", optimal_threshold$threshold)))
    
    
  }
  
  # Plot ROC curve
  plot <- ggplot(roc_data, aes(x = FPR, y = TPR, group = group, color = group)) +
    geom_line(linewidth = 1) +
    geom_point(data = optimal_points, aes(x = FPR, y = TPR)) +
    geom_text(data = optimal_points, aes(x = FPR, y = TPR, label = treshold), vjust = -0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "ROC Curve", x = "FPR (1 - Specificity)", y = "TPR (Sensitivity)") +
    standard_plot_settings() +
    theme(text = element_text(size = 10))
  
  return(plot)
}





plot.groups_boxplots <- function(data, title = "", xlabel = "", LEGEND = TRUE) {
  
  groups_labels <- unique(data$Group)
  
  names(m_colors) <- models
  
  data <- data %>% 
    mutate(Group = factor(Group, levels = groups_labels)) %>% 
    pivot_longer(cols = -Group) %>%
    mutate(name = factor(name, levels = models))
  
  plot <- ggplot(data = data, aes(x = Group, y = value, fill = name)) + 
    geom_boxplot(linewidth = 0.2, outlier.size = 0.5 ) +
    labs(x = xlabel, y = "", title = title) +
    theme_bw() +
    scale_fill_manual(values = m_colors, labels = m_names) +
    theme(
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      text = element_text(size=8),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.spacing.x = unit(0, 'cm'),
      legend.margin = margin(l = 1, unit = "cm"),
      legend.text = element_text(margin = margin(l = 0.2, r = 0.3, unit = "cm")),
      panel.grid = element_line(linewidth = 0.1),
      panel.border = element_rect(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.2)
    )
  
  if(LEGEND){
    plot <- plot +
      guides(fill=guide_legend(nrow=1,
                               byrow=TRUE))
  } else{
    plot <- plot +
      guides(fill="none")
  }
  
  
  if(length(unique(groups_labels)) > 1){
    plot <- plot +
      geom_vline(xintercept=seq(1.5, length(unique(groups_labels))-0.5, 1),
                 lwd=0.2, colour="grey")
  }
  
  return(plot)
}

