# Plots ----
# ||||||||||

## Settings ----
## |||||||||||||

fields_plot_settings <- function(){
  theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          text = element_text(size = 12),
          plot.title = element_text(
            color = "black",
            face = "bold",
            size = 14,
            hjust = 0.5,
            vjust = 1))
} 


## Plot fields ----
## ||||||||||||||||

# Plots a list of fields in a grid plot

# plot.fields <- function(nodes, fields, ncol, titles = NULL, type = "tile",
#                         limits = NULL, colormap = "D") {
#   
#   if(is.null(titles)){
#     titles <- rep("", length(fields))
#   }
#   
#   plots <- list()
#   for (i in 1:length(fields)) {
#     if(type == "tile"){
#       plots[[i]] <- plot.field_tile(nodes, fields[[i]], 
#                                     limits = limits, colormap = colormap)
#     } else{
#       plots[[i]] <- plot.field_points(nodes, fields[[i]],
#                                       limits = limits, colormap = colormap)
#     }
#     plots[[i]] <- plots[[i]] + ggtitle(titles[i]) + fields_plot_settings()
#   }
#   
#   grid_plot <- arrangeGrob(grobs = plots, ncol = ncol)
#   return(grid_plot)
# }


plot.fields <- function(fields, nodes, titles = "", ncol = 1,
                        legend = TRUE, size = 14) {
  
  # Settings
  plot_settings <- ggplot() +
    scale_fill_viridis() +
    theme_minimal() +
    coord_fixed() +
    labs(fill = "") + xlab("") + ylab("") +
    theme(
      text = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        color = "black", face = "bold", size = size,
        hjust = 0.5, vjust = 1
      ),
      legend.position = "bottom"
    )
  
  plots <- list()
  for (i in 1:length(fields)) {
    plots[[i]] <- plot_settings +
      geom_raster(data = data.frame(x = nodes[, 1],
                                    y = nodes[, 2],
                                    z = fields[[i]]), 
                  aes(x = x, y = y, fill = z))
    if(!legend){
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
    if(!is.logical(titles[[i]])){
      plots[[i]] <-  plots[[i]] + ggtitle(titles[[i]])
    }
    
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = ncol)
  return(grid_plot)
}


## Plot results comparison ----
## ||||||||||||||||||||||||||||

# plot.groups_boxplots <- function(errors, name,
#                                  LEGEND = TRUE, OUTLIERS = TRUE) {
#   
#   groups_labels <- unique(errors$Group)
#   
#   names(m_colors) <- models
#   
#   data <- errors %>% 
#     mutate(Group = factor(Group, levels = groups_labels)) %>% 
#     pivot_longer(cols = -Group) %>%
#     mutate(name = factor(name, levels = models))
#   
#   if(!OUTLIERS){
#     outlier.shape <- NA
#   }
#   else{
#     outlier.shape <- 19
#   }
#   
#   plot <- ggplot(data = data, aes(x = Group, y = value, fill = name)) + 
#     geom_boxplot(outlier.shape = outlier.shape) +
#     labs(x = "", y = "", title = name) +
#     theme_bw() +
#     scale_fill_manual(values = m_colors) +
#     theme(
#       plot.title = element_text(
#         color = "black",
#         face = "bold",
#         size = 14,
#         hjust = 0.5,
#         vjust = 1
#       ),
#       legend.title = element_blank(),
#       legend.position = "right",
#       panel.grid.major.x = element_blank()
#     )
#   
#   if(!OUTLIERS){
#     y_max <- 0
#     for(group in groups_labels){
#       for(model in models){
#         y_max <- max(y_max, quantile(data[data$Group == group & 
#                                           data$name == model,]$value,
#                                      na.rm = TRUE,
#                                      c(0.75)),
#                      na.rm = TRUE)
#       }
#     }
#     plot <- plot +
#      coord_cartesian(ylim = c(0, y_max*1.5))
#   }
#   
#   if(LEGEND){
#     plot <- plot +
#       guides(fill=guide_legend())
#   } else {
#     plot <- plot + guides(fill = "none")
#   }
#   
#   
#   if(length(unique(groups_labels)) > 1){
#     plot <- plot +
#       geom_vline(xintercept=seq(1.5, length(unique(groups_labels))-0.5, 1), lwd=0.2, colour="grey")
#   }
#   
#   return(plot)
#   
# }

plot.groups_boxplots <- function(data, title = "", xlabel = "",
                                 LEGEND = TRUE, OUTLIERS = TRUE) {
  
  groups_labels <- unique(data$Group)
  
  names(m_colors) <- models
  
  data <- data %>% 
    mutate(Group = factor(Group, levels = groups_labels)) %>% 
    pivot_longer(cols = -Group) %>%
    mutate(name = factor(name, levels = models))
  
  if(!OUTLIERS){
    outlier.shape <- NA
  }
  else{
    outlier.shape <- 19
  }
  
  plot <- ggplot(data = data, aes(x = Group, y = value, fill = name)) + 
    geom_boxplot(linewidth = 0.2, outlier.size = 0.5) + # outlier.shape = outlier.shape, 
    labs(x = xlabel, y = "", title = title) +
    theme_bw() +
    scale_fill_manual(values = m_colors, labels = m_names) +
    theme(
      plot.title = element_text(
        color = "black",
        # face = "bold",
        size = 10,
        hjust = 0.5,
        vjust = 1
      ),
      text = element_text(size=8),
      legend.title = element_blank(),
      legend.text.align = 0,
      # legend.position = "bottom",
      # legend.spacing.x = unit(0, 'cm'),
      legend.margin = margin(l = 1.5, r = 1.5, unit = "cm"),
      # legend.text = element_text(margin = margin(l = 0, r = 0.3, unit = "cm")),
      panel.grid = element_line(linewidth = 0.1),
      panel.border = element_rect(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.2)
    )
  
  if(!OUTLIERS){
    y_max <- 0
    for(group in groups_labels){
      for(model in models){
        y_max <- max(y_max, quantile(data[data$Group == group &
                                            data$name == model,]$value,
                                     na.rm = TRUE,
                                     c(0.75)),
                     na.rm = TRUE)
      }
    }
    plot <- plot +
      coord_cartesian(ylim = c(0, y_max*1.5))
  }
  
  if(LEGEND){
    plot <- plot # + guides(fill=guide_legend(nrow=1, byrow=TRUE))
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

plot.results_comparison <- function(errors, times, OUTLIERS = TRUE) {
  titles <- c(expression(RMSE[" Y"]), 
              expression(RMSE[" X"]),
              expression(RMSE[" B"]), 
              expression("Execution times"))
  
  p1 <- plot.groups_boxplots(errors$errors_Y, titles[[1]], "", FALSE, OUTLIERS)
  p2 <- plot.groups_boxplots(errors$errors_X, titles[[2]], "", FALSE, OUTLIERS)
  p3 <- plot.groups_boxplots(errors$errors_B, titles[[3]], "", FALSE, OUTLIERS)
  p4 <- plot.groups_boxplots(times, titles[[4]], "", TRUE, OUTLIERS)
  grid_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 1)
  return(grid_plot)
}


## plot Y scatterplots ----
## ||||||||||||||||||||||||

plot.Y_scatterplot <- function(Y1, Y2, title) {
  plot <-
    ggplot(data = data.frame(x = Y1, y = Y2), aes(x = x, y = y)) +
    ggtitle(title) +
    xlab("") + ylab("") +
    coord_fixed() +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      axis.text = element_blank(),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      )
    ) +
    geom_abline(intercept = 0,
                slope = 1,
                color = "red") +
    geom_point(aes(x = x, y = y), color = "black")
  
  return(plot)
}

plot.Y_scatterplots <- function(Y_clean, Y, Y_hat) {
  p <- list()
  p[[1]] <- plot.Y_scatterplot(Y_clean, Y, "Noisy data")
  for(h in 1:length(Y_hat)){
    p[[1+h]] <- plot.Y_scatterplot(Y_clean, Y_hat[[h]], paste(h,"comp."))
  }
  grid_plot <- arrangeGrob(grobs = p, ncol = 1)
  return(grid_plot)
}
