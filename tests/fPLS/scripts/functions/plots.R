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

plot.fields <- function(nodes, fields, ncol, titles = NULL, type = "tile",
                        limits = NULL, colormap = "D") {
  
  if(is.null(titles)){
    titles <- rep("", length(fields))
  }
  
  plots <- list()
  for (i in 1:length(fields)) {
    if(type == "tile"){
      plots[[i]] <- plot.field_tile(nodes, fields[[i]], 
                                    limits = limits, colormap = colormap)
    } else{
      plots[[i]] <- plot.field_points(nodes, fields[[i]],
                                      limits = limits, colormap = colormap)
    }
    plots[[i]] <- plots[[i]] + ggtitle(titles[i]) + fields_plot_settings()
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = ncol)
  return(grid_plot)
}


## Plot results comparison ----
## ||||||||||||||||||||||||||||

plot.groups_boxplots <- function(errors, name, LEGEND) {
  
  groups_labels <- unique(errors$Group)
  
  data <- errors %>% 
    mutate(Group = factor(Group, levels = groups_labels)) %>% 
    pivot_longer(cols = -Group)
  
  plot <- ggplot(data = data, aes(x = Group, y = value, fill = name)) + 
    geom_boxplot() +
    labs(x = "", y = "", title = name) +
    theme_bw() +
    scale_fill_manual(values = m_colors, labels = models) +
    theme(
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid.major.x = element_blank()
    )
  
  if(LEGEND){
    plot <- plot +
      guides(fill=guide_legend())
  } else {
    plot <- plot + guides(fill = "none")
  }
  
  
  if(length(unique(groups_labels)) > 1){
    plot <- plot +
      geom_vline(xintercept=seq(1.5, length(unique(groups_labels))-0.5, 1), lwd=0.2, colour="grey")
  }
  
  return(plot)
  
}

plot.results_comparison <- function(errors, times) {
  p1 <- plot.groups_boxplots(errors$errors_Y, "RMSE Y", FALSE)
  p2 <- plot.groups_boxplots(errors$errors_X, "RMSE X", FALSE)
  p3 <- plot.groups_boxplots(errors$errors_B, "RMSE B", FALSE)
  p4 <- plot.groups_boxplots(times, "Execution times", TRUE)
  grid_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 1)
  return(grid_plot)
}


## plot Y scatterplots ----
## ||||||||||||||||||||||||

plot.Y_scatterplot <- function(Y1, Y2, names) {
  plot <-
    ggplot(data = data.frame(x = Y1, y = Y2), aes(x = x, y = y)) +
    ggtitle(paste(names[1], "vs", names[2])) +
    xlab("") + ylab("") +
    coord_fixed() +
    theme_bw() +
    theme(
      text = element_text(size = 12),
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
  p[[1]] <- plot.Y_scatterplot(Y_clean, Y, c("Y_clean", "Y"))
  for(h in 1:length(Y_hat)){
    p[[1+h]] <- plot.Y_scatterplot(Y_clean, Y_hat[[h]], c("Y_clean", paste("Y_hat (",h,"comp)", sep = "")))
  }
  grid_plot <- arrangeGrob(grobs = p, ncol = 2)
  return(grid_plot)
}
