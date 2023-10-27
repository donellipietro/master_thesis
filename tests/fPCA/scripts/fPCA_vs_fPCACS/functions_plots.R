# Plots ----
# ||||||||||

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


## Fields ----
## |||||||||||

# Creates a grid of raster plots taking as input a list of fields, nodes, 
# and optional titles as input. The function generates individual plots for 
# each field and arranges them in a gridof specified number of columns.

plot.fields <- function(fields, nodes, titles = "", ncol = 1) {
  
  # Settings
  plot_settings <- ggplot() +
    scale_fill_viridis() +
    theme_bw() +
    coord_fixed() +
    labs(fill = "") + xlab("") + ylab("") +
    theme(
      text = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(
        color = "black", face = "bold", size = 14,
        hjust = 0.5, vjust = 1
      ),
      legend.position = "bottom"
    )
  
  plots <- list()
  for (i in 1:length(fields)) {
    plots[[i]] <- plot_settings +
      ggtitle(titles[i]) +
      geom_raster(data = data.frame(x = nodes[, 1],
                                    y = nodes[, 2],
                                    z = fields[[i]]), 
                  aes(x = x, y = y, fill = z))
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = ncol)
  return(grid_plot)
}


## Scatterplots ----
## |||||||||||||||||

# Creates a scatterplot to compare two sets of data, Y1 and Y2, the main 
# diagonal line (Y1 = Y2) is displayed as reference. 

plot.Y_scatterplot <- function(Y1, Y2, names) {
  
  plot <-
    ggplot(data = data.frame(x = Y1, y = Y2), aes(x = x, y = y)) +
    ggtitle(paste(names[1], "vs", names[2])) +
    # xlab(names[1]) + ylab(names[2]) +
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
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    geom_abline(intercept = 0,
                slope = 1,
                color = "red") +
    geom_point(aes(x = x, y = y), color = "black")
  
  return(plot)
}

# Generates a grid of scatterplots by calling theplot.Y_scatterplot function
# for each pair of true and estimated scores (one for each component).

plot.score_scatterplots <- function(S_true, S_hat) {
  
  plots <- list()
  names <- c("s_true", "s_hat")
  
  for(i in 1:nComp){
    plots[[i]] <- plot.Y_scatterplot(S_true[,i], S_hat[,i],
                                     paste(names, i))
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = nComp)
  return(grid_plot)
}


## Plot groups boxplots ----
## |||||||||||||||||||||||||

# Creates grouped boxplots taking as input a data frame data that contains
# a grouping variable Group. If LEGEND is set to TRUE, it includes a legend for
# the data categories, and if set to FALSE, it removes the legend.

plot.groups_boxplots <- function(data, title = "", xlabel = "", LEGEND = true) {
  
  groups_labels <- unique(data$Group)
  
  data <- data %>% 
    mutate(Group = factor(Group, levels = groups_labels)) %>% 
    pivot_longer(cols = -Group)
  
  plot <- ggplot(data = data, aes(x = Group, y = value, fill = name)) + 
    geom_boxplot() +
    labs(x = xlabel, y = "", title = title) +
    theme_bw() +
    scale_fill_manual(values = m_colors) +
    theme(
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      legend.title = element_blank(),
      legend.position = "bottom"
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

# Generates a grid of groups_boxplots for the different values of K.
# Tha argument 'data' is a list in K and each element of the list is a
# dataframe that must contain the 'Component' variable.
# The 'component' argument specifies which component should be displayed. 
# If the data does not contain any components its enought to fill these 
# varaibles with an empty string "".

plot.groups_boxplots_by_K <- function(data, component, title) {
  
  plots <- list()
  
  for(K in K_vect){
    
    name_K <- paste("K", K, sep = "")
    data_plot <- data[[name_K]][data[[name_K]]$Component == component, c(-2)]
    
    if(K == max(K_vect)){
      LEGEND = TRUE
    } else {
      LEGEND = FALSE
    }
    
    plots[[name_K]] <- plot.groups_boxplots(data_plot, 
                                            paste("K = ", K, sep = ""),
                                            "N", LEGEND)
    
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = 1,
                           top = grid::textGrob(title,
                                                gp=grid::gpar(fontsize=18,
                                                              fontface="bold")))
  return(grid_plot)
  
}


## Plot results ----
## |||||||||||||||||

# Creates and arranges two types of plots for the results of the analysis.
# In particular:
# - comparison between the "F_true" and "F_hat" components. 
# - comparison "S_true" and "S_hat" scores across multiple components.

plot.results <- function(data, results) {
  n <- min(ncol(data$F_true), ncol(results$F_hat))
  data_plot <- split(t(cbind(data$F_true[,1:n],
                             results$F_hat[,1:n],
                             data$F_true[,1:n] - results$F_hat[,1:n])), seq_len(3*n))
  plot_F <- plot.fields(data_plot, data$nodes,
                      c(paste("fPC_true", 1:n),
                        paste("fPC_hat", 1:n),
                        paste("fPC_true - fPC_hat", 1:n)),
                      n)
  plot_S <- plot.score_scatterplots(data$S_true, results$S_hat)
  plot <- arrangeGrob(plot_F, plot_S, ncol = 1, heights = c(3.5,1))
  
  return(plot)
}

# Creates and arranges a set of grouped boxplots (by components)
# visualizing error metrics and execution times.

plot.results_by_components <- function(errors, times) {
  plots <- list()
  plots[[1]] <- plot.groups_boxplots(errors$errors_X, "RMSE X", "", FALSE)
  plots[[2]] <- plot.groups_boxplots(errors$errors_F, "RMSE F", "", FALSE)
  plots[[3]] <- plot.groups_boxplots(errors$errors_S, "RMSE S", "", FALSE)
  plots[[4]] <- plot.groups_boxplots(times, "Times", "", FALSE)
  legend <- plot.groups_boxplots(times, "Times", "", TRUE)
  mylegend<-g_legend(legend)
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = 2)
  grid_plot <- arrangeGrob(grid_plot, mylegend, nrow=2, heights=c(10, 1))
}

# Generates the accuracy plot for:
# - Data reconstruction
# - Scores (component by component)
# - Loadings (component by component)

plot.accuracy <- function(errors){
  
  # Errors
  errors_X <- errors$errors_X
  errors_F <- errors$errors_F
  errors_S <- errors$errors_S
  
  # X
  component <- paste("Comp.", nComp)
  plot <- plot.groups_boxplots_by_K(errors_X, component, "RMSE X")
  ggsave(paste(images.directory_path, "RMSE_X.jpg", sep = ""),
         plot = plot, width = 10, height = 3*length(K_vect) + 1, dpi = 200)
  
  # F
  for(h in 1:nComp){ 
    component <- paste("Comp.", h)
    plot <- plot.groups_boxplots_by_K(errors_F, component, paste("RMSE F", h, sep = ""))
    ggsave(paste(images.directory_path, "RMSE_F",h,".jpg", sep = ""),
           plot = plot, width = 10, height = 3*length(K_vect) + 1, dpi = 200)
  }
  
  # S
  for(h in 1:nComp){ 
    component <- paste("Comp.", h)
    plot <- plot.groups_boxplots_by_K(errors_S, component, paste("RMSE S", h, sep = ""))
    ggsave(paste(images.directory_path, "RMSE_S",h,".jpg", sep = ""),
           plot = plot, width = 10, height = 3*length(K_vect) + 1, dpi = 200)
  }
  
}


## Plot time complexity analysis in N ----
## |||||||||||||||||||||||||||||||||||||||

plot.N_time_analysis_lines <- function(times, title, LOG, LEGEND) {
  
  norm <- times$value[1]
  
  plot <- ggplot(times, aes(x = N, y = value, color = Model)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = m_colors, labels = models) +
    labs(title = title,
         x = "N",
         y = "Time") +
    theme_bw() +
    theme(
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
  
  if(LOG){
    plot <- plot +
      scale_x_continuous(trans = "log10", breaks = N_vect, limits = c(N_vect[1], max(N_vect) + 500)) + 
      scale_y_continuous(trans = "log10") +
      geom_line(data = data.frame(N = N_vect, y = (N_vect/N_vect[1])^(1/2)*norm), aes(x = N, y = y), linetype = "dashed", color = "grey") +
      geom_line(data = data.frame(N = N_vect, y = N_vect/N_vect[1]*norm), aes(x = N, y = y), linetype = "dashed", color = "grey") +
      geom_line(data = data.frame(N = N_vect, y = (N_vect/N_vect[1])^2*norm), aes(x = N, y = y), linetype = "dashed", color = "grey") +
      geom_text(aes(x = max(N_vect) + 50, y = (max(N_vect)/N_vect[1])^(1/2)*norm, label = "N^1/2"), color = "grey", hjust = 0) +
      geom_text(aes(x = max(N_vect) + 50, y = (max(N_vect)/N_vect[1])*norm, label = "N^1"), color = "grey", hjust = 0) +
      geom_text(aes(x = max(N_vect) + 50, y = (max(N_vect)/N_vect[1])^2*norm, label = "N^2"), color = "grey", hjust = 0)
  } else {
    plot <- plot +
      scale_x_continuous(breaks = N_vect, limits = c(N_vect[1], max(N_vect)))
  }
  
  if(LEGEND){
    plot <- plot +
      guides(color=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(color = "none")
  }
}

plot.N_time_analysis <- function(times, LOG) {
  plots <- list()
  for(K in K_vect){
    name_K <- paste("K", K, sep = "") 
    norm <- times[[name_K]]$value[1]
    plots[[name_K]] <- plot.N_time_analysis_lines(times[[name_K]], paste("K = ",K), LOG, FALSE)
  }
  legend <- g_legend(plot.N_time_analysis_lines(times[[name_K]], paste("K = ",K), LOG, TRUE))
  grid_plot <- arrangeGrob(grobs = plots, ncol = 2)
  plot <- arrangeGrob(grid_plot, legend, ncol = 1,  heights=c(10, 2))
}


## Plot time complexity analysis in K ----
## |||||||||||||||||||||||||||||||||||||||

plot.K_time_analysis_lines <- function(times, title, color) {
  norm <- times$value[1]
  
  plot <- ggplot(times, aes(x = K, y = value, color = N)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = rep(color, length(N_vect)), labels = N_vect) +
    labs(title = title,
         x = "K",
         y = "Time") +
    scale_x_continuous(trans = "log10", breaks = K_vect, limits = c(K_vect[1], max(K_vect) + 150)) + 
    scale_y_continuous(trans = "log10") + 
    theme_bw() +
    theme(
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      legend.title = element_blank(),
      legend.position = "bottom"
    ) + guides(color = "none")  +
    geom_line(data = data.frame(K = K_vect, y = (K_vect/K_vect[1])^1*norm), aes(x = K, y = y), linetype = "dashed", color = "grey") +
    geom_line(data = data.frame(K = K_vect, y = (K_vect/K_vect[1])^2*norm), aes(x = K, y = y), linetype = "dashed", color = "grey") +
    geom_line(data = data.frame(K = K_vect, y = (K_vect/K_vect[1])^3*norm), aes(x = K, y = y), linetype = "dashed", color = "grey") +
    geom_text(aes(x = max(K_vect) + 30, y = (max(K_vect)/K_vect[1])^1*norm, label = "K^1"), color = "grey", hjust = 0) +
    geom_text(aes(x = max(K_vect) + 30, y = (max(K_vect)/K_vect[1])^2*norm, label = "K^2"), color = "grey", hjust = 0) +
    geom_text(aes(x = max(K_vect) + 30, y = (max(K_vect)/K_vect[1])^3*norm, label = "K^3"), color = "grey", hjust = 0)
}


plot.K_time_analysis <- function(times) {
  plots <- list()
  for(i in 1:length(models)){
    
    temp_times <- times[times$Model == models[i],]
    plots[[i]] <- plot.K_time_analysis_lines(temp_times, models[i], m_colors[i])
    
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = 2)
  
  return(grid_plot)
}