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

plot.fields <- function(fields, nodes, lim, titles = "", ncol = 1,
                        legend = TRUE) {
  
  lim_base = range(c(fields[[1]], fields[[3]], fields[[4]]))
  
  # Settings
  plot_settings <- ggplot() +
    theme_minimal() +
    coord_fixed() +
    labs(fill = "") + xlab("") + ylab("") +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        color = "black", face = "bold", size = 10,
        hjust = 0.5, vjust = -2
      ),
      legend.position = "bottom"
    )
  
  plots <- list()
  for (i in 1:length(fields)) {
    plots[[i]] <- plot_settings +
      ggtitle(titles[[i]]) +
      geom_raster(data = data.frame(x = nodes[, 1],
                                    y = nodes[, 2],
                                    z = fields[[i]]), 
                  aes(x = x, y = y, fill = z))
    if(!legend){
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
    
    if(i == 2){
      plots[[i]] <- plots[[i]] + scale_fill_viridis(lim = lim) 
    } else {
      plots[[i]] <- plots[[i]] + scale_fill_viridis(lim = lim_base) 
    }
      
  }
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = ncol)
  return(grid_plot)
}


## Plot groups boxplots ----
## |||||||||||||||||||||||||

# Creates grouped boxplots taking as input a data frame data that contains
# a grouping variable Group. If LEGEND is set to TRUE, it includes a legend for
# the data categories, and if set to FALSE, it removes the legend.

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
      text = element_text(size=10),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.spacing.x = unit(0, 'cm'),
      legend.margin = margin(l = 1, t = 0, unit = "cm"),
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

