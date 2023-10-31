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