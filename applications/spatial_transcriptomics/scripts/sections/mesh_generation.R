# Mesh generation ----
# ||||||||||||||||||||

set.seed(0)
locations <- SpatialPoints(locations)

## Domain ----
## |||||||||||

# Lattice
# h <- 10
# bbox <- NULL
# seed_point <- NULL
# type <- "hexagonal"
lattice <- generate_lattice(locations, h, bbox, seed_point, type)

# plot(lattice$domain)

# Domain simplification
# simplification <- 0.09
lattice_simplified <- simplify_domain(lattice, simplification)

# plot(lattice_simplified$domain)

# Plot
if(PLOT){
  plot <- plot.discretized_domain(lattice_simplified)
  plot <- plot + ggtitle("Discretized domain") 
  print(plot)
}

## Discarded locations ----
## ||||||||||||||||||||||||

indexes.discarded_locations <- is.na(over(locations, lattice_simplified$domain))
names.locations <- names.locations[!is.na(over(locations, lattice_simplified$domain))]
locations.final <- locations.initial[names.locations,]

if(PLOT){
  plot <- plot.final_locations(locations, SpatialPoints(locations.final), lattice_simplified)
  plot <- plot + ggtitle("Final locations")
  print(plot)
}


## Mesh ----
## |||||||||

# maximum_area <- 1000
mesh <- generate_mesh(lattice_simplified, maximum_area)

# Generated mesh
if(PLOT){
  plot <- plot.fdaPDE_mesh(mesh)
  plot <- plot + ggtitle("Mesh")
  print(plot)
}


## Save mesh ----
## ||||||||||||||

lattice <- lattice_simplified

save(locations.final,
     mesh,
     lattice,
     # Saving options
     file = paste(directory.processed_data, name.dataset, "_mesh", ".RData", sep = ""))


## Clean ----
## ||||||||||

rm(lattice, lattice_simplified, mesh,
   locations, locations.final,
   indexes.discarded_locations, names.locations, 
   bbox, h, seed_point, type, simplification, maximum_area,
   plot)



