# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Spatial Transcriptomics %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

setwd("~/master_thesis/applications/spatial_transcriptomics")
load("../../utils/functions/cat_utilities.RData")

cat.script_title("Spatial Transcriptomics")


# ||||||||||||||
# Libraries ----
# ||||||||||||||

cat.section_title("Libraries")

source("scripts/sections/libraries.R")


# ||||||||||||||
# Functions ----
# ||||||||||||||

cat.section_title("Functions")

load("../../utils/functions/plot_utilities.RData")

source("scripts/functions/plots.R")
source("scripts/functions/mesh.R")


# |||||||||||||||||||||
# Global variables ----
# |||||||||||||||||||||

cat.section_title("Global variables")

PLOT = TRUE
VERBOSE = TRUE

directory.initial_data <- "data/HER2/"
directory.processed_data <- "data/HER2/processed_data/"
name.dataset <- "HER2"

RUN <- list()
RUN[["Pre-Processing"]] <- FALSE
RUN[["Mesh Generation"]] <- TRUE


# |||||||||||||||||||
# Initialization ----
# |||||||||||||||||||

cat.section_title("Initialization")

source("scripts/sections/initialization.R")


# |||||||||
# Data ----
# |||||||||

cat.section_title("Data")

load(paste(directory.initial_data, name.dataset, ".RData", sep = ""))

# Data content
# - locations:      (data.frame)    #locations x 2
# - counts:         (dgCMatrix)     #genes x #locations
# - true_labels:    (data.frame)    #locations x 1  (NULL if not available)

# Save initial data
locations.initial <- locations
counts.initial <- counts
names.locations.initial <- rownames(locations)
names.genes.initial <- rownames(counts)


# |||||||||||||||||||
# Pre-Processing ----
# |||||||||||||||||||

sparkversion <- "spark" # "sparkX"
numCores_spark <- 1
number.genes <- NULL

# Pre-processing
if(RUN[["Pre-Processing"]]){
  tic()
  source("scripts/sections/preprocessing.R")
  elapsed <- toc(quiet = TRUE)
  # Execution time
  times[["Pre-Processing"]] <- as.numeric(elapsed$toc - elapsed$tic)
}

# Load processed data
load(paste(directory.processed_data, name.dataset, "_processed.RData", sep = ""))

# Processed data
locations <- locations.significant
counts <- counts.significant
counts.normalized <- counts.normalized
names.locations <- rownames(locations)
names.genes <- rownames(counts)

# Stats
cat("\nStats")
cat(paste("\n- Initial number of locations:", nrow(locations.initial)))
cat(paste("\n- Final number of locations:", nrow(locations)))
cat(paste("\n- Initial number of genes: ", nrow(counts.initial)))
cat(paste("\n- Final number of genes: ", nrow(counts)))

# Clean
rm(locations.significant, counts.significant)


# ||||||||||||||||||||
# Mesh Generation ----
# ||||||||||||||||||||

cat.section_title("Mesh Generation")

# Data exploration
if(PLOT){
  ggplot() +
    standard_plot_settings() + 
    xlab("x") + ylab("y") + ggtitle("Initial Locations") +
    geom_sf(data = st_as_sf(SpatialPoints(locations)), color = "black", size = 1)
}

# Hyper-parameters
h <- 1
bbox <- NULL
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))
type <- "square"
simplification <- 0.15
maximum_area <- 0.3

# Mesh generation
if(RUN[["Mesh Generation"]]){
  tic()
  source("scripts/sections/mesh_generation.R")
  elapsed <- toc(quiet = TRUE)
  # Execution time
  times[["Mesh Generation"]] <- as.numeric(elapsed$toc - elapsed$tic)
}

# Load generated data
load(paste(directory.processed_data, name.dataset, "_mesh.RData", sep = ""))

# Processed data
locations <- locations.final
names.locations <- rownames(locations)
mesh <- mesh
lattice <- lattice

# Stats
cat("\nStats")
cat(paste("\n- Number of locations:", nrow(locations)))
cat(paste("\n- Number of elements: ", nrow(mesh$triangles)))
cat(paste("\n- Number of nodes: ", nrow(mesh$nodes)))

# Clean
rm(locations.final)