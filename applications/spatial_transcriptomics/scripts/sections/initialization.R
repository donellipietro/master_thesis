# Initialization ----
# |||||||||||||||||||

## Variables ----
## ||||||||||||||

# Room for times
times <- list()


## Directories ----
## ||||||||||||||||

# Data
if (!file.exists(directory.processed_data)){
  dir.create(directory.processed_data)
  if(VERBOSE){
    cat("\n- Directory for pocessed data created!\n")
  }
}

# Images
if (!file.exists(directory.images)){
  dir.create(directory.images)
  if(VERBOSE){
    cat("\n- Directory for images created!\n")
  }
}
