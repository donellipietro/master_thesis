# Initialization ----
# |||||||||||||||||||

## Variables ----
## ||||||||||||||

# Room for times
times <- list()


## Directories ----
## ||||||||||||||||

if (!file.exists(directory.processed_data)){
  dir.create(directory.processed_data)
  if(VERBOSE){
    cat("\n- Directory for pocessed data created!\n")
  }
}
