# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% SpatialPCA installation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load("utils/functions/cat_utilities.RData")
cat.script_title("SpatialPCA installation")


# |||||||||||||||||||||||||
# Install dependencies ----
# |||||||||||||||||||||||||

cat.section_title("Install dependencies")

# devtools
if (!require(devtools)) {
  install.packages("devtools")
}


# |||||||||||||||||||||||||||||||
# Install SpatialPCA library ----
# |||||||||||||||||||||||||||||||

cat.section_title("Install SpatialPCA library ")

devtools::install_github("shangll123/SpatialPCA")