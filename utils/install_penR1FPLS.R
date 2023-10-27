# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% penR1FPLS installation %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%


load("utils/functions/cat_utilities.RData")
cat.script_title("penR1FPLS installation")


# |||||||||||||||||||||||||
# Install dependencies ----
# |||||||||||||||||||||||||

cat.section_title("Install dependencies")

# devtools
if (!require(devtools)) {
  install.packages("devtools")
}


# ||||||||||||||||||||||||||||||
# Install penR1FPLS library ----
# ||||||||||||||||||||||||||||||

cat.section_title("Install penR1FPLS library ")

devtools::install_github("hhroig/penR1FPLS", dependencies = TRUE)