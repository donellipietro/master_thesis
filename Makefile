# Define the Rscript command
RSCRIPT := Rscript

# Directories
UTILS_FUNCTIONS_DIR := utils/functions
UTILS_DIR := utils/

# Get a list of all R scripts
UTILS_FUNCTIONS_SCRIPTS := $(wildcard $(UTILS_FUNCTIONS_DIR)/*.R)

# List of the calls
UTILS_FUNCTIONS := $(patsubst %.R,%.RData,$(UTILS_FUNCTIONS_SCRIPTS))

# Targets
.PHONY: all build clean clone_fdaPDE install_fdaPDE2 install_penR1FPLS install_SpatialPCA install_libraries

# Default target (all): Run all the main tasks
all: install_libraries 

# Target: build
# Build the functions and directories needed
build: $(UTILS_FUNCTIONS)
	mkdir -p libraries

# Target: clean
# Clean intermediate files
clean:
	$(RM) $(UTILS_FUNCTIONS_DIR)/*.RData

# Target: distclean
# Clean all installation files and directories
distclean: clean
	$(RM) libraries/ -r

# Target: clone_fdaPDE
# Install fdaPDE core
DIR_TO_CHECK_fdaPDE := libraries/fdaPDE
clone_fdaPDE: build
	@if [ -d $(DIR_TO_CHECK_fdaPDE) ]; then \
    	echo "fdaPDE repository already exists."; \
    else \
		echo "Cloning fdaPDE repository."; \
        cd libraries && git clone https://github.com/donellipietro/fdaPDE.git; \
    fi
	cd libraries/fdaPDE && git checkout develop_donelli

# Target: install_libraries
# Install fdaPDE2, SpatialPCA and pR1FPLS
install_fdaPDE2: clone_fdaPDE utils/install_fdaPDE2.RExec
install_SpatialPCA: utils/install_SpatialPCA.RExec
install_penR1FPLS: utils/install_penR1FPLS.RExec
install_libraries: build install_fdaPDE2 install_SpatialPCA install_penR1FPLS

# Rule for R functions
%.RData: $(patsubst %.RData,%.R, $@)
	$(RSCRIPT) $(patsubst %.RData,%.R, $@)

# Rule for R scripts	
%.RExec: $(patsubst %.RExec,%.R, $@)
	$(RSCRIPT) $(patsubst %.RExec,%.R, $@)