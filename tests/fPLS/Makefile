# Define the Rscript command
RSCRIPT := Rscript

# Directories
FUNCTIONS_DIR := scripts/functions

# Get a list of all R scripts
FUNCTIONS_SCRIPTS := $(wildcard $(FUNCTIONS_DIR)/*.R)

# List of the calls
FUNCTIONS := $(patsubst %.R,%.RData,$(FUNCTIONS_SCRIPTS))


# Targets
.PHONY: all build clean distclean

# Default target (all): Run all the main tasks
all: build

# Target: build
# Build the functions and directories needed
build: $(FUNCTIONS)
	make -C ./../../. build
	mkdir -p images results data
	
# Target: run
# Run tests
run_data_generation_problems: build
	Rscript scripts/main_data_generation_problems.R
run_models_comparison: build
	Rscript scripts/main_models_comparison.R
run: run_data_generation_problems run_models_comparison
	
# Target: clean
# Clean intermediate files
clean:
	@$(RM) $(FUNCTIONS_DIR)/*.RData
	@$(RM) *.aux *.log *.pdf *.txt
	@$(RM) .Rhistory
	
# Target: distclean
# Clean all generated files and directories
distclean: clean
	@$(RM) data/ -r
	@$(RM) images/ -r
	@$(RM) results/ -r
	
# Rule for R functions
%.RData: $(patsubst %.RData,%.R, $@)
	$(RSCRIPT) $(patsubst %.RData,%.R, $@)

# Rule for R scripts	
%.RExec: $(patsubst %.RExec,%.R, $@)
	$(RSCRIPT) $(patsubst %.RExec,%.R, $@)
	

