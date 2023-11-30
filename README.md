# Master Thesis

This code was developed for the MSc. thesis in Mathematical Engineering at Politecnico di Milano titled **New Functional Data Analysis Methods, with Applications to Spatial Transcriptomics and Neuroimaging**.

## Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)
-   [Tests](#tests)
-   [Applications](#applications)
-   [Author](#author)
-   [Acknowledgment](#acknowledgment)

## Installation

To get started with this project, follow these steps:

1. Clone the Repository:

    ```bash
    git clone https://github.com/donellipietro/master_thesis.git
    ```

2. Navigate to the project directory:

    ```bash
    cd master_thesis
    ```

3. Launch the installer:
    ```bash
    make all
    ```

If the automatic installation fails you can perform it manually following these steps:

1. Clone the fdaPDE library into the `libraries` directory and move the the correct branch:
    ```bash
    cd libraries
    git clone https://github.com/donellipietro/master_thesis.git
    cd libraries/fdaPDE
    git checkout develop_donelli
    cd ../..
    ```
2. Install the `fdaPDE2` library using the `R` script [install_fdaPDE2.R](utils/install_fdaPDE2.R)
3. Install the `SpatialPCA` library using the `R` script [install_SpatialPCA.R](utils/install_SpatialPCA.R)
4. Install the `penR1FPLS` library using the `R` script [install_penR1FPLS.R](utils/install_penR1FPLS.R)

## Usage

Once the installation is done all tests and applications are ready to be used. The general procedure is to enter the desired test directory (the one containing the `Makefile`) and run the `make build` command from the terminal, this will initialize what is needed. If you are using `RStudio` the working directory must be set in this same directory. The following `make` targets are also available.

-   `make clean`: removes intermediate files
-   `make distclean`: removes all generated files and directories
-   `make run`: runs all the test
-   `make run_{test_name}`: runs the test `main_{test_name}`

For example if one wishes to run the `models_comparison` test for fPLS the procedure to be followed will be:

```bash
# Assuming the current directory to be master_thesis
cd tests/fPLS
make build
make run_models_comparison
```

## Tests

-   [fPCA](tests/fPCA/)
-   [FR-PDE](tests/FRPDE/)
-   [fPLS](tests/fPLS/)

## Applications

-   [Spatial Transcriptomics](applications/spatial_transcriptomics/)
-   [Brains](applications/brain/)

## Author

Pietro Donelli [@donellipietro](https://github.com/donellipietro)

### Acknowledgments

-   Alessandro Palummo [@AlePalu](https://github.com/AlePalu)
-   Harold A. Hdez. [@hhroig](https://github.com/hhroig)
-   Drew Burns [@Drew4495](https://github.com/Drew4495)
