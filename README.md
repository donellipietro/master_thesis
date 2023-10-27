# Master Thesis

This code was developed for the MSc. thesis in Mathematical Engineering at Politecnico di Milano titled **Insert Title**.

## Description

Provide a brief description of your project. Explain what it does, why it exists, and any other relevant information.

## Table of Contents

-   [Installation](#installation)
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

## Tests

-   [fPCA](#fCPA)
-   [FR-PDE](#FR-PDE)
-   [fPLS](#fPLS)

### fCPA

### FR-PDE

### fPLS

## Applications

-   [Spatial Transcriptomics](#spatial_transcriptomics)
-   [Brains](#brains)
-   [Pipe](#pipe)

### Spatial Transcriptomics

### Brains

### Pipe

## Author

Pietro Donelli [@donellipietro](https://github.com/donellipietro)

### Acknowledgment

-   Alessandro Palummo [@AlePalu](https://github.com/AlePalu)
