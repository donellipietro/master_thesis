# Fast and Accurate Spatial Domain Detection for Spatial transcriptomics, a Closed-Form Functiona-PCA Approach

This code was developed for the ...

## Description

Spatial domain detection aims to define regions in biological tissue based on genetic expression patterns using technologies such as Spatial Transcriptomics, which measures gene expression at the RNA level while retaining expression location. In the past few years, several computational methods have been developed to infer these spatial domains for purposes like identifying tumor microenvironments. Yet, their primary limitation in achieving widespread routine application is the slow algorithmic runtime, relying on iterative optimization instead of a closed-form solution. Further, many of these methods exhibit limited flexibility, often restricted to detecting only predefined spatial patterns, and leave potential for improved accuracy in spatial domain detection. Here, we introduce a 10-fold faster, compared to a current leading method, dimension reduction algorithm that has a 1) closed-form solution, 2) is more flexible in identifying spatial patterns, and 3) demonstrates improved accuracy for the datasets we analyzed. Our method synergizes several mathematical techniques such as 1) functional-PCA, 2) Finite Element Analysis, and 3) Graph Regularized Singular Value Decomposition, which allows us to create a continuous functional representation of the genetic expression data while taking into account the local neighborhood of each measurement. We illustrate that the improvements of our closed-form functional-PCA approach does not sacrifice other features of previous techniques, such as high resolution domain reconstruction, and while not shown here, we claim that our formulation also lends itself to 3D (multiple tissue slices) spatial domain detection and spatial trajectory inference.

## Table of Contents

-   [Installation](#installation)
-   [Downloads](#downloads)
-   [Tests](#tests)
-   [Authors](#authors)
-   [Acknowledgment](#acknowledgment)

## Installation

...

## Downloads

All the data used in the analysis are publicly available and can be downloaded from [here](#).

## Tests

-   [HER2.R](scripts/main_HER2.R)
-   [DLPFC](#)
-   [SlideseqCerebellum](#)
-   [SlideseqV2Hippocampus](#)

## Authors

Pietro Donelli [@donellipietro](https://github.com/donellipietro)

### Acknowledgment

-   Alessandro Palummo [@AlePalu](https://github.com/AlePalu)
-   Drew Burns [@Drew4495](https://github.com/Drew4495)
