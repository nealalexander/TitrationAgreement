# Analysis of agreement of titration assays in R

This repository contains R functions and example scripts related to the following paper:

> Alexander N, Schmidt WP. (2022). *Agreement and error of titration assays*. J Immunol Methods. 2022;502:113210.

Please cite this paper in any work resulting from the use of this repository.

## Overview

This repository provides:

- R functions implementing the paper's methods for analysis of titration data.
- Example datasets.
- Worked examples demonstrating how to use the functions and reproduce figures.


## Repository structure

```
R/
    R functions

examples/
    Example scripts demonstrating the functions

data/
    Example datasets

figures/
    Figures included with the repository


```

## Required packages

- Required packages:
  - ExtDist
  - mvtnorm

Install any missing packages with, for example:

```r
install.packages("ExtDist")
```

## Getting started

1. Download or clone this repository.
2. Open the `examples` folder.
3. Run one of the example scripts.

Each example script loads the required functions from the `R` folder and reads the example data from the `data` folder.

## Status

This repository contains a version of the code used for the published analyses in the accompanying paper. Minor bug fixes or documentation improvements may be made over time.

## License

Released under the MIT License (see `LICENSE`).

## Contact

Neal Alexander

Email: neal.alexander@lshtm.ac.uk