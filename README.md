# tRNA <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/tRNA/tRNA.png" height="200" align="right">

<!-- badges: start -->
[![R-CMD-check](https://github.com/FelixErnst/tRNA/workflows/R-CMD-check/badge.svg)](https://github.com/FelixErnst/tRNA/actions/)
[![Build Status](https://travis-ci.com/FelixErnst/tRNA.svg?branch=master)](https://travis-ci.com/FelixErnst/tRNA)
[![BioC Build](https://bioconductor.org/shields/build/release/bioc/tRNA.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/tRNA/)
[![codecov](https://codecov.io/gh/FelixErnst/tRNA/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/tRNA)
[![BioC Years](https://bioconductor.org/shields/years-in-bioc/tRNA.svg)](https://doi.org/doi:10.18129/B9.bioc.tRNA)
<!-- badges: end -->

The tRNA package allows feature information of tRNAs to be accessed and list of 
tRNA to be subset based on these features. The main purpose is to unify 
overlapping functions from the
[`tRNAscanImport`](https://doi.org/doi:10.18129/B9.bioc.tRNAscanImport) and 
[`tRNAdbImport`](https://github.com/FelixErnst/tRNAdbImport) packages.

The functionality is currently under development and may change. The package 
expects a `GRanges` object with certain columns as input. The following columns
are a requirement: `tRNA_length`, `tRNA_type`, `tRNA_anticodon`, `tRNA_seq`,
`tRNA_str`, `tRNA_CCA.end`. Outputs of `tRNAscanImport` and `tRNAdbImport` meet
these requirements.

## Installation

The current version of the `tRNA` package is available from Bioconductor.
 
```{r}
BiocManager::install("tRNA")
# Load and attach thepackage
library("tRNA")
```

## Functions

Have a look at the vignette for an overview of the functionality. Additional
functions are planned to be added in the future.
