# tRNA

The tRNA package allows featrue information of tRNAs to be accessed and list of 
tRNA to be subset based on these features. The main purpuse is to unify 
overlapping functions from the
[`tRNAscanImport`](https://doi.org/doi:10.18129/B9.bioc.tRNAscanImport) and 
[`tRNAdbImport`](https://github.com/FelixErnst/tRNAdbImport) packages.

The functionallity is currently under development and may change. The package 
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

Have a look at the vignette for an overview of the functionallity. Additional
functions are planned to be added in the future.
