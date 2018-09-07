#' @title
#' tRNA: analyzing tRNA sequences and structures
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' title
#'
#' @docType package
#' @name tRNA
NULL

#' @import methods
#' @import GenomicRanges
#' @import assertive
NULL
requireNamespace("assertive")

# constants tRNA ---------------------------------------------------------------

TRNA_FEATURES <- c(
  "tRNA_length",
  "tRNA_type",
  "tRNA_anticodon",
  "tRNA_seq",
  "tRNA_str",
  "tRNA_CCA.end"
)

TRNA_STRUCTURES <- c(
  "anticodonloop",
  "Dloop",
  "Tloop",
  "acceptorStem",
  "anticodonStem",
  "DStem",
  "TStem",
  "variableLoop",
  "discriminator"
)

# data -------------------------------------------------------------------------

#' @name tRNA-data
#' @title tRNA example data
#' @description Example data for using the tRNA package
#' @docType data
#' @usage tRNA
#' @format object of class \code{GRanges}
#' @keywords datasets tRNA
"gr"
