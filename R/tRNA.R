#' @title tRNA: analyzing tRNA sequences and structures
#'
#' @author Felix G M Ernst [aut]
#'
#' @description The tRNA package allows feature information of tRNAs to be
#' accessed and list of tRNA to be subset based on these features. The main
#' purpose is to unify overlapping functions from the tRNAscanImport and
#' tRNAdbImport packages. The functionality is currently under development and
#' may change. The package expects a GRanges object with certain columns as
#' input. The following columns are a requirement: \code{tRNA_length},
#' \code{tRNA_type}, \code{tRNA_anticodon}, \code{tRNA_seq}, \code{tRNA_str},
#' \code{tRNA_CCA.end}. Outputs of tRNAscanImport and tRNAdbImport meet these
#' requirements.
#'
#' Have a look at the vignette for an overview of the functionality. Additional
#' functions are planned to be added in the future.
#'
#' @docType package
#' @name tRNA
NULL

#' @import methods
#' @import GenomicRanges
#' @import Biostrings
#' @import Structstrings
#' @import Modstrings
NULL

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
  "acceptorStem",
  "Dprime5",
  "DStem",
  "Dloop",
  "Dprime3",
  "anticodonStem",
  "anticodonLoop",
  "variableLoop",
  "TStem",
  "Tloop",
  "discriminator"
)

TRNA_STRUCTURES_LOOP <- c(
  "Dprime5",
  "Dloop",
  "Dprime3",
  "anticodonLoop",
  "variableLoop",
  "Tloop",
  "discriminator"
)

TRNA_STRUCTURES_PAIRED <- c(
  "acceptorStem",
  "DStem",
  "anticodonStem",
  "variableLoop",
  "TStem"
)

tRNAStructureFunctionList <- list(
  acceptorStem = ".getAcceptorStem",
  Dprime5 = ".getDprime5",
  DStem = ".getDstem",
  Dloop = ".getDloop",
  Dprime3 = ".getDprime3",
  anticodonStem = ".getAnticodonStem",
  anticodonLoop = ".getAnticodonLoop",
  variableLoop = ".getVariableLoop",
  TStem = ".getTstem",
  Tloop = ".getTloop",
  discriminator = ".getDiscriminator")


TRNA_STRUCTURE_ORDER <- c("acceptorStem.prime5",
                          "Dprime5",
                          "DStem.prime5",
                          "Dloop",
                          "DStem.prime3",
                          "Dprime3",
                          "anticodonStem.prime5",
                          "anticodonLoop",
                          "anticodonStem.prime3",
                          "variableLoop",
                          "TStem.prime5",
                          "Tloop",
                          "TStem.prime3",
                          "acceptorStem.prime3",
                          "discriminator")

TRNA_COLOUR_PALETTE <- "Set1"
TRNA_COLOUR_YES <- "green"
TRNA_COLOUR_NO <- "red"
# options
.onLoad <- function(libname,pkgname){
  options("tRNA_colour_palette" = TRNA_COLOUR_PALETTE)
  options("tRNA_colour_yes" = TRNA_COLOUR_YES)
  options("tRNA_colour_no" = TRNA_COLOUR_NO)
}
# get a qualitative color palette
.get_colour <- function(ident){
  .checkValueValidity(ident,c("palette","yes","no"))
  # construct complete ident
  ident <- paste0("tRNA_colour_",ident)
  colour_palette <- getOption(ident)
  if(!.is_a_string(colour_palette)){
    identVar <- toupper(ident)
    colour_palette <- get(identVar)
    warning("The option '",
            ident,
            "' is not a valid palette ", 
            " identifier for RColorBrewer. ",
            "Please set '",
             ident,
             "' to a valid ",
             " palette identifier for a RColorBrewer.",
            call. = FALSE)
  }
  colour_palette
}

# data -------------------------------------------------------------------------

#' @name tRNA-data
#' @title tRNA example data
#' @description Example data for using the tRNA package
#' @docType data
#' 
#' @format object of class \code{GRanges}
#' @keywords datasets tRNA
#' 
#' @usage data(gr)
"gr"
#' @name tRNA-data
#' @usage data(gr_human)
"gr_human"
#' @name tRNA-data
#' @usage data(gr_human2)
"gr_human2"
#' @name tRNA-data
#' @usage data(gr_eco)
"gr_eco"
