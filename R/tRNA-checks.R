#' @include tRNA.R
NULL

#' @name checktRNAGRanges
#' @aliases checktRNAGRanges
#'
#' @title tRNA compatibility check
#'
#' @description
#' \code{checktRNAGRanges} checks whether a GRanges object contains the
#' information expected for a tRNA result.
#'
#' @param gr the \code{GRanges} object to test
#'
#' @return a logical value
#'
#' @examples
#' data("gr", package = "tRNA", envir = environment())
#' checktRNAGRanges(gr)
NULL
#' @rdname checktRNAGRanges
#' @export
setMethod(
  f = "checktRNAGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr) .check_trna_granges(gr,
                                                TRNA_FEATURES))

# checks whether a GRanges object is trnascan compatible
.check_trna_granges <- function(gr,features){
  if(class(gr) != "GRanges"){
    stop("Input is not a GRanges object.",
         call. = FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. Please refer to the vignette of tRNAscanImport for ",
         "an exmaple on what information is expected.",
         call. = FALSE)
  }
  return(TRUE)
}
