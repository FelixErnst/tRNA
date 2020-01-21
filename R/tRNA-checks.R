#' @include tRNA.R
NULL

#' @name istRNAGRanges
#' @aliases istRNAGRanges
#'
#' @title tRNA compatibility check
#'
#' @description
#' \code{istRNAGRanges} checks whether a GRanges object contains the
#' information expected for a tRNA result. This is used internally to ensure the
#' the required data is present in the input.
#'
#' @param x the \code{GRanges} object to test for compatibility.
#'
#' @return a logical value
#'
#' @examples
#' data("gr", package = "tRNA")
#' istRNAGRanges(gr)
NULL
#' @rdname istRNAGRanges
#' @export
setMethod(
  f = "istRNAGRanges",
  signature = signature(x = "GRanges"),
  definition = function(x) .check_trna_granges(x, TRNA_FEATURES))

# checks whether a GRanges object is tRNA compatible
.check_trna_granges <- function(gr,features){
  if(!is(gr,"GRanges")){
    warning("Input is not a GRanges object.", call. = FALSE)
    return(FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    warning("Input GRanges object does not meet the requirements of the ",
            "function. The following columns are expected:\n'",
            paste(features, collapse = "', '"),
            "'.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}
