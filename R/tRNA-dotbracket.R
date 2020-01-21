#' @include tRNA.R
NULL

#' @name gettRNABasePairing
#' @aliases gettRNALoopIDs
#'
#' @title Accessing Dot Bracket annotation of tRNAs
#'
#' @description
#' \code{gettRNABasePairing} converts the dot bracket annotation into a
#' \code{DotBracketDataFrame}. Base pairing is indicated by cosrresponding 
#' numbers in the forward and reverse columns. For more detail have a look at
#' \code{\link[Structstrings:getBasePairing]{getBasePairing}}.
#' 
#' \code{gettRNALoopIDs} converts the dot bracket annotation into a 
#' \code{LoopIDList}. For more details have a look at 
#' \code{\link[Structstrings:getBasePairing]{getLoopIDList}}.
#'
#' @return
#' \code{gettRNABasePairing}: 
#' The result is a \code{DotBracketDataFrame} with following columns: pos, 
#' forward, reverse, character and base. If a position is unpaired, forward and 
#' reverse will be \code{0}, otherwise it will match the base paired positions.
#' 
#' \code{gettRNALoopIDs}: return a list of list of loop ids.
#'
#' @param x a GRanges object created by \code{import.tRNAscanAsGRanges} or
#' GRanges with equivalent information. The \code{tRNA_str} and \code{tRNA_seq} 
#' columns will be used to construct a StructuredXStringSet and used for input 
#' into \code{getBasePairing}.
#' @param with.nucleotides a single logical value: should the nucleotides be 
#' saved alongside the base pairing information in the 'base' column?
#' 
#' @examples 
#' data("gr", package = "tRNA")
#' gettRNABasePairing(gr[1])
#' gettRNALoopIDs(gr[1])
NULL

#' @rdname gettRNABasePairing
#' @export
setMethod(
  f = "gettRNABasePairing",
  signature = signature(x = "GRanges"),
  definition = function(x, with.nucleotides = FALSE) {
    .check_trna_granges(x, TRNA_FEATURES)
    seq <- x$tRNA_seq
    str <- x$tRNA_str
    if(!is(seq,"RNAStringSet") && !is(seq,"ModRNAStringSet")){
      seq <- as(seq,"RNAStringSet")
    }
    if(!is(str,"DotBracketStringSet")){
      str <- as(str,"DotBracketStringSet")
    }
    if(with.nucleotides){
      strseq <- do.call(paste0("Structured",class(seq)),
                        list(x = seq,
                             structure = str))
    } else {
      strseq <- str
    }
    Structstrings::getBasePairing(strseq, return.sequence = with.nucleotides)
  }
)

#' @rdname gettRNABasePairing
#' @export
setMethod(
  f = "gettRNALoopIDs",
  signature = signature(x = "GRanges"),
  definition = function(x) {
    .check_trna_granges(x, TRNA_FEATURES)
    str <- x$tRNA_str
    if(!is(str,"DotBracketStringSet")){
      str <- as(str,"DotBracketStringSet")
    }
    Structstrings::getLoopIndices(str)
  }
)
