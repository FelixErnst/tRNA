#' @include tRNA.R
NULL

#' @name tRNA-subset
#'
#' @title Subsetting tRNA
#'
#' @description
#' The functions \code{has*} can be used to subset the GRanges object containing
#' information about tRNAs.
#'
#' Please not that the settings \code{mismatches} and \code{bulged} take
#' precedence before \code{unpaired} or \code{paired}. This means that by
#' setting either \code{mismatches} or \code{bulged} to either \code{TRUE} or
#' \code{FALSE}, \code{unpaired = TRUE} or \code{paired = TRUE} are
#' automatically set to allow specific subsetting. If this removes elements from
#' the results, please consider constructing a logical vectors with two calls as
#' suggested in the examples.
#'
#' @param x a GRanges object from a tRNAscan import or with equivalent
#' information
#' @param length the length as integer
#' @param unpaired logical: has unpaired nucleotides
#' @param paired logical: has paired nucleotides (only used for loops)
#' @param mismatches logical: has mismatched nucleotides
#' @param bulged logical: has mismatched nucleotides of different length
#' creating a bulge
#'
#' @return a logical vector of the length or input GRanges object
#'
#' @examples
#' data("gr", package = "tRNA", envir = environment())
#' hasTStem(gr, length = 5, mismatches = TRUE)
#' gr[hasTStem(gr, length = 5, mismatches = TRUE)]
#' gr[hasDStem(gr, unpaired = FALSE) & hasDStem(gr, mismatches = FALSE)]
NULL

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasTStem",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("TStem",
                                                  x,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasDStem",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("DStem",
                                                  x,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasAcceptorStem",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("acceptorStem",
                                                  x,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasAnticodonStem",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("anticodonStem",
                                                  x,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

.subset_tRNA_stem <- function(ident,
                              gr,
                              length,
                              unpaired,
                              mismatches,
                              bulged){
  ans <- rep(TRUE,length(gr))
  # input check
  .check_trna_granges(gr, TRNA_FEATURES)
  if(!(ident %in% names(tRNAStructureFunctionList))){
    stop("Unknown identifier '",ident,"'.")
  }
  if(!is.na(mismatches)){
    assertive::assert_is_a_bool(mismatches)
    unpaired <- TRUE
  }
  if(!is.na(bulged)){
    assertive::assert_is_a_bool(bulged)
    unpaired <- TRUE
  }
  if(!is.na(unpaired) &&
     is.na(mismatches) &&
     is.na(bulged)){
    assertive::assert_is_a_bool(unpaired)
    if(assertive::is_true(unpaired)){
      mismatches <- TRUE
      bulged <- TRUE
    } else if(assertive::is_false(unpaired)) {
      mismatches <- FALSE
      bulged <- FALSE
    }
  }
  # get structure information, only one information is returned
  # browser()
  strList <- getBasePairing(gr$tRNA_str)
  str <- .get_tRNA_structures(ident,
                              gr,
                              strList)
  strList <- .get_ident_structures(ident,
                                   gr,
                                   strList,
                                   str)[[1]]
  # check if structure was found
  structureFound <- lapply(lapply(unlist(str),width), 
                           function(x){ x != 0 })
  ans <- ans & Reduce("&",structureFound)
  # apply structure subsetting
  if(!is.na(unpaired)){
    ansMismatches <- ans
    ansBulged <- ans
    if(assertive::is_false(mismatches)){
      ansMismatches <- .get_mismatches(strList = strList)[["F"]]
    } else if(assertive::is_true(mismatches)) {
      ansMismatches <- .get_mismatches(strList = strList)[["T"]]
    }
    if(assertive::is_false(bulged)){
      ansBulged <- .get_bulges(strList = strList)[["F"]]
    } else if(assertive::is_true(bulged)) {
      ansBulged <- .get_bulges(strList = strList)[["T"]]
    }
    if(assertive::is_true(mismatches) &&
       assertive::is_true(bulged)){
      ans <- ans & (ansMismatches | ansBulged)
    } else {
      ans <- ans & ansMismatches & ansBulged
    }
  }
  # apply length subsetting
  if(!is.na(length)){
    assertive::assert_all_are_whole_numbers(length)
    isLength <- unlist(.get_structures_length(str))
    ansLength <- isLength == length
    # if no length is found, NA is returned. Set these to FALSE
    ansLength[is.na(ansLength)] <- FALSE
    ans <- ans & ansLength
  }
  return(unname(ans))
}

# returns the length of the structure elements
.get_structures_length <- function(str){
  ans <- lapply(str,
                function(s){
                  if(is.list(s)){
                    z <- sapply(s, .get_widths)
                    z <- rowMeans(z)
                  } else {
                    z <- .get_widths(s)
                  }
                  z
                })
  ans
}


# loop subsetting --------------------------------------------------------------

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasTloop",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length) .subset_tRNA_loop("Tloop",
                                                  x,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasDloop",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length) .subset_tRNA_loop("Dloop",
                                                  x,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasAnticodonLoop",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length) .subset_tRNA_loop("anticodonLoop",
                                                  x,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasVariableLoop",
  signature = signature(x = "GRanges"),
  definition = function(x,
                        length,
                        paired,
                        mismatches,
                        bulged) .subset_tRNA_loop("variableLoop",
                                                  x,
                                                  length,
                                                  paired,
                                                  mismatches,
                                                  bulged))

.subset_tRNA_loop <- function(ident,
                              gr,
                              length,
                              paired,
                              mismatches,
                              bulged){
  ans <- rep(TRUE,length(gr))
  # input check
  .check_trna_granges(gr, TRNA_FEATURES)
  if(!(ident %in% names(tRNAStructureFunctionList))){
    stop("Unknown identifier '",ident,"'.")
  }
  if(!is.na(mismatches)){
    assertive::assert_is_a_bool(mismatches)
    paired <- TRUE
  }
  if(!is.na(bulged)){
    assertive::assert_is_a_bool(bulged)
    paired <- TRUE
  }
  if(!is.na(paired)){
    assertive::assert_is_a_bool(paired)
  }
  # get structure information, only one information is returned
  # browser()
  strList <- getBasePairing(gr$tRNA_str)
  str <- .get_tRNA_structures(ident,
                              gr,
                              strList)
  strList <- .get_ident_structures(ident,
                                   gr,
                                   strList,
                                   str)[[1]]
  # check if structure was found
  structureFound <- lapply(lapply(unlist(str),width), 
                           function(x){ x != 0 })
  ans <- ans & Reduce("&",structureFound)
  # apply paired subsetting
  if(!is.na(paired)){
    pairedAns <- vapply(
      strList,
      function(str){
        any(str$forward > 0)
      },
      logical(1))
    # if any paired nucleotides are found
    if(any(pairedAns) &&
       (!is.na(mismatches) ||
        !is.na(bulged))){
      ansMismatches <- ans
      ansBulged <- ans
      # preparation since the assumption is that equal length is the same on 5'
      # and 3'end
      strList[pairedAns] <- lapply(strList[pairedAns],
                                     function(str){
                                       str[seq_len(floor(nrow(str)/2)),]
                                     })
      #
      if(assertive::is_false(mismatches)){
        ansMismatches <- .get_mismatches(strList = strList)[["F"]]
      } else if(assertive::is_true(mismatches)) {
        ansMismatches <- .get_mismatches(strList = strList)[["T"]]
      }
      if(assertive::is_false(bulged)){
        ansBulged <- .get_bulges(strList = strList)[["F"]]
      } else if(assertive::is_true(bulged)) {
        ansBulged <- .get_bulges(strList = strList)[["T"]]
      }
      if(assertive::is_true(mismatches) &&
         assertive::is_true(bulged)){
        ans <- ans & (ansMismatches | ansBulged)
      } else {
        ans <- ans & ansMismatches & ansBulged
      }
    }
    if(assertive::is_true(paired)){
      ans <- ans &
        pairedAns
    } else {
      ans <- ans &
        !pairedAns
    }
  }
  # apply length subsetting
  if(!is.na(length)){
    assertive::assert_all_are_whole_numbers(length)
    isLength <- unlist(.get_structures_length(str))
    ansLength <- isLength == length
    # if no length is found, NA is returned. Set these to FALSE
    ansLength[is.na(ansLength)] <- FALSE
    ans <- ans & ansLength
  }
  return(unname(ans))
}
