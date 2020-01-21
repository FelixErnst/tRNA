#' @include tRNA.R
NULL

#' @name gettRNASummary
#' @aliases gettRNASummary
#' 
#' @title Summary of tRNA features
#' 
#' @description
#' \code{gettRNASummary} prepares a DataFrame with the aggregated features
#' of tRNAs from a GRanges object. Logical values are converted to numeric
#' values.
#' 
#' @param x a GRanges or a GRangesList object. All elements have to pass the
#' \code{istRNAGRanges} test.
#'
#' @return a DataFrame object 
#' 
#' @export
#' @examples
#' data("gr", package = "tRNA")
#' gettRNASummary(gr)
NULL

TRNA_SUMMARY_FEATURES <- c("gc" = ".get_gc_content",
                           "width" = ".get_widths",
                           "length" = ".get_seq_lengths",
                           "cca" = ".get_cca_ends")


TRNA_SUMMARY_FEATURES_LIST <- c("feature_valid" = ".get_features_valid",
                                "feature_dstem" = ".get_features_dstem",
                                "feature_tstem" = ".get_features_tstem")
TRNA_SUMMARY_FEATURES_LIST_SUBSET <- 
  c("feature_length" = ".get_features_length")
TRNA_SUMMARY_FEATURES_LIST_SUBSET_PAIRED <- 
  c("unpaired" = ".get_features_unpaired",
    "mismatches" = ".get_features_mismatches",
    "bulged" = ".get_features_bulged")

TRNA_SUMMARY_FEATURES_OPTIONAL <- c(
  "score" = ".get_score",
  "tRNAscan_potential.pseudogene" = ".get_potential_pseudogene",
  "tRNAscan_intron.start" = ".get_introns",
  "tRNAscan_score" = ".get_scan_score",
  "tRNAscan_hmm.score" = ".get_hmm_score",
  "tRNAscan_sec.str.score" = ".get_secondary_structure_score",
  "tRNAscan_infernal" = ".get_infernal_score")

TRNA_SUMMARY_FEATURES_RENAMED <- c("tRNAscan_intron.start" = "tRNAscan_intron")
#' @rdname gettRNASummary
#' @export
setMethod(
  f = "gettRNASummary",
  signature = signature(x = "GRangesList"),
  definition = function(x) {
    ans <- lapply(x, gettRNASummary)
    names(ans) <- names(x)
    return(ans)
  }
)

#' @rdname gettRNASummary
#' @export
setMethod(
  f = "gettRNASummary",
  signature = signature(x = "GRanges"),
  definition = function(x) {
    # Input check
    istRNAGRanges(x)
    # get default features
    df <-  S4Vectors::DataFrame(lapply(TRNA_SUMMARY_FEATURES,
                                       function(f){
                                         do.call(f,list(x))
                                       })) 
    # get list features from structure informations
    strList <- getBasePairing(x$tRNA_str)
    str <- .get_tRNA_structures(TRNA_STRUCTURES,
                                x,
                                strList)
    strListSubsetAll <- .get_ident_structures(ident = TRNA_STRUCTURES,
                                              gr = x,
                                              strList = strList,
                                              str = str)
    strListSubsetPaired <- strListSubsetAll[TRNA_STRUCTURES_PAIRED]
    data <- lapply(TRNA_SUMMARY_FEATURES_LIST,
                   function(f){
                     do.call(f,list(x, str))
                   })
    dataSubset <- lapply(TRNA_SUMMARY_FEATURES_LIST_SUBSET,
                         function(f){
                           do.call(f,list(x, strListSubsetAll))
                         })
    dataSubsetPaired <- lapply(TRNA_SUMMARY_FEATURES_LIST_SUBSET_PAIRED,
                   function(f){
                     do.call(f,list(x, strListSubsetPaired))
                   })
    dfAdd <- do.call(cbind,
                     c(unname(lapply(data, S4Vectors::DataFrame)),
                       unname(lapply(dataSubset, S4Vectors::DataFrame)),
                       unname(lapply(dataSubsetPaired, S4Vectors::DataFrame))))
    colnames(dfAdd) <- unname(c(unlist(lapply(data,names)),
                                unlist(lapply(dataSubset,names)),
                                unlist(lapply(dataSubsetPaired,names))))
    df <- cbind(df,dfAdd)
    # get optional features
    dataFeatures <- colnames(S4Vectors::mcols(x))
    dataFeatures <- dataFeatures[dataFeatures %in% 
                                   names(TRNA_SUMMARY_FEATURES_OPTIONAL)]
    if(length(dataFeatures) > 0){
      dfAdd <- S4Vectors::DataFrame(
        lapply(TRNA_SUMMARY_FEATURES_OPTIONAL[dataFeatures],
               function(f){
                 do.call(f,list(x))
               })) 
      df <- cbind(df, dfAdd)
    }
    # rename some colums
    cnames <- colnames(df)
    f <- match(names(TRNA_SUMMARY_FEATURES_RENAMED), cnames)
    cnames[f] <- TRNA_SUMMARY_FEATURES_RENAMED
    colnames(df) <- cnames
    # convert columns to fix types
    df$width <- as.integer(df$width)
    df$length <- as.integer(df$length)
    df$cca <- as.logical(df$cca)
    df$features_all_valid <- as.logical(df$features_all_valid)
    df$features_Dstem_found <- as.logical(df$features_Dstem_found)
    df$features_Tstem_found <- as.logical(df$features_Tstem_found)
    df$acceptorStem_length <- as.integer(df$acceptorStem_length)
    df$Dprime5_length <- as.integer(df$Dprime5_length)
    df$DStem_length <- as.integer(df$DStem_length)
    df$Dloop_length <- as.integer(df$Dloop_length)
    df$Dprime3_length <- as.integer(df$Dprime3_length)
    df$anticodonStem_length <- as.integer(df$anticodonStem_length)
    df$anticodonLoop_length <- as.integer(df$anticodonLoop_length)
    df$variableLoop_length <- as.integer(df$variableLoop_length)
    df$TStem_length <- as.integer(df$TStem_length)
    df$Tloop_length <- as.integer(df$Tloop_length)
    df$discriminator_length <- as.integer(df$discriminator_length)
    df$acceptorStem_unpaired <- as.logical(df$acceptorStem_unpaired)
    df$DStem_unpaired <- as.logical(df$DStem_unpaired)
    df$anticodonStem_unpaired <- as.logical(df$anticodonStem_unpaired)
    df$variableLoop_unpaired <- as.logical(df$variableLoop_unpaired)
    df$TStem_unpaired <- as.logical(df$TStem_unpaired)
    df$acceptorStem_mismatches <- as.logical(df$acceptorStem_mismatches)
    df$DStem_mismatches <- as.logical(df$DStem_mismatches)
    df$anticodonStem_mismatches <- as.logical(df$anticodonStem_mismatches)
    df$variableLoop_mismatches <- as.logical(df$variableLoop_mismatches)
    df$TStem_mismatches <- as.logical(df$TStem_mismatches)
    df$acceptorStem_bulges <- as.logical(df$acceptorStem_bulges)
    df$DStem_bulges <- as.logical(df$DStem_bulges)
    df$anticodonStem_bulges <- as.logical(df$anticodonStem_bulges)
    df$variableLoop_bulges <- as.logical(df$variableLoop_bulges)
    df$TStem_bulges <- as.logical(df$TStem_bulges)
    df$tRNAscan_potential.pseudogene <- as.logical(df$tRNAscan_potential.pseudogene)
    df$tRNAscan_intron <- as.logical(df$tRNAscan_intron)
    #
    return(df)
  }
)

# default features -------------------------------------------------------------

# returns the length of tRNAs
.get_widths <- function(gr){
  as.numeric(BiocGenerics::width(gr))
}
.get_seq_lengths <- function(gr){
  as.numeric(S4Vectors::mcols(gr)$tRNA_length)
}

# get GC content
.get_gc_content <- function(gr){
  freq <- Biostrings::alphabetFrequency(S4Vectors::mcols(gr)$tRNA_seq, 
                                        baseOnly=TRUE)
  gc <- vapply(seq_len(nrow(freq)), function(i){
    sum(freq[i,c("G","C")]) / sum(freq[i,])
  },numeric(1))
  gc
}

# fractions of tRNA with encoded CCA ends
.get_cca_ends <- function(gr){
  as.numeric(S4Vectors::mcols(gr)$tRNA_CCA.end)
}

# Summarize features
.get_features_length <- function(gr, strList){
  lengths <- lapply(TRNA_STRUCTURES,
                    function(ident){
                      l <- .get_str_lengths(ident = ident,
                                            strList = strList[[ident]])
                      # if no length is found, NA is returned. Set these to 0
                      l[is.na(l)] <- 0
                      l
                    })
  names(lengths) <- paste0(TRNA_STRUCTURES,"_length")
  return(lengths)
}
.get_features_valid <- function(gr,
                                str){
  ans <- .get_validity_of_structures(str)
  ans <- lapply(ans,as.numeric)
  names(ans) <- "features_all_valid"
  return(ans)
}
.get_features_dstem <- function(gr,
                                str){
  ans <- .get_validity_of_structures(str["DStem"])
  ans <- lapply(ans,as.numeric)
  names(ans) <- "features_Dstem_found"
  return(ans)
}
.get_features_tstem <- function(gr,
                                str){
  ans <- .get_validity_of_structures(str["TStem"])
  ans <- lapply(ans,as.numeric)
  names(ans) <- "features_Tstem_found"
  return(ans)
}
#
.get_features_unpaired <- function(gr,
                                   strList){
  unpaired <- lapply(TRNA_STRUCTURES_PAIRED,
                    function(ident){
                      if(ident == "variableLoop"){
                        strList[[ident]] <- lapply(strList[[ident]],
                                          function(str){
                                            str[seq_len(floor(nrow(str)/2)),]
                                          })
                      }
                      ansMismatches <- 
                        .get_mismatches(strList = strList[[ident]])[["T"]]
                      ansBulged <- 
                        .get_bulges(strList = strList[[ident]])[["T"]]
                      as.numeric(ansMismatches | ansBulged)
                    })
  names(unpaired) <- paste0(TRNA_STRUCTURES_PAIRED,"_unpaired")
  return(unpaired)
}
.get_features_mismatches <- function(gr,
                                     strList){
  mismatches <- lapply(
    TRNA_STRUCTURES_PAIRED,
    function(ident){
      if(ident == "variableLoop"){
        strList[[ident]] <- lapply(strList[[ident]],
                                   function(str){
                                     str[seq_len(floor(nrow(str)/2)),]
                                   })
      }
      as.numeric(.get_mismatches(strList = strList[[ident]])[["T"]])
    })
  names(mismatches) <- paste0(TRNA_STRUCTURES_PAIRED,"_mismatches")
  return(mismatches)
}
.get_features_bulged <- function(gr,
                                 strList){
  bulges <- lapply(
    TRNA_STRUCTURES_PAIRED,
    function(ident){
      if(ident == "variableLoop"){
        strList[[ident]] <- lapply(strList[[ident]],
                                   function(str){
                                     str[seq_len(floor(nrow(str)/2)),]
                                   })
      }
      as.numeric(.get_bulges(strList = strList[[ident]])[["T"]])
    })
  names(bulges) <- paste0(TRNA_STRUCTURES_PAIRED,"_bulges")
  return(bulges)
}

# optional features ------------------------------------------------------------

# fractions of tRNA with pseudogene
.get_potential_pseudogene <- function(gr){
  as.numeric(S4Vectors::mcols(gr)$tRNAscan_potential.pseudogene)
}

# fractions of tRNA with introns
.get_introns <- function(gr){
  introns <- S4Vectors::mcols(gr)$tRNAscan_intron.start
  introns[is.na(introns)] <- 0
  introns[introns > 0] <- 1
  introns
}

# aggregates the scores
.get_score <- function(gr){
  as.numeric(S4Vectors::mcols(gr)[,c("score")])
}
.get_scan_score <- function(gr){
  as.numeric(S4Vectors::mcols(gr)[,c("tRNAscan_score")])
}
.get_hmm_score <- function(gr){
  as.numeric(S4Vectors::mcols(gr)[,c("tRNAscan_hmm.score")])
}
.get_secondary_structure_score <- function(gr){
  as.numeric(S4Vectors::mcols(gr)[,c("tRNAscan_sec.str.score")])
}
.get_infernal_score <- function(gr){
  as.numeric(S4Vectors::mcols(gr)[,c("tRNAscan_infernal")])
}

# Feature utility functions ----------------------------------------------------
# also used for subsetting functions

.get_validity_of_structures <- function(str,
                                     combine = TRUE){
  widths <- lapply(str,
                   function(s){
                     if(is.list(s)){
                       return(lapply(s,width))
                     }
                     width(s)
                   })
  widths <- c(widths[!vapply(widths,is.list,logical(1))],
              unlist(widths[vapply(widths,is.list,logical(1))],
                     recursive = FALSE))
  ans <- lapply(widths,"!=",0)
  if(combine){
    ans <- list(Reduce("&",ans))
  }
  ans
}

.get_ident_structures <- function(ident, gr,  strList, str){
  if(missing(strList) || missing(str)){
    strList <- getBasePairing(gr$tRNA_str)
    str <- .get_tRNA_structures(ident, gr, strList)
  }
  # Since we will subset the DotBracketDataFrameList it is best to convert them
  # to a normal CompressedDataFrameList, since not all positions remain paired
  # resulting in an invalid DotBracketDataFrame
  strList <- lapply(strList,as.data.frame)
  strList <- lapply(ident,
                    function(id){
                      if(is.list(str[[id]])){
                        strL <- mapply(
                          .subset_structure,
                          split(c(str[[id]]$prime5,
                                  str[[id]]$prime3),
                                c(seq_along(str[[id]]$prime5),
                                  seq_along(str[[id]]$prime3))),
                          strList,
                          SIMPLIFY = FALSE)
                      } else {
                        strL <- mapply(
                          .subset_structure,
                          split(str[[id]],seq_along(str[[id]])),
                          strList,
                          MoreArgs = list(pairedOnly = FALSE),
                          SIMPLIFY = FALSE)
                      }
                      strL
                    })
  names(strList) <- ident
  return(strList)
}

.get_structures_details <- function(strList){
  forwardContinously <- vapply(
    strList,
    function(str){
      if(length(str[str$forward != 0,]$forward) == 0) return(TRUE)
      ans <- .is_continous_evenly_spaced(str[str$forward != 0,]$forward)
    },
    logical(1))
  reverseContinously <- vapply(
    strList,
    function(str){
      if(length(str[str$reverse != 0,]$reverse) == 0) return(TRUE)
      .is_continous_evenly_spaced(str[str$reverse != 0,]$reverse)
    },
    logical(1))
  forwardLength <- vapply(
    strList,
    function(str){
      ifelse(length(str[str$forward != 0,]$forward) > 0,
             max(str[str$forward != 0,]$forward) -
               min(str[str$forward != 0,]$forward) + 1,
             0)
    },
    double(1))
  reverseLength <- vapply(
    strList,
    function(str){
      ifelse(length(str[str$reverse != 0,]$reverse) > 0,
             max(str[str$reverse != 0,]$reverse) -
               min(str[str$reverse != 0,]$reverse) + 1,
             0)
    },
    double(1))
  return(list(forwardContinously = unname(forwardContinously),
              reverseContinously = unname(reverseContinously),
              forwardLength = unname(forwardLength),
              reverseLength = unname(reverseLength)))
}

.get_mismatches <- function(gr, ident, strList){
  if(!missing(gr) && 
     !missing(ident) &&
     missing(strList)){
    strList <- .get_ident_structures(ident, gr)[[1]]
  }
  strListDetails <- .get_structures_details(strList)
  ansMismatchesFalse <- strListDetails[["forwardLength"]] == 
    strListDetails[["reverseLength"]] &
    strListDetails[["forwardContinously"]] &
    strListDetails[["reverseContinously"]]
  ansMismatchesTrue <- strListDetails[["forwardLength"]] == 
    strListDetails[["reverseLength"]] &
    ((strListDetails[["forwardContinously"]] & 
        !strListDetails[["reverseContinously"]]) |
       (!strListDetails[["forwardContinously"]] & 
          strListDetails[["reverseContinously"]]) |
       (!strListDetails[["forwardContinously"]] & 
          !strListDetails[["reverseContinously"]]))
  return(list("T" = unname(ansMismatchesTrue),
              "F" = unname(ansMismatchesFalse)))
}

.get_bulges <- function(gr,
                        ident,
                        strList){
  if(!missing(gr) && 
     !missing(ident) &&
     missing(strList)){
    strList <- .get_ident_structures(ident,
                                     gr)[[1]]
  }
  strListDetails <- .get_structures_details(strList)
  ansBulgedFalse <- strListDetails[["forwardLength"]] == 
    strListDetails[["reverseLength"]] &
    ((strListDetails[["forwardContinously"]] & 
        strListDetails[["reverseContinously"]]) |
       (!strListDetails[["forwardContinously"]] & 
          !strListDetails[["reverseContinously"]]))
  ansBulgedTrue <- strListDetails[["forwardLength"]] != 
    strListDetails[["reverseLength"]] &
    ((strListDetails[["forwardContinously"]] & 
        !strListDetails[["reverseContinously"]]) |
       (!strListDetails[["forwardContinously"]] & 
          strListDetails[["reverseContinously"]]) |
       (!strListDetails[["forwardContinously"]] & 
          !strListDetails[["reverseContinously"]]))
  return(list("T" = unname(ansBulgedTrue),
              "F" = unname(ansBulgedFalse)))
}

.get_str_lengths <- function(gr, ident, strList){
  if(missing(strList)){
    strList <- .get_ident_structures(ident, gr)[[1]]
  }
  if(!(ident %in% TRNA_STRUCTURES_LOOP)){
    length <- lapply(strList,
                     function(str){
                       if(nrow(str) == 0) return(NA)
                       max(max(str[str$forward != 0,]$forward) -
                             min(str[str$forward != 0,]$forward),
                           max(str[str$reverse != 0,]$reverse) -
                             min(str[str$reverse != 0,]$reverse)) + 1
                     })
  } else {
    length <- lapply(strList,
           function(str){
             if(nrow(str) == 0) return(NA)
             max(str$pos) - 
               min(str$pos) + 1
           })
  }
  length <- unlist(length)
  return(unname(length))
}
