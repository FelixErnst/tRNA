#' @include tRNA.R
NULL

#' @name getBasePairing
#' @aliases getBasePairing gettRNABasePairing gettRNALoopIDs getLoopIDs
#'
#' @title Accessing Dot Bracket annotation
#'
#' @description
#' \code{getBasePairing} converts a dot bracket annotation into a
#' \code{data.frame}. Base pairing is indicated by corresponding numbers
#' in the forward and reverse columns.
#' 
#' \code{getHairpinLoops} and \code{getStems} are helper functions for getting
#' coordinates for stems and hairpin loops as \code{IRanges} objects.
#'
#' @return
#' \code{getBasePairing}: 
#' The result is a data.frame with following columns: pos, forward, reverse, chr
#' and base (if \code{sequence} was provided or a \code{GRanges} objects was 
#' used). If a position is unpaired, forward and reverse will be \code{0}, 
#' otherwise it will match the base paired positions.
#' \code{gettRNALoopIDs}, \code{getLoopIDs}: 
#' return a list of list of loop ids.
#'
#' @param x a GRanges object created by \code{import.tRNAscanAsGRanges} or
#' GRanges with equivalent information. The \code{tRNA_str} and \code{tRNA_seq} 
#' columns will be used for input into \code{getBasePairing}.
#'
#' @importFrom stringr str_locate_all str_locate
#' 
#' @examples 
#' data("gr", package = "tRNA", envir = environment())
#' gettRNABasePairing(gr[1])
#' getBasePairing(gr[1]$tRNA_str)
#' gettRNALoopIDs(gr[1])
#' getLoopIDs(gr[1]$tRNA_str)
NULL

STRUCTURE_OPEN_CHR <- c("<","\\[","\\(","\\{")
STRUCTURE_CLOSE_CHR <- c(">","\\]","\\)","\\}")

#' @rdname getBasePairing
#'
#' @export
setMethod(
  f = "gettRNABasePairing",
  signature = signature(x = "GRanges"),
  definition = function(x) {
    .check_trna_granges(x, TRNA_FEATURES)
    .get_base_pairing(x$tRNA_str,
                      as.character(x$tRNA_seq))
  }
)

#' @rdname getBasePairing
#'
#' @param dotBracket character vectors describing a nucleotide sequence
#' structure in the dot bracket annotations. Valid characters are:
#' \code{.(\\{[><]\\})}
#' @param sequence optional: character vectors describing a nucleotide sequence.
#' The same number of sequences with the same length as the dot bracket string 
#' have to be used. Each nucleotide sequence has to be a character vector.
#' The identity of the nucleotides are not control, so in theory all letters 
#' can be used.
#'
#' @export
getBasePairing <- function(dotBracket,
                           sequence){
  assertive::assert_is_non_empty(dotBracket)
  dotBracket <- unlist(dotBracket)
  .check_dot_bracket(dotBracket)
  assertive::assert_all_are_non_missing_nor_empty_character(dotBracket)
  if(!missing(sequence)){
    assertive::assert_is_non_empty(sequence)
    sequence <- unlist(sequence)
    assertive::assert_all_are_non_missing_nor_empty_character(sequence)
    # check that the number and the length of sequences matches those of the
    # dot bracket annotation
    if(length(sequence) != length(dotBracket)){
      stop("Number of dot bracket strings and sequences do not match.",
           call. = FALSE)
    }
    check <- unlist(mapply(
      function(d,s){
        nchar(d) != nchar(s)
      },
      dotBracket,sequence))
    if(any(check)){
      stop("Length of dot bracket annotation and sequence do not match for ",
           "the following pairs:\n'",
           paste(which(check),collapse = "', '"),
           "'.",
           call. = FALSE)
    }
  } else {
    sequence <- rep(NA,length(dotBracket))
  }
  ans <- .get_base_pairing(dotBracket,
                           sequence)
  return(ans)
}

# convert dot bracket annotation in ct like format
.get_base_pairing <- function(x,
                              seq){
  # special case for <>. This is the usually used orientation (eg. by ViennaRNA)
  # , but tRNAscan files use a different orientation. If the first occurance is 
  # < switch out the orientation.
  f <- which(stringr::str_locate(x, ">")[,"start"] < 
               stringr::str_locate(x, "<")[,"start"])
  if(length(f) > 0){
    tmp <- gsub("<","a",x[f])
    tmp <- gsub(">","b",tmp)
    tmp <- gsub("a",">",tmp)
    tmp <- gsub("b","<",tmp)
    x[f] <- tmp
  }
  open <- lapply(STRUCTURE_OPEN_CHR,
                 function(chr){
                   stringr::str_locate_all(x, chr)
                 })
  close <- lapply(STRUCTURE_CLOSE_CHR,
                  function(chr){
                    stringr::str_locate_all(x, chr)
                  })
  lengthOpen <- lapply(open,function(z){lapply(z,length)})
  lengthClose <- lapply(close,function(z){lapply(z,length)})
  lengthMatch <- lapply(
    seq_along(lengthOpen),
    function(i){
      which(unlist(lengthOpen[[i]]) != unlist(lengthClose[[i]]))
    })
  # check for unmatched positions
  if(any(unlist(lapply(lengthMatch,length)) != 0)){
    stop("Following structures are invalid: \n'",
         paste(unique(unlist(lengthMatch)),
               collapse = "', '"),
         "'.\nThey contain unmatched positions.",
         call. = FALSE)
  }
  structure <- mapply(.get_base_pairing_data_frame,
                      open,
                      close,
                      STRUCTURE_OPEN_CHR,
                      STRUCTURE_CLOSE_CHR)
  structure <- split(structure,seq_len(length(x)))
  # check if any positions are unmatched due to orientation
  check <- vapply(structure,
                  function(z){
                    any(vapply(z,
                               function(zz){
                                 any(is.na(zz$forward))
                               },
                               logical(1)))
                  },
                  logical(1))
  if(any(check)){
    stop("Following structures are invalid: \n'",
         paste(unique(unlist(which(check))),
               collapse = "', '"),
         "'.\nThe order of the opening and closing characters is wrong.",
         call. = FALSE)
  }
  #
  structure <- mapply(.complete_base_pairing_data_frame,
                      structure,
                      lapply(x,nchar),
                      seq,
                      SIMPLIFY = FALSE)
  return(structure)
}
# assembles base pairing data.frame
.get_base_pairing_data_frame <- function(open,
                                         close,
                                         opchr,
                                         clchr){
  ident <- gsub("\\\\","",paste0(opchr,clchr))
  ans <- mapply(
    function(op,cl,id){
      op <- op[,"start"]
      cl <- cl[,"start"]
      if(length(cl) > 0){
        forward <- rep(0, length(cl))
        for(j in seq_along(cl)){
          forward[j] <- rev(op[op < cl[j]])[1]
          op <- op[op != forward[j]]
        }
        ans <- data.frame(forward = forward,
                          reverse = cl,
                          chr = id)
        return(ans)
      }
      return(NULL)
    },
    open,
    close,
    MoreArgs = list(ident),
    SIMPLIFY = FALSE)
  ans
}
# add missing value to base pairing data.frame
.complete_base_pairing_data_frame <- function(z,
                                              n,
                                              s){
  z <- do.call(rbind,z[!vapply(z,is.null,logical(1))])
  z2 <- z[,c("reverse","forward","chr")]
  colnames(z2) <- colnames(z)
  z <- rbind(z,z2)
  z$pos <- z$forward
  missing <- seq_len(n)
  missing <- missing[!(missing %in% z$forward)]
  if(length(missing) > 0){
    z <- rbind(z,
               data.frame(pos = missing,
                          forward = rep(0,length(missing)),
                          reverse = rep(0,length(missing)),
                          chr = rep(".",length(missing))))
  }
  z <- z[order(z$pos),c("pos","forward","reverse","chr")]
  # add sequence if not NA
  if(!is.na(s)){
    z$base <- strsplit(s,"")[[1]]
  }
  #
  rownames(z) <- NULL
  return(z)
}

.get_hairpin_loops <- function(strList){
  tmp <- lapply(strList,
                function(z){
                  ans <- z[z$forward != 0,]
                  if(nrow(ans) == 0) return(list(NA,NA))
                  f <- c(ans$reverse,NA)
                  f <- f[2:length(f)]
                  ans <- ans[ans$forward == f & !is.na(f),]
                  if(nrow(ans) == 0) return(list(NA,NA))
                  return(list(start = ans$forward + 1,
                              end = ans$reverse - 1))
                })
  ans <- lapply(tmp,
                function(z){
                  .get_IRanges2(z$start,
                                z$end,
                                seq_along(z$start))
                })
  return(ans)
}

#' @rdname getBasePairing
#' @export
setMethod(
  f = "gettRNALoopIDs",
  signature = signature(x = "GRanges"),
  definition = function(x) {
    .check_trna_granges(x, TRNA_FEATURES)
    .get_ids_of_loops(getBasePairing(x$tRNA_str))
  }
)

#' @rdname getBasePairing
#' @export
getLoopIDs <- function(dotBracket){
  strList <- getBasePairing(dotBracket)
  .get_ids_of_loops(strList)
}


.get_ids_of_loops <- function(strList){
  ans <- mapply(.get_loop_ids,
                strList,
                SIMPLIFY = FALSE)
  return(ans)
}

# z = the base pairing table with forward, reverse and chr column
.get_loop_ids_c <- function(z){
  len <- nrow(z)
  nl <- 0
  l <- 0
  hx <- 1
  stack <- c()
  loop <- c()
  for(i in seq_len(len)){
    # opening
    if(z[i,"forward"] != 0 && 
       i < z[i,"reverse"] &&
       z[i,"chr"] != "."){
      nl <- nl + 1
      l <- nl
      stack[hx] <- i
      hx <- hx + 1
    }
    loop[i] <- l
    # closing
    if(z[i,"forward"] != 0 && 
       i > z[i,"reverse"] &&
       z[i,"chr"] != "."){
      hx <- hx - 1
      if(hx > 1){
        l <- loop[stack[hx - 1]]
      } else {
        l <- 0
      }
      if(hx < 1) return(NA)
    }
  }
  return(loop)
}
.get_loop_ids <- compiler::cmpfun(.get_loop_ids_c)