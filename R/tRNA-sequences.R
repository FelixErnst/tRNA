#' @include tRNA.R
NULL

#' @rdname gettRNAstructureSeqs
#'
#' @importFrom stringr str_locate
#' @importFrom BiocGenerics width
#' @importFrom IRanges reverse
#' @importFrom XVector subseq
#' @importFrom Biostrings DNAStringSet
#'
#' @export
setMethod(
  f = "gettRNAstructureSeqs",
  signature = signature(x = "GRanges"),
  definition = function(x, structure, joinCompletely = FALSE,
                        joinFeatures = FALSE, padSequences = TRUE) {
    # input check
    .check_trna_granges(x, TRNA_FEATURES)
    .check_trna_structure_ident(structure)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    if(!.is_a_bool(joinCompletely)){
      stop("'joinCompletely' must TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(joinFeatures)){
      stop("'joinFeatures' must TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(padSequences)){
      stop("'padSequences' must TRUE or FALSE.", call. = FALSE)
    }
    if(joinCompletely == TRUE && joinCompletely == joinFeatures){
      warning("Both 'joinCompletely' and 'joinFeatures' are set to TRUE.
              'joinCompletely' takes precedence.", call. = FALSE)
    }
    # join completly or get splitup sequences
    if(joinCompletely){
      # get Ranges
      strList <- getBasePairing(x$tRNA_str)
      res <- .get_tRNA_structures(names(tRNAStructureFunctionList),
                                  x,
                                  strList)
      # get sequences
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(x,
                                     joinFeatures = FALSE,
                                     padSequences = TRUE,
                                     strList))
      # assemble boundaries IRanges
      z <- unlist(seqs)[TRNA_STRUCTURE_ORDER]
      ir <- lapply(lapply(unlist(seqs)[TRNA_STRUCTURE_ORDER],
                          BiocGenerics::width),
                   unique)
      start <- c(1,unlist(ir[seq_len((length(ir) - 1))]))
      end <- unlist(ir)
      ir <- IRanges(start = unlist(lapply(seq_along(start),
                                          function(i){
                                            sum(start[seq_len(i)])
                                          })),
                    end = unlist(lapply(seq_along(end),
                                        function(i){
                                          sum(end[seq_len(i)])
                                        })))
      names(ir) <- names(end)
      # concat sequences
      seqs <- do.call(Biostrings::xscat,
                      z)
      # store boundaries as metadata
      S4Vectors::metadata(seqs) <- list("tRNA_structures" = ir)
    } else {
      # get Ranges
      strList <- gettRNABasePairing(x)
      res <- .get_tRNA_structures(structure, x, strList)
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(x,
                                     joinFeatures,
                                     padSequences,
                                     strList))
    }
    return(seqs)
    }
)

################################################################################
# join sequences of tRNA with correct padding
.assemble_sequences <- function(ir, name, gr, joinFeatures, padSequences,
                                strList){
  seqs <- gr$tRNA_seq
  ans_class_fun <- match.fun(class(seqs))
  ##############################################################################
  # if it is a stem and features should be not be joined nor padded
  ##############################################################################
  if(is.list(ir) && !joinFeatures && !padSequences){
    ans <- mapply(.assemble_sequences,
                  ir,
                  MoreArgs = list(name,
                                  gr,
                                  joinFeatures,
                                  padSequences,
                                  strList))
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should be joined, but not padded
  ##############################################################################
  if(is.list(ir) && joinFeatures && !padSequences){
    prime5 <- XVector::subseq(seqs, ir$prime5)
    prime3 <- XVector::subseq(seqs, ir$prime3)
    ans <- Biostrings::xscat(
      prime5,
      prime3
    )
    names(ans) <- names(ir$prime5)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should not be joined, but padded
  ##############################################################################
  if(is.list(ir) && !joinFeatures && padSequences){
    x <- .pad_unpaired_in_stem_region(seqs, ir, strList, ans_class_fun)
    prime5 <- .pad_right(x$prime5, ans_class_fun)
    prime3 <- .pad_left(x$prime3, ans_class_fun)
    ans <- list(prime5,
                prime3)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should be joined and padded
  ##############################################################################
  if(is.list(ir) && joinFeatures){
    return(.join_list(seqs, ir, strList, ans_class_fun))
  }
  ##############################################################################
  # special cases below
  ##############################################################################
  ##############################################################################
  # if padding should be done in the center
  ##############################################################################
  if(name == "center" && padSequences){
    ans <- .pad_center(XVector::subseq(seqs, ir),ans_class_fun)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if padding should be done left or right
  ##############################################################################
  if(name == "left" && padSequences){
    ans <- .pad_left(XVector::subseq(seqs, ir), ans_class_fun)
    names(ans) <- names(ir)
    return(ans)
  }
  if(name == "right" && padSequences){
    ans <- .pad_right(XVector::subseq(seqs, ir), ans_class_fun)
    names(ans) <- names(ir)
    return(ans)
  }
  #############################################
  # special rule for Dloop
  #############################################
  if(name == "Dloop" && padSequences){
    return(.pad_Dloop(seqs, ir, gr, strList))
  }
  #############################################
  # special rule for variable loop
  #############################################
  if(name == "variableLoop" && padSequences){
    return(.pad_variableLoop(seqs, ir, gr, strList, ans_class_fun))
  }
  ##############################################################################
  # if it is a loop and should not be padded
  ##############################################################################
  if(!is.list(ir) && !padSequences){
    ans <- XVector::subseq(seqs, ir)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a loop and should be padded
  ##############################################################################
  ans <- .pad_outside(XVector::subseq(seqs, ir),ans_class_fun)
  names(ans) <- names(ir)
  return(ans)
}

#
.join_list <- function(seqs, ir, strList, ans_class_fun){
  x <- .pad_unpaired_in_stem_region(seqs, ir, strList, ans_class_fun)
  prime5 <- .pad_right(x$prime5, ans_class_fun)
  prime3 <- .pad_left(x$prime3, ans_class_fun)
  maxWidth <- max(BiocGenerics::width(prime5) + BiocGenerics::width(prime3))
  addN <- maxWidth - 
    (BiocGenerics::width(prime5) + 
       BiocGenerics::width(prime3))
  addString <- ans_class_fun(rep("--",length(seqs)))
  addString <- Biostrings::xscat(
    addString,
    ans_class_fun(unlist(lapply(addN,
                                function(n){
                                  do.call(paste0,as.list(rep("-",(n + 1))))
                                }))
    )
  )
  ans <- Biostrings::xscat(
    prime5,
    addString,
    prime3
  )
  names(ans) <- names(ir[[1]])
  return(ans)
}

#
.pad_Dloop <- function(seqs, ir, gr, strList){
  # three at 5'-end and one 3'-end
  # search for GGG, GG, GC, AG, AC, AT, TT, CT
  # if found put in der middle and the rest at the 3'-end
  # if not put everthing at the 5'-end
  GGfix <- c("GGG","GG","GC","AG","AC","AT","TT","CT")
  #
  getGGPos <- function(seqs, ir, i, searchString){
    if(is.null(searchString[i])) return(NULL)
    x <- stringr::str_locate(reverse(XVector::subseq(seqs, ir)),
                             reverse(searchString[i]))
    x <- BiocGenerics::start(ir) +
      width(reverse(XVector::subseq(seqs, ir))) -
      x[,"end"]
    names <-  rep(searchString[i], length(ir))
    names[is.na(x)] <- ""
    if(any(is.na(x))){
      y <- getGGPos(seqs[is.na(x)],
                    ir[is.na(x)],
                    (i + 1),
                    searchString)
      if(!is.null(y)){
        names[is.na(x)] <- names(y)
        x[is.na(x)] <- y
      }
    }
    names(x) <- names
    return(x)
  }
  #
  ir5 <- ir
  BiocGenerics::end(ir5) <- BiocGenerics::start(ir5) + 1
  ir3 <- ir
  BiocGenerics::start(ir3) <- BiocGenerics::end(ir3)
  irM <- ir
  BiocGenerics::start(irM) <- BiocGenerics::start(irM) + 2
  BiocGenerics::end(irM) <- BiocGenerics::end(irM) - 1
  #
  GGpos <- getGGPos(seqs, ir, 1, GGfix)
  BiocGenerics::start(irM[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)]
  BiocGenerics::end(ir5[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)] - 1
  # split in the middle if no GG like found
  if(any(is.na(GGpos))){
    BiocGenerics::start(irM[is.na(GGpos)]) <-
      BiocGenerics::start(irM[is.na(GGpos)]) +
      floor(BiocGenerics::width(irM[is.na(GGpos)])/2)
    BiocGenerics::end(ir5[is.na(GGpos)]) <-
      BiocGenerics::end(ir5[is.na(GGpos)]) +
      floor(BiocGenerics::width(irM[is.na(GGpos)])/2)
  }
  #
  ans <- Biostrings::xscat(
    .assemble_sequences(ir5,
                        "right",
                        gr,
                        joinFeatures = FALSE,
                        padSequences = TRUE,
                        strList),
    .assemble_sequences(irM,
                        "right",
                        gr,
                        joinFeatures = FALSE,
                        padSequences = TRUE,
                        strList),
    .assemble_sequences(ir3,
                        "left",
                        gr,
                        joinFeatures = FALSE,
                        padSequences = TRUE,
                        strList)
  )
  #
  names(ans) <- names(ir)
  return(ans)
}

#
.pad_variableLoop <- function(seqs, ir, gr, strList, ans_class_fun){
  strs <- gr$tRNA_str
  # circumvent working on empty variable loops
  zero_width <- BiocGenerics::width(ir) == 0L
  if(any(zero_width)){
    ir_orig <- ir
    seqs_orig <- seqs
    ir <- ir[!zero_width]
    seqs <- seqs[!zero_width]
  }
  # prepad loops with width == 1L
  one_width <- BiocGenerics::width(ir) == 1L
  if(any(one_width)){
    subseq(seqs[one_width], ir[one_width]) <- 
      Biostrings::xscat(subseq(seqs[one_width], ir[one_width]),
                        ans_class_fun(rep("-",length(which(one_width)))))
    end(ir[one_width]) <- end(ir[one_width]) + 1L
  }
  # keep one pos at 5'- and 3'-end
  ir5 <- ir
  BiocGenerics::end(ir5) <- BiocGenerics::start(ir5)
  prime5 <- XVector::subseq(seqs, ir5)
  ir3 <- ir
  BiocGenerics::start(ir3) <- BiocGenerics::end(ir3)
  prime3 <- XVector::subseq(seqs, ir3)
  # get the middle sequnce
  irM <- ir
  BiocGenerics::start(irM) <- BiocGenerics::start(irM) + 1L
  BiocGenerics::end(irM) <- BiocGenerics::end(irM) - 1L
  middle <- XVector::subseq(seqs, irM)
  # if longer variable loops exist search for paired regions
  facLong <- BiocGenerics::width(ir) > 4
  facBasePaired <- rep(FALSE, length(ir))
  if(any(facLong)){
    startPaired5 <- BiocGenerics::start(irM) - 1 +
      stringr::str_locate(substr(strs[!zero_width],
                                 BiocGenerics::start(irM),
                                 BiocGenerics::end(irM)),
                          ">")[,"start"]
    facBasePaired <- !is.na(startPaired5)
    # if paired regions exist pad them accordingly
    if(any(facBasePaired)){
      startPaired5 <- BiocGenerics::start(irM[facBasePaired])
      # test for different loop types
      endPaired5 <- BiocGenerics::start(irM[facBasePaired]) - 1 +
        stringr::str_locate(substr(strs[!zero_width][facBasePaired],
                                   BiocGenerics::start(irM[facBasePaired]),
                                   BiocGenerics::end(irM[facBasePaired])),
                            "<\\.\\.|<\\.>|<>")[,"start"]
      startPaired3 <- BiocGenerics::start(irM[facBasePaired]) + 1 +
        stringr::str_locate(substr(strs[!zero_width][facBasePaired],
                                   BiocGenerics::start(irM[facBasePaired]),
                                   BiocGenerics::end(irM[facBasePaired])),
                            "\\.\\.>|<<>|\\.<>")[,"start"]
      endPaired3 <- BiocGenerics::end(irM[facBasePaired])
      # assemble middle sequences
      stem <- list(prime5  = IRanges::IRanges(start = startPaired5,
                                              end = endPaired5),
                   prime3  = IRanges::IRanges(start = startPaired3,
                                              end = endPaired3))
      stem <- .assemble_sequences(stem,
                                  "variableloopstem",
                                  gr[!zero_width][facBasePaired],
                                  joinFeatures = FALSE,
                                  padSequences = TRUE,
                                  strList[facBasePaired])
      m <- XVector::subseq(seqs[facBasePaired], irM[facBasePaired])
      # proceed with sequences which have a loop
      f <- endPaired5 < startPaired3 - 1
      m[f] <- .assemble_sequences(IRanges::IRanges(start = endPaired5[f] + 1,
                                                   end = startPaired3[f] - 1),
                                  "right",
                                  gr[!zero_width][facBasePaired][f],
                                  joinFeatures = FALSE,
                                  padSequences = TRUE,
                                  strList[facBasePaired][f])
      # add spacer for missing loops
      addWidth <- rep(max(BiocGenerics::width(m[f])),length(m[!f]))
      m[!f] <- ans_class_fun(unlist(lapply(addWidth,
                                           function(n){
                                             do.call(paste0,as.list(rep("-",n)))
                                           }))
      )
      # combien everything
      middle[facBasePaired] <- Biostrings::xscat(stem$prime5,m,stem$prime3)
    }
  }
  # pad non paired region sequences in the middle
  maxWidth <- max(BiocGenerics::width(middle))
  addNMiddle <- maxWidth - BiocGenerics::width(middle)
  addMiddle <- ans_class_fun(unlist(lapply(addNMiddle,
                                           function(n){
                                             do.call(paste0,as.list(rep("-",n)))
                                           }))
  )
  middle[addNMiddle > 0] <- Biostrings::xscat(
    XVector::subseq(middle[addNMiddle > 0],
                    1,
                    ceiling(width(middle[addNMiddle > 0])/2)),
    addMiddle,
    XVector::subseq(middle[addNMiddle > 0],
                    (ceiling(width(middle[addNMiddle > 0])/2) + 1),
                    width(middle[addNMiddle > 0]))
  )
  # assemble left right and middle sequences
  ans <- Biostrings::xscat(
    .assemble_sequences(ir5,
                        "right",
                        gr[!zero_width],
                        joinFeatures = FALSE,
                        padSequences = TRUE,
                        strList),
    middle,
    .assemble_sequences(ir3,
                        "left",
                        gr[!zero_width],
                        joinFeatures = FALSE,
                        padSequences = TRUE,
                        strList)
  )
  if(any(zero_width)){
    ir_orig[!zero_width] <- ir
    seqs_orig[!zero_width] <- ans
    width <- max(BiocGenerics::width(ans))
    seqs_orig[zero_width] <- 
      ans_class_fun(vapply(seq_along(which(zero_width)),
                           function(i){
                             paste(rep("-",width), collapse = "")
                           },
                           character(1)))
    ir <- ir_orig
    ans <- seqs_orig
  }
  names(ans) <- names(ir)
  return(ans)
}

#
.pad_left <- function(ans, ans_class_fun){
  maxWidth <- max(BiocGenerics::width(ans))
  addNLeft <- maxWidth - BiocGenerics::width(ans)
  if(any(addNLeft > 0)){
    addLeft <- ans_class_fun(unlist(lapply(addNLeft,
                                           function(n){
                                             do.call(paste0,as.list(rep("-",n)))
                                           }))
    )
    ans[addNLeft > 0] <- Biostrings::xscat(
      addLeft,
      ans[addNLeft > 0]
    )
  }
  return(ans)
}

#
.pad_right <- function(ans, ans_class_fun){
  maxWidth <- max(BiocGenerics::width(ans))
  addNRight <- maxWidth - BiocGenerics::width(ans)
  if(any(addNRight > 0)){
    addRight <- 
      ans_class_fun(unlist(lapply(addNRight,
                                  function(n){
                                    do.call(paste0,as.list(rep("-",n)))
                                  }))
                          
    )
    ans[addNRight > 0] <- Biostrings::xscat(
      ans[addNRight > 0],
      addRight
    )
  }
  return(ans)
}

#
.pad_center <- function(ans, ans_class_fun){
  maxWidth <- max(BiocGenerics::width(ans))
  addNMiddle <- maxWidth - BiocGenerics::width(ans)
  if(any(addNMiddle > 0)){
    addMiddle <- 
      ans_class_fun(unlist(lapply(addNMiddle,
                                  function(n){
                                    do.call(paste0,as.list(rep("-",n)))
                                  }))
                           
    )
    ans[addNMiddle > 0] <- Biostrings::xscat(
      XVector::subseq(ans[addNMiddle > 0],
                      1,
                      ceiling(width(ans[addNMiddle > 0])/2)),
      addMiddle,
      XVector::subseq(ans[addNMiddle > 0],
                      (ceiling(width(ans[addNMiddle > 0])/2) + 1),
                      width(ans[addNMiddle > 0]))
    )
  }
  return(ans)
}

#
.pad_outside <- function(ans, ans_class_fun){
  maxWidth <- max(BiocGenerics::width(ans))
  missingWidth <- maxWidth - BiocGenerics::width(ans)
  # if sequences should be padded not in the center but outside
  addNLeft <- floor(missingWidth / 2)
  addNRight <- ceiling(missingWidth / 2)
  f <- addNRight > 0 & (missingWidth %% 2) == 1
  addNLeft[f] <- addNLeft[f] + 1
  addNRight[f] <- addNRight[f] - 1
  if(any(addNLeft > 0)){
    addLeft <- ans_class_fun(unlist(lapply(addNLeft,
                                           function(n){
                                             do.call(paste0,as.list(rep("-",n)))
                                           }))
    )
    ans[addNLeft > 0] <- Biostrings::xscat(
      addLeft,
      ans[addNLeft > 0]
    )
  }
  if(any(addNRight > 0)){
    addRight <- 
      ans_class_fun(unlist(lapply(addNRight,
                                  function(n){
                                    do.call(paste0,as.list(rep("-",n)))
                                  }))
    )
    ans[addNRight > 0] <- Biostrings::xscat(
      ans[addNRight > 0],
      addRight
    )
  }
  return(ans)
}

# detect bulges in stem region and pad on opposite site
.pad_unpaired_in_stem_region <- function(seqs, ir, strList, ans_class_fun){
  prime5 <- XVector::subseq(seqs, ir$prime5)
  prime3 <- XVector::subseq(seqs, ir$prime3)
  if(length(unique(BiocGenerics::width(prime5))) > 1){
    f <- max(BiocGenerics::width(prime5)) > BiocGenerics::width(prime5)
    prime5[f] <- .add_padding_unpaired(prime5[f],
                                       ir$prime5[f],
                                       strList[f],
                                       ans_class_fun)
  }
  if(length(unique(BiocGenerics::width(prime3))) > 1){
    f <- max(BiocGenerics::width(prime3)) > BiocGenerics::width(prime3)
    prime3[f] <- .add_padding_unpaired(prime3[f],
                                       ir$prime3[f],
                                       strList[f],
                                       ans_class_fun)
  }
  return(list(prime5 = prime5, prime3 = prime3))
}

#
.add_padding_unpaired <- function(seqs, ir, strList, ans_class_fun){
  dims <- lapply(
    seq_along(strList),
    function(i){
      # subset to relevant structure
      z <- strList[[i]][strList[[i]]$pos %in% 
                          BiocGenerics::start(ir[i]):BiocGenerics::end(ir[i]),]
      z <- z[z$reverse != 0,]
      # check for bulge
      zz <- vapply(
        seq_len(nrow(z)),
        function(j){
          # skip last
          if(j >= nrow(z)) return(FALSE)
          # missing pos on reverse
          z[j,]$reverse - 1 > z[j + 1,]$reverse &
            # not the same bulge on the other side
            (z[j,]$reverse - 1 - z[j + 1,]$reverse) != 
            (z[j + 1,]$forward - 1 - z[j,]$forward)
        },
        logical(1))
      # if not unpaired position can be detected it is just a shorter
      # stem
      w_zz <- which(zz)
      if(length(w_zz) == 0){
        return(NULL)
      }
      length <- z[zz,]$reverse - z[w_zz + 1,]$reverse - 1L -
        (z[w_zz + 1,]$forward - z[w_zz,]$forward - 1)
      # if bulge is shorter in the other side, skip
      length_neg <- length <= 0L
      if(all(length_neg)){
        return(NULL)
      }
      start <- z[zz,]$forward + 1 - BiocGenerics::start(ir[i])
      stop <- z[w_zz + 1,]$forward + 1L - BiocGenerics::start(ir[i]) -
        (z[w_zz + 1,]$forward - z[w_zz,]$forward - 1)
      return(list(start = start[!length_neg], 
                  stop = stop[!length_neg], 
                  length = length[!length_neg]))
    })
  need_editing <- !vapply(dims,is.null,logical(1))
  if(!any(need_editing)){
    return(seqs)
  }
  dims <- dims[need_editing]
  ans_class_fun_element <- match.fun(S4Vectors::elementType(seqs))
  seqs[need_editing] <- mapply(.add_padding_to_pos,
                               seqs[need_editing],
                               lapply(dims,"[[","start"),
                               lapply(dims,"[[","stop"),
                               lapply(dims,"[[","length"),
                               MoreArgs = list(ans_class_fun_element))
  return(seqs)
}

#
.add_padding_to_pos <- function(seq, start, stop, length, ans_class_fun){
  length_add <- c(0L,length[seq_len(length(length) - 1L)])
  start <- start + length_add
  stop <- stop + length_add
  add <- lapply(length,
                function(n){
                  ans_class_fun(paste(rep("-",n),collapse = ""))
                })
  # for-loop needed since multiple addition per sequence can occur
  for(i in seq_along(add)){
    seq <- Biostrings::xscat(
      XVector::subseq(seq, start = 1, end = start[i]),
      add[[i]],
      XVector::subseq(seq, start = stop[i], end = length(seq))
    )
  }
  seq
}
