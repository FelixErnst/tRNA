#' @include tRNA.R
NULL

#' @name gettRNAstructureSeqs
#' @aliases gettRNAstructureSeqs gettRNAstructureGRanges
#'
#' @title tRNA structures and sequences
#'
#' @description
#' \code{gettRNAstructureGRanges} returns a list of GRanges describing the
#' boundaries of tRNA structures as extracted from the dot bracket annotation.
#' The dot bracket annotation is parsed using \code{gettRNABasePairing}, which
#' internally uses \code{getBasePairing}.
#'
#' \code{gettRNAstructureSeq} returns split or partial tRNA sequences based on
#' the structure information of tRNAscan. Variances in length of certain
#' structure features can be padded. If sequences are joined by setting
#' \code{joinCompletely = FALSE}, the boundaries of the tRNA structure are
#' stored in the result as metadata. They can be accessesed as an IRanges
#' object by using \code{metadata()[["tRNA_structures"]]}.
#'
#' @param gr a GRanges object with tRNA information. It has to pass the
#' \code{istRNAGRanges} function.
#' @param structure optional parameter for returning just partial sequences.
#' The following values are accepted:
#' \code{anticodonStem}, \code{Dprime5}, \code{DStem}, \code{Dloop},
#' \code{Dprime3}, \code{acceptorStem}, \code{anticodonloop},
#' \code{variableLoop}, \code{TStem}, \code{Tloop}, \code{discriminator}.
#' (default: \code{structure = ""})
#' @param joinCompletely Should the sequence parts, which are to be returned, be
#' joined into one sequence? (default: \code{joinCompletely = FALSE}))
#' Setting this to TRUE excludes \code{joinFeatures} be set to TRUE as well. In
#' addition, \code{joinCompletely = TRUE} uses automatically all sequence
#' structures.
#' @param joinFeatures Should the sequence parts, which are to be returned and
#' are from the same structure type, be joined into one sequence?
#' (default: \code{joinCompletely = FALSE})) Setting this to TRUE excludes
#' \code{joinCompletely} be set to TRUE as well. \code{joinCompletely} takes
#' precedence.
#' @param padSequences parameter whether sequences of the same type
#' should be returned with the same length. For stems missing positions will be
#' filled up in the middle, for loops at the ends.
#' (default: \code{padSequences = TRUE}). If \code{joinCompletely == TRUE} this
#' is set to TRUE automatically.
#'
#' @return a list of \code{GRanges} or \code{DNAStringSet} objects. In case
#' joinCompletly is set to TRUE a single \code{DNAStringSet} is returned.
#'
#' @export
#' @examples
#' data("gr", package = "tRNA", envir = environment())
#' gettRNAstructureGRanges(gr, structure = "anticodonLoop")
#' gettRNAstructureSeqs(gr, structure = "anticodonLoop")
#' gettRNABasePairing(gr[1:10])
#' getBasePairing(gr[1:10]$tRNA_str)
NULL

#' @rdname gettRNAstructureSeqs
#'
#' @importFrom IRanges IRanges
#'
#' @export
setMethod(
  f = "gettRNAstructureGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        structure) {
    # input check
    .check_trna_granges(gr, TRNA_FEATURES)
    .check_trna_structure_ident(structure)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    #
    strList <- gettRNABasePairing(gr)
    res <- .get_tRNA_structures(structure, gr, strList)
    return(res)
  }
)

# positions for each requested feature are returned. If one of the features is
# not recognized, NA is returned for one or both start/end. In this case the
# width is zero, which will result in downstream recognition of the missing
# feature
.get_tRNA_structures <- function(structure,
                                 gr,
                                 strList){
  # get start data
  ans <- list()
  loopPositions <- .get_loop_positions(strList)
  tRNAStructureTypes <- .get_tRNA_structure_type(loopPositions)
  # discriminator stem
  discriminator <- .getDiscriminator(gr,
                                     strList,
                                     tRNAStructureTypes)
  ans[["discriminator"]] <- .get_IRanges(discriminator,
                                         gr$tRNA_anticodon)
  if(all(structure %in% names(ans))){
    return(ans[structure])
  }
  # acceptor stem
  acceptorStem <- .getAcceptorStem(gr,
                                   strList,
                                   tRNAStructureTypes,
                                   loopPositions,
                                   discriminator)
  # T stem
  TStem <- .getTstem(gr,
                     strList,
                     tRNAStructureTypes,
                     loopPositions,
                     acceptorStem)
  # D stem
  DStem <- .getDstem(gr,
                     strList,
                     tRNAStructureTypes,
                     loopPositions,
                     acceptorStem)
  # anticodon stem
  anticodonStem <- .getAnticodonStem(gr,
                                     strList,
                                     tRNAStructureTypes,
                                     loopPositions)
  ##############################################################################
  # Intermediate control step
  ##############################################################################
  # if Tstem does not immediately follow acceptor stem, adjust coordinates of
  # acceptor stem
  fBase <- TStem$prime3$end + 1 != acceptorStem$prime3$start &
    !is.na(acceptorStem$prime3$start) &
    !is.na(TStem$prime3$end) &
    !is.na(anticodonStem$prime3$end) &
    !is.na(DStem$prime5$start)
  # 1 pos and only acceptor stem short
  fA <- which(fBase & 
                (TStem$prime3$end - TStem$prime3$start) == 4 &
                (acceptorStem$prime3$end - acceptorStem$prime3$start) < 6)
  if(length(fA) > 0){
    diff <- acceptorStem$prime3$start[fA] - TStem$prime3$end[fA] - 1
    acceptorStem$prime3$start[fA] <- TStem$prime3$end[fA] + 1
    acceptorStem$prime5$end[fA] <- acceptorStem$prime5$end[fA] + diff
  }
  fDcheckOverlap <- which(acceptorStem$prime5$end >= DStem$prime5$start)
  if(length(fDcheckOverlap) > 0){
    acceptorStem$prime5$end[fDcheckOverlap] <- 
      DStem$prime5$start[fDcheckOverlap] - 1
  }
  # 1 pos and only t stem short
  fT <- which(fBase & 
                (TStem$prime3$end - TStem$prime3$start) < 4 &
                (acceptorStem$prime3$end - acceptorStem$prime3$start) == 6)
  if(length(fT) > 0){
    diff <- acceptorStem$prime3$start[fT] - TStem$prime3$end[fT] - 1
    TStem$prime3$end[fT] <- acceptorStem$prime3$start[fT] - 1
    TStem$prime5$start[fT] <- TStem$prime5$start[fT] - diff 
  }
  fACcheckOverlap <- which(TStem$prime5$start <= anticodonStem$prime3$end)
  if(length(fACcheckOverlap) > 0){
    TStem$prime5$start[fACcheckOverlap] <- anticodonStem$prime3$end + 1
  }
  # additional space: add to acceptr 3'-end
  fRest <- which(TStem$prime3$end + 1 != acceptorStem$prime3$start &
                   !is.na(acceptorStem$prime3$start) &
                   !is.na(TStem$prime3$end) &
                   !is.na(anticodonStem$prime3$end) &
                   !is.na(DStem$prime5$start))
  if(length(fRest) > 0){
    diff <- acceptorStem$prime3$start[fRest] - TStem$prime3$end[fRest] - 1
    acceptorStem$prime3$start[fRest] <- acceptorStem$prime3$start[fRest] - diff
  }
  # dstem
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ans[["acceptorStem"]] <- .get_IRanges(acceptorStem,
                                        gr$tRNA_anticodon)
  ans[["TStem"]] <- .get_IRanges(TStem,
                                 gr$tRNA_anticodon)
  ans[["DStem"]] <- .get_IRanges(DStem,
                                 gr$tRNA_anticodon)
  ans[["anticodonStem"]] <- .get_IRanges(anticodonStem,
                                         gr$tRNA_anticodon)
  if(all(structure %in% names(ans))){
    return(ans[structure])
  }
  ##############################################################################
  if(any(structure %in% c("Dloop","Tloop","anticodonLoop"))){
    # D loop, T loop and anticodon loop
    Dloop <- .getDloop(gr,
                       strList,
                       tRNAStructureTypes,
                       loopPositions,
                       DStem)
    Tloop <- .getTloop(gr,
                       strList,
                       tRNAStructureTypes,
                       loopPositions,
                       TStem)
    anticodonLoop <- .getAnticodonLoop(gr,
                                       strList,
                                       tRNAStructureTypes,
                                       loopPositions,
                                       anticodonStem)
    ans[["Dloop"]] <- .get_IRanges(Dloop,
                                   gr$tRNA_anticodon)
    ans[["Tloop"]] <- .get_IRanges(Tloop,
                                   gr$tRNA_anticodon)
    ans[["anticodonLoop"]] <- .get_IRanges(anticodonLoop,
                                           gr$tRNA_anticodon)
  }
  if(all(structure %in% names(ans))){
    return(ans[structure])
  }
  # Dprime5, Dprime3 and vairable loop
  Dprime5 <- .getDprime5(gr,
                         strList,
                         acceptorStem,
                         DStem)
  Dprime3 <- .getDprime3(gr,
                         strList,
                         DStem,
                         anticodonStem)
  variableLoop <- .getVariableLoop(gr,
                                   strList,
                                   TStem,
                                   anticodonStem)
  ans[["Dprime5"]] <- .get_IRanges(Dprime5,
                                   gr$tRNA_anticodon)
  ans[["Dprime3"]] <- .get_IRanges(Dprime3,
                                   gr$tRNA_anticodon)
  ans[["variableLoop"]] <- .get_IRanges(variableLoop,
                                        gr$tRNA_anticodon)
  if(!all(structure %in% names(ans))){
    stop("Something went wrong.")
  }
  ans
}

# returns the boundary positions of the center of the tRNA
.get_loop_positions <- function(strList){
  loopids <- .get_loop_ids(strList)
  min <- .local_min(loopids)
  ans <- mapply(
    function(z,zz){
      m <- unlist(lapply(zz,
                         function(zzz){
                           # not length but max number of distinct occurances
                           zzzz <- diff(which(z == z[zzz]))
                           length(zzzz[zzzz > 1]) + 1
                         }))
      f <- z[zz[which(m == max(m))]]
      # no acceptor found, return empty IRanges
      if(0 %in% f) return(IRanges::IRanges())
      #
      f <- z[zz[which(m == max(m))[1]]]
      f <- which(z == f)
      forw <- c(2,diff(f))
      rev <- rev(c(-2,diff(rev(f))))
      IRanges::IRanges(start = f[which(forw > 1)],
                       end = f[which(rev < -1)])
    },
    loopids,
    min,
    SIMPLIFY = FALSE)
  return(ans)
}

.local_max <- function(list){
  .local_min_max(list)
}
.local_min <- function(list){
  .local_min_max(list, type = "min")
}

.local_min_max <- function(list,
                           type = "max"){
  z <- -.Machine$integer.max
  if(type == "min"){
    z <- .Machine$integer.max
  }
  ans <- lapply(
    list,
    function(x){
      y <- diff(c(z, x))
      y <- diff(c(z, x)) > 0L
      y <- cumsum(rle(y)$lengths)
      y <- y[seq.int(1L, length(y), 2L)]
      y <- y[y != 1 & y != max(x)]
      y
    })
  ans
}

# returns the type of tRNA
# either canonical, armless_dloop, armless_tloop, armless_ac, 
# armless_d_t, centerless, unrecognized
.get_tRNA_structure_type <- function(loopPositions){
  ans <- rep("unrecognized", length(loopPositions))
  noLoops <- lapply(loopPositions,length)
  # if four our five loops are found, it is a canonical tRNA
  ans[noLoops == 4 | noLoops == 5] <- "canonical"
  # one arm is missing
  one_arm_missing <- lapply(loopPositions[noLoops == 3],
                            function(z){
                              w <- BiocGenerics::width(z)
                              which(w == max(w))[1]
                            })
  tmp <- rep("unrecognized",length(one_arm_missing))
  tmp[one_arm_missing == 1] <- "armless_dloop"
  tmp[one_arm_missing == 2] <- "armless_ac"
  tmp[one_arm_missing == 3] <- "armless_tloop"
  ans[noLoops == 3] <- tmp
  # two arms are missing
  ans[noLoops == 2] <- "armless_d_t"
  # only the acceptor stem exists
  ans[noLoops == 1] <- "centerless"
  ans
}


# calculate the width of the ranges.
# In case of invalid structures NA are returned, which require the width
# to be set. For this case width = 0 is returned
.get_width <- function(start,
                       end){
  width <- rep(0,length(start))
  f <- which(assertive::is_negative(as.numeric(start)))
  width[f] <- end[f] - 1
  f <- which(assertive::is_negative(as.numeric(end)))
  width[f] <- 0
  f <- which(!is.na(start) & 
               !is.na(end))
  width[f] <- end[f] - start[f] + 1
  width
}
# set positions from a zero width element to 1
.get_start_end <- function(pos,
                           width,
                           value){
  f <- width == 0
  pos[f] <- value
  pos
}
# get an IRanges object. Make sure that invalid entries have width = 0 by 
# applying correct start and end coordinates
.get_IRanges <- function(pos,
                         names){
  if(!is.list(pos)){
    return(.get_IRanges2(pos,
                         pos,
                         names))
  }
  if(!is.list(pos[[1]])){
    return(.get_IRanges2(pos$start,
                         pos$end,
                         names))
  }
  return(lapply(pos,function(p){
    .get_IRanges2(p$start,
                  p$end,
                  names)
  }))
}
.get_IRanges2 <- function(start,
                          end,
                          names){
  width <- .get_width(start,
                      end)
  IRanges::IRanges(start = .get_start_end(start, width, NA),
                   end = .get_start_end(end, width, 0),
                   width = width,
                   names = names)
}

################################################################################
################################################################################
################################################################################

# checks if a invalid char is in the seq
.has_invalid_char <- function(seq,
                              start,
                              end,
                              char){
  seq <- substr(seq,
                start,
                end)
  unname(!is.na(stringr::str_locate(seq,
                                    char)[,"start"]))
}


################################################################################
################################################################################
################################################################################
.getAcceptorStem <- function(x,
                             strList,
                             tRNAStructureTypes,
                             loopPositions,
                             discriminator){
  # 3'end
  # remove discriminator
  end <- discriminator - 1
  # 5'start is opposite of the discriminator 
  start5 <- mapply(
    function(z,zend){
      if(is.na(zend)) return(NA)
      diff <- zend - max(z$reverse)
      ans <- z[z$reverse == max(z$reverse),]$forward - diff
      if(ans < 0) return(NA)
      ans
    },
    strList,
    end)
  # use center positions to get search parameters
  tmp <- sapply(loopPositions,
                     function(l){
                       if(length(l) == 0) return(c(NA,NA))
                       list(start(l)[1],
                            end(l)[length(l)])
                     })
  maxStart <- tmp[1,]
  minEnd <- tmp[2,]
  # get the last base paired position, which is before the D stem
  # returns both end5 and start3
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$forward <= zend & 
                 z$reverse >= zend & 
                 z$forward != 0, ]
      ans <- ans[ans$forward <= zstart,]
      if(nrow(ans) == 0) return(c(NA,NA))
      return(c(max(ans$forward),min(ans$reverse)))
    },
    strList,
    maxStart,
    minEnd)
  end5 <- tmp[1,]
  start3 <- tmp[2,]
  coord <- list("prime5" = list(start = unname(unlist(start5)),
                                end = unname(unlist(end5))),
                "prime3" = list(start = unname(unlist(start3)),
                                # discriminator base is remove
                                end = end))
  return(coord)
}

################################################################################
################################################################################
################################################################################
.getDprime5 <- function(x,
                        strList,
                        acceptor,
                        dstem){
  start <- acceptor$prime5$end + 1
  end <- dstem$prime5$start - 1
  coord <- list(start = unname(start),
                end = unname(end))
  # if width is negative return NA
  f <- which(is.na(coord$start) | 
               is.na(coord$end) |
               coord$start > coord$end)
  if(length(f) > 0){
    coord$start[f] <- NA
    coord$end[f] <- NA
  }
  return(coord)
}

################################################################################
.getDstem <- function(x,
                      strList,
                      tRNAStructureTypes,
                      loopPositions,
                      acceptor){
  # for start search use acceptor
  width <- .get_width(acceptor$prime5$start,
                      acceptor$prime5$end)
  minStart <- .get_start_end(acceptor$prime5$end, width, NA) + 1
  # for end search use center positions
  maxEnd <- mapply(
    function(l, type){
      if(type %in% c("armless_dloop",
                     "armless_d_t",
                     "unrecognized")) return(NA)
      return(start(l)[2])
    },
    loopPositions,
    tRNAStructureTypes)
  # get the first basepair which is completly in first half of tRNA
  # make sure if the start position is a mismatch to extend the end for this as
  # well
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos < zend & 
                 z$reverse < zend & 
                 z$forward != 0, ]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- ans[1,]
      return(c(ans$forward,ans$reverse))
    },
    strList,
    minStart,
    maxEnd)
  start5 <- tmp[1,]
  end3 <- tmp[2,]
  # remove all unpaired position between the start and end. the last row in the
  # middle contain both end5 and start3
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos >= zstart &
              z$pos <= zend &
              z$forward != 0,]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- ans[max(1,floor(nrow(ans) / 2)),]
      return(c(ans$forward,ans$reverse))
    },
    strList,
    start5,
    end3)
  end5 <- tmp[1,]
  start3 <- tmp[2,]
  #
  coord <- list("prime5" = list(start = unname(start5),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  return(coord)
}

################################################################################
.getDloop <- function(x,
                      strList,
                      tRNAStructureTypes,
                      loopPositions,
                      dstem){
  coord <- list(start = (dstem$prime5$end + 1),
                end = (dstem$prime3$start - 1))
  # If no D stem was found, assign D loop from center positions
  f <- which(is.na(dstem$prime5$start) |
               is.na(dstem$prime5$end) |
               is.na(dstem$prime3$start) |
               is.na(dstem$prime3$end))
  if(length(f) > 0) {
    tmp <- mapply(
      function(l, type){
        start <- NA
        end <- NA
        if(type %in% c("armless_dloop",
                       "armless_d_t")){
          start <- start(l[1]) + 1
          end <- end(l[1])
        }
        return(c(start,
                 end))
      },
      loopPositions[f],
      tRNAStructureTypes[f])
    coord$start[f] <- tmp[1,]
    coord$end[f] <- tmp[2,]
  }
  return(coord)
}

################################################################################
.getDprime3 <- function(x,
                        strList,
                        dstem,
                        anticodon){
  start <- dstem$prime3$end + 1
  end <- anticodon$prime5$start - 1
  coord <- list(start = unname(start),
                end = unname(end))
  # if width is negative return NA
  f <- which(is.na(coord$start) | 
               is.na(coord$end) |
               coord$start > coord$end)
  if(length(f) > 0){
    coord$start[f] <- NA
    coord$end[f] <- NA
  }
  return(coord)
}

################################################################################
################################################################################
################################################################################
.getAnticodonStem <- function(x,
                              strList,
                              tRNAStructureTypes,
                              loopPositions){
  # for search parameters use center positions 
  tmp <- mapply(
    function(l, type){
      start <- NA
      end <- NA
      if(type == "canonical") {
        start <- end(l)[2]
        end <- start(l)[3]
      }
      if(type %in% c("armless_dloop")){
        start <- end(l)[1] 
        end <- start(l)[2]
      }
      if(type %in% c("armless_tloop")){
        start <- end(l)[length(l)-1]
        end <- start(l)[length(l)]
      }
      if(type %in% c("armless_d_t")){
        start <- end(l)[1]
        end <- start(l)[2]
      }
      return(c(start,
               end))
    },
    loopPositions,
    tRNAStructureTypes)
  minStart <- tmp[1,]
  maxEnd <- tmp[2,]
  # 
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos >= zstart & 
                 z$pos <= zend & 
                 z$forward != 0, ]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- ans[1,]
      return(c(ans$forward,ans$reverse))
    },
    strList,
    minStart,
    maxEnd)
  start5 <- tmp[1,]
  end3 <- tmp[2,]
  # get the last base paired position, which is before the D stem
  # returns both end5 and start3
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos >= zstart & 
                 z$pos <= zend &
                 z$forward != 0,]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- ans[max(1,floor(nrow(ans) / 2)),]
      return(c(ans$forward,ans$reverse))
    },
    strList,
    start5,
    end3)
  end5 <- tmp[1,]
  start3 <- tmp[2,]
  # adjust coordinated for tRNA with non canonical base pairing at end of
  # anticodon stem
  f <- which((start3 - 9) > end5 & 
               !is.na(start3) & 
               !is.na(end5))
  if(length(f) > 0){
    end5[f] <- end5[f] + 1
    start3[f] <- start3[f] - 1
  }
  #
  coord <- list("prime5" = list(start = unname(start5),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  return(coord)
}

################################################################################
.getAnticodonLoop <- function(x,
                              strList,
                              tRNAStructureTypes,
                              loopPositions,
                              anticodon){
  coord <- list(start = (anticodon$prime5$end + 1),
                end = (anticodon$prime3$start - 1))
  # if anticodon stem not found, set d and t stem coordinates
  f <- which(is.na(coord$start) | is.na(coord$end))
  if(length(f) == 0){
    return(coord)
  }
  # if some anticodon stems are not present, use center positions
  tmp <- mapply(
    function(l, type){
      start <- NA
      end <- NA
      if(type == "centerless"){
        start <- start(l[1]) + 1
        end <- end(l[1]) - 1
      }
      if(type == "armless_ac"){
        start <- start(l[2])
        end <- end(l[2])
      }
      return(c(start,
               end))
    },
    loopPositions[f],
    tRNAStructureTypes[f])
  coord$start[f] <- tmp[1,]
  coord$end[f] <- tmp[2,]
  return(coord)
}

################################################################################
################################################################################
################################################################################
.getVariableLoop <- function(x,
                             strList,
                             tstem,
                             anticodon){
  coord <- list(start = (anticodon$prime3$end + 1),
                end = (tstem$prime5$start - 1))
  # if no T stem found or directly adjacent, no variable loop
  f <- which(is.na(coord$start) | 
               is.na(coord$end) |
               coord$start > coord$end)
  if(length(f) > 0){
    coord$start[f] <- NA
    coord$end[f] <- NA
  }
  return(coord)
}


################################################################################
################################################################################
################################################################################
.getTstem <- function(x,
                      strList,
                      tRNAStructureTypes,
                      loopPositions,
                      acceptor){
  # get maximal end from acceptor
  width <- .get_width(acceptor$prime3$start,
                      acceptor$prime3$end)
  maxEnd <- .get_start_end(acceptor$prime3$start, width, NA) - 1
  # for start search use center positions
  minStart <- mapply(
    function(l, type){
      if(type %in% c("armless_tloop",
                     "armless_d_t",
                     "unrecognized")) return(NA)
      return(end(l)[length(l)-1])
    },
    loopPositions,
    tRNAStructureTypes)
  # get the last basepair which is completly in second half of tRNA and
  # make sure to avoid the variable loop potential basepairs
  # make sure if the start position is a mismatch to extend the end for this as
  # well
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos > zstart & 
                 z$reverse > zstart & 
                 z$forward != 0, ]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- z[z$forward == max(ans$forward),]
      return(c(ans$reverse,ans$forward))
    },
    strList,
    minStart,
    maxEnd)
  start5 <- tmp[1,]
  end3 <- tmp[2,]
  # remove all unpaired position between the start and end. the last row in the
  # middle contain both end5 and start3
  tmp <- mapply(
    function(z,zstart,zend){
      if(is.na(zstart) | is.na(zend)) return(c(NA,NA))
      ans <- z[z$pos >= zstart &
                 z$pos <= zend &
                 z$forward != 0,]
      if(nrow(ans) == 0) return(c(NA,NA))
      ans <- ans[max(1,floor(nrow(ans) / 2)),]
      return(c(ans$forward,ans$reverse))
    },
    strList,
    start5,
    end3)
  end5 <- tmp[1,]
  start3 <- tmp[2,]
  #
  coord <- list("prime5" = list(start = unname(start5),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  return(coord)
}

################################################################################
.getTloop <- function(x,
                      strList,
                      tRNAStructureTypes,
                      loopPositions,
                      tstem){
  coord <- list(start = (tstem$prime5$end + 1),
                end = (tstem$prime3$start - 1))
  # If no T stem was found, assign loop from center positions
  f <- which(is.na(tstem$prime5$start) |
               is.na(tstem$prime5$end) |
               is.na(tstem$prime3$start) |
               is.na(tstem$prime3$end))
  if(length(f) > 0) {
    tmp <- mapply(
      function(l, type){
        start <- NA
        end <- NA
        if(type %in% c("armless_tloop",
                       "armless_d_t")){
          start <- start(l[length(l)])
          end <- end(l[length(l)]) - 1
        }
        return(c(start,
                 end))
      },
      loopPositions[f],
      tRNAStructureTypes[f])
    coord$start[f] <- tmp[1,]
    coord$end[f] <- tmp[2,]
  }
  return(coord)
}

################################################################################
################################################################################
################################################################################
.getDiscriminator <- function(x,
                              strList,
                              tRNAStructureTypes){
  ans <- as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length(x)-3,
                           .get_tRNA_length(x)))
  ans[tRNAStructureTypes == "unrecognized"] <- NA
  return(ans)
}
