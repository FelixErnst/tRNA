#' @include tRNA.R
NULL

#input type checks
.checkValueValidity <- function(value, checkValues,
                                .xvalue = .get_name_in_parent(value)){
  if(!all(value %in% checkValues)){
    stop("'",gsub("\"","",.xvalue),
         "' must be one of the following values: '",
         paste(checkValues, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}

# checks whether a string is a valid tRNA structure
.check_trna_structure_ident <-
  function(value, .xvalue = .get_name_in_parent(value)){
  # check input
  checkValues <- c("",TRNA_STRUCTURES)
  if(!(value %in% checkValues)){
    stop("'",gsub("\"","",.xvalue),
         "' must be one of the following values: '",
         paste(checkValues, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}

# get the tRNA length without the intron
.get_tRNA_length <- function(x){
  nchar(as.character(x$tRNA_seq))
}

.is_continous_evenly_spaced_c <- function(n){
  if(length(n) < 2) return(FALSE)
  n <- n[order(n)]
  n <- n - min(n)
  step <- n[2] - n[1]
  test <- seq(from = min(n), to = max(n), by = step)
  if(length(n) == length(test) &&
     all(as.character(n) == as.character(test))){
    return(TRUE)
  }
  return(FALSE)
}
.is_continous_evenly_spaced <- compiler::cmpfun(.is_continous_evenly_spaced_c)

# subset a structure data.frame using one or two iRanges
.subset_structure_c <- function(ir, str, pairedOnly = TRUE){
  if(pairedOnly){
    if(length(ir) == 2){
      f <- str$forward >= BiocGenerics::start(ir)[1] &
        str$forward <= BiocGenerics::end(ir)[1] &
        str$reverse >= BiocGenerics::start(ir)[2] &
        str$reverse <= BiocGenerics::end(ir)[2]
    } else {
      f <- str$forward >= BiocGenerics::start(ir) &
        str$forward <= BiocGenerics::end(ir)
    }
  } else {
    if(length(ir) == 2){
      f <- (str$pos >= BiocGenerics::start(ir)[1] &
              str$pos <= BiocGenerics::end(ir)[1]) |
        (str$pos >= BiocGenerics::start(ir)[2] &
           str$pos <= BiocGenerics::end(ir)[2])
    } else {
      f <- str$pos >= BiocGenerics::start(ir) &
        str$pos <= BiocGenerics::end(ir)
    }
  }
  str[f,]
}
.subset_structure <- compiler::cmpfun(.subset_structure_c)
