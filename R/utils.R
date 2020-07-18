
.is_a_bool <- function(x){
    is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
    is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
    .is_non_empty_character(x) && length(x) == 1L
}

.is_a_string <- function(x){
    is.character(x) && length(x) == 1L
}

.are_whole_numbers <- function(x){
    tol <- 100 * .Machine$double.eps
    abs(x - round(x)) <= tol && !is.infinite(x)
}

.get_name_in_parent <- function(x) {
    .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
    paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}
