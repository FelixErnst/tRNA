#' @include tRNA.R
NULL

################################################################################
# tRNA
################################################################################

#' @rdname istRNAGRanges
#' @export
setGeneric (
  name = "istRNAGRanges",
  def = function(x) standardGeneric("istRNAGRanges")
)

# Structures and Sequences -----------------------------------------------------

#' @rdname gettRNAstructureSeqs
#' @export
setGeneric (
  name = "gettRNAstructureGRanges",
  def = function(x,
                 structure = "") standardGeneric("gettRNAstructureGRanges")
)
#' @rdname gettRNAstructureSeqs
#' @export
setGeneric (
  name = "gettRNAstructureSeqs",
  def = function(x, structure = "", joinCompletely = FALSE,joinFeatures = FALSE,
                 padSequences = TRUE) standardGeneric("gettRNAstructureSeqs")
)
#' @rdname gettRNABasePairing
#' @export
setGeneric (
  name = "gettRNABasePairing",
  def = function(x,
                 with.nucleotides = FALSE) standardGeneric("gettRNABasePairing")
)
#' @rdname gettRNABasePairing
#' @export
setGeneric (
  name = "gettRNALoopIDs",
  def = function(x) standardGeneric("gettRNALoopIDs")
)

# Features ---------------------------------------------------------------------

#' @rdname gettRNASummary
#' @export
setGeneric (
  name = "gettRNASummary",
  def = function(x) standardGeneric("gettRNASummary")
)

# Subsetting -------------------------------------------------------------------

#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasTStem",
  def = function(x,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasTStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasDStem",
  def = function(x,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasDStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAcceptorStem",
  def = function(x,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasAcceptorStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAnticodonStem",
  def = function(x,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasAnticodonStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasTloop",
  def = function(x,
                 length = NA) standardGeneric("hasTloop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasDloop",
  def = function(x,
                 length = NA) standardGeneric("hasDloop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAnticodonLoop",
  def = function(x,
                 length = NA) standardGeneric("hasAnticodonLoop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasVariableLoop",
  def = function(x,
                 length = NA,
                 paired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasVariableLoop")
)

# Visualization ----------------------------------------------------------------

#' @rdname gettRNAFeaturePlots
#' @export
setGeneric (
  name = "gettRNAFeaturePlots",
  def = function(x,
                 plotScores = FALSE,
                 scores = NA,
                 scoreLabel = "Score") standardGeneric("gettRNAFeaturePlots")
)
