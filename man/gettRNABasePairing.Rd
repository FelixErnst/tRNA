% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/tRNA-dotbracket.R
\name{gettRNABasePairing}
\alias{gettRNABasePairing}
\alias{gettRNALoopIDs}
\alias{gettRNABasePairing,GRanges-method}
\alias{gettRNALoopIDs,GRanges-method}
\title{Accessing Dot Bracket annotation of tRNAs}
\usage{
gettRNABasePairing(x, with.nucleotides = FALSE)

gettRNALoopIDs(x)

\S4method{gettRNABasePairing}{GRanges}(x, with.nucleotides = FALSE)

\S4method{gettRNALoopIDs}{GRanges}(x)
}
\arguments{
\item{x}{a GRanges object created by \code{import.tRNAscanAsGRanges} or
GRanges with equivalent information. The \code{tRNA_str} and \code{tRNA_seq} 
columns will be used to construct a StructuredXStringSet and used for input 
into \code{getBasePairing}.}

\item{with.nucleotides}{a single logical value: should the nucleotides be 
saved alongside the base pairing information in the 'base' column?}
}
\value{
\code{gettRNABasePairing}: 
The result is a \code{DotBracketDataFrame} with following columns: pos, 
forward, reverse, character and base. If a position is unpaired, forward and 
reverse will be \code{0}, otherwise it will match the base paired positions.

\code{gettRNALoopIDs}: return a list of list of loop ids.
}
\description{
\code{gettRNABasePairing} converts the dot bracket annotation into a
\code{DotBracketDataFrame}. Base pairing is indicated by cosrresponding 
numbers in the forward and reverse columns. For more detail have a look at
\code{\link[Structstrings:getBasePairing]{getBasePairing}}.

\code{gettRNALoopIDs} converts the dot bracket annotation into a 
\code{LoopIDList}. For more details have a look at 
\code{\link[Structstrings:getBasePairing]{getLoopIDList}}.
}
\examples{
data("gr", package = "tRNA")
gettRNABasePairing(gr[1])
gettRNALoopIDs(gr[1])
}
