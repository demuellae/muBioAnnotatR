% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RegionSetDB-class.R
\docType{methods}
\name{regionSetGr,RegionSetDB-method}
\alias{regionSetGr,RegionSetDB-method}
\alias{regionSetGr}
\alias{regionSetGr,regionSetGrDB-method}
\title{regionSetGr-methods}
\usage{
\S4method{regionSetGr}{RegionSetDB}(.object, name, col, gen)
}
\arguments{
\item{.object}{\code{\linkS4class{regionSetGrDB}} object}

\item{name}{name of the region set to be retrieved (character)}

\item{col}{name of the collection the region set is in (character)}

\item{gen}{genome assembly for the region set (character)}
}
\value{
\code{GRanges} object containing the region set or \code{NULL} if it does not exist
}
\description{
Return region set \code{GRanges} object for a given name, collection and genome
}
\author{
Fabian Mueller
}
