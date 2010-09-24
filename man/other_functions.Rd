\name{other_functions}
\title{Other functions}
\alias{getOTUs}
\alias{getSamples}
\alias{o_diversity}
\alias{o_estimateR}
\alias{otuseqplot}
\alias{otusize}
\alias{otuspersample}
\alias{rseqplot}
\alias{seqspersample}
\alias{sharedotus}

\description{
These are other functions available. Caution is advised when using them. Some are still in development and others only work on specific objects (OTUset or TAXset).

}

\usage{
    getOTUs(object, colnum, value, exact)
    getSamples(object, colnum, value, exact)
    o_diversity(object, ...)
    o_estimateR(object, ...)
}

\arguments{
    \item{object}{An OTUset or a TAXset object.}
    \item{colnum}{The column of the \code{sampleData} or \code{assignmentData} dataframe that contains the value.}
    \item{value}{The desired value.}
    \item{exact}{If exact=T \code{value} must match perfectly. If exact=F \code{value} will grep instead of match. }
    \item{...}{Other arguments. Often these are passed to \code{abundance}}
}

\details{
    \itemize{
        \item{getOTUs}{ Returns OTU names that match given values in the assignmentData dataframe.}

	\item{getSamples}{ Returns sample names that match given values in the sampleData dataframe.}

	\item{o_diversity}{ Wrapper for vegan's diversity function.}

	\item{o_estimateR}{ Wrapper for vegan's estimateR function.}

	\item{otuseqplot}{ Plots the samples acording to number of OTUs and number of sequences.}

	\item{otusize}{ Returns the size of each OTU.}

	\item{otuspersample}{ Lists the number of OTUs in each sample.}

	\item{rseqplot}{ Plots the samples by estimated richness and number of sequences.}

	\item{seqspersample}{ Returns the number of sequences in each sample.}

	\item{sharedotus}{ Returns the number of OTUs shared between samples.}
    }
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

getSamples(soginOTU, colnum="Site", value="Labrador", exact=FALSE)

o_estimateR(soginOTU)

}

