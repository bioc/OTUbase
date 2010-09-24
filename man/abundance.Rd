\name{abundance}
\title{abundance}
\alias{abundance,.TAXset-method}
\alias{abundance,.OTUset-method}
\alias{abundance}

\description{
    \code{abundance} generates an abundance table. This table can be either weighted or unweighted.
}
\usage{
abundance(object, ...)
}

\arguments{

    \item{object}{An OTUset or a TAXset object}
    \item{...}{Additional arguments. These will depend on if the object is an OTUset or a TAXset object.}

}

\details{
    These are other arguments passed to \code{abundance}
    \itemize{
        \item{taxCol}{ If generating the abundance from a TAXset object, \code{taxCol} selects the column of the \code{tax} dataframe from which to calculate the abundance.}
        \item{assignmentCol}{ If generating the abundance from an OTUset object \code{assignmentCol} will select a column of the \code{assignmentData} dataframe to use when calculating abundance. This will override the default of creating an abundance table of the OTUs and instead create an abundance table of a column in the \code{assignmentData} dataframe.}
	\item{sampleCol}{ \code{sampleCol} generates the abundance table using a column in the sampleData data fram instead of the default of using the sampleID.}
        \item{collab}{ An optional parameter that selets a column of the \code{sampleData} dataframe to use when labeling the columns of the abundance table.}
        \item{weighted}{ By default this is FALSE. When set to TRUE \code{abundance} will return proportional abundances.}

    }
}

\value{
    The returned value will be a data.frame.
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## calculate abundance
abundance(soginOTU, collab="Site")
}


\keyword{abundance}
