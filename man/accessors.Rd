\name{accessors}
\title{Accessor functions for OTUset and TAXset objects}
\alias{otuID}
\alias{otuID,.OTUset-method}
\alias{sampleID}
\alias{sampleID,.OTUbase-method}
\alias{tax}
\alias{tax,.TAXset-method}
\alias{tax<-}
\alias{tax<-,.TAXset-method}
\alias{tax<-,.TAXset,data.frame-method}

\description{
    These functions provide access to some of the slots of OTUset and TAXset objects. \code{otuID} returns the \code{otuID} slot of OTUset objects. \code{sampleID} returns the \code{sampleID} slot of both OTUset and TAXset objects. \code{tax} and \code{tax<-} return and replace the \code{tax} slot of TAXset objects.
}

\usage{
    sampleID(object, ...)
    otuID(object, ...)
    tax(object, ...)
    tax(object)<-value
}


\arguments{
    \item{object}{An OTUset or a TAXset object}
    \item{value}{The replacement value for \code{tax}}
    \item{...}{Added for completeness. Enables the passing of arguments.}
}

\value{
    \code{sampleID} and \code{otuID} return a character.
    \code{tax} returns a data.frame.
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## get the sampleID slot
sampleID(soginOTU) 

## get the otuID slot
otuID(soginOTU)

}

\seealso{
\code{\link{ShortRead}}
}


