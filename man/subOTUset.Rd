\name{subOTUset}
\title{subOTUset}
\alias{subOTUset}

\description{
Function to get a subset of an OTUset object.
}

\usage{

subOTUset(object, samples, otus)

}

\arguments{
    \item{object}{An OTUset object}
    \item{samples}{A list of sample names}
    \item{otus}{A list of OTU names}
}

\value{
    \code{subOTUset} returns an OTUset
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## get subset of soginOTU

subOTUset(soginOTU, samples=getSamples(soginOTU, colnum="Site", value="Labrador", exact=FALSE))

}
