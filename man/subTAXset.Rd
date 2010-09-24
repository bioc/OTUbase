\name{subTAXset}
\title{subTAXset}
\alias{subTAXset}

\description{
Function to get a subset of an TAXset object.
}

\usage{

subTAXset(object, samples)

}

\arguments{
    \item{object}{An TAXset object}
    \item{samples}{A list of sample names}
}

\value{
    \code{subTAXset} returns an TAXset
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into TAXset object
soginTAX <- readTAXset(dirPath=dirPath, samplefile="sogin.groups", fastafile="sogin.fasta", taxfile="sogin.unique.fix.tax",namefile="sogin.names", sampleADF="sample_metadata.txt")

## get subset of soginTAX

subTAXset(soginTAX, samples=getSamples(soginTAX, colnum="Site", value="Labrador", exact=FALSE))

}
