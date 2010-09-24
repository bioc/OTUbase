\name{otherGenerics}
\title{Other Generics}
\alias{notus}
\alias{notus,.OTUset-method}
\alias{nsamples}
\alias{nsamples,.OTUbase-method}
\alias{seqnames}
\alias{seqnames,.OTUbase-method}

\description{
Various functions. \code{notus} returns the number of OTUs in an OTUset object. \code{nsamples} returns the number of samples in either an OTUset or a TAXset object. \code{seqnames} returns the sequence names of the OTUset or TAXset object without the extra information commonly present with the id.
}

\usage{
    notus(object, ...)
    nsamples(object, ...)
    seqnames(object, ...)
}

\arguments{
    \item{object}{An OTUset or a TAXset object.}
    \item{...}{Other arguments. These are currently nonfunctional.}
}

\examples{

dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## get the number of OTUs
notus(soginOTU)

## get the number of samples
nsamples(soginOTU)

}
