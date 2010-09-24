\name{clusterSamples}
\title{clusterSamples}
\alias{clusterSamples-generic}
\alias{clusterSamples-method}
\alias{clusterSamples,.OTUset-method}
\alias{clusterSamples,.TAXset-method}
\alias{clusterSamples}

\description{
This function is a wrapper for the vegan function vegedist and hclust. It allows the user to cluster samples using a number of different distance measure and clustering methods. Please see the documentation for vegedist and hclust for a more indepth explanation.
}

\usage{
    clusterSamples(object, ...)
}

\arguments{

    \item{object}{An OTUset or a TAXset object}
    \item{...}{Additional arguments. These will depend on if the object is an OTUset or a TAXset object.}

}

\details{
    These are other arguments passed to \code{clusterSamples}. For further information on specific arguments, please see \code{\link{abundance}}, \code{\link{vegdist}}, or \code{\link{hclust}}.
    \itemize{
        \item{taxCol}{ Column of the \code{tax} slot dataframe on which to cluster (unique to TAXset objects). Passed to the \code{abundance} function.}
        \item{assignmentCol}{ Column of the \code{assignmentData} dataframe used to classify sequences for clustering. This overrides the default of using the OTUs to cluster samples. This is passed to the \code{abundance} function.}
        \item{collab}{ Specifies a column of the \code{sampleData} dataframe that will provide the sample lables for the cluster analysis. This is passed to the \code{abundance} function.}
        \item{distmethod}{ The distance method to be used. This value is passed to the \code{vegedist} function. The default is the Bray-Curtis distance.}
        \item{clustermethod}{ The clustering method to be used. This value is passed to the \code{hclust} function. The default is \code{complete} clustering.}
    }

}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## cluster samples
clusterSamples(soginOTU, collab="Site", distmethod="jaccard")

}
