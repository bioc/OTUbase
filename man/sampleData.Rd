\name{sampleData}
\title{sampleData}
\alias{sData}
\alias{sData,.OTUbase-method}
\alias{sData<-}
\alias{sData<-,.OTUbase-method}
\alias{sData<-,.OTUbase,data.frame-method}
\alias{sampleData}
\alias{sampleData,.OTUbase-method}
\alias{sampleData<-}
\alias{sampleData<-,.OTUbase-method}
\alias{sampleData<-,.OTUbase,AnnotatedDataFrame-method}
\alias{sampleLabels}
\alias{sampleLabels,.OTUbase-method}
\alias{sampleLabels<-}
\alias{sampleLabels<-,.OTUbase-method}
\alias{sampleLabels<-,.OTUbase,character-method}




\description{
These functions access and replace the \code{sampleData} slot of OTUbase objects. \code{sampleData} and \code{sampleData<-} access and replace the AnnotatedDataFrame \code{sampleData}. \code{sampleLabels} and \code{sampleLabels<-} access and replace the labels of this AnnotatedDataFrame. \code{sData} and \code{sData<-} access and replace the dataframe component of the AnnotatedDataFrame.

}


\usage{

    sData(object,...)
    sData(object)<-value
    sampleData(object,...)
    sampleData(object)<-value
    sampleLabels(object,...)
    sampleLabels(object)<-value

}

\arguments{
    \item{object}{An OTUset or a TAXset object}
    \item{value}{The replacement value for \code{sampleData} or \code{sampleLabels}}
    \item{...}{Added for completeness. Enables the passing of arguments.}
}

\value{
    \code{sData} returns a dataframe.
    \code{sampleData} returns an AnnotatedDataFrame.
    \code{sampleLabels} returns a character.
    \code{assignmentNames} returns a character.
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## get the sData dataframe
sData(soginOTU)

## get the sampleData slot
sampleData(soginOTU)

}
