\name{assignmentData}
\title{assignmentData}
\alias{aData}
\alias{aData,.OTUbase-method}
\alias{aData<-}
\alias{aData<-,.OTUbase-method}
\alias{aData<-,.OTUbase,data.frame-method}
\alias{assignmentData}
\alias{assignmentData,.OTUbase-method}
\alias{assignmentData<-}
\alias{assignmentData<-,.OTUbase-method}
\alias{assignmentData<-,.OTUbase,AnnotatedDataFrame-method}
\alias{assignmentLabels}
\alias{assignmentLabels,.OTUbase-method}
\alias{assignmentLabels<-}
\alias{assignmentLabels<-,.OTUbase-method}
\alias{assignmentLabels<-,.OTUbase,character-method}
\alias{assignmentNames,.OTUbase-method}
\alias{assignmentNames}

\description{
These accessors access and replace the assignmentData slot of OTUbase objects. \code{assignmentData} is an AnnotatedDataFrame. \code{assignmentData} and \code{assignmentData<-} access and replace this AnnotatedDataFrame. \code{assignmentLabels} and \code{assignmentLabels<-} access and replace the labels of this AnnotatedDataFrame. \code{aData} and \code{aData<-} access and replace the dataframe component of the AnnotatedDataFrame. 

\code{assignmentNames} returns the assignment names present in the \code{assignmentData} slot.
}


\usage{

    aData(object,...)
    aData(object)<-value
    assignmentData(object,...)
    assignmentData(object)<-value
    assignmentLabels(object,...)
    assignmentLabels(object)<-value
    assignmentNames(object,...)

}

\arguments{
    \item{object}{An OTUset or a TAXset object}
    \item{value}{The replacement value for \code{assignmentData} or \code{assignmentLabels}}
    \item{...}{Added for completeness. Enables the passing of arguments.}
}

\value{
    \code{aData} returns a dataframe.
    \code{assignmentData} returns an AnnotatedDataFrame.
    \code{assignmentLabels} returns a character.
    \code{assignmentNames} returns a character.
}

\examples{

## locate directory with data
dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

## get the aData dataframe
aData(soginOTU)

## get the assignmentData slot
assignmentData(soginOTU)

}







