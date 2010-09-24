\name{readTAXset}
\title{readTAXset}
\alias{readTAXset}
\description{
Function to read in data and create a TAXset object
}

\usage{
    readTAXset(dirPath, taxfile, namefile, fastafile, qualfile, samplefile, sampleADF, assignmentADF, sADF.names, aADF.names, type, ...)
}

\arguments{
    \item{dirPath}{The directory path were the datafile are located. This is the current directory by default.}
    \item{taxfile}{The classification file. The default format is RDP's fixed format.}
    \item{namefile}{A names file in the Mothur format. This is used to add removed unique sequences back into the dataset.}
    \item{fastafile}{The fasta file. This is read in by ShortRead.}
    \item{qualfile}{The quality file. This is read in by ShortRead.}
    \item{samplefile}{The sample file. Currently this must be in Mothur format (.groups).}
    \item{sampleADF}{The sample meta data file. This is in AnnotatedDataFrame format.}
    \item{assignmentADF}{The assignment meta data file (the OTU meta data) This is in AnnotatedDataFrame format.}
    \item{sADF.names}{The column of the sampleADF file that has the sample names.}
    \item{aADF.names}{The column of the assignmentADF file that has the assignment names.}
    \item{type}{This is the type of taxfile. By default this is the RDP fixed format. However, if \code{type} is changed to anything else the \code{read.table} function is used to read in the \code{taxfile}. In this case the first column of the taxfile must be the sequence names.}
    \item{...}{Additional arguments passed to \code{read.table} to read in the taxfile.}

}

\examples{

dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into TAXset object
soginTAX <- readTAXset(dirPath=dirPath, samplefile="sogin.groups", fastafile="sogin.fasta", taxfile="sogin.unique.fix.tax", namefile="sogin.names", sampleADF="sample_metadata.txt")

soginTAX

}
