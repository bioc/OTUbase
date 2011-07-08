\name{readOTUset}
\title{readOTUset}
\alias{readOTUset}

\description{
This function reads in data and creates an OTUset object
}


\usage{
    readOTUset(dirPath, otufile, level, fastafile, qualfile, samplefile, sampleADF, assignmentADF, sADF.names, aADF.names, rdp=F, otufiletype)
}

\arguments{
    \item{dirPath}{The directory path were the datafile are located. This is the current directory by default.}
    \item{otufile}{The OTU file. The only format currently supported is the Mothur format for .list files.}
    \item{level}{The OTU clustering level. By default this is 0.03. This level must correspond to levels present in the otufile.}
    \item{fastafile}{The fasta file. This is read in by ShortRead.}
    \item{qualfile}{The quality file. This is read in by ShortRead.}
    \item{samplefile}{The sample file. Currently this must be in Mothur format (.groups).}
    \item{sampleADF}{The sample meta data file. This is in AnnotatedDataFrame format.}
    \item{assignmentADF}{The assignment meta data file (the OTU meta data). This is generally in AnnotatedDataFrame format although it is also possible to read in an RDP classification file if there is only one read classification for each cluster and \code{rdp}=TRUE.}
    \item{sADF.names}{The column of the sampleADF file that has the sample names.}
    \item{aADF.names}{The column of the assignmentADF file that has the assignment names.}
    \item{rdp}{By default this is FALSE. Change to TRUE if assignmentADF is an RDP classification file. The RDP file must be in the fixed format.}
    \item{otufiletype}{The type of OTU file. Takes values "mothur" and "cdhit". Defaults to "mothur".}
}

\examples{

dirPath <- system.file("extdata/Sogin_2006", package="OTUbase")

## read in data into OTUset object
soginOTU <- readOTUset(dirPath=dirPath, level="0.03", samplefile="sogin.groups", fastafile="sogin.fasta", otufile="sogin.unique.filter.fn.list", sampleADF="sample_metadata.txt")

soginOTU

}
