\name{.OTUset-class}
\alias{OTUset}
\alias{OTUsetB-class}
\alias{OTUsetF-class}
\alias{OTUsetQ-class}
\alias{OTUsetB}
\alias{OTUsetF}
\alias{OTUsetQ}



\docType{class}
\title{"OTUset" class for OTU data}
\description{
	This class provides a way to store and manipulate operational taxonomic unit data. ".OTUset" is inherited by "OTUsetQ", "OTUsetF", and "OTUsetB". The user will want to use "OTUsetQ" when quality data is available, "OTUsetF" when sequence data (without quality data) is available, and "OTUsetB" when only OTU and sample data are available.
}

\section{Slots}{
  OTUsetB includes Slots \code{id}, \code{sampleID}, \code{otuID}, \code{sampleData}, \code{assignmentData}.

  OTUsetF includes Slots \code{id} \code{sampleID}, \code{otuID}, \code{sampleData}, \code{assignmentData}, \code{sread}.

  OTUsetQ includes Slots \code{id} \code{sampleID}, \code{otuID}, \code{sampleData}, \code{assignmentData}, \code{sread}, \code{quality}.


}

\section{Methods}{
    
    Methods include:

    \describe{
        
	\item{id}{provides access to the id slot of \code{object}}
        \item{sampleID}{provides access to the sampleID slot of \code{object}}
        \item{otuID}{provides access the otuID slot of \code{object}}
        \item{sampleData}{provides accesss the sampleData slot of \code{object}}
        \item{assignmentData}{provides access the assignmentData slot of \code{object}}
	\item{sread}{provides access to the sread slot of \code{object}}
	\item{quality}{provides access to the quality slot of \code{object}}
        
        \item{seqnames}{returns the first word of the id line. Intended to extract the sequence name from other sequence information.}
        \item{nsamples}{returns the number of samples in an OTUset object}
        \item{notus}{returns the number of OTUs in an OTUset object}

	\item{show}{\code{signature(object=".OTUset")}: provides a brief summary of the object, including its class, number of sequences, number of samples, and number of OTUs.}
    }
}

\examples{
showClass(".OTUset")
showMethods(class=".OTUset")
showClass("OTUsetQ")
}
\keyword{classes}



