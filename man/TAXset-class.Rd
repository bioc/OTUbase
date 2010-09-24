\name{.TAXset-class}
\alias{TAXset}
\alias{TAXsetB-class}
\alias{TAXsetF-class}
\alias{TAXsetQ-class}
\alias{TAXsetB}
\alias{TAXsetF}
\alias{TAXsetQ}



\docType{class}
\title{"TAXset" class for TAX data}
\description{
	This class provides a way to store and manipulate read-classification data. ".TAXset" is inherited by "TAXsetQ", "TAXsetF", and "TAXsetB". The user will want to use "TAXsetQ" when quality data is available, "TAXsetF" when sequence data (without quality data) is available, and "TAXsetB" when only classification and sample data are available.
}

\section{Slots}{
  TAXsetB includes Slots \code{id}, \code{sampleID}, \code{tax}, \code{sampleData}, \code{assignmentData}.

  TAXsetF includes Slots \code{id} \code{sampleID}, \code{tax}, \code{sampleData}, \code{assignmentData}, \code{sread}.

  TAXsetQ includes Slots \code{id} \code{sampleID}, \code{tax}, \code{sampleData}, \code{assignmentData}, \code{sread}, \code{quality}.


}

\section{Methods}{
    
    Methods include:

    \describe{
        
	\item{id}{provides access to the id slot of \code{object}}
        \item{sampleID}{provides access to the sampleID slot of \code{object}}
        \item{tax}{provides access the tax slot of \code{object}}
        \item{sampleData}{provides accesss the sampleData slot of \code{object}}
        \item{assignmentData}{provides access the assignmentData slot of \code{object}}
	\item{sread}{provides access to the sread slot of \code{object}}
	\item{quality}{provides access to the quality slot of \code{object}}
        
        \item{seqnames}{returns the first word of the id line. Intended to extract the sequence name from other sequence information.}
        \item{nsamples}{returns the number of samples in an TAXset object}

	\item{show}{\code{signature(object=".TAXset")}: provides a brief summary of the object, including its class, number of sequences, and number of samples.}
    }
}

\examples{
showClass(".TAXset")
showMethods(class=".TAXset")
showClass("TAXsetQ")
}
\keyword{classes}



