\name{examples_data}
\alias{examples_data}
\docType{data}
\title{Data sets to run example scripts}
\description{
Contains source data sets to run all the scripts provided in the \emph{Examples} sections in the documentation 
for the package \code{PMXStan}.
}
\usage{data("examples_data")}
\details{
All data sets in \code{examples_data} were directly imported from literature, or slightly modified to fit the
illustration purpose in a specific example included in this package.

\code{d1_nm_poppk}: this dataset is an expanded version of the original theophylline population PK data.
                    Covariates \code{AGE} and \code{GENDER} were randomly assigned.             

\code{d2_nm_poppkpd}: this dataset was simulated by using theophylline population PK data and an indirect
                      response PD model.
                      Covariates \code{BMK1} and \code{BMK2} were randomly assigned.             

\code{d3_nm_idvpkpd}: this dataset is similarly simulated as \code{d2_nm_poppkpd}, but with only a single
                      subject.            
}
\author{
Yuan Xiong and Wenping Wang
}
\keyword{datasets}
