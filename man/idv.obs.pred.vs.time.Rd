% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PMXStanFit.R
\name{idv.obs.pred.vs.time}
\alias{idv.obs.pred.vs.time}
\title{Comparison of the time profiles of observations with the predictions for each individual.}
\usage{
idv.obs.pred.vs.time(fit)
}
\arguments{
\item{fit}{a \code{PMXStanFit} object.}
}
\description{
plots time profiles of the prediction medians and 95% predictive intervals on top of the time 
profile of observations for individual patients. This function is a generic version of the method 
\code{plot.gof.idv()} for the \code{\link{PMXStanFit}} class.
}
\value{No return value, called for side effects}
\examples{
\donttest{
m1 <- PMXStanModel(path = tempfile("pk_m1"), pk.struct = "1-cmpt", compile=TRUE)

data("examples_data")
dat <- prepareInputData(data.source = d1_nm_poppk, model = m1)

fit <- PMXStanFit(m1, dat, iter=100, chains=1)

fit$plot.gof.idv()
#idv.obs.pred.vs.time(fit)
}
}
\seealso{
\code{\link{PMXStanFit}} for the methods \code{plot.gof.idv()}.
}

