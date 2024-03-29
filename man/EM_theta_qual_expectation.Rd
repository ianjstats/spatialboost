% Generated by roxygen2 (4.0.2): do not edit by hand
\name{EM_theta_qual_expectation}
\alias{EM_theta_qual_expectation}
\title{EM E-step on theta}
\usage{
EM_theta_qual_expectation(KAPPA, GENEBOOST, BETA, SIG2_0)
}
\arguments{
\item{KAPPA}{the value of kappa in the prior on beta}

\item{GENEBOOST}{the negative gene boost of each marker}

\item{BETA}{the current value of beta}

\item{SIG2_0}{the current value of sigma squared}
}
\description{
This function performs the E-step of our EM algorithm for the
spatial boost model and computes the expected value of each
indicator variable in our model (the theta's)
}
\examples{
EM_theta_qual_expectation()
}
\keyword{Expected}
\keyword{of}
\keyword{theta}
\keyword{value}

