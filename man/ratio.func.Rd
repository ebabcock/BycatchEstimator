% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bycatchFunctions.R
\name{ratio.func}
\alias{ratio.func}
\title{Basic ratio estimator with variance (Cochran)}
\usage{
ratio.func(x, y, g, X, N, G)
}
\arguments{
\item{x}{x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.}

\item{y}{x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.}

\item{g}{x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.}

\item{X}{X is the total effort/catch by stratum}

\item{N}{N is the total number of sample units by stratum, if available, otherwise total effort}

\item{G}{Value}
}
\description{
Output is mean and standard error of bycatch by stratum and the total bycatch with SE. Assumes unobserved strata have zero catch
}
\keyword{internal}
