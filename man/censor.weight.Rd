\name{censor.weight}
\alias{censor.weight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates censoring probability for weighting
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; calculates censoring probability needed for inverse probability of censoring weighting
}
\usage{
censor.weight(data.x, data.delta, t, weight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.x}{
%%     ~~Describe \code{data.x} here~~
numeric vector, the observed event time: X = min(T, C) where T is the time of the primary outcome, C is the censoring time
}
  \item{data.delta}{
%%     ~~Describe \code{data.delta} here~~
numeric vector of 0/1, the censoring indicator: D = I(T<C) where T is the time of the primary outcome, C is the censoring time
}
  \item{t}{
%%     ~~Describe \code{t} here~~
number, the time of interest
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
a numeric vector or matrix of weights used for perturbation-resampling, default is null.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Computes the Kaplan Meier estimate of survival for the censoring random variable at the specified time
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Kaplan Meier estimate of survival for censoring at time t

}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(ExampleData)
censor.weight(data.x = ExampleData$x1, data.delta = ExampleData$delta1, t=5)

}

\keyword{nonparametric}
\keyword{survival}
\keyword{internal}