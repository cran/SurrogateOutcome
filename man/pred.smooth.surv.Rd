\name{pred.smooth.surv}
\alias{pred.smooth.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the conditional probability of survival for control group values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Helper function; calculates the estimated probability of survival for control 
control group surrogate values using treatment group surrogate and 
outcome information.
}
\usage{
pred.smooth.surv(xone.f, deltaone.f, sone.f, szero.one, myt, bw = NULL, 
weight, transform, extrapolate = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xone.f}{
%%     ~~Describe \code{xone.f} here~~
observed event times in the treatment group.
}
  \item{deltaone.f}{
%%     ~~Describe \code{deltaone.f} here~~
censoring indicators in the treatment group.
}
  \item{sone.f}{
%%     ~~Describe \code{sone.f} here~~
surrogate marker values in the treatment group.
}
  \item{szero.one}{
%%     ~~Describe \code{szero.one} here~~
surrogate marker values in the control group.
}
  \item{myt}{
%%     ~~Describe \code{myt} here~~
time of interest.
}
  \item{bw}{
%%     
bandwidth.
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
weight used for perturbation resampling.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether a transformation should be used, default is FALSE.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used, default is TRUE.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned:                
 \item{Phat.ss}{conditional probability of survival past t for control group}
 \item{warn.flag}{warning flag equal to 1 if extrapolation was used; 0 otherwise}
}

\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{survival}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{smooth}
\keyword{robust}
\keyword{internal}