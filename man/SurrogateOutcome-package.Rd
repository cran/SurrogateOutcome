\name{SurrogateOutcome-package}
\alias{SurrogateOutcome-package}
\alias{SurrogateOutcome}
\docType{package}
\title{
\packageTitle{SurrogateOutcome}
}
\description{
\packageDescription{SurrogateOutcome}
}
\details{

This package implements all methods proposed in "Assessing the Value of a Censored Surrogate Outcome" by Parast L, Tian L, and Cai T which is currently in press at Lifetime Data Analysis. The main functions are (1) R.q.event() which calculates the proportion of the treatment effect (the difference in restricted mean survival time at time t) explained by surrogate outcome information observed up to a selected landmark time, (2) R.t.estimate() which calculates the proportion of the treatment effect explained by primary outcome information only observed up to a selected landmark time, and (3) IV.event() which calculates the incremental value of the surrogate outcome information.
}
\author{
\packageAuthor{SurrogateOutcome}

Maintainer: \packageMaintainer{SurrogateOutcome}
}
\references{
Parast L, Tian L, and Cai T (2019). Assessing the Value of a Censored Surrogate Outcome. Lifetime Data Analysis, In press.}

\keyword{ package }

\examples{
data(ExampleData)
names(ExampleData)

\donttest{
R.q.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np")
R.t.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, landmark=2)
IV.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np")
R.q.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np", std = TRUE, conf.int = TRUE)
R.t.estimate(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, t = 5, landmark=2, std = TRUE, conf.int = TRUE)
IV.event(xone = ExampleData$x1, xzero = ExampleData$x0, deltaone = ExampleData$delta1, 
deltazero = ExampleData$delta0, sone = ExampleData$s1, szero = ExampleData$s0, t = 5, 
landmark=2, type = "np", std = TRUE, conf.int = TRUE)
}

}
