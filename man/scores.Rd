\name{Scoring functions}
\alias{extremal_score}
\alias{apl_score}
\alias{ase_score}
\title{Scoring functions}
\description{
Implementations of some scoring functions discussed in the paper.
}
\usage{
extremal_score(x, y, theta, functional = "expectile", alpha = 0.5)

apl_score(x, y, alpha = 0.5)
ase_score(x, y, alpha = 0.5)

}
\arguments{
  \item{x}{Numeric vector of forecasts}
  \item{y}{Numeric vector of realizations (same length as \code{x})}
  \item{theta}{Threshold parameter for extremal score (must be a numeric scalar)}
  \item{functional}{String, either "expectile" or "quantile"}
  \item{alpha}{Level of the quantile or expectile, must be a numeric scalar in the (0,1) interval}
}
\value{
All functions return a vector of scores (same length as \code{x} and \code{y}). Smaller scores correspond to better forecasts. 

\code{extremal_score} is the scoring function defined in Equations (10) and (12) of Ehm et al (2016). \code{apl_score} is the asymmetric piecewise scoring function for quantiles, see Equation (6) in Ehm et al (2016). \code{ase_score} is the asymmetric squared error for expectiles, see Equation (8) in Ehm et al (2016).

}
\references{
Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2016): Of Quantiles and Expectiles: Consistent Scoring Functions, Choquet Representations, and Forecast Rankings. \emph{Journal of the Royal Statistical Society (Series B)} 78, 1-29. \doi{10.1111/rssb.12154} (open access).
}
\author{
Fabian Krueger
}

\keyword{scoring functions}

