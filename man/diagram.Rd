\name{Murphy Diagrams}
\alias{murphydiagram}
\alias{murphydiagram_diff}
\title{Murphy diagrams to visualize forecast comparisons}
\description{
Visual comparisons of two forecasting methods, allowing to study whether the ranking is robust across the class of elementary or extremal
scoring functions. See Ehm et al (2015, esp. Sections 3 and 4) for details. 
}
\usage{
murphydiagram(f1, f2, y, functional = "expectile", 
			        alpha = 0.5, labels = c("Method 1", "Method 2"), 
			        colors =  NULL, equally_spaced = FALSE)

murphydiagram_diff(f1, f2, y, functional = "expectile", 
				           alpha = 0.5, equally_spaced = FALSE, 
				           lag_truncate = 0, conf_level = 0.95)
}
\arguments{
  \item{f1, f2}{Vectors of point forecasts}
  \item{y}{Vector of realizing observations.}
  \item{functional}{Either "expectile" (the default) or "quantile". Note that the probability of a binary event is an expectile at level \code{alpha} = 0.5 (see below).}
  \item{alpha}{Level of the expectile or quantile, must be between 0 and 1. Defaults to 0.5, which is the mean (if functional is set to "expecile") or median (if functional is set to "quantile"). }
  \item{labels}{Method labels for \link{murphydiagram} to be used in plot legend. Character vector of length two, or \code{NULL} (in order to omit labels).}
  \item{colors}{Colors used. Defaults to NULL, such that the colors are as in Ehm et al (2015). Alternative colors can be specified as a character vector of length two.}
  \item{equally_spaced}{Method for choosing the grid of values on the horizontal axis. If set to FALSE (the default), the set of points that is relevant for dominance (c.f. Section 3.4 of the paper) is chosen. This can be somewhat time consuming for large data sets. If set to TRUE, an auxiliary grid of equally spaced points is used.}
  \item{lag_truncate}{Largest order of autocorrelation that is accounted for in the variance estimator for \link{murphydiagram_diff} (defaults to zero).}
  \item{conf_level}{Level of the confidence bands plotted in \link{murphydiagram_diff}, defaults to 0.95.}
}
\value{
None, used for the effect of creating a plot.
}
\references{
Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2015): Of Quantiles and Expectiles: Consistent Scoring Functions, Choquet Representations, and Forecast Rankings. Preprint available at \url{http://arxiv.org/abs/1503.08195}.
}
\author{
Fabian Krueger
}

\keyword{forecasting methods}

\examples{

\dontrun{
# Comparison of Inflation Forecasts: Survey of Professional Forecasters (SPF) 
# versus Michigan Survey of Consumers

data(inflation_mean)
murphydiagram(inflation_mean$spf, inflation_mean$michigan, inflation_mean$rlz, 
              labels = c("SPF", "Michigan"))
murphydiagram_diff(inflation_mean$spf, inflation_mean$michigan, inflation_mean$rlz)
}

}
