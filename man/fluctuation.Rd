\name{Fluctuation Test}
\alias{fluctuation_test}
\title{Fluctuation test}
\description{
Test to analyze whether the ranking of two forecasts is stable over time. The variant implemented here has been proposed in Proposition 1 of Giacomini and Rossi (2010); the critical values are tabulated in their Table 1. The null hypothesis of the test is that both forecasting methods perform equally well (same expected score) at all time points. The alternative is that their performance differs in at least one time point. 
}
\usage{
fluctuation_test(loss1, loss2, mu = 0.5, dmv_fullsample = TRUE,
                 lag_truncate = 0, time_labels = NULL,
                 conf_level = 0.05)
}
\arguments{
  \item{loss1, loss2}{Vectors of losses corresponding to two forecast methods (smaller losses correspond to better forecasts).}
  \item{mu}{Size of the rolling window (relative to evaluation sample). Must be in 0.1, 0.2, ..., 0.9.}
  \item{dmv_fullsample}{Logical; if \code{TRUE} (the default), the full sample is used to estimate the variance of the Diebold-Mariano type
  statistic employed in the test. See page 14/footnote 16 in the working paper version of Rossi (2013).}
  \item{lag_truncate}{Truncation lag used when estimating the variance of the Diebold-Mariano type test statistic.}
  \item{time_labels}{Vector of labels to be used for the time axis. If \code{NULL} (the default), integer labels are used.}
  \item{conf_level}{Confidence level, either \code{0.05} or \code{0.1} (two-sided test).}
}
\value{
List with two elements: 1) Data frame containing the time path of the test statistic, and 2) the relevant critical values. In addition, the function draws a plot which illustrates the test.
}
\references{
Giacomini, R. and Rossi, B. (2010): Forecast Comparisons in Unstable Environments. Journal of Applied Econometrics 25, 595-620.  

Rossi, B. (2013): Advances in Forecasting under Model Instability. In: Handbook of Economic Forecasting, vol. 2, Graham Elliott and Alan Timmermann (eds), pp. 1203-1324. 
}
\author{
Fabian Krueger
}

\keyword{forecasting methods}

\examples{

# Comparison of Inflation Forecasts: 
# Survey of  Professional Forecasters (SPF) 
# versus Michigan Survey of Consumers

data(inflation_mean)

# Compute extremal scores of SPF/Michigan (theta = 3)
score_spf <- extremal_score(x = inflation_mean$spf, 
                            y = inflation_mean$rlz, theta = 3)
score_michigan <- extremal_score(x = inflation_mean$michigan, 
                                 y = inflation_mean$rlz, theta = 3)

# Make simplified label for time axis
tml <- as.numeric(substr(inflation_mean$dt, 1, 4))

# Fluctuation test
fluct_test <- fluctuation_test(score_spf, score_michigan, 
                               time_labels = tml, lag_truncate = 4)

}
