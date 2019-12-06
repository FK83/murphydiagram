\name{Data sets}
\alias{inflation_mean}
\alias{recession_probability}
\docType{data}
\title{
Data sets with forecasts and realizations
}
\description{
Data sets with forecasts and corresponding realizations, as used in the paper by Ehm et al (2016). In the \link{inflation_mean} data, the outcome variable is continuous; in the \link{recession_probability}
data, the outcome is binary.}
\usage{
data(inflation_mean)
data(recession_probability)
}
\format{
Both data sets are data frames, with the following layout: First column contains the quarterly date, in string format (e.g. "1998Q4" for the fourth quarter of 1998). The second and third columns contain forecasts by two alternative methods. The fourth column contains realizations.}

\source{
Forecasts are generated as described in Section 4 of Ehm et al (2016).

Data sources: Inflation - ``spf'' forecasts and realizations based on data from the Federal Reserve Bank of Philadelphia, \url{http://www.phil.frb.org/research-and-data/real-time-center/} 
(individual-level CPI forecasts, and real-time data for CPI realizations). ``michigan'' forecasts based on data from the Michigan Survey of Consumers, \url{https://data.sca.isr.umich.edu/tables.php}, Table 32. Recessions - ``spf'' forecasts and realizations based on data from the 
Federal Reserve Bank of Philadelphia, \url{http://www.phil.frb.org/research-and-data/real-time-center/} (``anxious index'' and real-time data for real GDP growth). The Probit forecasts uses the same real-time data on GDP growth, as
well as interest rate data from the Federal Reserve Bank of St. Louis, \url{http://research.stlouisfed.org/fred2/} (series TB3MS and GS10).

\emph{Disclaimer: The providers of the raw data take no responsibility for the accuracy of the forecast and realization data sets posted here. Furthermore, the raw data may be revised over time, and the websites linked above 
should be consulted for the official, most recent versions.}

Code and raw data to construct the two data sets can be found at \url{https://sites.google.com/site/fk83research/code}.

}

\examples{

\dontrun{

# Load inflation forecasts
data(inflation_mean)

# Make numeric time axis
tm <- as.numeric(substr(inflation_mean$dt, 1, 4)) + 
      0.25*(as.numeric(substr(inflation_mean$dt, 6, 6))-1)

# Plot
matplot(x = tm, y = inflation_mean[,2:4], type = "l", bty = "n",
        xlab = "Time", ylab= "Inflation (percent)", col = 3:1)
legend("topright", legend = c("SPF", "Michigan", "Actual"), fill = 3:1, bty = "n")

}

}

\references{
Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2016): Of Quantiles and Expectiles: Consistent Scoring Functions, Choquet Representations, and Forecast Rankings. \emph{Journal of the Royal Statistical Society (Series B)} 78, 1-29. \doi{10.1111/rssb.12154} (open access).

}

\keyword{datasets}