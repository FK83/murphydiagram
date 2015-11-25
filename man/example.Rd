\name{Synthetic Example}
\alias{expected_score_mean}
\alias{expected_score_quantile}
\title{Analytical Expressions from the Synthetic Example in Section 3.3 and Appendix B}
\description{
Functions to compute the analytical expressions in Table 3 of the paper by Ehm et al (2015). These expressions yield the expected score of various forecasters, given the synthetic setup 
studied in Section 3.3 and Appendix B of the paper. The expressions can be used to replicate Figure 2 in the paper.
}
\usage{
expected_score_mean(theta, forecaster = "P")
expected_score_quantile(theta, alpha, forecaster = "P")
}
\arguments{
  \item{theta}{Value of the parameter $theta$, indexing the extremal score}
  \item{alpha}{Quantile level, between zero and one}
  \item{forecaster}{ID of the forecaster, string of length one. Either "P" (perfect forecaster), "C" (climatological forecaster), "U" (unfocused forecaster), or "SR" (sign-reversed forecaster).}
}
\value{
Expected value of the extremal score, given the synthetic setup described in Section 3.3 of Ehm et al (2015).
}
\references{
Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2015): Of Quantiles and Expectiles: Consistent Scoring Functions, Choquet Representations, and Forecast Rankings. Preprint available at \url{http://arxiv.org/abs/1503.08195}.
}
\author{
Alexander Jordan, Fabian Krueger
}

\keyword{synthetic example}

\examples{

\dontrun{
# Color palette, obtained from http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73")
cbbPalette <- cbbPalette[c(1, 4, 2, 3)]

# Labeling stuff
forecasters <- c("P", "C", "U", "SR")
names <- c("Perfect", "Climatological", "Unfocused", "Sign-Reversed")
x_label <- expression(paste("Parameter ", theta))

# Figure 2, top left

# Grid for theta
theta_grid1 <- seq(-3, 3, 0.01)
# Expected scores for all forecasters
scores1 <- sapply(forecasters, expected_score_mean, theta = theta_grid1)
# Plot
matplot(x = theta_grid1, y = scores1[, 4:1], type = "l", lty = 1, col = cbbPalette[4:1], 
        lwd = 2, bty = "n", xlab = x_label, ylab = expression("Expected Score"))
legend("topright", names, col = cbbPalette, lwd = 2, bty = "n")
}

}
