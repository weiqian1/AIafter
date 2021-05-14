\name{AIafter}
\alias{AIafter}
\title{Forecast combination by AI-AFTER.}
\description{
This function generates forecast combination results using the AI-AFTER algorithm.
}
\usage{
AIafter(X, y, Xnew, ynew = NULL, ai.test = TRUE, safe = TRUE, sig.level = 0.1, cfi.nskip = 35, after.nskip = 5)
}
\arguments{
	\item{X}{Matrix of candidate forecasts for observed series.}
	\item{y}{Observed series.}
	\item{Xnew}{Matrix of candidate forecasts for unobserved series.}
	\item{ynew}{(Optional) True values of unobserved series, usually supplied for testing forecast accurary. If \code{ynew} is given, a forecast from \code{Xnew} can use previous observations in \code{ynew} as in any real forecasting scenarios.}
	\item{ai.test}{If \code{TRUE}, use AI-AFTER test to determine forecast direction.}
	\item{safe}{If \code{TRUE}, the prediction safeguard feature is activated.}
	\item{sig.level}{Significance level.}
	\item{cfi.nskip}{Initial training sample size for creating potentially stronger forecasts.}
	\item{after.nskip}{Number of ``burn-in'' observations for AFTER.} 
}
\details{This function first creates the potentially stronger forecasts from the CFI family, then tests whether the goal of CFA or CFI is more appropriate (optional) and uses AFTER to combine the forecasts according to the more favorable goal.
}
\value{It returns a list that contains combined forecasts (both fitted values from \code{X} and prediction from \code{Xnew}):
	\item{fcst}{a forecast matrix, with first column being forecasts from L1-AI-AFTER and second column from L2-AI-AFTER.}
	\item{pval_dm}{p-value.}
	}
\author{Wei Qian, Craig A. Rolling, Gang Cheng and Yuhong Yang\cr
Maintainer: Wei Qian  <weiqian@udel.edu>}
\references{
Qian, W., Rolling, C. A., Cheng, G., and Yang, Y. (2021), ``Combining Forecasts for Universally Optimal Performance'', International Journal of Forecasting, to appear. \cr

}

\seealso{\code{\link{AFTER1}}.}
\examples{
library(AIafter)
data(rgdp)

# observed series and its candidate forecasts
neval <- 20
n <- nrow(rgdp)
y <- rgdp[1:(n-neval),1]
X <- rgdp[1:(n-neval),-1]
# candidate forecasts for unobserved series
Xnew <- rgdp[(n-neval+1):n,-1]

# generate forecasts by AI-AFTER (both fitted and predicted values). The first cfi.nskip fitted values should be NA under CFI.
fcst <- AIafter(X,y,Xnew) 

# true values of unobserved series
ynew <- rgdp[(n-neval+1):n,1]
# generate forecasts by AI-AFTER (usually used for evaluation of forecast accuracy)
fcst <- AIafter(X,y,Xnew,ynew)

}
\keyword{combining for adaptation}
\keyword{CFA}
\keyword{combining for improvement}
\keyword{CFI}




