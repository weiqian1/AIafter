\name{AFTER1}
\alias{AFTER1}
\title{Forecast combination by AFTER.}
\description{
This function generates forecast combination results using the AFTER algorithm (Yang, 2004).
}
\usage{
AFTER1(X, y, Xnew, ynew = NULL, n0 = 5)
}
\arguments{
	\item{X}{Matrix of candidate forecasts for observed series.}
	\item{y}{Observed series.}
	\item{Xnew}{Matrix of candidate forecasts for unobserved series.}
	\item{ynew}{(Optional) True values of unobserved series, usually supplied for testing forecast accurary. If \code{ynew} is given, a forecast from \code{Xnew} can use previous observations in \code{ynew} as in any real forecasting scenarios.}
	\item{n0}{Number of ``burn-in'' observations for AFTER.} 
}
\details{This function generates forecasts from the AFTER algorithm. For detailed updating scheme, see Yang (2004).}
\value{It returns a list that contains combined forecasts (both fitted values from \code{X} and prediction from \code{Xnew}):
	\item{fcst}{a forecast matrix, with first column being forecasts from L1-AFTER and second column from L2-AFTER.}
	\item{L1Wt}{weight matrix used for L1-AFTER.}
	\item{L2Wt}{weight matrix used for L2-AFTER.}	
	}
\author{Wei Qian, Craig A. Rolling, Gang Cheng and Yuhong Yang\cr
Maintainer: Wei Qian  <weiqian@udel.edu>}
\references{
Yang, Y. (2004), ``Combining Forecasting Procedures: Some Theoretical Results,'' \emph{Econometric Theory}, 20, 175-222.\cr

Qian, W., Rolling, C. A., Cheng, G., and Yang, Y. (2021), ``Combining Forecasts for Universally Optimal Performance'', International Journal of Forecasting, to appear. \cr
}

\seealso{\code{\link{AIafter}}.}
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

# generate forecasts by AFTER (both fitted and predicted values)
fcst <- AFTER1(X,y,Xnew)

# true values of unobserved series
ynew <- rgdp[(n-neval+1):n,1]
# generate forecasts by AFTER (usually used for evaluation of forecast accuracy)
fcst <- AFTER1(X,y,Xnew,ynew)

}
\keyword{combining for adaptation}
\keyword{CFA}
