\name{CLR1}
\alias{CLR1}
\title{Forecast combination by constrained linear regression.}
\description{
This function uses constrained linear regression (CLR) to find the optimal combination of all forecasts.
}
\usage{
CLR1(X, y, Xnew, ynew = NULL, n0 = 5)
}
\arguments{
	\item{X}{Matrix of candidate forecasts for observed series.}
	\item{y}{Observed series.}
	\item{Xnew}{Matrix of candidate forecasts for unobserved series.}
	\item{ynew}{(Optional) True values of unobserved series, usually supplied for testing forecast accurary. If \code{ynew} is given, a forecast from \code{Xnew} can use previous observations in \code{ynew} as in any real forecasting scenarios.} 
	\item{n0}{Number of ``burn-in'' observations for CLR.} 

}
\details{This function generates forecasts from CLR, with the constraint that all weights are non-negative and sum to unity.}
\value{A vector with combined forecasts (both fitted values from \code{X} and prediction from \code{Xnew}).	
	}
\author{Wei Qian, Craig A. Rolling, Gang Cheng and Yuhong Yang\cr
Maintainer: Wei Qian  <weiqian@udel.edu>}
\references{
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

# generate forecasts by CLR (both fitted and predicted values)
fcst <- CLR1(X,y,Xnew)

# true values of unobserved series
ynew <- rgdp[(n-neval+1):n,1]
# generate forecasts by CLR (usually used for evaluation of forecast accuracy)
fcst <- CLR1(X,y,Xnew,ynew)

}

