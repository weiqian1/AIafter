\name{BG1}
\alias{BG1}
\title{Forecast combination by the Bates-Granger method.}
\description{
This function generates forecast combination results using Bates-Granger (BG) method, without considering the correlations between candiate forecasts.
}
\usage{
BG1(X, y, Xnew, ynew = NULL, n0 = 5, rhos = c(1,0.9))
}
\arguments{
	\item{X}{Matrix of candidate forecasts for observed series.}
	\item{y}{Observed series.}
	\item{Xnew}{Matrix of candidate forecasts for unobserved series.}
	\item{ynew}{(Optional) True values of unobserved series, usually supplied for testing forecast accurary. If \code{ynew} is given, a forecast from \code{Xnew} can use previous observations in \code{ynew} as in any real forecasting scenarios.} 
	\item{n0}{Number of ``burn-in'' observations for BG.} 
	\item{rhos}{A vector of discount parameters.} 


}
\details{This function generates forecasts from the BG method.}
\value{A matrix with combined forecasts (both fitted values from \code{X} and prediction from \code{Xnew})	
	}
\author{Wei Qian, Craig A. Rolling, Gang Cheng and Yuhong Yang\cr
Maintainer: Wei Qian  <weiqian@udel.edu>}
\references{
Qian, W., Rolling, C. A., Cheng, G., and Yang, Y. (2020), ``Combining Forecasts for Universally Optimal Performance'', submitted. \cr

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
# candidate forecasts for unobserved observations
Xnew <- rgdp[(n-neval+1):n,-1]

# generate forecasts by BG (both fitted and predicted values)
fcst <- BG1(X,y,Xnew)

# true values of unobserved series
ynew <- rgdp[(n-neval+1):n,1]
# generate forecasts by BG (usually used for evaluation of forecast accuracy)
fcst <- BG1(X,y,Xnew,ynew)

}
