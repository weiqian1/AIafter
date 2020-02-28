\name{rgdp}
\alias{rgdp}

\title{Real GDP growth with forecasts.}
\usage{
data(rgdp)
}

\description{
Real GDP growth data for Canada and its two-year ahead forecasts.   
}

\details{This data set contains 91 real GDP growth time points and 42 candidate forecasts, as described in Qian et al. (2020).
}

\value{
The first column gives the observed time series, and the remaining 42 columns give the candidate forecasts.
}

\references{
 Stock, J. H. and Watson, M. W. (2003), ``Forecasting Output and Inflation: The Role
of Asset Prices'', \emph{Journal of Economic Literature}, XLI, 788-829.\cr
	
Qian, W., Rolling, C. A., Cheng, G., and Yang, Y. (2020), ``Combining Forecasts for Universally Optimal Performance'', submitted. \cr

}



\examples{
library(AIafter)
data(rgdp)

# observed series
y <- rgdp[,1]
# candidate forecasts
X <- rgdp[,-1]

}


\keyword{datasets}
