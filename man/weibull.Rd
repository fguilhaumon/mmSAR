\name{weibull}
\Rdversion{1.1}
\alias{weibull}
\docType{data}
\title{The cumulative weibull distribution SAR regression model (3 parameters)}
\description{S == c(1 - exp(-z * A^f))
}
\usage{data(weibull)}
\format{
  The format is:
List of 8
 $ name       : "Cumulative Weibull"
 $ formula    : expression(S == c(1 - exp(-z * A^f)))
 $ paramnumber: 3
 $ paramnames : "c" "z" "f"
 $ parLim     : "Rplus" "Rplus" "Rplus"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The cumulative Weibull distribution (Weibull, 1951) is sigmoid asymptotic. Parameter C is the upper asymptote.It goes through the origin of axes.
}
\source{Weibull (1951), Tjørve (2003).
}
\references{
Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.

Weibull, W. (1951) A statistical distribution function of wide applicability. Journal of Applied Mathematics, 18, 293–296.

}
\examples{
#loading the weibull model object before using it
data(weibull)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
