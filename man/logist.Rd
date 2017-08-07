\name{logist}
\Rdversion{1.1}
\alias{logist}
\docType{data}
\title{The logistic SAR regression model}
\description{S == over(c, (1 + exp(-z * A + f)))}
\usage{data(logist)}
\format{
  The format is:
List of 8
 $ name       : "Logistic function"
 $ formula    : expression(S == over(c, (1 + exp(-z * A + f))))
 $ paramnumber: 3
 $ paramnames : "c" "z" "f"
 $ parLim     : "Rplus" "Rplus" "R"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{This model is sigmoid asymptotic. Parameter c is the upper asymptote. The function is symmetrical about the inflection point. The curve has a lower asymptote of zero, and the intersection with the y-axis depends on all three parameters. Thus it never goes through the origin.  For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Ratkowsky (1990), Tjørve (2003).
}
\references{

Ratkowsky, D.A. (1990) Handbook of nonlinear regression models. Marcel Dekker, New York.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the logist model object before using it
data(logist)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
