\name{power}

\alias{power}
\docType{data}
\title{The power law SAR model}
\description{\deqn{S = cA^{z}}{%
s = c * a^z}}
\usage{data(power)}
\format{
  The format is:
List of 9
 $ name       : "Power"
 $ formula    : expression(s == c * a^z)
 $ paramnumber: 2
 $ parLim     : "Rplus" "unif"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The power SAR mdoel is convex non asymptotic. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Arrhenius (1921), Tjørve (2003).
}
\references{
Arrhenius, O. 1921. Species and area. Journal of Ecology. 9: 95-99.

Tjørve, E. 2003. Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the power model object before using it
data(power)
}
\keyword{mmSAR model object}
\keyword{non linear regression}

