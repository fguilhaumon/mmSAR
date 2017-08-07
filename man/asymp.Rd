\name{asymp}
\Rdversion{0.1}
\alias{asymp}
\docType{data}
\title{The asymptotic SAR regression model}
\description{}
\usage{data(asymp)}
\format{
  The format is:
List of 8
 $ name       : "asymptotic regression"
 $ formula    : expression(S == c - z * f^A)
 $ paramnumber: 3
 $ paramnames : "c" "z" "f"
 $ parLim     : "Rplus" "Rplus" "unif"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation from Ratkowsky (1983) p178
 $ form       : a formula object for further calculations
}
\details{This model is convex asymptotic. Parameter c is the upper asymptote. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Tjorve (2003)}
\references{
Ratkowsky, D.A. (1983) Nonlinear regression modelling: a unified approach. Marcel Dekker, New York.

Tjorve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.

}
\seealso{
\code{\link{rssoptim}} ; \code{\link{multiSAR}} ; \code{\link{power}}
}
\examples{
#loading the asymp model object before using it
data(asymp)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
