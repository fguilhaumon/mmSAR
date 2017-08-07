\name{ratio}
\Rdversion{1.1}
\alias{ratio}
\docType{data}
\title{The rational function SAR regression model}
\description{S == (c + z * A) / (1 + f * A)}
\usage{data(ratio)}
\format{
  The format is:
List of 8
 $ name       : "Rational function"
 $ formula    : expression(S == over((c + z * A), (1 + f * A)))
 $ paramnumber: 3
 $ paramnames : "c" "z" "f"
 $ parLim     : "R" "R" "R"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The rational function model (Ratkowsky, 1990) is convex asymptotic. Each parameter in the model affects several curve aspects (shape, upper asymptote and y-axis intersection).
}
\source{Ratkowsky (1990), Tjørve (2003).
}
\references{
Ratkowsky, D.A. (1990) Handbook of nonlinear regression models. Marcel Dekker, New York.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the ratio model object before using it
data(ratio)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
