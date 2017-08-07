\name{lomolino}
\Rdversion{1.1}
\alias{lomolino}
\docType{data}
\title{The lomolino SAR regression model}
\description{S == c / 1 + z^log(f/A)}
\usage{data(lomolino)}
\format{
  The format is:
List of 8
 $ name       : "Lomolino"
 $ formula    : expression(S == over(c, 1 + z^log(over(f, A))))
 $ paramnumber: 3
 $ paramnames : "c" "z" "f"
 $ parLim     : "Rplus" "Rplus" "Rplus"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The lomolino model is sigmoid asymptotic. Parameter c is the upper asymptote. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Lomolino, M.V. (2000), Tjørve (2003).
}
\references{
Lomolino, M.V. (2000) Ecology’s most general, yet protean pattern: the species–area relationship. Journal of Biogeography, 27, 17–26.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the lomolino model object before using it
data(lomolino)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
