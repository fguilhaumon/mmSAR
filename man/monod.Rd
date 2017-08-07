\name{monod}
\Rdversion{1.1}
\alias{monod}
\docType{data}
\title{The Monod function SAR regression model}
\description{s == c * a /z + a}
\usage{data(monod)}
\format{
  The format is:
List of 8
 $ name       : "Monod"
 $ formula    : expression(s == over(c * a, z + a))
 $ paramnumber: 2
 $ paramnames : "c" "z"
 $ parLim     : "Rplus" "Rplus"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The Monod function SAR model is convex asymptotic. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Monod (1950), Tjørve (2003).
}
\references{
Monod, J. (1950) La technique de culture continue, the´orie et applications. Annales de l’Institut Pasteur, 79, 390–410.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the monod model object before using it
data(monod)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
