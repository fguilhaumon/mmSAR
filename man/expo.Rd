\name{expo}
\Rdversion{1.1}
\alias{expo}
\docType{data}
\title{The exponetial SAR regression model}
\description{s == z * log(a) + c}
\usage{data(expo)}
\format{
  The format is:
List of 8
 $ name       : "Exponential"
 $ formula    : expression(s == z * log(a) + c)
 $ paramnumber: 2
 $ paramnames : "c" "z"
 $ parLim     : "R" "Rplus"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{This model is convex non asymptotic. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Gleason (1922, 1925), Fisher et al. (1943), Tjørve (2003).
}
\references{
Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943) The relation between the number of species and the number of individuals in a random sample of an animal population. Journal of Animal Ecology, 12, 42–58.

Gleason, H.A. (1922) On the relation between species and area. Ecology, 3, 158–162.

Gleason, H.A. (1925) Species and area. Ecology, 6, 66–74.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the expo model object before using it
data(expo)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
