\name{negexpo}
\Rdversion{1.1}
\alias{negexpo}
\docType{data}
\title{The negative exponential SAR regression model}
\description{s == c * (1 - exp(-z * a))}
\usage{data(negexpo)}
\format{
  The format is:
List of 8
 $ name       : "Negative exponential"
 $ formula    : expression(s == c * (1 - exp(-z * a)))
 $ paramnumber: 8
 $ paramnames : "c" "z"
 $ parLim     : "Rplus" "unif"
 $ fun        : model function
 $ rssfun     : Residual Sum of Squares function
 $ init       : initial values calculation
 $ form       : a formula object for further calculations
}
\details{The negative exponential model is convex asymptotic. For further details and documentation please see the online documentation at http://mmsar.r-forge.r-project.org/
}
\source{Holdridge et al. (1971), Miller & Wiegert (1989), Ratkowsky (1990), Tjørve (2003).
}
\references{
Holdridge, L.R., Grenke, W.C., Hatheway, W.H., Liang, T. & Tosi, J.A. (1971) Forest environments in tropical life zones. Pergamon Press, Oxford.

Miller, R.I. & Wiegert, R.G. (1989) Documenting completeness, species–area relations, and the species abundance distribution of a regional flora. Ecology, 70, 16–22.

Ratkowsky, D.A. (1990) Handbook of nonlinear regression models. Marcel Dekker, New York.

Tjørve, E. (2003) Shapes and functions of species–area curves: a review of possible models. Journal of Biogeography, 30, 827–835.
}
\examples{
#loading the negexpo model object before using it
data(negexpo)
}
\keyword{mmSAR model object}
\keyword{non linear regression}
