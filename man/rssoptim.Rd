\name{rssoptim}
\Rdversion{1.1}
\alias{rssoptim}
\title{The Residual Sum of Squares optimisation function}
\description{The rssoptim function is used to fit SAR models (as mmSAR model objects) to dataset (as mmSAR data objects). This function allow the optimisation of the model parameters, fitted values calculation, information criteria (AIC, AICc, BIC) and R² calculations and fit evaluation with residuals normality and homoscedasticity tests.
}
\usage{
rssoptim(model, data, norTest = "lillie", verb = TRUE)
}
\arguments{
  \item{model}{A mmSAR model object representing the model to fit.
}
  \item{data}{A mmSAR data object representing the data to fit the model on.
}
  \item{norTest}{A character string from "lillie" or "shapiro" for respectively test for the normality of the residuals with a Lilliefors or Shapiro test. 
}
  \item{verb}{Bolean : should the function reports informations while running ?
}
}
\details{
The rssoptim function is basically a wrapper for the optim R built-in function. rssoptim obtains model parameters estimates by minimizing the residual sum of squares using the unconstrained Nelder–Mead optimization algorithm. Assuming normality of the observations, this approach produces optimal maximum likelihood estimates of model parameters. rssoptim provides two tests for the normality of the residuals : the Lilliefors extension of the Kolmogorov normality test, which is advocated when sample size is large or when the data show a substantial variability (e.g. continental scale studies) and the Shapiro-Wilk test for normality, which focuses on skewness and kurtosis of the empirical distribution of the residuals and is useful for small sample size or when data results from small scale sampling. mmSAR tests for homoscedasticity by evaluating the correlation between residuals magnitude and areas (Pearson’s product moment correlation coefficient).
}
\value{
A list :
\enumerate{
\item $par : the parameters estimates.

\item $value : the value of the Residual Sum of Squares (RSS).

\item $counts : number of iterations for the convergence of the fitting algorithm

\item $convergence : optim convergence code (only 0 is OK).

\item $message : optim convergence message (only NULL is OK).

\item $sigConf : a p x 6 matrix with columns for the estimated coefficient, its standard error, t-statistic, corresponding (two-sided) p-value and lower & upper bounds of the 95% confidence interval.

\item $pearson : Pearson’s product moment correlation coefficient value.

\item $pearpval : Pearson’s product moment correlation coefficient p.value.

\item $normaTest : residuals normality test used.

\item $normaStat : residuals normality test statistic value.

\item $normaPval : residuals normality test statistic p.value.

\item $startvalues : starting values for the fitting algorithm.

\item $data : the data.

\item $model : the mmSAR model object.

\item $calculated : vector of fitted values.

\item $residuals : vector of residuals.

\item $AIC : AIC value.

\item $AICc : AICc value.

\item $BIC : BIC value.

\item $R2a : R² value.

\item $data.name : the name of the dataset (from the mmSAR data object).
}

}

\author{Francois Guilhaumon
}

\seealso{\code{\link{multiSAR}}
}

\examples{

\dontrun{
#FITTING THE POWER MODEL TO THE GALAPAGOS DATASET

#loading the data
data(data.galap)

#loading the model
data(power)

#fitting the model with verbosity
rssoptim(power,data.galap,verb=TRUE)
}


}
