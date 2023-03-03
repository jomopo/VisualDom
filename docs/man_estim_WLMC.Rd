\name{estim_WLMC}
\alias{estim_WLMC}
\alias{estimation_WLMC}
\title{Estimates the wavelet local multiple correlation 
}
\description{The \code{\link{estim_WLMC}} function estimates the wavelet local multiple correlation (WLMC) for a set of multivariate time series. The function is based mainly on the work of Fernández-Macho (2018) and to lesser extent in the work of Polanco-Martínez et al. (2020). 
}
\usage{
estim_WLMC(inputdata, wf="la8", J, window, M, Ymaxr=NULL)
}
\arguments{
\item{inputdata}{ 
A matrix of N columns by P rows: the first column is time (regular/evenly spaced) and the other columns are the variables under study. 
}
\item{wf}{
Name of the wavelet filter used in the wavelet transform (MODWT) decomposition. There are several wavelet filters to be used, but we use by default the Daubechies orthonormal compactly supported wavelet of length L=8, that is, ``la8''. 
}
\item{J}{
The  maximum level of the MODWT decomposition. It is recommended to use \eqn{log2(N) - 3}, where N is the number of rows or elements of \code{inputdata} (Fernández-Macho 2018, Polanco-Martínez et al. 2020). 
}
\item{window}{
Weight or window function, by the default is the Gaussian window (\emph{gaussian}) but other five window functions can be used, e.g., uniform, Bartlett's triangular, Cleveland's tricube, Wendland's truncated power or Epanechnikov's parabolic. Please look at the function \code{wave.local.multiple.correlation} from the \R package wavemulcor (Fernandez-Macho 2018).
}
\item{M}{ 
 The length of the weight/window function, it is recommended to use \eqn{N/8}, where N is the number of rows or columns (Fernández-Macho 2018, Polanco-Martínez et al. 2020).
}
\item{Ymaxr}{ 
This parameter is used to indicate which variable will be used to maximize the multiple correlation for each wavelet scale, by default is 'NULL', that is, we do not define a priori an specific variable but instead let the WLMC select one (Fernández-Macho 2018, Polanco-Martínez et al. 2020). 
}
}
\details{
The \code{\link{estim_WLMC}} function estimates the WLMC for multivariate time series including the correlation coefficients and their statistical significance. The \code{\link{estim_WLMC}} function uses the function \code{wave.local.multiple.correlation} (package:wavemulcor) to estimate the wavelet local multiple correlation, and the functions \code{modwt} and \code{brick.wall} (package:waveslim) to carry out the wavelet decomposition of the time series under study.
} 
\value{
Outputs: 
A list named \code{LISTvals} that contains four elements: \code{CORCOEF} that contains the correlation coefficients, \code{CIlo} and \code{CIup} are the lower and upper confidence intervals (CI), and \code{YmaxR} contains the indices (numbers from 1 to the number of columns or variables) of the corresponding variables whose correlation is calculated against a linear combination of the rest. 
} 
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Excellence Unit GECOS, IME, Universidad de Salamanca, Salamanca, SPAIN.\cr
BC3 - Basque Centre for Climate Change, Leioa, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue-Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} \cr
Acknowledgement:\cr 
We acknowledge to the Excellence Unit GECOS (reference number CLU-2019-03), Universidad de Salamanca for its funding support. 
}  
\references{
Fernández-Macho, J. (2018). Time-localized wavelet multiple regression and correlation. Physica A: Statistical Mechanics and its Applications, 492, 1226-1238. <URL: \doi{10.1016/j.physa.2017.11.050}>. \cr 

Polanco-Martínez, J. M., Fernández-Macho, J., & Medina-Elizalde, M. (2020). Dynamic wavelet correlation analysis for multivariate climate time series. Scientific Reports, 10(1), 1-11. <URL: \doi{10.1038/s41598-020-77767-8}>. \cr 

Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. Journal of Geophysical Research: Atmospheres, 105(D11), 14941-14962. <URL: \doi{10.1029/2000JD900110}>.  
}
\examples{
\donttest{
# We reproduce Figure 3 left presented in Polanco-Martínez et al. (2020). 
inputdata  <- rdata_climate
N          <- nrow(inputdata)
wf         <- "la8"               
window     <- "gaussian"          
J          <- 7 
M          <- 168 

output_WLMC <- estim_WLMC(inputdata, wf=wf, J=J, window=window, M=M, Ymaxr=NULL) 
}
}
\keyword{estimate_the_WLMC}
\keyword{estimate_the_Wavelet_Local_Multiple_Correlation}
