\name{plot_dominant_WLMC}
\alias{plot_dominant_WLMC}
\alias{plot_dominante_WLMC}
\title{Plot as a heat map the 'dominant' variable/s contained in the output of the function 'estim_WLMC'
}
\description{The \code{\link{plot_dominant_WLMC}} function plots as a heat map the 'dominant' (that is, the variable/s that maximizes the multiple correlation through time and scale) variable/s of the time series analysed. This function also plots the time series under analysis and discern the correlation coefficients that are not statistically significant, which are plotted as blanks.
}
\usage{
plot_dominant_WLMC(inputdata, LISTvals, J, fac=1, FLAG=TRUE, FLAGNA=1, 
                   COLS=c(1:5), LTY=c(rep(1,5)), LWD=c(rep(1.2,5)), 
                   DIST=c(seq(0,10,2)))
}
\arguments{
\item{inputdata}{ 
A matrix of N columns by P rows: the first column is time (regular/evenly spaced) and the other columns are the variables under study. 
}
\item{LISTvals}{
Output of the function \code{\link{estim_WLMC}}. This is a list that contains: (1) CORCOEF (correlation coefficients), (2) CIlo (CI lower bounds), (3) CIup (CI upper bounds), and (4) YmaxR (the index numbers, from 1 to number of variables, of the variable whose correlation is calculated against a linear combination of the rest).
}
\item{J}{
The  maximum level of the MODWT decomposition. It is recommended to use \eqn{log2(N) - 3}, where N is the number of rows or elements of \code{inputdata} (Fernández-Macho 2018, Polanco-Martínez et al. 2020). 
}
\item{fac}{
This factor is used to scale the wavelet time-scales or 'periods' when the time scale is not the unit, by the default is 1. 
}
\item{FLAG}{
This 'flag' is used to plot the Y axis of the multivariate time series if the number of these series is less than four, by default is TRUE. 
}
\item{FLAGNA}{
This is used to plot (by the default is 1) or not (please us 0) the correlation coefficients that are/not statistically significant. 
}
\item{COLS}{
The colors used to plot the multivariate time series. 
}
\item{LTY}{
The type of lines used to plot the multivariate time series.
}
\item{LWD}{
The tick sizes used to plot the multivariate time series. 
}
\item{DIST}{
This parameter is used to define the distances between the Y axis when the multivariate time series are plotted. 
}
}
\details{
The \code{\link{plot_dominant_WLMC}} function plots as a heat map the 'dominant' variable/s, that is, the variable/s that maximizes the multiple correlation through time and scale of the variables (time series) analysed. This function also plot the time series under analysis and discern the correlation coefficients that are not statistically significant, we use blanks to plot these coefficients.
}
\value{
A plot of the time series under analysis and a heat map (a multi-plot via screen) of the dominant variable/s whose correlation coefficients are statistically significant. This multi-plot can be saved in your preferred format.
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

Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. Journal of Geophysical Research: Atmospheres, 105(D11), 14941-14962. <URL: \doi{10.1029/2000JD900110}>. \cr 

Polanco-Martínez, J.M. (2023, in press), A computational and graphical approach to analyze the dynamic wavelet correlation among components of a nonlinear dynamical system, Journal of Applied Nonlinear Dynamics, 1-13.  
}
\examples{
# We reproduce Figure 3 below presented in Polanco-Martínez (2023). 
#\donttest{
inputdata  <- rdata_Lorenz 
N          <- nrow(inputdata)
wf         <- "la8"               
window     <- "gaussian"          
J          <- 6 
M          <- 20 

LISTvals <- estim_WLMC(inputdata, wf=wf, J=J, window=window, M=M, Ymaxr=NULL) 

plot_dominant_WLMC(inputdata, LISTvals=LISTvals, J=J, fac=0.05, 
  FLAG=TRUE, FLAGNA=1, COLS=c("blue", "green", "red"), 
  LTY=c(rep(1,5)), LWD=c(rep(1.2,5)), DIST=c(seq(0, 10, 2.75)))
#}
}
\keyword{plot_dominant_Wavelet_Local_Multiple_Correlation}
\keyword{plot_dominante_Wavelet_Local_Multiple_Correlation}
 
