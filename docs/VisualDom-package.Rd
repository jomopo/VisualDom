\name{VisualDom-package}
\alias{VisualDom-package}
\alias{VisualDom}
\docType{package}
\title{Visualize Dominant Variables in Wavelet Multiple Correlation
}
\description{
'VisualDom' estimates and plots as a heat map the correlation coefficients obtained via the wavelet local multiple correlation 'WLMC' (Fernández-Macho 2018). We improve the graphical outputs of WLMC proposing a didactic and useful way to visualize the 'dominant' variable(s) that the maximizes the multiple correlation through time for a set of time series (Polanco-Martínez et al. 2020; Polanco-Martínez 2023). The WLMC was designed for financial time series, but other kinds of data (e.g., climatic, ecological, etc.) can be used. 'VisualDom' contains three functions: (1) 'estim_WLMC': estimates the wavelet local multiple correlation; (2) 'plot_estim_WLMC': plots the time series under analysis and a heat map the correlation coefficients obtained via 'estim_WLMC'; (3) 'plot_dominant_WLMC': this function also plot the time series under study and a heat maps of the dominant variable/s using the outputs of 'estim_WLMC'. The functions contained in 'VisualDom' are highly flexible since these contains several parameters to personalize the time series under analysis and the heat maps. In addition, we have also included two data sets (named 'rdata_climate' and 'rdata_Lorenz') to exemplify the use of the functions contained in 'VisualDom'. 
}
\details{
\tabular{ll}{
Package: \tab VisualDom\cr
Type: \tab Package\cr
Version: \tab 0.8\cr
Date: \tab 2022-12-21\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
VisualDom package contains three functions: (1) \code{\link{estim_WLMC}} that estimates the wavelet local multiple correlation (WLMC); (2) \code{\link{plot_estim_WLMC}} that plots the time series under analysis and the correlation coefficients as a heat map obtained via \code{\link{estim_WLMC}}; and (3) \code{\link{plot_dominant_WLMC}} that plots the time series under study and the dominant variable/s as a heat maps, this function also use the output of \code{\link{estim_WLMC}}. }
\note{
 Dependencies: \emph{wavemulcor}, \emph{waveslim}, \emph{plot3D}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Excellence Unit GECOS, IME, Universidad de Salamanca, Salamanca, SPAIN.\cr
BC3 - Basque Centre for Climate Change, Leioa, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue-Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} \cr
Acknowledgement:\cr 
We acknowledge to the Excellence Unit GECOS (reference number CLU-2019-03), Universidad de Salamanca for its funding support. }  
\references{
Fernández-Macho, J. (2018). Time-localized wavelet multiple regression and correlation. Physica A: Statistical Mechanics and its Applications, 492, 1226-1238. <URL: \doi{10.1016/j.physa.2017.11.050}>. \cr 

Polanco-Martínez, J. M., Fernández-Macho, J., & Medina-Elizalde, M. (2020). Dynamic wavelet correlation analysis for multivariate climate time series. Scientific Reports, 10(1), 1-11. <URL: \doi{10.1038/s41598-020-77767-8}>. \cr 

Polanco-Martínez, J.M. (2023, in press). A computational and gaphical approach to analyze the dynamic wavelet correlation among   components of a nonlinear dynamical system. Journal of Applied Nonlinear Dynamics, 1-13. \cr

Whitcher, B., Guttorp, P., & Percival, D. B. (2000). Wavelet analysis of covariance with application to atmospheric time series. Journal of Geophysical Research: Atmospheres, 105(D11), 14941-14962. <URL: \doi{10.1029/2000JD900110}>.  
}
