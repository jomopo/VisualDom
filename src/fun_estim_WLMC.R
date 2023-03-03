########################################################################
# "estim_WLMC" function from the R package VisualDom                #
########################################################################
# The function "estim_WLMC" estimates the Wavelet Local Multiple    #
# Correlation (WLMC) that is developed by Fernández-Macho (2018), Time-#
# localized wavelet multiple regression and correlation. Physica A,    # 
# 492, 1226-1238. It is based on the "wave.local.multiple.correlation" #
# function from the wavemulcor R package (Fernández-Macho (2018), but  # 
# some improvements have been added. The wavelet transform (MODWT) is  #
# carried out through the R package waveslim, Whitcher, B., Guttorp,   #
# P. & Percival, D. B (2000), Wavelet analysis of covariance with      #
# application to atmospheric time series, J. Geophys. Res. Atmos. 105, #
# 14941-14962. In addition to these key references, you can also look: #
# Polanco-Martínez, J.M., Fernández-Macho, J. & Medina-Elizalde,       #
# M. (2020), Dynamic wavelet correlation analysis for multivariate     #
# climate time series, Scientific Reports, 10(1), 1-11.  	       #
########################################################################

########################################################################
#   Copyright (C) 2022 by Josué M. Polanco-Martínez                    #
#   This file/code is part of the R package VisualDom                  #
########################################################################
#								     
#   VisualDom is free software:  
#   you can redistribute it and/or modify it under the terms of the GNU 
#   General Public License as published by the Free Software 
#   Foundation, either version 3 of the License, or (at your option) 
#   any later version.
#
#   VisualDom is distributed in the hope that it will be 
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with VisualDom If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

 estim_WLMC <- function(inputdata, wf="la8", J, window, M, Ymaxr=NULL) 

 { # Open the main function 

 #---------------------------------------------------------------------#
 #                Description of arguments (INPUTS):                   #
 #---------------------------------------------------------------------#
 # 1. inputdata:  Matrix of N columns: time and the N-1 variables      # 
 # 							               # 
 # 2. wf: name of the wavelet filter used in the decomposition. By     #
 #    default we use "la8," i.e., the Daubechies orthonormal compactly #
 #    supported wavelet of length L=8, but other wavelets can be used. # 
 #    Please look at the manual of modwt in R package waveslim.        # 
 #								       #
 # 3. J: maximum level of the MODWT decomposition. It is recommended   # 
 #    to use "trunc(log2(N)) - 3", where N is the number of rows       # 
 #    (elements) of inputdata (Fernández-Macho 2018, Polanco-Martínez  # 
 #    et al. 2020).  						       #
 #								       # 
 # 4. window: Weight (window) function, by the default is the Gaussian # 
 #    window, but other five windows are available: uniform,  	       #
 #    Bartlett’s triangular, Cleveland’s tricube, Wendland’s truncated #
 #    power & Epanechnikov’s parabolic), please look at the function   #
 #    wave.local.multiple.correlation from the R package wavemulcor    # 
 #    (Fernandez-Macho 2018).  	                                       #
 #							               #
 # 5. M: the length of the weight/window function, it is recommended   # 
 #    to use "trunc(N/2^3)", and N has been previously defined         #
 #    (Fernández-Macho 2018, Polanco-Martínez et al. 2020).   	       # 
 #							               #
 # 6. Ymaxr: it is used to maximize the multiple correlation for each  #
 #    wavelet scale, by default is NULL, that is, we do not define     #
 #    a priori specific variable but instead let the WLMC select one.  # 
 #---------------------------------------------------------------------#


 ########################################################################
 # Step 1: Check the input data  
 ########################################################################
 # "Check" 1: time steps MUST be regular/evenly - no gaps! 
 warning("\n W A R N I N G: The input data must be regular (evenly spaced 
  in time), please verify this condition. Otherwise, please, consider 
  to address this drawback before using VisualDom. We recommend our 
  BINCOR package and method (also in CRAN: 
  https://cran.r-project.org/package=BINCOR), but other packages and 
  methods can be used. \n")
 
 # Determine the number columns & rows, first column MUST be the time 
 NR  <- dim(inputdata)[1]  
 NC  <- dim(inputdata)[2]  

 ########################################################################
 # Step 2: Estimating the MODWT (wavelet transform) via waveslim package 
 ########################################################################
 datin.modwt.LST <- list()

 for (n in 2:(NC)) {  
  datin.modwt.LST[[n-1]] <- waveslim::modwt(inputdata[,n], wf, J)  
  datin.modwt.LST[[n-1]] <- waveslim::brick.wall(datin.modwt.LST[[n-1]], wf)
 } 

 ########################################################################
 # Step 3: Estimating the WLMC via wavemulcor package  
 ########################################################################
 wlmc.LST <- wavemulcor::wave.local.multiple.correlation(datin.modwt.LST, M=M, 
              window=window, ymaxr=Ymaxr)

 corcof <- as.matrix(wlmc.LST$val)   # Correlation coefficients
 CIlo   <- as.matrix(wlmc.LST$lo)    # CI lower bounds  
 CIup   <- as.matrix(wlmc.LST$up)    # CI upper bounds
 YmaxR  <- as.matrix(wlmc.LST$YmaxR) # The index numbers of the 
 #                        variable whose correlation is calculated
 #                        against a linear combination of the rest

 # Output: 
 namesLST        <- c("CORCOEF", "CIlo", "CIup", "YmaxR") 
 LISTvals        <- list(corcof, CIlo, CIup, YmaxR) 
 names(LISTvals) <- namesLST
 return(LISTvals)  
 
 } # Close the main function  

 
