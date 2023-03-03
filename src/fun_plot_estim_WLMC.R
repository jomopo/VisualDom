########################################################################
# "plot_estim_WLMC" function from the R package VisualDom              #
########################################################################
# The function "plot_estim_WLMC" plots as a heat map the output of the #
# function "estim_WLMC" from the R package VisualDom. One of the       # 
# features of our "estim_heatmap function is that this discern the     # 
# correlation coefficients that are not statistically significant.     # 
# The Wavelet Local Multiple Correlation (WLMC) was developed by       #
# Fernández-Macho (2018), Time-localized wavelet multiple regression   #
# and correlation. Physica A, 492, 1226-1238. The wavelet transform    #
# (MODWT) is carried out through the R package waveslim, Whitcher, B., #
# Guttorp, P. & Percival, D. B (2000), Wavelet analysis of covariance  #
# with application to atmospheric time series, J. Geophys. Res. Atmos. #
# 105, 14941-14962. In addition to these key references, you can also  #
# take a look at: Polanco-Martínez, J.M., Fernández-Macho, J. &        #
# Medina-Elizalde, M. (2020), Dynamic wavelet correlation analysis for #
# multivariate  climate time series, Scientific Reports, 10(1), 1-11.  #
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


 plot_estim_WLMC <- function(inputdata, LISTvals, J, fac=1, FLAG=TRUE, 
   		     FLAGNA=1, COLS=c(1:5), LTY=c(rep(1,5)), 
                     LWD=c(rep(1.2,5)), DIST=c(seq(0,10,2))) 

 {  # Open the main function 

 #---------------------------------------------------------------------#
 #                Description of arguments (INPUTS):                   #
 #---------------------------------------------------------------------#
 # 1. inputdata:  Matrix of N columns: time and the N-1 variables      # 
 # 							               # 
 # 2. LISTvals: this is the output of the function estim_heatmap, it   #
 #    is a list that contains: (1) CORCOEF (correlation coefficients), #
 #    (2) CIlo (CI lower bounds ), (3) CIup (CI upper bounds), and (4) #
 #    YmaxR (the index numbers, from 1 to number of variables, of the  #
 #    variable whose correlation is calculated against a linear        # 
 #    combination of the rest). 			               #
 # 							               # 
 # 3. J: maximum level of the MODWT decomposition. It is recommended   # 
 #    to use "J=trunc(log2(N)) - 3", where N is the number of rows     # 
 #    (elements) of inputdata (Fernández-Macho 2018, Polanco-Martínez  # 
 #    et al. 2020).  						       #
 # 							               # 
 # 4. fac: this factor is used to scale the wavelet time-scales or     #
 #    "periods" when the basic time scale is not the unit, by the      #
 #    default is 1.                                                    # 
 # 							               # 
 # 5. FLAG: this "flag" is used to plot the Y axis of the multivariate # 
 #    time series if the number of these series is less than four, by  # 
 #    default is TRUE.                                                 # 
 # 6. FLAGNA: this is used to plot (by the default is 1) or not (0)    #
 #    the correlation coefficients that are/not statistically          # 
 #    significant.                                                     # 
 # 							               # 
 # 7. COLS: the colors used to plot the multivariate time series.      #  
 # 							               # 
 # 8. LTY: the type of lines used to plot the multivariate time series.#  
 # 							               # 
 # 9. LWD: the tick size used to plot the multivariate time series.    #
 # 							               # 
 # 10 DIST: this parameter is used to define the distances between the # 
 #    Y axis. 							       #
 #---------------------------------------------------------------------#


 ########################################################################
 # Step 1: Check the input data (multivariate time series) 
 ########################################################################
 # "Check" 1: time steps MUST be regular/evenly - no gaps! 
 ########################################################################
 
 ########################################################################
 # Step 2: Settings the parameters to plot the time series analysed  
 ########################################################################
 # Determine the number columns & rows, first column MUST be the time 
 NR     <- dim(inputdata)[1]  
 NC     <- dim(inputdata)[2]  
 LABELS <- colnames(inputdata)[2:NC]

 ########################################################################
 # "Check" 2: check the number of time series to plot 
 # if these are > 4 then don't plot left y axis
 if ( (NC-1) > 4 ) { 
  FLAG <- FALSE  
 }  

 ########################################################################
 # Step 3: Settings parameters for WLMC heat map 
 ########################################################################
 # This is used to build the time-scales or periods 
 jscales <- 1:J  
 twopowj <- rep(NA, J+1)  
 
 for (j in 1:(J+1)) {
  twopowj[j] <- 2^j  
 }

 ########################################################################
 # "fac" is used to scale the wavelet time-scales or "periods" when the 
 # basic time scale is not the unit.
 twopowj <- fac*twopowj 

 PAR_left    <- rep("(", J) 
 PAR_left[1] <- "[" 
 PAR_righ    <- rep("]", J) 
 scale.names <- paste(PAR_left, twopowj[1:(J)], "-", twopowj[2:(J+1)], 
                 PAR_righ, sep="")

 scale.names <- c(scale.names[1:J],"Smooth")
 xlab <- "Time"
 ylab <- ""

 ########################################################################
 # Step 4: Settings parameters for WLMC heat map 
 ########################################################################
 corcof <- as.matrix(LISTvals$CORCOEF) # Correlation coefficients
 CIlo   <- as.matrix(LISTvals$CIlo)    # CI lower bounds  
 CIup   <- as.matrix(LISTvals$CIup)    # CI upper bounds
 YmaxR  <- as.matrix(LISTvals$YmaxR)   # The index numbers of the 
 #                        variable whose correlation is calculated
 #                        against a linear combination of the rest

 ########################################################################
 # Step 5: Plot correlation coefficients that are statistically significant 
 ########################################################################
 # This piece of code is used to plot blank marks to indicate that these 
 # points (correlation coefficients) are not statistically significant. 
if(FLAGNA){
 new_corcof <- corcof
 for (l in 1:dim(corcof)[2]) { 
  distan <- CIup[,l] - CIlo[,l] 
  id.0   <- which(CIup[,l] >= 0 & CIlo[,l] <= 0) 
  corcof[id.0,l] <- NA  
 }  
}

 ########################################################################
 # Step 6: Setting the plotting parameters 
 ########################################################################
 oldpar <- par(no.readonly = TRUE)
 on.exit(par(oldpar)) 

 layout(matrix(c(1,2), 2, 1, byrow = FALSE), heights=c(2.15, 3.35))
 par(oma=c(0, 0, 0, 1), mar=c(4.1, 4.2, 2.2, 3.5)+0.1, 
  mai = c(0.95, 1.5, 0.6, 0.05))

 ########################################################################
 # Step 7: Plot the multivariate time series  
 ########################################################################
 # This trick is used to pass the labels argument to "main" when N time
 # series are plotted. 
 #
 LABELS.back <- LABELS 
 NL          <- length(LABELS) 
 labels      <- list()
 for (l in 1:NL) { 
  labels[[l]] <- paste(LABELS[l], ",", sep="")
 } 
 RABELS <- paste(LABELS, collapse=", ")

 #  
 plot(inputdata[,1], scale(inputdata[,2]), t="l", xlab="Time", ylab="", 
  lwd=1.2, col=COLS[1], yaxt="n", xaxs="i", cex.lab=1.65, cex.axis=1.55, 
  cex.main=1.75, main=RABELS)
 axis(side=2, at=pretty(scale(inputdata[,2])), 
  labels=pretty(scale(inputdata[,2])), las=1, col.axis=COLS[1])
 mtext(2, line=DIST[1]+1.5, text=LABELS[1], col=COLS[1])
 for (p in 3:NC) {
  par(new=T)
  plot(inputdata[,1], scale(inputdata[,p]), t="l", lty=LTY[p-1], 
   col=COLS[p-1], lwd=LWD[p-1], xlab="", ylab="", xaxt="n", 
   yaxt="n", xaxs="i")
  if (FLAG==TRUE) { 
   axis(side=2, at=pretty(scale(inputdata[,p])), 
    labels=pretty(scale(inputdata[,p])), line=DIST[p-1], las=1, 
    col.axis=COLS[p-1])   
   mtext(2, line=DIST[p-1]+1.5, text=LABELS[p-1], col=COLS[p-1])
  } 
 }

 # This is used to avoid to plot the smooth scale! 
 corcof.back      <- corcof
 corcof           <- corcof[,1:J]
 scale.names.back <- scale.names 
 scale.names      <- scale.names[1:J]
 YmaxR.back       <- YmaxR
 YmaxR    	  <- YmaxR[,1:J]

 ########################################################################
 # Step 8: Plot the heat map  
 ########################################################################
 plot3D::image2D(z=corcof, x=inputdata[,1], y=1:(ncol(corcof)), cex=1.55, 
  cex.lab=1.55, main="", colkey=list(cex.axis=1.5, side=1), #sub=sub,
  xlab="", ylab=ylab, axes=FALSE, clab=expression(varphi), 
  rasterImage=FALSE) 
 axis(side=2, at=1:ncol(corcof),labels=scale.names, las=1, cex.axis=1.65)
 mtext(2, text="Periods", line=5.95, cex=1.75)
 axis(side=3, at=pretty(inputdata[,1]), labels=pretty(inputdata[,1]), 
  cex.axis=1.65) 

 } # Close the main function 

