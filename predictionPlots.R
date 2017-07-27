################################################################################
#
#   M. Karnauskas, Jul 21, 2017                                     
#   Part II of code - run analysis.R first
#
#   Compares model predictions from GAMMs to GLMMs
#   Creates regularly spaced prediction grid, extracts depth measurements for 
#   prediction points, calculated predictions and plots predictions on map
#
################################################################################
rm(list=ls())
load("model_parameters.RData")
load("model_data.RData")  
source("plotSAmap.r")                                                           # plotting code  

################################  libraries  ###################################
if (!"lme4" %in% installed.packages()) install.packages("lme4", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
if (!"AICcmodavg" %in% installed.packages()) install.packages("AICcmodavg", repos='http://cran.us.r-project.org')
if (!"mgcv" %in% installed.packages()) install.packages("mgcv", repos='http://cran.us.r-project.org')
if (!"splancs" %in% installed.packages()) install.packages("splancs", repos='http://cran.us.r-project.org')
if (!"PBSmapping" %in% installed.packages()) install.packages("PBSmapping", repos='http://cran.us.r-project.org')
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
if (!"matlab" %in% installed.packages()) install.packages("matlab", repos='http://cran.us.r-project.org')
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
library(lme4)
library(maps)
library(AICcmodavg)
library(mgcv)
library(splancs)
library(PBSmapping)
library(sp)
library(matlab)
library(ncdf4)          
################################################################################

#######################   CREATE PREDICTION POLYGON   ##########################
                                                              
nc <- nc_open('GEBCO_2014_2D_-84_24_-74_36.nc')                                 # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)

z[which(z>0)] <- 0
#image(x,y,z)                                                                    # plot bathymetry
abline(h=min(d$lat)); abline(h=max(d$lat))                                      # check extent with respect to data limits

c15<- contourLines(x,y,z, levels=c(-15))                                        # isobaths for 0, 15, 85m
  for (i in 1:1) {  lines(c15[[i]]$x, c15[[i]]$y)  }
c85 <- contourLines(x,y,z, levels=c(-85))
  for (i in 1:1) {  lines(c85[[i]]$x, c85[[i]]$y)  }
c0 <- contourLines(x,y,z, levels=c(0))
  for (i in 1:1) {  lines(c0[[i]]$x, c0[[i]]$y)  }
  
###############   define shallow and deep polygon boundaries    ################
sh <- cbind(c15[[1]]$x[seq(100,5000,5)], c15[[1]]$y[seq(100,5000,5)])           # use every 5 points to approximate isobath
dp <- cbind(c85[[1]]$x[seq(1,2600,3)], c85[[1]]$y[seq(1,2600,3)])
sh <- sh[which(sh[,2] < max(d$lat) & sh[,2] > min(d$lat)),]
dp <- dp[which(dp[,2] < max(d$lat) & dp[,2] > min(d$lat)),]
points(sh, pch=20, col=2)                                                       # check extraction of points
points(dp, pch=20, col=2)
polygon(rbind(sh, dp), col=51)

############  convert data points and polygon from LL to UTM (km) ##############
p <- as.data.frame(rbind(sh, dp))                                               # convert point list to polygon
p.utm <- p
names(p.utm) = c("X","Y")
attr(p.utm, "projection") = "LL"
attr(p.utm, "zone") = 17                                                        # Use UTM zone 17 based on http://www.dmap.co.uk/utmworld.htm
UTMPts_pol = convUL(p.utm, km=TRUE)                                             # convert coordinates to km so that prediction grid is regularly spaced
names(UTMPts_pol) = c("Easting","Northing")

########################   CREATE PREDICTION GRID   ############################
UTMPts_samp <- gridpts(as.matrix(UTMPts_pol), xs=4, ys=4)                       # prediction grid  DEFINE RESOLUTION HERE  (now 4km grid)
UTMPts_samp <- as.data.frame(UTMPts_samp)
names(UTMPts_samp) <- c("X", "Y")
attr(UTMPts_samp, "projection") = "UTM"
attr(UTMPts_samp, "zone") = 17                                                  # Use UTM zone 17 based on http://www.dmap.co.uk/utmworld.htm
samp <- convUL(UTMPts_samp, km=TRUE)                                            # convert back to coordinates

points(samp$X, samp$Y, col=5, pch=19, cex=0.2)                                  # check point selection 
plot(UTMPts_samp$X, UTMPts_samp$Y, col=5, pch=19, cex=0.1)                      # map points in prediction grid
map('usa', xlim=c(-82, -75), ylim=c(26, 36))
points(samp$X, samp$Y, col=5, pch=19, cex=0.1)

####################   add depths to prediction grid   #########################
for (i in 1:nrow(samp)) {  samp$dep[i] <- -z[which.min(abs(x - samp$X[i])), which.min(abs(y - samp$Y[i]))] }
head(samp)
samp$lat <- samp$Y
cols <- rainbow(104)
plot(samp$X, samp$Y, col=cols[samp$dep], pch=15, cex=0.5)                       # check assignment of depths

#################   EXTRAPOLATE MODEL PREDICTIONS TO NEW GRID   ################
gamPAfin
gamNPfin

#  add other required factors to prediction grid 
samp$doy <- mean(d$doy)                                                         # add other factors to prediction grid
samp$lunar <- mean(d$lunar)                                                     
samp$temp <- max(d$temp, na.rm=T)
samp$year <- 2012

##########################   predict on new grid   #############################

comb.var <- function(A, Ase, P, Pse) { (P^2 * Ase^2 + A^2 * Pse^2 - Ase^2 * Pse^2)  }   # combined var function

predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)       #  predlogit$fit     = predicted occurrences in probability space; predlogit$se.fit  = predicted SE in probability space
predpos <- predict.gam(gamNPfin, samp, type="response", se.fit=T)     #  predpos           = predicted abundance when present neg binomial 
#        co <- as.numeric(cor(predlogit$fit, predpos$fit, method="pearson"))       # pearson correlation between prop positive and abundance when present

samp$se <- comb.var(predpos$fit, predpos$se.fit, predlogit$fit, predlogit$se.fit)     # combined variance - (pred N when present, neg bin SE, prob of occ, logistic SE)
samp$N <-  predlogit$fit * predpos$fit                                                  # estimated abundance is prob. of occurrence * estimated abundance when present

hist(predlogit$fit)
hist(predpos$fit)
hist(predlogit$se.fit)
hist(predpos$se.fit)

################### compare GLMM vs GAMM predictions  ##########################
par(mfrow=c(1,3)) 

plotSAmap(predlogit$fit, samp$X, samp$Y, pchnum=15, cexnum=0.5)
text(-76.8, 27.6, "spawning female\nprobability of occurrence")

plotSAmap(predpos$fit , samp$X, samp$Y, pchnum=15, cexnum=0.5)
text(-76.8, 27.6, "spawning female\nabundance when present")

plotSAmap(samp$N, samp$X, samp$Y, pchnum=15, cexnum=0.5)
text(-76.8, 27.6, "spawning female\nindex")
  

###################    plot GAMM predictions and SE   ##########################
par(mfrow=c(1,2)) 

plotSAmap(predlogit$fit, samp$X, samp$Y, pchnum=15, cexnum=0.35)
text(-76.8, 27.6, "spawning female\nprobability of occurrence")

plotSAmap(predlogit$se.fit, samp$X, samp$Y, pchnum=15, cexnum=0.35)
text(-76.8, 27.6, "model S.E.")


################  plot spawning activity by day of year  #######################
par(mfrow=c(3,4), mex=0.6) 

    yloc <- seq(28, 32, length.out=100)
    cols <- rainbow(100, start=0.01, end=0.7)[100:1]
    
for (i in seq(min(d$doy), max(d$doy), length.out=12))  {
  samp$doy <- i                                                                   # loop through days of year
  predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)      

  map("state", interior = TRUE, xlim=c(-81.75, -75), ylim=c(26.8, 35.2)); axis(1); axis(2); box(); 
  mtext(side=3, paste("RS spawning activity\nday", round(mean(samp$doy*365)), "of year"), cex=0.8)
  points(samp$X, samp$Y, col=cols[round(predlogit$fit*300)+1], pch=15, cex=0.5)
  for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
  text(x=-76.1, y=yloc[seq(0,100,10)]+0.2, seq(0,0.9,0.1), pos=1)
  text(-77.2, 27.6, "spawning female\nprobability of occurrence")  } 
  
################  plot spawning activity by lunar phase  #######################
boxplot(d$lunar ~ d$lun4)

par(mfrow=c(2,2), mex=0.5) 
samp$doy <- mean(d$doy)
n <- c("New", "Waxing", "Full", "Waning")

for (i in 1:4)  {
  samp$lunar <- seq(0.5, 5, length.out=4)[i]                                    # loop through lunar phases
  predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)    

  map("state", interior = TRUE, xlim=c(-81.75, -75), ylim=c(26.8, 35.2)); axis(1); axis(2); box(); 
  mtext(side=3, paste("RS spawning activity\n", n[i], "Moon"), cex=0.8)
  points(samp$X, samp$Y, col=cols[round(predlogit$fit*300)+1], pch=15, cex=0.5)
  for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
  text(x=-76.1, y=yloc[seq(0,100,10)]+0.2, seq(0,0.9,0.1), pos=1)  
  text(-76.8, 27.6, "spawning female\nprobability of occurrence") } 

####################################  END  #####################################




                                                                 