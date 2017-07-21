
rm(list=ls())
load("model_parameters.RData") 

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

#######################   CREATE PREDICTION POLYGON   ###########################

nc <- nc_open('GEBCO_2014_2D_-84_24_-74_36.nc')                                 # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)

z[which(z>0)] <- 0
image(x,y,z)                                                                    # plot bathymetry
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

#################   extrapolate model predictions to new grid   ################
gamf
glmf

####################   add depths to prediction grid   #########################
for (i in 1:nrow(samp)) {  samp$dep[i] <- -z[which.min(abs(x - samp$X[i])), which.min(abs(y - samp$Y[i]))] }
head(samp)
samp$lat <- samp$Y
cols <- rainbow(104)
plot(samp$X, samp$Y, col=cols[samp$dep], pch=15, cex=0.5)                       # check assignment of depths

samp$doy <- mean(d$doy)                                                         # add other factors to prediction grid
samp$lunar <- mean(d$lunar)                                                     
samp$temp <- mean(d$temp)
samp$year <- 2014

samp$dep2 <- samp$dep
samp$dep2[which(samp$dep2 > 85)] <- 84                                          # replace deeper values so that prediction bins match
samp$depbins <- cut(samp$dep2, breaks=c(10, 25, 30, 35, 40, 50, 60, 85))
samp$latbins <- cut(samp$lat, breaks=seq(27.1, 35.1, 1))
samp$mon <- as.factor(8)
samp$lun4 <- "Full"
samp$tempbins <- "(22,24]"

predglm <- predictSE  (glmf, samp, type="response", se.fit=TRUE) 
predgam <- predict.gam(gamf, samp, type="response", se.fit=TRUE) 

samp$Nglm <- predglm$fit
samp$Nglmse <- predglm$se.fit
samp$Ngam <- predgam$fit
samp$Ngamse <- predgam$se.fit

cols=rainbow(100, start=0.01, end=0.8)[100:1];  plot(1:100, col=cols, pch=20)   # choose colormap

par(mfrow=c(1,2))                                                               # compare results between GAM and GLM
map('usa', xlim=c(-82, -75), ylim=c(26.5, 36), main="RS spawning activity"); axis(1); axis(2); box()
points(samp$X, samp$Y, col=cols[round(samp$Ngam*100)], pch=15, cex=0.5)

yloc <- seq(28, 32, length.out=100)
for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
text(x=-76.2, y=yloc[seq(0,100,10)], seq(0,0.9,0.1), pos=1)

map('usa', xlim=c(-82, -75), ylim=c(26.5, 36), main="RS spawning activity"); axis(1); axis(2); box()
points(samp$X, samp$Y, col=cols[round(samp$Nglm*100)], pch=15, cex=0.5)

yloc <- seq(28, 32, length.out=100)
for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
text(x=-76.2, y=yloc[seq(0,100,10)], seq(0,0.9,0.1), pos=1)

                                                                 