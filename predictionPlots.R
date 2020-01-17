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
if (!"lunar" %in% installed.packages()) install.packages("lunar", repos='http://cran.us.r-project.org')
library(lme4)
library(maps)
library(AICcmodavg)
library(mgcv)
library(splancs)
library(PBSmapping)
library(sp)
library(matlab)
library(ncdf4)          
library(lunar) 
################################################################################

#######################   CREATE PREDICTION POLYGON   ##########################
                                                              
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

########################   CREATE PREDICTION GRID   ############################
UTMPts_samp <- gridpts(as.matrix(UTMPts_pol), xs=10, ys=10)                       # prediction grid  DEFINE RESOLUTION HERE  (was 4km grid - changed to 10km Sep2019 to match Gulf map resolution)
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
cols <- rainbow(104)
plot(samp$X, samp$Y, col=cols[samp$dep], pch=15, cex=0.5)                       # check assignment of depths

save(samp, file="SApredictionGrid.RData")                                            # save prediction grid

#################   EXTRAPOLATE MODEL PREDICTIONS TO NEW GRID   ################
gamPAfin
gamNfin

#  add other required factors to prediction grid 
samp$doy <- mean(d$doy)                                                         # add other factors to prediction grid
samp$lunar <- mean(d$lunar)                                                     
samp$temp <- mean(d$temp, na.rm=T)
samp$year <- 2009           
samp$lon <- samp$X
samp$lat <- samp$Y

######################################  functions   #############################################
comb.var   <- function(A, Ase, P, Pse, p) { (P^2 * Ase^2 + A^2 * Pse^2 + 2 * p * A * P * Ase * Pse)  }   # combined variance function
lnorm.mean <- function(x1, x1e) {  exp(x1 + 0.5 * x1e^2)   }
lnorm.se   <- function(x1, x1e) {  ((exp(x1e^2)-1)*exp(2 * x1 + x1e^2))^0.5  } 
#################################################################################################

predmat <- c()
predsemat <- c()
for (i in seq(min(d$lunar), max(d$lunar), length.out=4))  {
   for (j in seq(min(d$temp, na.rm=T), max(d$temp, na.rm=T), length.out=10))  {
      for (k in unique(d$year)[2:11])  {
        samp$lunar <- i                                                          # add other factors to prediction grid
        samp$temp <- j
        samp$year <- k

        predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)           # predict occurrences
        predposlog <- predict(gamNfin, samp, type="response", se.fit=T)           # predict eggs when present
        predpos   <- lnorm.mean(predposlog$fit, predposlog$se.fit)                # convert lognormal mean and SE to normal space
        predposse <- lnorm.se(predposlog$fit, predposlog$se.fit)
        co <- as.numeric(cor(predlogit$fit, predpos, method="pearson"))           # calculate covariance
        predse <- (comb.var(predpos, predposse, predlogit$fit, predlogit$se.fit, co))^0.5    # calculate combined variance
        predind <-  predlogit$fit * predpos                                       # calculate combined index
    predmat   <- cbind(predmat, predind)
    predsemat <- cbind(predsemat, predse)   } }
        cat(i, "\n")    }
    
tab <- cor(predmat)
min(cor(predmat))
min(cor(round(predmat/100)))
min(cor(round(predmat/10000)))
which.min(tab)
min(tab)
tab[3,400]
cor(predmat[,3], predmat[,400])
plot(predmat[,3], predmat[,400])

par(mfrow=c(1,2))
plotSAmap(round(predmat[,3]/800), samp$X, samp$lat, pchnum=15, cexnum=0.8)
plotSAmap(round(predmat[,400]/200000), samp$X, samp$lat, pchnum=15, cexnum=0.8)
cor(round(predmat[,3]/800), round(predmat[,400]/200000))
plot(round(predmat[,3]/800), round(predmat[,400]/200000))
      
      
mp <- rowMeans(predmat)
mv <- rowMeans(predsemat)        

plot(data.frame(predmat[,1:12]))
plot(data.frame(predmat[, seq(1, 60, 5)]))
cor(data.frame(predmat[,1:12]))
cor(data.frame(predmat[, seq(1:60)]))

#  add other required factors to prediction grid 
samp$doy <- mean(d$doy)                                                         # add other factors to prediction grid
samp$lunar <- mean(d$lunar)                                                     
samp$temp <- mean(d$temp, na.rm=T)
samp$year <- 2010

##########################   predict on new grid   #############################
predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)                 # predict occurrences 
predposlog <- predict(gamNfin, samp, type="response", se.fit=T)                 # predict eggs when present   
predpos   <- lnorm.mean(predposlog$fit, predposlog$se.fit)                      # convert lognormal mean and SE to normal space
predposse <- lnorm.se(predposlog$fit, predposlog$se.fit)
co <- as.numeric(cor(predlogit$fit, predpos, method="pearson"))                 # calculate covariance 
predse <- (comb.var(predpos, predposse, predlogit$fit, predlogit$se.fit, co))^0.5    # calculate combined variance
predind <-  predlogit$fit * predpos                                             # calculate combined index
        
samp$se <- predse
samp$N <-  predind   
par(mfrow=c(2,2))                                             
hist(predind)
hist(predse)

plot(samp$N, mp)
cor(samp$N, mp)
plot(samp$se, mv)
cor(samp$se, mv)

################### compare GLMM vs GAMM predictions  ##########################
par(mfrow=c(1,2))

plotSAmap(samp$N/100000, samp$X, samp$Y, pchnum=15, cexnum=0.8)
text(-76.8, 27.6, "spawning female\nprobability of occurrence")

plotSAmap(mp/100000, samp$X, samp$Y, pchnum=15, cexnum=0.8)
text(-76.8, 27.6, "spawning female\nprobability of occurrence")
  
###################    plot GAMM predictions and SE   ##########################
par(mfrow=c(1,2))
 
plotSAmap(mp/100000, samp$X, samp$Y, pchnum=15, cexnum=0.8)
text(-76.8, 27.6, "total red snapper\negg production")

plotSAmap(mv/100000, samp$X, samp$Y, pchnum=15, cexnum=0.8)
text(-76.8, 27.6, "model S.E.")

################  plot spawning activity by day of year  #######################
windows()
par(mfrow=c(3,4), mex=0.6) 

    yloc <- seq(28, 32, length.out=100)
    cols <- rainbow(100, start=0.01, end=0.7)[100:1]
    
for (i in seq(min(d$doy), max(d$doy), length.out=12))  {
  samp$doy <- i                                                                 # loop through days of year
predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)                 # predict occurrences 
predposlog <- predict(gamNfin, samp, type="response", se.fit=T)                 # predict eggs when present   
predpos   <- lnorm.mean(predposlog$fit, predposlog$se.fit)                      # convert lognormal mean and SE to normal space
predposse <- lnorm.se(predposlog$fit, predposlog$se.fit)
co <- as.numeric(cor(predlogit$fit, predpos, method="pearson"))                 # calculate covariance 
predse <- (comb.var(predpos, predposse, predlogit$fit, predlogit$se.fit, co))^0.5    # calculate combined variance
predind <-  predlogit$fit * predpos                                             # calculate combined index

  map("state", interior = TRUE, xlim=c(-81.75, -75), ylim=c(26.8, 35.2)); axis(1); axis(2); box(); 
  mtext(side=3, paste("RS spawning activity\nday", round(mean(samp$doy*365)), "of year"), cex=0.8)
  points(samp$X, samp$Y, col=cols[round(predlogit$fit*300)+1], pch=15, cex=0.5)
  for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
  text(x=-76.1, y=yloc[seq(0,100,10)]+0.2, seq(0,0.9,0.1), pos=1)
  text(-77.2, 27.6, "spawning female\nprobability of occurrence")  } 
  
################  plot spawning activity by lunar phase  #######################
windows()
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


################  plot spawning activity by year  #######################

windows()
par(mfrow=c(3,4), mex=0.5) 
samp$doy <- mean(d$doy)
yr <- 2004:2014

for (i in 1:11)  {
  samp$lunar <- 2.75                                    
  samp$year <- yr[i]           # loop through year
  predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)    
  
  map("state", interior = TRUE, xlim=c(-81.75, -75), ylim=c(26.8, 35.2)); axis(1); axis(2); box(); 
  mtext(side=3, yr[i], cex=0.8)
  points(samp$lon, samp$lat, col=cols[round(predlogit$fit*250)+1], pch=15, cex=0.75)
  for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
  text(x=-76.1, y=yloc[seq(0,100,10)]+0.2, seq(0,0.9,0.1), pos=1)  
  text(-76.8, 27.6, "spawning female\nprobability of occurrence") } 

###############################  END PLOTS  ####################################

########################   groundtruthing with Sue dataset  ####################

d1 <- d[which(d$year == 2012),]       # use year specific to Sue data b/c year effects strong
# rerun GAM for 2012 only
gam9 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d1, method="REML")
plot(gam9)

su <- read.table("mandy_rs.csv", sep=",", header=T)        # Sue data
                                                           # add necessary factors
su$doy <- NA                                                                    #  for GAM, can use continuous day of year instead of month
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(su))  {  su$doy[i] <- (sum(dinmon[1:su$month[i]]) + su$day[i]) / 365  }
su$lunar <- lunar.phase(as.Date(paste(su$year,"-",su$mon,"-",su$day,sep="")), name=F)  # also can use continuous lunar phase                                   
out <- gam(temp ~ s(doy), data=d1)
plot(out)
su$temp <- predict(out, su, se.fit=F)
points(su$doy, scale(su$temp))
su$dep <- su$depth
su$lat <- su$Latitude

predsu <- predict(gam9, su, type="response", se.fit=T) 
su$pred <- predsu$fit        

su$Reprophase <- as.factor(su$Reprophase)
su$mature <- as.factor(su$mature)
su$activespawn <- as.factor(su$activespawn)

labs <- c("immature", "developing", "spawning\ncapable", "actively\nspawning", "regressing", "regenerating")
boxplot(su$pred ~ su$Reprophase, names=labs, las=1, tck=0)
summary(lm(su$pred ~ su$Reprophase))

plot(as.numeric(su$Reprophase), su$pred, axes=F); axis(1, at=1:6, lab=labs); axis(2); box()
abline(h=mean(range(su$pred)), col=2, lwd=2)

su$phase <- "nospawn"
su$phase[which(su$Reprophase==3 | su$Reprophase==4)] <- "spawn"
                               
boxplot(su$pred ~ su$phase)
t.test(su$pred ~ su$phase)
wilcox.test(su$pred ~ su$phase)

boxplot(su$pred ~ su$mature)
t.test(su$pred ~ su$mature)
wilcox.test(su$pred ~ su$mature)

par(mgp=c(3,1,0))
boxplot(su$pred ~ su$activespawn, axes=F, ylab="predicted index of spawning activity", notch=T)
axis(2, las=2)
par(mgp=c(0,2,0))
axis(1, at=1:2, lab=c("actively \nspawning (1)", "not actively \nspawning (2)")); box()
t.test(su$pred ~ su$activespawn)
wilcox.test(su$pred ~ su$activespawn)
text(1.5, 0.43, "difference \nbetween means: \nP < 0.001")

# Repro phases are: 
# 1=immature 
# 2=developing (just beginning to develop for the season, not yet capable of spawning) 
# 3=spawning capable
# 4=actively spawning
# 5=regressing (shutting down for the season) 
# 6=regenerating

# If a fish is mature it is given a 2, if it is actively spawning (within 2 h) it is given a 1

lat <- tapply(su$Latitude, su$Reference, mean)
lon <- tapply(su$Longitude, su$Reference, mean)

f <- data.frame(cbind(table(su$Reference, su$phase), lon, lat))

f$prspawn <- f$spawn/(f$nospawn + f$spawn)
f$dep <- tapply(su$dep, su$Reference, mean)
f$doy <- tapply(su$doy, su$Reference, mean)
f$lunar <- tapply(su$lunar, su$Reference, mean)
f$N <- f$nospawn + f$spawn

points(f$lon, f$lat, cex=f$prspawn)

su$bin <- as.numeric(as.factor(su$phase)) - 1

gamsu <- gam(cbind(spawn, N) ~ s(dep) + s(lat) + s(doy) + s(lunar), family=binomial, data=f, method="REML")
plot(gamsu)


####################################  END  #####################################

#  for CMS file:  
#  samp - predict on dep, lat
#  loop through doy; use average temp for that month
#  use average lunar phase (little effect)
#  calculate for all years and average

f <- gam (temp ~ s(doy), data=d)
plot(f)
doy <- seq(min(d$doy), max(d$doy), length.out=32)
diff(doy*365)
temps <- predict(f, data.frame(doy))
plot(doy, temps)

sampold <- samp
#  add other required factors to prediction grid 
samp <- data.frame(cbind(sampold$lon, sampold$lat, samp$dep))
names(samp) <- c("lon", "lat", "dep")
samp$lunar <- mean(d$lunar)                                                     

windows()
par(mfrow=c(6,6), mex=0.5)
for (j in 1:length(doy))  {
        predmat <- c()
        predsemat <- c()                                                                                         
        samp$temp <- temps[j]
        samp$doy <- doy[j]
   for (k in unique(d$year))  {
        samp$year <- k    
        pred <- predict(gamPAfin, samp, type="response", se.fit=F)       #  predlogit$fit     = predicted occurrences in probability space; predlogit$se.fit  = predicted SE in probability space
        predmat   <- cbind(predmat, pred)
#        predsemat <- cbind(predsemat, pred$se.fit)   
        } 
  mp <- rowMeans(predmat)
#  mv <- rowMeans(predsemat) 
  plotSAmap(mp*10, samp$lon, samp$lat, pchnum=15, cexnum=0.35)
  mtext(side=3, paste("RS spawning activity\nday", round(mean(doy[j]*365)), "of year"), cex=0.8)
    } 

