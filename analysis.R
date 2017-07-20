
rm(list=ls())

setwd("C:/Users/mkarnauskas/Desktop/RSmap_SA")           
source("Xvalidate.r")

##########  libraries  ##############

if (!"chron" %in% installed.packages()) install.packages("chron", repos='http://cran.us.r-project.org')
if (!"lme4" %in% installed.packages()) install.packages("lme4", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
if (!"MASS" %in% installed.packages()) install.packages("MASS", repos='http://cran.us.r-project.org')
if (!"AICcmodavg" %in% installed.packages()) install.packages("AICcmodavg", repos='http://cran.us.r-project.org')
if (!"mgcv" %in% installed.packages()) install.packages("mgcv", repos='http://cran.us.r-project.org')
if (!"gamm4" %in% installed.packages()) install.packages("gamm4", repos='http://cran.us.r-project.org')
library(mgcv)
library(gamm4)
library(AICcmodavg)
library(chron)
library(lme4)
library(maps)
library(MASS)
library(ncdf4)
library(matlab)
library(grDevices)
library(sp)
library(PBSmapping)
library(splancs)

################################################################################
dat <- read.table("RedSnapperMatData.csv", sep=",", header=T, na.strings = c("NA"))             # read in data

length(unique(dat$PCG))
length(unique(dat$Latitude))

table(dat$Species)
table(dat$Gear)
barplot(table(dat$Day))
barplot(table(dat$Month))
map('usa')
points(dat$Longitude, dat$Latitude)
hist(dat$StationDepth)
hist(dat$Temp)
hist(dat$TL)
hist(dat$FL)
hist(dat$SL)
hist(dat$WholeWt)
table(dat$Sex)
table(dat$Mat)
table(dat$Inc)
table(dat$Mat, dat$Year)
 
matcodes <- c(3, 7, "B", "C", "D", "G", "H") 
dat$fem <- "M"
dat$fem[which(dat$Sex==2)] <- "NF"
dat$fem[which(dat$Mat %in% matcodes & dat$Sex==2)] <- "SF"
table(dat$fem, useNA="always")

d <- dat
names(d)
names(d)[4:11] <- c("day", "mon", "year", "Date", "lat", "lon", "dep", "temp")
 
d$fem <- as.factor(d$fem)
plot(d$lon, d$lat, col=as.numeric(d$fem)-1, pch=20)

d$lunim <- lunar.illumination(as.Date(paste(d$year, "-", d$mon, "-", d$day, sep="")))
d$lun4 <- lunar.phase(as.Date(paste(d$year, "-", d$mon, "-", d$day, sep="")), name=T)

dim(d); table(d$fem)
d <- d[which(d$fem!="M"), ]                                                         # take all samples which have at least 1 female
dim(d); table(d$fem)
table(d$fem)
table(as.numeric(d$fem)-2)
d$fem <- (as.numeric(d$fem)-2) 

#################  variable to define position across shelf  ###################
a.x <- max(d$lon)+0.1
a.y <- min(d$lat)-0.1
d$ang <- atan((d$lat-a.y)/(d$lon-a.x))*180/pi
cols <- rainbow(100, start=0.1)
plot(d$lon, d$lat, col=cols[round(d$ang+90)])
################################################################################
 
d$mon <- as.factor(d$mon) 
d$year <- as.factor(d$year) 

# bin variables as finely as possible while maintaining adequate number of samples per bin

d$angbins <- cut(d$ang, breaks=c(-90, -79, -60, -45, -30, -17, 0))       # -79, 
map('usa', xlim=c(-82, -75), ylim=c(26, 36))
points(d$lon, d$lat, col=d$angbins)
 
d$depbins <- cut(d$dep, breaks=c(10, 25, 30, 35, 40, 50, 60, 85))
d$tempbins <- cut(d$temp, breaks=c(10, 20, 22, 24, 30))

table(d$angbins, useNA="always")
table(d$depbins, useNA="always")
table(d$tempbins, useNA="always")
table(d$mon, useNA="always")

windows()
barplot(table(d$angbins, useNA="always"))
barplot(table(d$depbins, useNA="always"))
barplot(table(d$tempbins, useNA="always"))
barplot(table(d$mon, useNA="always"))

barplot(tapply(d$fem, d$angbins, mean, na.rm=T))
barplot(tapply(d$fem, d$depbins, mean, na.rm=T))
barplot(tapply(d$fem, d$tempbins, mean, na.rm=T))
barplot(tapply(d$fem, d$mon, mean, na.rm=T))
barplot(tapply(d$fem, d$lun4, mean, na.rm=T))

table(d$angbins, d$mon)
table(d$angbins, d$depbins)
table(d$angbins, d$tempbins)
table(d$depbins, d$tempbins)
                  
round(tapply(d$fem, list(d$angbins, d$mon), mean, na.rm=T), 3)  # temperature * latitude interaction would be nice, but not enough data
round(tapply(d$fem, list(d$angbins, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction possible!
barplot(tapply(d$fem, list(d$angbins, d$mon), mean, na.rm=T), beside=T, col=1:6)
barplot(tapply(d$fem, list(d$angbins, d$depbins), mean, na.rm=T), beside=T, col=1:6)
matplot(tapply(d$fem, list(d$angbins, d$mon), mean, na.rm=T), type="l")
matplot(tapply(d$fem, list(d$angbins, d$depbins), mean, na.rm=T), type="l")

barplot(tapply(d$fem, d$year, mean, na.rm=T))
               
# Month and temp are highly correlated - can't use both
  
yrs0 <- as.numeric(names(which(tapply(d$fem, d$year, mean, na.rm=T)>0)))   # need to take out 1999
dim(d)
d <- d[which(d$year %in% yrs0),]; dim(d)   
  
# fixed effects model
outfixed <- glm(fem ~ depbins + angbins + mon + tempbins + lun4 + year,  family="binomial", data=d)
summary(outfixed) 
stepAIC(outfixed) 

# mixed effects model - year as random effect
out1 <- glmer(fem ~ depbins + angbins + mon + lun4 + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out1)                  
extractAIC(out1) 

out2 <- glmer(fem ~ depbins + angbins + mon + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out2)                  
extractAIC(out2)

out3 <- glmer(fem ~ depbins + angbins + mon + lun4 + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out3)                  
extractAIC(out3)

out4 <- glmer(fem ~ depbins + angbins + mon + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out4)                  
extractAIC(out4) 

extractAIC(out1)
extractAIC(out2)
extractAIC(out3)
extractAIC(out4) 

################################################################################

outnull <- glm(fem ~ 1, family=binomial(logit), data=d)
(deviance(outnull)-deviance(out3))/deviance(outnull)      # deviance explained by all factors combined

outnorand <- glmer(fem ~ 1 + (1|year), family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
(deviance(outnorand)-deviance(out4))/deviance(outnorand)  # deviance explained by fixed factors combined

############################   CROSS-VALIDATION  ###############################

# x <- xvalid(out4, d, kfold=5)
# colMeans(x)
#                                                                               Area.Under.Curve    False.Positive.Rate False.Negative.Rate 

#############################   GAM MODEL   ####################################

d$doy <- NA
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(d))  {  d$doy[i] <- (sum(dinmon[1:d$mon[i]]) + d$day[i]) / 365  }
d$lunar <- lunar.phase(as.Date(paste(d$year, "-", d$mon, "-", d$day, sep="")), name=F)

#  FACTORS:  year   mon     depbins    tempbins    lunar     angbins           
#                   doy     dep        temp        lunim     ang
#                                                            lat
                                                  
gam1 <- gam(fem ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam1)
par(mfrow=c(5,6), mex=0.5)
plot(gam1)

gam2 <- gam(fem ~ s(dep) + s(ang) + s(doy) + s(lunim) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam2)
plot(gam2)

gam3 <- gam(fem ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam3)
plot(gam3); plot.new()

gam4 <- gam(fem ~ s(dep) + s(ang) + s(doy) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam4)
plot(gam4); plot.new()

gam5 <- gam(fem ~ s(dep) + s(ang) + s(doy) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam5)
plot(gam5)

gam6 <- gam(fem ~ te(dep, ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam6)
windows()
plot(gam6)

gam7 <- gam(fem ~ s(dep) + te(doy, ang) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam7)
plot(gam7)

gam8 <- gam(fem ~ s(dep) + te(ang, temp) + s(doy) + s(lunar) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam8)
plot(gam8)

gam9 <- gam(fem ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam9)
plot(gam9)

gam10 <- gam(fem ~ te(dep, lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML" )
summary(gam10)
plot(gam10)

extractAIC(gam1)
extractAIC(gam2)
extractAIC(gam3)
extractAIC(gam4)
extractAIC(gam5)
extractAIC(gam6)
extractAIC(gam7)
extractAIC(gam8)
extractAIC(gam9)
extractAIC(gam10)

min(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10))[2,])
plot(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10))[2,], ylab="")

#  compare to GAMM package 
gamm1 <- gamm(fem ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(temp), family=binomial, random=list(year=~1), data=d)
summary(gam1)
summary(gamm1$gam)
par(mfrow=c(2,6), mex=0.5)
plot(gam1)
plot(gamm1$gam)

############################   CROSS-VALIDATION  ###############################

d <- d[which(d$year!=2006),]

dlast5 <- d[which(d$year==2010 | d$year==2011 | d$year==2012 | d$year==2013 | d$year==2014),]

x <- xvalid(gam1, d, kfold=5)
x; colMeans(x)     

x <- xvalid(gam9, d, kfold=5)
x; colMeans(x)              

x <- xvalid(gam9, dlast5, kfold=5)
x; colMeans(x)    

x <- xvalid(out1, d, kfold=5)
colMeans(x, na.rm=T)

#############################  END MODELING  ###################################

############################   PREDICTION PLOTS  ###############################

d$latbins <- cut(d$lat, breaks=seq(27.1, 35.1, 1))
out1 <- glmer(fem ~ depbins + latbins + mon + lun4 + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out1)

barplot(summary(out1)$coefficients[grep('dep', rownames(summary(out1)$coefficients)),1], las=2)
barplot(summary(out1)$coefficients[grep('lat', rownames(summary(out1)$coefficients)),1], las=2)
barplot(summary(out1)$coefficients[grep('mon', rownames(summary(out1)$coefficients)),1], las=2)
barplot(summary(out1)$coefficients[grep('lun', rownames(summary(out1)$coefficients)),1], las=2)
barplot(summary(out1)$coefficients[grep('temp', rownames(summary(out1)$coefficients)),1], las=2)
#barplot(summary(outfrand)$coefficients[grep('C', rownames(summary(outfrand)$coefficients)),1], las=2)

samp <- expand.grid(d$year[1], unique(d$latbins), unique(d$depbins), unique(d$mon), unique(d$lun4), unique(d$tempbins)) 
names(samp) <- c("year", "latbins", "depbins", "mon", "lun4", "tempbins")

dim(samp)
out1
pred <- predictSE(out1, samp, type="response", se.fit=T)       # takes a long time

head(samp)
samp$pred <- pred$fit
samp$glmse <- pred$se.fit

windows()
par(mfrow=c(2, 6), mex=0.56)

plot(gam9)

barplot(tapply(samp$pred, samp$depbins, mean, na.rm=T), xlab="depth (bins)")
barplot(tapply(samp$pred, samp$latbins, mean, na.rm=T), xlab="lat (bins)")
barplot(tapply(samp$pred, samp$mon, mean, na.rm=T), xlab="month")
barplot(tapply(samp$pred, samp$lun4, mean, na.rm=T), xlab="lunar phase")
barplot(tapply(samp$pred, samp$temp, mean, na.rm=T), xlab="temp")

#######################   CREATE PREDICTION POLYGON   ###########################

nc <- nc_open('C://Users/mkarnauskas/Desktop/red_snapper/Mar2016/GEBCO_2014_2D_-84_24_-74_36.nc')        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)

z[which(z>0)] <- 0
image(x,y,z)
abline(h=min(d$lat))
abline(h=max(d$lat))

c15<- contourLines(x,y,z, levels=c(-15))                  # isobaths for 0, 15, 85m
  for (i in 1:1) {  lines(c15[[i]]$x, c15[[i]]$y)  }
c85 <- contourLines(x,y,z, levels=c(-85))
  for (i in 1:1) {  lines(c85[[i]]$x, c85[[i]]$y)  }
c0 <- contourLines(x,y,z, levels=c(0))
  for (i in 1:1) {  lines(c0[[i]]$x, c0[[i]]$y)  }
  
##### define shallow and deep polygon boundaries
sh <- cbind(c15[[1]]$x[seq(100,5000,5)], c15[[1]]$y[seq(100,5000,5)])
dp <- cbind(c85[[1]]$x[seq(1,2600,3)], c85[[1]]$y[seq(1,2600,3)])

sh <- sh[which(sh[,2] < max(d$lat) & sh[,2] > min(d$lat)),]
dp <- dp[which(dp[,2] < max(d$lat) & dp[,2] > min(d$lat)),]

points(sh, pch=20, col=2)
points(dp, pch=20, col=2)

polygon(rbind(sh, dp), col=51)

############  convert data points and polygon from LL to UTM (km) ##############

p <- as.data.frame(rbind(sh, dp))
p.utm <- p
names(p.utm) = c("X","Y")
attr(p.utm, "projection") = "LL"
attr(p.utm, "zone") = 17
UTMPts_pol = convUL(p.utm, km=TRUE)
names(UTMPts_pol) = c("Easting","Northing")#
##########  create regular sample grid for kriging
UTMPts_samp <- gridpts(as.matrix(UTMPts_pol), xs=4, ys=4)       # kriging prediction grid  DEFINE RESOLUTION HERE
UTMPts_samp <- as.data.frame(UTMPts_samp)
names(UTMPts_samp) <- c("X", "Y")
attr(UTMPts_samp, "projection") = "UTM"
attr(UTMPts_samp, "zone") = 17      # Use UTM zone 16 based on http://www.dmap.co.uk/utmworld.htm
samp <- convUL(UTMPts_samp, km=TRUE)

points(samp$X, samp$Y, col=5, pch=19, cex=0.2)
plot(UTMPts_samp$X, UTMPts_samp$Y, col=5, pch=19, cex=0.1)
map('usa', xlim=c(-82, -75), ylim=c(26, 36))
points(samp$X, samp$Y, col=5, pch=19, cex=0.1)

#################   extrapolate model predictions to new grid   ################
gam9
out1

#######  add depth
for (i in 1:nrow(samp)) {  samp$dep[i] <- -z[which.min(abs(x - samp$X[i])), which.min(abs(y - samp$Y[i]))] }
head(samp)
samp$lat <- samp$Y

cols <- rainbow(104)
plot(samp$X, samp$Y, col=cols[samp$dep], pch=15, cex=0.5)

samp$doy <- mean(d$doy)
samp$lunar <- mean(d$lunar)
samp$temp <- mean(d$temp)
samp$year <- 2014

samp$dep2 <- samp$dep
samp$dep2[which(samp$dep2 > 85)] <- 84
samp$depbins <- cut(samp$dep2, breaks=c(10, 25, 30, 35, 40, 50, 60, 85))
samp$latbins <- cut(samp$lat, breaks=seq(27.1, 35.1, 1))
samp$mon <- as.factor(5)
samp$lun4 <- "New"
samp$tempbins <- "(24,30]"

predglm <- predict(out1, samp, type="response", se.fit=TRUE) 
pred <- predict.gam(gam9, samp, type="response", se.fit=TRUE) 

samp$N <- pred$fit
samp$Nse <- pred$se.fit
samp$Nglm <- predglm

cols=rainbow(100, start=0.01, end=0.8)[100:1];  plot(1:100, col=cols, pch=20)

par(mfrow=c(1,2))
map('usa', xlim=c(-82, -75), ylim=c(26.5, 36), main="RS spawning activity"); axis(1); axis(2); box()
points(samp$X, samp$Y, col=cols[round(samp$N*100)], pch=15, cex=0.5)

yloc <- seq(28, 32, length.out=100)
for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
text(x=-76.2, y=yloc[seq(0,100,10)], seq(0,0.9,0.1), pos=1)

map('usa', xlim=c(-82, -75), ylim=c(26.5, 36), main="RS spawning activity"); axis(1); axis(2); box()
points(samp$X, samp$Y, col=cols[round(samp$Nglm*100)], pch=15, cex=0.5)

yloc <- seq(28, 32, length.out=100)
for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
text(x=-76.2, y=yloc[seq(0,100,10)], seq(0,0.9,0.1), pos=1)

