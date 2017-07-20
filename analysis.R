
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
(deviance(outnorand)-deviance(out3))/deviance(outnorand)  # deviance explained by fixed factors combined

############################   CROSS-VALIDATION  ###############################

x <- xvalid(out3, d, rand=T, kfold=5)
colMeans(x)
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

gamm1 <- gamm(fem ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(temp), family=binomial, random=list(year=~1), data=d)
summary(gam1)
summary(gamm1$gam)
par(mfrow=c(2,6), mex=0.5)
plot(gam1)
plot(gamm1$gam)

x <- xvalid(gam1, d, rand=T, kfold=5)
colMeans(x)
#               


#############################  END MODELING  ###################################

############################   PREDICTION PLOTS  ###############################

windows()
par(mar=c(10,1,1,1), mfrow=c(3,2))
barplot(summary(outfrand)$coefficients[grep('ang', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('dep', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('mon', rownames(summary(outfrand)$coefficients)),1], las=2)
#barplot(summary(outfrand)$coefficients[grep('temp', rownames(summary(outfrand)$coefficients)),1], las=2)
#barplot(summary(outfrand)$coefficients[grep('C', rownames(summary(outfrand)$coefficients)),1], las=2)

levels(d$mon)[3]
samp <- expand.grid(d$year[1], unique(d$angbins), unique(d$depbins), d$mon[3]) 
names(samp) <- c("year", "angbins", "depbins", "mon")

dim(samp)
outfrand
pred <- predictSE(outfrand, samp, type="response", se.fit=T) 

head(samp)
samp$pred <- pred$fit
samp$glmse <- pred$se.fit

par(mfrow=c(3,2))
barplot(tapply(samp$pred, samp$angbins, mean, na.rm=T), xlab="lat (bins)")
barplot(tapply(samp$pred, samp$depbins, mean, na.rm=T), xlab="depth (bins)")
barplot(tapply(samp$pred, samp$year, mean, na.rm=T), xlab="year")
barplot(tapply(samp$pred, samp$mon, mean, na.rm=T), xlab="bathymetry")

samp2 <- samp       
tab <- tapply(samp2$pred, list(samp2$depbins, samp2$angbins), mean, na.rm=T); tab
tabv <- tapply(samp2$glmse, list(samp2$depbins, samp2$angbins), mean, na.rm=T); tabv
r <- nrow(tab); w <- ncol(tab)
cols <- hcl(h= seq(20,300,length=r), c=100); cols2 <- hcl(seq(20,300,length=r), alpha=0.3)
matplot(t(tab), type="l", lty=1, axes=F, ylim=c(0, 0.8), col=cols, lwd=2, xlab=" ", ylab="probability of spawning condition female")
axis(1, las=2, at=1:ncol(tab), lab=round(tapply(d$lat, d$angbins, mean), 2)) 
axis(2, las=2); box()
up <- t(tab) + t(tabv)
lo <- t(tab) - t(tabv)
for (i in 1:r)  {  polygon(c(1:w,w:1), c(up[1:w,i], lo[w:1,i]), col=cols2[i], border=NA) } 
legend("topright", paste(substr(rownames(tab),2,3), "-", substr(rownames(tab),5,6), "m", sep=""), col=cols, lty=1, bty="n", lwd=2)



########### stopped here  




###############  redo below #############
db <- read.dbf("SPAGgrid270_320N_fwcmmPredict.dbf")                             # read in prediction grid
head(db)

db$latbins <- cut(db$Lat, breaks=c(27, 31.2, 32.2, 35))               # create factor breaks on prediction grid to match model output
db$depbins <- cut(-db$MAXCDEPTH, breaks=c(10, 30, 40, 50, 70))                  #######  check this - right depth to be using?  
db$bathbin <- cut(db$meancFSBPI, breaks=c(-7.07, -0.0001, 0, 0.25, 0.667, 11.4))

unique(d$latbins)
unique(db$latbins)
unique(d$depbins)
unique(db$depbins)
unique(d$MEANCFSBPI)
unique(db$bathbin)

table(db$latbins, useNA="always")                                               # why are there NAs?  I'm not sure.
table(db$depbins, useNA="always")
table(db$bathbin, useNA="always")

db$RSpred <- NA                                                                # column for new predictions
db$RSpredSE <- NA                                                              # column for SE of predictions
for (i in 1:nrow(samp)) {                                                       # merge predictions!  takes a long time.
  db$RSpred[which(db$latbins == samp$latbins[i] & db$depbins==samp$depbins[i])] <- samp$pred[i]                        #& db$bathbin == samp$MEANCFSBPI[i]
  db$RSpredSE[which(db$latbins == samp$latbins[i] & db$depbins==samp$depbins[i])] <- samp$glmse[i] }                    # & db$bathbin == samp$MEANCFSBPI[i]

db <- db[names(db)!="latbins"]                                                  # remove the new columns that are no longer needed
db <- db[names(db)!="depbins"]
db <- db[names(db)!="bathbin"]

cols <- rainbow
map('usa', xlim=c(-8
