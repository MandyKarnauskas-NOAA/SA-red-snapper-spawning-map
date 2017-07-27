################################################################################
#
#   M. Karnauskas, Jul 21, 2017                                     
#   Code to update Farmer et al. 2017 analysis using more recent MARFIN data
#   Results may differ substantially because sample sizes greatly increased 
#       (515 valid sets increased from 158 in paper)
#   Also explored use of GAMMs and compared model performance to GLMMs
#   This chunk of code imports and cleans data, experiments with different model 
#   formulations, selected best GLMM and GAMM models via cross-validation, and 
#   compares parameter estimates of GLMM versus GAMM
#
################################################################################

rm(list=ls())

setwd("C:/Users/mkarnauskas/Desktop/RSmap_SA")           
source("Xvalidate.r")                                                           # cross-validation code

################################  libraries  ###################################
if (!"chron" %in% installed.packages()) install.packages("chron", repos='http://cran.us.r-project.org')
if (!"lme4" %in% installed.packages()) install.packages("lme4", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
if (!"MASS" %in% installed.packages()) install.packages("MASS", repos='http://cran.us.r-project.org')
if (!"AICcmodavg" %in% installed.packages()) install.packages("AICcmodavg", repos='http://cran.us.r-project.org')
if (!"mgcv" %in% installed.packages()) install.packages("mgcv", repos='http://cran.us.r-project.org')
if (!"grDevices" %in% installed.packages()) install.packages("grDevices", repos='http://cran.us.r-project.org')
if (!"lunar" %in% installed.packages()) install.packages("lunar", repos='http://cran.us.r-project.org')
if (!"gamm4" %in% installed.packages()) install.packages("gamm4", repos='http://cran.us.r-project.org')
library(chron)
library(lme4)
library(maps)
library(MASS)
library(AICcmodavg)
library(mgcv)
library(grDevices)
library(lunar)
library(gamm4)
################################################################################

dat <- read.table("RedSnapperMatData.csv", sep=",", header=T, na.strings = c("NA"))       # read in data
p <- read.table("RedSnapperSamplingandAbundance.csv", sep=",", header=T)

head(dat)
head(p)
table(dat$Year)
table(p$Year)

table(dat$Year, dat$Mat)
p <- p[which(p$Year < 2015),]
dat <- dat[which(dat$Year < 2015),]
table(dat$Year)
table(p$Year)
names(table(dat$Year)) == names(table(p$Year))
table(dat$Year) / table(p$Year)
barplot(table(dat$Year) / table(p$Year))                                        # remove years where data are missing

p <- p[which(!is.na(p$Latitude)),]                                              # remove NA locations
table((p$StationDepth <= 140))                                                  # use only stations in <140m depth (limit for RS)
p <- p[which(p$StationDepth <= 140),]

d1 <- as.data.frame(table(dat$PCG))
d2 <- merge(p, d1, by.x="PCG", by.y="Var1")
summary(d2$Abundance==d2$Freq)                                                  # ensure reference site file matches up with data file
d2[which(d2$Abundance!=d2$Freq),]                                               # two sites are off by N of 1; recheck when updated data are available.

length(unique(dat$PCG))                                                         # getting to know data
length(unique(dat$Latitude))

table(dat$Species)
table(dat$Gear)
barplot(table(dat$Day))
barplot(table(dat$Month))
map('usa')                                                                      # map of data
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
 
matcodes <- c(3, 7, "B", "C", "D", "G", "H")                                    # considered "spawning females" -- see Excel tab 2 in data sheet
dat$fem <- "M"
dat$fem[which(dat$Sex==2)] <- "NF"
dat$fem[which(dat$Mat %in% matcodes & dat$Sex==2)] <- "SF"
table(dat$fem, useNA="always")
table(dat$fem, dat$Mat, dat$Sex, useNA="always")
dim(dat)
#  dat <- dat[which(dat$fem=="SF"), ]  
#  dim(dat)

mis <- which(is.na(dat$TL))
out <- lm(dat$TL ~ dat$FL)
for (i in mis)  {  dat$TL[i] <- round(out$coef[1] + out$coef[2] * dat$FL[i])   }

# dat$lenbins <- cut(dat$TL, breaks=c(300, 500, 600, 700, 905))  
# table(dat$lenbins, useNA="always")

########  explore relationship between TL and distribution and maturity  #######

plot(dat$StationDepth, dat$TL)
g <- gam(TL ~ s(StationDepth), data=dat)
summary(g)
plot(g)                                # fish get larger with depth

plot(dat$Latitude, dat$TL)
g <- gam(TL ~ s(Latitude), data=dat)
summary(g)                             # fish get larger at high latitudes
plot(g)

boxplot(dat$TL ~ as.factor(dat$fem))
g <- lm(TL ~ as.factor(dat$fem), data=dat)     
summary(g)                                  #  but, spawning fish are larger
#plot(g)

datf <- dat[which(dat$fem=="SF"),]
plot(datf$StationDepth, datf$TL)
g <- gam(TL ~ s(StationDepth), data=datf)
summary(g)                              # significance much lower
plot(g)                                 # little trend of size with depth for mature fish

plot(datf$Latitude, datf$TL)
g <- gam(TL ~ s(Latitude), data=datf)
summary(g)
plot(g)                                # trend more important with latitude


#################  merge maturity categories with site database  ###############
tab <- as.data.frame.matrix(table(dat$PCG, dat$fem))
tab$PCG <- rownames(tab)

dim(p)
dim(tab)   
dmerge <- merge(p, tab, by="PCG", all.x = T)
dim(dmerge)
names(dmerge)

table((dmerge$M + dmerge$NF + dmerge$SF)==dmerge$Abundance)
dmerge[which((dmerge$M + dmerge$NF + dmerge$SF)!=dmerge$Abundance),]

dmerge <- dmerge[-which(names(dmerge)=="Species" | names(dmerge)=="Gear")]
                                                                        # rename to match existing code
names(dmerge)[4:13] <- c("day", "mon", "year", "Date", "lat", "lon", "dep", "type", "abundance", "temp")
names(dmerge)
d <- dmerge

d$M[which(is.na(d$M))] <- 0
d$NF[which(is.na(d$NF))] <- 0
d$SF[which(is.na(d$SF))] <- 0

# ll <- as.numeric(substr((strsplit(names(dmerge)[grep('0', names(dmerge))], ",")), 5, 7))        # get simpler labels for length classes
# ul <- as.numeric(substr((strsplit(names(dmerge)[grep('0', names(dmerge))], ",")), 12, 14))
# names(dmerge)[grep('0', names(dmerge))] <- paste("L",(ll+ul)/2, sep="")
# names(dmerge)
# datcol <- grep('L', names(dmerge))
# d <- c()                                                                      # for reformatting data matrix by length class
# for (i in datcol)  {  dsub <- cbind(dmerge[1:(min(datcol-1))], dmerge[i], names(dmerge[i]))  
#                       names(dsub)[(ncol(dsub)-1): ncol(dsub)] <- c("N", "lenclass") 
#                       d <- rbind(d, dsub)   }  
# dim(dmerge) * length(unique(d$lenclass))
# dim(d)
# d$N[which(is.na(d$N))] <- 0
  
map('usa', xlim=c(-82, -75), ylim=c(26, 36))                                            
points(d$lon, d$lat, cex=(d$abundance+1)/3)                                         # view locations of mature females
#points(d$lon[!is.na(d$M)], d$lat[!is.na(d$M)], pch=20, col="#FF000012")   
#points(d$lon[!is.na(d$NF)], d$lat[!is.na(d$NF)], pch=20, col="#00FF0012")   
#points(d$lon[!is.na(d$SF)], d$lat[!is.na(d$SF)], pch=20, col="#0000FF12")   

###############  extract lunar phase data using lunar package  #################
d$lunim <- lunar.illumination(as.Date(paste(d$year, "-", d$mon, "-", d$day, sep="")))   # lunar illumination
d$lun4 <- lunar.phase(as.Date(paste(d$year, "-", d$mon, "-", d$day, sep="")), name=T)   # lunar phase - 4-name format

###  variable to define position across shelf - use as alternate to latitude ###
a.x <- max(d$lon)+0.1
a.y <- min(d$lat)-0.1
d$ang <- atan((d$lat-a.y)/(d$lon-a.x))*180/pi
cols <- rainbow(100, start=0.1)
plot(d$lon, d$lat, col=cols[round(d$ang+90)])                                   # more orthogonal to depth; potentially better for analysis
################################################################################
 
d$mon <- as.factor(d$mon)                                                       # convert month and year to factors
d$year <- as.factor(d$year) 

# bin variables as finely as possible while maintaining adequate number of samples per bin

d$angbins <- cut(d$ang, breaks=c(-90, -79, -60, -45, -30, -17, 0))              # for fitting models, with latitude*depth interaction, remove -79
d$angbin2 <- cut(d$ang, breaks=c(-90, -60, -45, -30, 0))              

map('usa', xlim=c(-82, -75), ylim=c(26, 36))
points(d$lon, d$lat, col=d$angbins)
 
d$depbins <- cut(d$dep, breaks=c(10, 25, 30, 35, 40, 50, 60, 140))
d$tempbins <- cut(d$temp, breaks=c(10, 20, 22, 24, 30))
d$latbins <- cut(d$lat, breaks=seq(27.1, 35.1, 1))

table(d$angbins, useNA="always")                                                # look at binning
table(d$angbin2, useNA="always")  
table(d$depbins, useNA="always")
table(d$tempbins, useNA="always")
table(d$latbins, useNA="always")
table(d$mon, useNA="always")

windows()
barplot(table(d$angbins, useNA="always"))
barplot(table(d$depbins, useNA="always"))
barplot(table(d$tempbins, useNA="always"))
barplot(table(d$latbins, useNA="always"))
barplot(table(d$mon, useNA="always"))
#barplot(table(d$lenclass, useNA="always"))

prop <- d$SF
barplot(tapply(prop, d$angbins, mean, na.rm=T))                                # look at percent spawning F by bin
barplot(tapply(prop, d$depbins, mean, na.rm=T))
barplot(tapply(prop, d$tempbins, mean, na.rm=T))
barplot(tapply(prop, d$latbins, mean, na.rm=T))
barplot(tapply(prop, d$mon, mean, na.rm=T))
barplot(tapply(prop, d$lun4, mean, na.rm=T))
#barplot(tapply(prop, d$lenclass, mean, na.rm=T))

table(d$angbins, d$mon)
d <- d[which(d$mon!=3),]
d$mon <- droplevels(d$mon)
table(d$angbins, d$mon)
table(d$angbins, d$depbins)
table(d$angbins, d$tempbins)
table(d$depbins, d$tempbins)
                  
prop <- d$SF
round(tapply(prop, list(d$angbins, d$mon), mean, na.rm=T), 3)  # month * latitude interaction can be done with fewer bins
round(tapply(prop, list(d$angbins, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction can be done with fewer bins
round(tapply(prop, list(d$angbin2, d$tempbins), mean, na.rm=T), 3)  # temp * latitude interaction can be done
round(tapply(prop, list(d$angbin2, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction can be done 
#round(tapply(prop, list(d$lenclass, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction can be done 
#round(tapply(prop, list(d$lenclass, d$angbins), mean, na.rm=T), 3)   # depth * latitude interaction can be done 
barplot(tapply(prop, list(d$angbins, d$tempbins), mean, na.rm=T), beside=T, col=1:6)     # looks to be important
barplot(tapply(prop, list(d$angbins, d$depbins), mean, na.rm=T), beside=T, col=1:6)
matplot(tapply(prop, list(d$angbins, d$tempbins), mean, na.rm=T), type="l")         # interactions appear to be potentially important
matplot(tapply(prop, list(d$depbins, d$angbins), mean, na.rm=T), type="l")
#matplot(tapply(prop, list(d$latbins, d$lenclass), mean, na.rm=T), type="l")         # interactions appear to be potentially important
#matplot(tapply(prop, list(d$depbins, d$lenclass), mean, na.rm=T), type="l")
#barplot(tapply(prop, list(d$lenclass, d$depbins), mean, na.rm=T), beside=T, col=1:6)

barplot(tapply(prop, d$year, mean, na.rm=T))

#yrs0 <- as.numeric(names(which(tapply(d$N, d$year, mean, na.rm=T)>0)))        # take out years with no female data
#dim(d)
#d <- d[which(d$year != 2006),];                                                 # take out 2006 - only one observation
#dim(d)
#d <- d[which(d$year %in% yrs0),]; dim(d)   

d$SF
d$pres <- d$SF                                                                  # model catch of spawning females in delta-GAM
d$pres[which(d$pres>1)] <- 1
d$SF[which(d$SF==0)] <- NA
d$SF

#############################   GAM MODEL   ####################################

d$doy <- NA                                                                     #  for GAM, can use continuous day of year instead of month
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(d))  {  d$doy[i] <- (sum(dinmon[1:d$mon[i]]) + d$day[i]) / 365  }
d$lunar <- lunar.phase(as.Date(paste(d$year,"-",d$mon,"-",d$day,sep="")), name=F)  # also can use continuous lunar phase

#  FACTORS:  year   mon     depbins    tempbins    lunar     angbins           
#                   doy     dep        temp        lunim     ang
#                                                            lat

# d <- d[which(d$year==2010 | d$year==2011 | d$year==2012 | d$year==2013 | d$year==2014),]

##########################    model presence/absence   #########################
gam1 <- gam(pres ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam1)
par(mfrow=c(5,6), mex=0.5)
plot(gam1)

gam2 <- gam(pres ~  s(dep) + s(ang) + s(doy) + s(lunim) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam2)                                                                   # lunar illumination not significant
plot(gam2)

gam3 <- gam(pres ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam3)                                                                   # lunar phase not significant here
plot(gam3); plot.new()

gam4 <- gam(pres ~ s(dep) + s(ang) + s(doy) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam4)
plot(gam4); plot.new()

gam5 <- gam(pres ~ s(dep) + s(ang) + s(doy) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam5)
plot(gam5)                                                                      # including temp appears to help doy fit more dome-shaped as would be expected

gam6 <- gam(pres ~ te(dep, ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam6)                                                                   # much lower deviance
# windows()
plot(gam6)

gam7 <- gam(pres ~ s(dep) + te(doy, ang) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam7)
plot(gam7)

gam8 <- gam(pres ~ s(dep) + te(ang, temp) + s(doy) + s(lunar) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam8)                                                                   # also lower deviance
plot(gam8)

gam9 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam9)                                                                   # highest deviance explained
plot(gam9)                                                                      # lunar not significant

gam9b <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam9b)                                                                   # highest deviance explained
plot(gam9b)

gam10 <- gam(pres ~ te(dep, lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam10)
plot(gam10)

gam11 <- gam(pres ~ ti(dep, lat) + ti(dep) + ti(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML")
summary(gam11)
plot(gam11)

extractAIC(gam1)
extractAIC(gam2)
extractAIC(gam3)
extractAIC(gam4)
extractAIC(gam5)
extractAIC(gam6)
extractAIC(gam7)
extractAIC(gam8)
extractAIC(gam9)                                                                #  gam9 is best model by AIC and deviance explained 
extractAIC(gam10)
extractAIC(gam11)

min(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10), extractAIC(gam11))[2,])
plot(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10), extractAIC(gam11))[2,], ylab="")

#  compare to GAMM package 
gamm1 <- gamm(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, random=list(year=~1), data=d)       # does not converge
summary(gam9)                                                                   
summary(gamm1$gam)
par(mfrow=c(2,6), mex=0.5)                                                      # parameter estimates are similar 
plot(gam9)
plot(gamm1$gam)

gamPAfin <- gam9

##############   optimize smoothing parameter for best GAM model  ##############
sp <- gam9$sp
tuning.scale <- c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5)
scale.exponent <- log10(tuning.scale)
n.tuning <- length(tuning.scale)
edf   <- rep(NA,n.tuning)
mn2ll <- rep(NA,n.tuning)
aic   <- rep(NA,n.tuning)
bic   <- rep(NA,n.tuning)

for (i in 1:n.tuning) {
gamobj <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML", 
          sp=tuning.scale[i]*sp)
mn2ll[i] <- -2*logLik(gamobj)
edf[i] <- sum(gamobj$edf) + 1
aic[i] <- AIC(gamobj)
bic[i] <- BIC(gamobj)   }
par(mfrow=c(2,2))
plot(scale.exponent, mn2ll, type="b", main="-2 log likelihood")
plot(scale.exponent, edf, ylim=c(0,70), type="b", main="effective number of parameters")
plot(scale.exponent, aic, type="b", main="AIC")
plot(scale.exponent, bic, type="b", main="BIC")
opt.sp <- tuning.scale[which.min(bic)] * sp

gamopt <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=binomial, data=d, method="REML", sp=opt.sp)
summary(gamopt)             # does not seem to improve fit
plot(gamopt)

############################   CROSS-VALIDATION  ###############################

drmna <- d[which(!is.na(d$temp)),]
dlast5 <- drmna[which(drmna$year==2010 | drmna$year==2011 | drmna$year==2012 | drmna$year==2013 | drmna$year==2014),]

par(mfrow=c(5,6), mex=0.5)
x <- xvalid(gam1, drmna, kfold=5)
x; colMeans(x)     

x <- xvalid(gam9, drmna, kfold=5)       # gam9 outperforms other gam models
x; colMeans(x)                      # total FPR + FNR = 0.5382341

x <- xvalid(gamopt, drmna, kfold=5)     # smoothing optimizer does not really improve performance
x; colMeans(x)                      # total FPR + FNR = 0.5480571

x <- xvalid(gam9, dlast5, kfold=5)  # using older data does not seem to reduce performance
x; colMeans(x)                        # total FPR + FNR = 0.5286315

x <- xvalid(out4, d, kfold=5)
colMeans(x, na.rm=T)                # GAM outperforms GLM


######################    model abundance when present    #####################                                                    
gam1 <- gam(SF ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), family=poisson, data=d, method="REML")
summary(gam1)
par(mfrow=c(6,6), mex=0.5)
plot(gam1)

gam2 <- gam(SF ~  s(dep) + s(lat) + s(doy) + s(lunim) + s(temp) + s(year, bs="re"), family="poisson", data=d, method="REML")
summary(gam2)                                                                   
plot(gam2)

gam3 <- gam(SF ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(year, bs="re"), family="poisson", data=d, method="REML")
summary(gam3)                                                                   
plot(gam3); plot.new()

gam4 <- gam(SF ~ s(dep) + s(lat) + s(doy) + s(temp) + s(year, bs="re"), family="poisson", data=d, method="REML")
summary(gam4)
extractAIC(gam4)
plot(gam4); plot.new()            

gam5 <- gam(SF ~ s(dep) + s(doy) + s(temp) + s(year, bs="re"), data=d, family=poisson, method="REML")
summary(gam5)
plot(gam5); plot.new(); plot.new()                                                                       

gam6 <-  gam(SF ~ s(doy) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam6)                                                                   
plot.new() ; plot(gam6)

gam7 <- gam(SF ~ s(dep) + te(doy, lat) + s(lunar) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam7)
plot(gam7)

gam8 <- gam(SF ~ s(dep) + te(lat, temp) + s(doy) + s(lunar) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam8)                                                                  
plot(gam8)

gam9 <- gam(SF ~ s(dep) + s(ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam9)                                                                   
plot(gam9)                                                                   

gam9b <- gam(SF ~ s(dep) + s(ang) + s(doy) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam9b)                                                                  
plot(gam9b)

gam10 <- gam(SF ~ te(dep, ang) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam10)
plot(gam10)

gam11 <- gam(SF ~ ti(dep, ang) + ti(dep) + ti(lat) + s(doy) + s(lunar) + s(temp) + s(year, bs="re"), data=d, family="poisson", method="REML")
summary(gam11)
plot(gam11)

extractAIC(gam1)
extractAIC(gam2)
extractAIC(gam3)
extractAIC(gam4)                                                                #  gam5 is best model by AIC and deviance explained 
extractAIC(gam5)
extractAIC(gam6)
extractAIC(gam7)
extractAIC(gam8)
extractAIC(gam9)                                                                
extractAIC(gam10)
extractAIC(gam11)

min(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10), extractAIC(gam11))[2,])
plot(cbind(extractAIC(gam1), extractAIC(gam2), extractAIC(gam3), extractAIC(gam4), extractAIC(gam5), extractAIC(gam6), extractAIC(gam7), extractAIC(gam8), extractAIC(gam9), extractAIC(gam10), extractAIC(gam11))[2,], ylab="")

#  compare to GAMM package 
gamm1 <- gamm(SF ~ s(dep) + s(doy) + s(temp), random=list(year=~1), data=d, family="poisson")     
summary(gam5)                                                                   
summary(gamm1$gam)
par(mfrow=c(2,4), mex=0.5)                                                      # parameter estimates are similar 
plot(gam5)
plot(gamm1$gam)                                                                 # 

gamNPfin <- gam5

#####################################  GLM  ####################################
# mixed effects model - year as random effect
out1 <- glmer(pres ~ depbins + latbins + mon + lun4 + tempbins +  (1|year),  family="binomial", data=d)
summary(out1)                  
extractAIC(out1) 

out2 <- glmer(fem ~ depbins + angbins + mon + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out2)                  
extractAIC(out2)

out3 <- glmer(fem ~ depbins + angbins + mon + lun4 + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out3)                  
extractAIC(out3)

out4 <- glmer(pres ~ depbins + angbins + mon + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out4)                  
extractAIC(out4) 

out5 <- glmer(fem ~ depbins + latbins + mon + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out5)                  
extractAIC(out5) 

out6 <- glmer(fem ~ depbins * angbin2 + mon + lun4 + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out6)                  # does not converge well 
extractAIC(out6) 

out7 <- glmer(fem ~ depbins + angbin2 * tempbins + lun4 + mon + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(out7)                  # does not converge well
extractAIC(out7) 

extractAIC(out1)                                                                # 1218.533
extractAIC(out2)                                                                # 1217.303
extractAIC(out3)                                                                # 1218.389
extractAIC(out4)                                                                # 1217.029    lowest AIC but not by much
extractAIC(out5)                                                                # 1218.478
extractAIC(out6)                                                                # 1222        interaction factors don't seem to help
extractAIC(out7)                                                                # 1290

################################################################################

outnull <- glm(fem ~ 1, family=binomial(logit), data=d)
(deviance(outnull)-deviance(out3))/deviance(outnull)      # deviance explained by all factors combined

outnorand <- glmer(fem ~ 1 + (1|year), family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
(deviance(outnorand)-deviance(out4))/deviance(outnorand)  # deviance explained by fixed factors combined

############################   CROSS-VALIDATION  ###############################

par(mfrow=c(5,6))
x1 <- xvalid(out1, d, kfold=5)
x2 <- xvalid(out2, d, kfold=5)
x3 <- xvalid(out3, d, kfold=5)
x4 <- xvalid(out4, d, kfold=5)                                                  # best model? 
x5 <- xvalid(out5, d, kfold=5)
#                                                                                 Area.Under.Curve    False.Positive.Rate False.Negative.Rate  total.FPR_FNR 
colMeans(x1, na.rm=T)                                                           #  0.7711004           0.3333786           0.2710542           0.6044328
colMeans(x2, na.rm=T)                                                           #  0.7657702           0.3201122           0.2773341           0.5974464
colMeans(x3, na.rm=T)                                                           #  0.7649905           0.3366355           0.2842509           0.6208864 
colMeans(x4, na.rm=T)                                                           #  0.7602100           0.2799078           0.3173104           0.5972182 
colMeans(x5, na.rm=T)                                                           #  0.7615463           0.3205856           0.2984228           0.6190084 

plot(out4)
sum(residuals(out4, type = "pearson")^2)
deviance(out4)
1 - pchisq(deviance(out4), df.residual(out4))      # not a great fit ?


############################  END GLM MODEL  ###################################




#############################  END MODELING  ###################################

#################   COMPARE GLM AND GAM PARAMETER ESTIMATES  ###################

gamf <- gam9
summary(gam9)
glmf <- glmer(fem ~ depbins + latbins + mon + lun4 + tempbins + (1|year),       # use same GLM model as comparison
family="binomial", data=d, control=glmerControl(optimizer="bobyqa")) 
summary(outf)

par(mfrow=c(3,2))
barplot(summary(glmf)$coefficients[grep('dep', rownames(summary(glmf)$coefficients)),1], las=2)
barplot(summary(glmf)$coefficients[grep('lat', rownames(summary(glmf)$coefficients)),1], las=2)
barplot(summary(glmf)$coefficients[grep('mon', rownames(summary(glmf)$coefficients)),1], las=2)
barplot(summary(glmf)$coefficients[grep('lun', rownames(summary(glmf)$coefficients)),1], las=2)
barplot(summary(glmf)$coefficients[grep('tem', rownames(summary(glmf)$coefficients)),1], las=2)

pred <- expand.grid(d$year[1], unique(d$latbins), unique(d$depbins), unique(d$mon), unique(d$lun4), unique(d$tempbins)) 
names(pred) <- c("year", "latbins", "depbins", "mon", "lun4", "tempbins")

dim(pred)
glmf
res <- predictSE(glmf, pred, type="response", se.fit=T)  # takes a long time   # predict for all combinations of factors

head(pred)
pred$pred <- res$fit
pred$glmse <- res$se.fit

windows()
par(mfrow=c(2, 6), mex=0.56)

##############  compare parameter estimates for gam versus glm  ################
plot(gam9)
barplot(tapply(pred$pred, pred$depbins, mean, na.rm=T), xlab="depth (bins)")
barplot(tapply(pred$pred, pred$latbins, mean, na.rm=T), xlab="lat (bins)")
barplot(tapply(pred$pred, pred$mon, mean, na.rm=T), xlab="month")
barplot(tapply(pred$pred, pred$lun4, mean, na.rm=T), xlab="lunar phase")
barplot(tapply(pred$pred, pred$temp, mean, na.rm=T), xlab="temp")               # results look very similar

##########################  SAVE MODEL OUTPUTS  ################################

save("gamNPfin", "gamPAfin", file="model_parameters.RData")                             # save final model results
# save("gamf", "glmf", file="model_parameters.RData")                             # save final model results
save("d", file="model_data.RData")   

##################################  END  #######################################

