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
p <- p[which(p$Year < 2015 & p$Year > 2003),]                                   # use years 2004 - 2014!!!!!!!!!!!!!!!
dat <- dat[which(dat$Year < 2015 & dat$Year > 2003),]                           # use years 2004 - 2014!!!!!!!!!!!!!!!
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

mis <- which(is.na(dat$TL))
out <- lm(dat$TL ~ dat$FL)
for (i in mis)  {  dat$TL[i] <- round(out$coef[1] + out$coef[2] * dat$FL[i])   }

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
  
map('usa', xlim=c(-82, -75), ylim=c(26, 36))                                            
points(d$lon, d$lat, cex=(d$abundance+1)/3)                                         # view locations of mature females

tapply(d$M, d$year, sum)
tapply(d$NF, d$year, sum)
tapply(d$SF, d$year, sum, na.rm=T)

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
d$depbins <- cut(d$dep, breaks=c(10, 25, 30, 35, 40, 50, 60, 140))
d$tempbins <- cut(d$temp, breaks=c(10, 20, 22, 24, 30))
d$latbins <- cut(d$lat, breaks=seq(27.1, 35.1, 1))

################  convert to presence - absence 
d$SF
d$pres <- d$SF                                                                  # model catch of spawning females in delta-GAM
d$pres[which(d$pres>1)] <- 1
d$SF[which(d$SF==0)] <- NA
#d$SF

#############################   GAM MODEL   ####################################

d$doy <- NA                                                                     #  for GAM, can use continuous day of year instead of month
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(d))  {  d$doy[i] <- (sum(dinmon[1:d$mon[i]]) + d$day[i]) / 365  }
d$lunar <- lunar.phase(as.Date(paste(d$year,"-",d$mon,"-",d$day,sep="")), name=F)  # also can use continuous lunar phase

#  FACTORS:  year   mon     depbins    tempbins    lunar     angbins           
#                   doy     dep        temp        lunim     ang
#                                                            lat

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
plot(gam9)      
acf(resid(gam9)) 
plot(residuals(gam9))   

gam9a <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), correlation = corAR1(form=~year), family=binomial, data=d, method="REML")
summary(gam9a)                                                                   # highest deviance explained
plot(gam9a)                                                                      # lunar not significant
acf(resid(gam9a)) 
plot(residuals(gam9a))   
                                                             
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

d04 <- d[which(d$year==2004),]
d08 <- d[which(d$year==2008),]
d09 <- d[which(d$year==2009),]
d10 <- d[which(d$year==2010),]
d11 <- d[which(d$year==2011),]
d12 <- d[which(d$year==2012),]
d13 <- d[which(d$year==2013),]
d14 <- d[which(d$year==2014),]

g04 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d04, method="REML")
g08 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d08, method="REML")
g09 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d09, method="REML")
g10 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d10, method="REML")

g11 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d12, method="REML")
g12 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d12, method="REML")
g13 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d13, method="REML")
g14 <- gam(pres ~ s(dep) + s(lat) + s(doy) + s(lunar) + s(temp), family=binomial, data=d14, method="REML")

par(mfrow=c(4,5))

summary(g04); plot(g04)      
summary(g08); plot(g08)      
summary(g09); plot(g09) 
summary(g10); plot(g10) 

summary(g11); plot(g11)      
summary(g12); plot(g12)      
summary(g13); plot(g13) 
summary(g14); plot(g14) 


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
x <- xvalid(gam1, drmna, kfold=5)      # total FPR + FNR = 0.4933012
x; colMeans(x)     

x <- xvalid(gam9, drmna, kfold=5)       # gam9 outperforms other gam models
x; colMeans(x)                      # total FPR + FNR = 0.4301723

x <- xvalid(gamopt, drmna, kfold=5)     # smoothing optimizer does not really improve performance
x; colMeans(x)                      # total FPR + FNR = 0.4622951

x <- xvalid(gam9, dlast5, kfold=5)  # using older data does not seem to reduce performance
x; colMeans(x)                        # total FPR + FNR = 0.4850022

# mixed effects model - year as random effect
out2 <- glmer(pres ~ depbins + latbins + mon + lun4 + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
out2 <- glm(pres ~ depbins + latbins + mon + lun4 + tempbins + year,  family="binomial", data=d)
summary(out2)                  
extractAIC(out2)

x <- xvalid(out2, drmna, kfold=5)    # total FPR + FNR = 0.5201849
colMeans(x, na.rm=T)                # GAM outperforms GLM

################################################################################

outnull <- gam(pres ~ 1, family=binomial, data=d, method="REML")
(deviance(outnull)-deviance(gam9))/deviance(outnull)      # deviance explained by all factors combined

outnorand <- gam(pres ~ 1 + s(year, bs="re"), family=binomial, data=d, method="REML")
(deviance(outnorand)-deviance(gam9))/deviance(outnorand)  # deviance explained by fixed factors combined

#############################  END MODELING  ###################################

##########################  SAVE MODEL OUTPUTS  ################################

save("gamPAfin", file="model_parameters.RData")                             # save final model results
save("d", file="model_data.RData")   

##################################  END  #######################################

