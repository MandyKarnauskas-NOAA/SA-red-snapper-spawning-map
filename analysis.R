
rm(list=ls())

setwd("C:/Users/mkarnauskas/Desktop/RSmap_SA")           
source("Xvalidate.r")

##########  libraries  ##############

if (!"chron" %in% installed.packages()) install.packages("chron", repos='http://cran.us.r-project.org')
if (!"lme4" %in% installed.packages()) install.packages("lme4", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
if (!"MASS" %in% installed.packages()) install.packages("MASS", repos='http://cran.us.r-project.org')
library(chron)
library(lme4)
library(maps)
library(MASS)
#library(Hmisc) 
#library(plyr)  
#library(arm)

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

tab <- table(dat$PCG, dat$fem)
d <- data.frame(cbind(rownames(tab), tab[,1:3]), row.names=NULL, stringsAsFactors = FALSE)
d$M <-  as.numeric(d$M)
d$NF <-  as.numeric(d$NF)
d$SF <-  as.numeric(d$SF)

d$TF <- d$NF + d$SF

d$day  <- tapply(dat$Day, dat$PCG, mean)
d$mon  <- tapply(dat$Month, dat$PCG, mean)
d$year <- tapply(dat$Year, dat$PCG, mean)
d$lat  <- tapply(dat$Latitude, dat$PCG, mean)
d$lon  <- tapply(dat$Longitude, dat$PCG, mean)
d$dep  <- tapply(dat$StationDepth, dat$PCG, mean)
d$temp <- tapply(dat$Temp, dat$PCG, mean)
d$avTL <- tapply(dat$TL, dat$PCG, mean)
d$avwt <- tapply(dat$WholeWt, dat$PCG, mean)
 
d$day <- as.numeric(d$day) 
d$mon <- as.numeric(d$mon) 
d$year <- as.numeric(d$year) 
d$lat <- as.numeric(d$lat)
d$lon <- as.numeric(d$lon)
d$dep <- as.numeric(d$dep)
d$temp <- as.numeric(d$temp)
d$avTL <- as.numeric(d$avTL)
d$avwt <- as.numeric(d$avwt)

f <- d$SF / d$TF
plot(d$lon, d$lat, cex=f*3)

# data source:  http://www.somacon.com/p570.php  
lun <- read.table("moon-phases-1990-2016-America_New_York.csv", header=T, sep=",")
lun$jln <- julian(lun$month, lun$day, lun$year)
d$jln <- julian(d$mon, d$day, d$year)
d$lunar <- NA
for (i in 1:nrow(d)) { d$lunar[i] <- as.character(lun$phase[which.min(abs(lun$jln - d$jln[i]))])  }

dim(d)
d <- d[which(d$TF>0), ]                                                         # take all samples which have at least 1 female
dim(d)

#################  variable to define position across shelf  ###################
a.x <- max(d$lon)+0.1
a.y <- min(d$lat)-0.1
b.x <- d$lon
b.y <- d$lat
d$ang <- atan((b.y-a.y)/(b.x-a.x))*180/pi
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

barplot(tapply(d$SF/d$TF, d$angbins, mean, na.rm=T))
barplot(tapply(d$SF/d$TF, d$depbins, mean, na.rm=T))
barplot(tapply(d$SF/d$TF, d$tempbins, mean, na.rm=T))
barplot(tapply(d$SF/d$TF, d$mon, mean, na.rm=T))

table(d$angbins, d$mon)
table(d$angbins, d$depbins)
table(d$angbins, d$tempbins)
table(d$depbins, d$tempbins)
                  
round(tapply(d$SF/d$TF, list(d$angbins, d$tempbins), mean, na.rm=T), 3)  # temperature * latitude interaction would be nice, but not enough data
round(tapply(d$SF/d$TF, list(d$angbins, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction possible!
barplot(tapply(d$SF/d$TF, list(d$angbins, d$tempbins), mean, na.rm=T), beside=T, col=1:6)
barplot(tapply(d$SF/d$TF, list(d$angbins, d$depbins), mean, na.rm=T), beside=T, col=1:6)

tapply(d$SF/d$TF, d$year, mean, na.rm=T)
               
# Month and temp are highly correlated - can't use both
  
yrs0 <- as.numeric(names(which(tapply(d$SF/d$TF, d$year, mean, na.rm=T)>0)))   # need to take out 1999
dim(d)
d <- d[which(d$year %in% yrs0),]; dim(d)   
  
# fixed effects model
outfixed <- glm(cbind(SF, NF) ~ depbins + angbins + mon + tempbins + lunar + year,  family="binomial", data=d)
summary(outfixed) 
stepAIC(outfixed) 

# mixed effects model - year as random effect
outrand <- glmer(cbind(SF, NF) ~ depbins + angbins + mon + lunar + tempbins + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(outrand)                  
extractAIC(outrand) 

outfrand <- glmer(cbind(SF, NF) ~ depbins + angbins + mon + (1|year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(outfrand)                  
extractAIC(outfrand) 

################################################################################

outnull <- glm(cbind(SF, NF) ~ 1, family=binomial(logit), data=d)
(deviance(outnull)-deviance(outfrand))/deviance(outnull)      # deviance explained by all factors combined

outnorand <- glmer(cbind(SF, NF) ~ 1 + (1|year), family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
(deviance(outnorand)-deviance(outfrand))/deviance(outnorand)  # deviance explained by fixed factors combined

############################   CROSS-VALIDATION  ###############################

x <- xvalid(outfrand, d, rand=T, kfold=10)
colMeans(x)
#                                                                               Area.Under.Curve    False.Positive.Rate False.Negative.Rate 

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
matplot(t(tab), type="l", lty=1, axes=F, ylim=c(0, 0.8), col=cols, lwd=2, xlab=" ", ylab="probability of spawning condition female", 
main= paste("Temperature = ", substr(samp2$tempbins[1], 2,3), "-", substr(samp2$tempbins[1], 5,6), "degrees;  full moon" ))
axis(1, las=2, at=1:ncol(tab), lab=paste(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(2,w*2,2)], "-", substr(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(1,w*2,2)], 2, 5), " N", sep="")) 
axis(2, las=2); box()
up <- t(tab) + t(tabv)
lo <- t(tab) - t(tabv)
for (i in 1:r)  {  polygon(c(1:w,w:1), c(up[1:w,i], lo[w:1,i]), col=cols2[i], border=NA) } 
legend("bottomleft", paste(substr(rownames(tab),2,3), "-", substr(rownames(tab),5,6), "m", sep=""), col=cols, lty=1, bty="n", lwd=2)



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
