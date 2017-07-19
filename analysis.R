
rm(list=ls())

setwd("C:/Users/mkarnauskas/Desktop/RSmap_SA")           
source("Xvalidate.r")

##########  libraries  ##############
#library(arm)
#library(lme4)
library(maps)
#library(MASS)
#library(Hmisc) 
#library(plyr)  
#library(AICcmodavg)

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
 
d$day <- as.factor(d$day) 
d$mon <- as.factor(d$mon) 
d$year <- as.factor(d$year) 
d$lat <- as.numeric(d$lat)
d$lon <- as.numeric(d$lon)
d$dep <- as.numeric(d$dep)
d$temp <- as.numeric(d$temp)
d$avTL <- as.numeric(d$avTL)
d$avwt <- as.numeric(d$avwt)

f <- d$SF / d$TF
plot(d$lon, d$lat, cex=f*3)

#dat$Lunar3 <- ifelse(dat$Lunar3==-1,1,dat$Lunar3)     #-1 and 1 are both Full Moon
#dat$Lunar2 <- ifelse(dat$Lunar2==-1,1,dat$Lunar2)     #-1 and 1 are both Full Moon                                  
#dat$Lunar <- as.factor(dat$Lunar3)

dim(d)
d <- d[which(d$TF>0), ]                                                         # take all samples which have at least 1 female
dim(d)


v <- which(names(d) == "MAXCASPECT") : which(names(d) == "MEANCSLOPE")      # columns for bathymetric variables
ptab <- matrix(NA, nrow=length(v), ncol=3)
colnames(ptab) <- c("lessthan0", "zero", "greaterthan0")
for (i in v)  {                                     # look at distributions of bathymetric variables
  vr <- d[,i]
  ptab[i-min(v)-1,1] <- round(length(which(vr<0))/length(vr),2)
  ptab[i-min(v)-1,2] <- round(length(which(vr==0))/length(vr),2)
  ptab[i-min(v)-1,3] <- round(length(which(vr>0))/length(vr),2)  }
cbind(names(d)[v], ptab)  
  
#  Worked with binned bathymetric variables.  Using a systematic binning method 
#  (quantiles) it ultimately led to better statistical support (lower AICs and 
#  better diagnostics). Also, outputs suggest the relationships are not linear.  
for (i in v)  {          # bin bathymetric variables 
  vr <- d[,i]          # negative numbers are a separate bin unless there are fewer than 100 negative obs 
  table(vr)              # zeros are a separate bin
  a <- min(vr) *1.01     # positive numbers are binned by quantile (thirds unless highly skewed distribution, then 50% quantile used)
  b <- max(vr) *1.01
  lev <- round(quantile(vr[which(vr>0)], probs=seq(0,1,0.5)),5) 
  if (lev[2] == lev[3]) {   lev <- round(quantile(vr[which(vr>0)], probs=seq(0,1,0.5)),5)  }
  lev[1] <- 0
  brks <- c(a, -0.0001, lev)
  if (length(which(vr<0)) < 10)  {   brks <- c(a, lev)}  
  if (a==0)  { brks <- c(-0.0001, lev) }
  lev[which.max(lev)] <- b
  print(names(d)[i])
  print(table(cut(d[,i], breaks = brks)))
  d[,i] <- cut(d[,i], breaks = brks)  }



# bin variables as finely as possible while maintaining adequate number of samples per bin
#  d$Month <- cut(d$Month, breaks=seq(2,10,2))

quantile(d$Latitude, probs=c(0.33, 0.66))
quantile(d$Samp_Depth)
quantile(d$Temp, na.rm=T)
d$lonbins <- cut(d$Longitude, breaks=seq(-81, -75, 2))
d$latbins <- cut(d$Latitude, breaks=c(27, 31.2, 32.2, 35))      # based on quantiles and also location of Charleston Bump
d$depbins <- cut(d$Samp_Depth, breaks=c(10, 30, 40, 50, 70))
d$lunarbins <- cut(d$Lunar2, breaks=seq(-1,1,0.5))
d$tempbins <- cut(d$Temp, breaks=c(15, 20, 22, 24, 28))
d$salbins <- cut(d$Sal, breaks=c(33, 34, 35, 36, 37, 40))

table(d$latbins, useNA="always")
table(d$depbins, useNA="always")
table(d$tempbins, useNA="always")
table(d$lunarbins, useNA="always")
table(d$Lunar, useNA="always")

map('usa', fill = 1, interior=F, col = gray(0.95), ylim=c(24,37), xlim=c(-90, -75)); axis(1); axis(2)
points(d$Lon, d$Lat, pch=19)
points(-78.8, 31.7, col=2, pch=19)
abline(h=c(27, 31.2, 32.2, 35), col=3)         # view latitude binning

windows()
plotfactors()                   # look at raw percentages by factor
for (i in v)  { barplot(tapply(d$Females/d$allF, d[i], mean, na.rm=T), main=names(d)[i]) }

table(d$latbins, d$Month)
table(d$latbins, d$depbins)
table(d$latbins, d$tempbins)
table(d$depbins, d$tempbins)
                  
round(tapply(d$Females/d$allF, list(d$latbins, d$tempbins), mean, na.rm=T), 3)  # temperature * latitude interaction would be nice, but not enough data
round(tapply(d$Females/d$allF, list(d$latbins, d$depbins), mean, na.rm=T), 3)   # depth * latitude interaction possible!
barplot(tapply(d$Females/d$allF, list(d$latbins, d$tempbins), mean, na.rm=T), beside=T)
barplot(tapply(d$Females/d$allF, list(d$depbins, d$latbins), mean, na.rm=T), beside=T)

tapply(d$Females/d$allF, d$Year, mean, na.rm=T)
               
# Month and temp are highly correlated - can't use both
# no interaction factors because would have to exclude data where highest spawning F encounter rates occur
  
# yrs0 <- as.numeric(names(which(tapply(d$Females/d$allF, d$Year, mean, na.rm=T)>0)))   # need to take out 1999
# d <- d[which(d$Year %in% yrs0),]; dim(d)   
  
# fixed effects model
outfixed <- glm(cbind(Females, NSFemales)  ~ lunarbins + latbins + depbins + tempbins + Year,  family="binomial", data=d)
summary(outfixed)  

# mixed effects model - year as random effect
outrand <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year),  family="binomial", data=d, control=glmerControl(optimizer="bobyqa"))
summary(outrand)                  
extractAIC(outrand) 

v <- which(names(d) == "MAXCASPECT") : which(names(d) == "MEANCSLOPE")
addn <- data.frame(matrix(NA, nrow=length(v), ncol=3))
# determine which bathymetric variable has best AIC support. 
# Can use random effects model by changing out commented lines.  Loop takes much longer and gives similar result. 
for (i in v)  {
     d$nvar <- d[,i]
      out1 <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + nvar, family="binomial", data=d, na.action=na.omit, control=glmerControl(optimizer="bobyqa"))
#      out1 <- glm(cbind(Females, NSFemales) ~ Lunar + latbins + depbins + tempbins + Year + nvar, family="binomial", data=d)
      addn[(i- v[1]+1),1] <- names(d[i]) 
      addn[(i- v[1]+1),2] <- extractAIC(out1)[2]         }
  names(addn) <- c("variable", "model AIC", "deltaAIC")
  addn$deltaAIC <- addn$model - extractAIC(outrand)[2]   
#  addn$deltaAIC <- addn$model - extractAIC(outfixed)[2]       # calculate difference in AIC with bathymetric variable included 
  addn      
  addn[order(addn$deltaAIC),]
  
########  final random effects model specified here -- check carefully  ########
outfrand <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + MEANCBSBPI, family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
summary(outfrand)
################################################################################

outnull <- glm(cbind(Females, NSFemales) ~ 1, family=binomial(logit), data=d)
(deviance(outnull)-deviance(outfrand))/deviance(outnull)      # deviance explained by all factors combined

outnorand <- glmer(cbind(Females, NSFemales) ~ 1 + (1|Year), family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
(deviance(outnorand)-deviance(outfrand))/deviance(outnorand)  # deviance explained by fixed factors combined

outnobath <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year), family=binomial(logit), data=d, control=glmerControl(optimizer="bobyqa"))
bathdev <- (deviance(outnobath)-deviance(outfrand))/deviance(outnobath); bathdev  # deviance explained by bathymetric variable

############################   RANDOM DATA TEST  ###############################
#######   this takes a really long time 
#######   for initial exploration, can convert to fixed factors which is much faster
reps <- 10
permtest <- rep(NA, reps)                                              # set up vector to contain results
for (k in 1:reps)  {                                                   # repeats over defined loop
  drand <- d;  addn$rand_AIC <- NA; addn$rand_deltaAIC <- NA           # clear random data simulation storage structures
    for (i in v)  {                                                    # loop through bathymetric variables
      drand[,i] <- sample(d[,i], nrow(d), replace=FALSE)               # resample each bathymetric variable
      drand$nvar <- drand[,i]                                          # assign to name "nvar" and fit in random effects model
        out2 <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + nvar, family="binomial", data=drand, na.action=na.omit, control=glmerControl(optimizer="bobyqa",  tol = 1e-1))
        addn$rand_AIC[(i-v[1]+1)] <- extractAIC(out2)[2]       }       # grab AIC from each random variable model
  addn$rand_deltaAIC <- addn$rand_AIC - extractAIC(outrand)[2]         # calculate difference in AIC with random bathymetric variable included 
  drand$nvar <- drand[,which.min(addn$rand_deltaAIC) + v[1]-1]         # find random variable of the group which minimized AIC, then refit model with that variable  
  outr <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + nvar, family="binomial", data=drand, na.action=na.omit, control=glmerControl(optimizer="bobyqa"))
  permtest[k] <- (deviance(outnobath)-deviance(outr))/deviance(outnobath)   
  cat(k)   }
  
  hist(permtest); abline(v=bathdev, col=2)
  abline(v=quantile(permtest, probs=c(0.95)), col=3)
  quantile(permtest, probs=seq(0,1,0.01)) 
  summary(permtest > bathdev)

#save("permtest", file="random_test_results.RData")  

cat("% of time that random data explain more deviance than actual data: ", length(which(permtest > bathdev))/(length(permtest))*100, "%\n")    # deviance explained
#  % of time that random data explain more deviance than actual data: 30 %

############################   CROSS-VALIDATION  ###############################

d1 <- d[which(!is.na(d$Temp)),]
outfrand2 <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + MEANCBSBPI, data=d1, family=binomial(logit), control=glmerControl(optimizer="bobyqa"))
x <- xvalid(outfrand2, d1, rand=T, kfold=5)
colMeans(x)
#                                                                               Area.Under.Curve    False.Positive.Rate False.Negative.Rate 
#  lunarbins + latbins + depbins + tempbins + (1|Year) + MEANCBSBPI             0.8093112           0.3347222           0.4834233
#  lunarbins + latbins + depbins + (1|Year) + MEANCBSBPI                        0.7733537           0.3141026           0.4686325
#  lunarbins + latbins + depbins + tempbins + (1|Year)                          0.7594643           0.3920860           0.4651462 
#  lunarbins + latbins + depbins + tempbins + (1|Month)+ (1|Year) + MEANCBSBPI  0.7199885           0.3633333           0.6510206 

################################################################################

outfrand <- glmer(cbind(Females, NSFemales) ~ lunarbins + latbins + depbins + tempbins + (1|Year) + MEANCBSBPI, data=d, family=binomial(logit), control=glmerControl(optimizer="bobyqa"))

save("outfrand", "d", file="final_RS_model.RData")

################################  END BSB  #####################################

############################   PREDICTION PLOTS  ###############################

windows()
par(mar=c(10,1,1,1), mfrow=c(3,2))
barplot(summary(outfrand)$coefficients[grep('lat', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('dep', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('lun', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('temp', rownames(summary(outfrand)$coefficients)),1], las=2)
barplot(summary(outfrand)$coefficients[grep('C', rownames(summary(outfrand)$coefficients)),1], las=2)

levels(d$lunarbins)[2]
levels(d$tempbins)[4]

samp <- expand.grid(levels(d$lunarbins)[2], levels(d$tempbins)[4], 
d$Year[1], unique(d$latbins), unique(d$depbins), unique(d$MEANCBSBPI)) 
names(samp) <- c("lunarbins", "tempbins", "Year", "latbins", "depbins", "MEANCBSBPI")

dim(samp)
outfrand

pred <- predictSE(outfrand, samp, type="response", se.fit=T) 

head(samp)
samp$pred <- pred$fit
samp$glmse <- pred$se.fit

par(mfrow=c(3,2))
barplot(tapply(samp$pred, samp$tempbins, mean, na.rm=T), xlab="Temp (bins)")
barplot(tapply(samp$pred, samp$latbins, mean, na.rm=T), xlab="lat (bins)")
barplot(tapply(samp$pred, samp$lunarbins, mean, na.rm=T), xlab="lunar phase")
barplot(tapply(samp$pred, samp$depbins, mean, na.rm=T), xlab="depth (bins)")
barplot(tapply(samp$pred, samp$Year, mean, na.rm=T), xlab="year")
barplot(tapply(samp$pred, samp$MEANCBSBPI, mean, na.rm=T), xlab="bathymetry")

windows()
par(mfrow=c(2,1))

samp2 <- samp[which(samp$MEANCBSBPI=="(0,1]"),]       #(-7,0.0152]           
tab <- tapply(samp2$pred, list(samp2$depbins, samp2$latbins), mean, na.rm=T); tab
tabv <- tapply(samp2$glmse, list(samp2$depbins, samp2$latbins), mean, na.rm=T); tabv
r <- nrow(tab); w <- ncol(tab)
cols <- hcl(h= seq(20,300,length=r), c=100); cols2 <- hcl(seq(20,300,length=r), alpha=0.3)
matplot(t(tab), type="l", lty=1, axes=F, ylim=c(0.7, 1), col=cols, lwd=2, xlab=" ", ylab="probability of spawning condition female", 
main= paste("Temperature = ", substr(samp2$tempbins[1], 2,3), "-", substr(samp2$tempbins[1], 5,6), "degrees;  full moon" ))
axis(1, las=2, at=1:ncol(tab), lab=paste(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(2,w*2,2)], "-", substr(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(1,w*2,2)], 2, 5), " N", sep="")) 
axis(2, las=2); box()
up <- t(tab) + t(tabv)
lo <- t(tab) - t(tabv)
for (i in 1:r)  {  polygon(c(1:w,w:1), c(up[1:w,i], lo[w:1,i]), col=cols2[i], border=NA) } 
legend("bottomleft", paste(substr(rownames(tab),2,3), "-", substr(rownames(tab),5,6), "m", sep=""), col=cols, lty=1, bty="n", lwd=2)

samp2 <- samp[which(samp$MEANCBSBPI=="(1,3]"),]                  
tab <- tapply(samp2$pred, list(samp2$depbins, samp2$latbins), mean, na.rm=T); tab
tabv <- tapply(samp2$glmse, list(samp2$depbins, samp2$latbins), mean, na.rm=T); tabv
r <- nrow(tab); w <- ncol(tab)
cols <- hcl(h= seq(20,300,length=r), c=100); cols2 <- hcl(seq(20,300,length=r), alpha=0.3)
matplot(t(tab), type="l", lty=1, axes=F, ylim=c(0.7, 1), col=cols, lwd=2, xlab=" ", ylab="probability of spawning condition female", 
main= paste("Temperature = ", substr(samp2$tempbins[1], 2,3), "-", substr(samp2$tempbins[1], 5,6), "degrees;  full moon" ))
axis(1, las=2, at=1:ncol(tab), lab=paste(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(2,w*2,2)], "-", substr(unlist(strsplit(unlist(strsplit(colnames(tab), "]")), ","))[seq(1,w*2,2)], 2, 5), " N", sep="")) 
axis(2, las=2); box()
up <- t(tab) + t(tabv)
lo <- t(tab) - t(tabv)
for (i in 1:r)  {  polygon(c(1:w,w:1), c(up[1:w,i], lo[w:1,i]), col=cols2[i], border=NA) } 
legend("bottomleft", paste(substr(rownames(tab),2,3), "-", substr(rownames(tab),5,6), "m", sep=""), col=cols, lty=1, bty="n", lwd=2)

#  plot with depth bins on x-axis                        
cols <- hcl(h= seq(30,300,length=w), alpha=1); cols2 <- hcl(seq(30,300,length=w), alpha=0.2)
matplot(tab, type="l", lty=1, axes=F, ylim=c(0.9, 1), col=cols, lwd=2, xlab="depths", ylab="probability of spawning female")
axis(1, at=1:r, lab= paste(substr(rownames(tab),2,3), "-", substr(rownames(tab),5,6), "m", sep=""))
axis(2, las=2); box()
up <- tab + tabv
lo <- tab - tabv
for (i in 1:w)  {  polygon(c(1:r,r:1), c(up[,i], lo[r:1,i]), col=cols2[i], border=NA) }
legend("bottom", names(table(d$latbins)), col=cols, lty=1, bty="n", lwd=2, horiz=F, ncol=2)



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
