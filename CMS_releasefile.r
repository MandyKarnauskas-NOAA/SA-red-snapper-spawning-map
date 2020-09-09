################################################################################
#
#   M. Karnauskas, Jul 21, 2017
#   Part III of code - run analysis.R and predictionPlots.R first
#
#   Creates release file for input to larval transport simulation
#
################################################################################
rm(list=ls())
setwd("C:/Users/mandy.karnauskas/Desktop/RSmap_SA")

#######################   libraries and functions   ############################
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
library(ncdf4)
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
library(sp)
if (!"splancs" %in% installed.packages()) install.packages("splancs", repos='http://cran.us.r-project.org')
library(splancs)
if (!"mgcv" %in% installed.packages()) install.packages("mgcv", repos='http://cran.us.r-project.org')
library(mgcv)
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(maps)

source("plotSAmap.r")                                                           # plotting code

#####################  import data  ############################################

load("model_parameters.RData")
load("SApredictionGrid.RData")
load("model_data.RData")

gamPAfin
gamNfin

names(samp)[1:2] <- c("lon", "lat")

#############################  functions   #####################################
comb.var   <- function(A, Ase, P, Pse, p) { (P^2 * Ase^2 + A^2 * Pse^2 + 2 * p * A * P * Ase * Pse)  }   # combined variance function
lnorm.mean <- function(x1, x1e) {  exp(x1 + 0.5 * x1e^2)   }
lnorm.se   <- function(x1, x1e) {  ((exp(x1e^2)-1)*exp(2 * x1 + x1e^2))^0.5  }
################################################################################

################   make predictions for release file by doy   ##################
    
plot(d$doy, d$temp)
g <- gam(temp ~ s(doy), data=d)
summary(g)
plot(g)
ref <- data.frame(unique(d$doy)); names(ref) <- "doy"
ref$temp <- predict(g, ref)
plot(ref$doy, ref$temp)
    
preddoy <- data.frame()

for (i in seq(min(d$doy), max(d$doy), 6))  {      ########  NOTE!!!!   #####   revisit this later.
        doy <- i                                      # inclusion of temp prediction based on doy causes bimodal peak in spawning season
        samp$doy <- i
        samp$lunar <- mean(d$lunar)                                             # use average lunar phase
        samp$year  <- 2014                                                      # use most recent year
        samp$temp  <- mean(d$temp, na.rm=T)   # predict(g, data.frame(doy))              # temp correlated with day of year;

        predlogit <- predict(gamPAfin, samp, type="response", se.fit=T)         # predict occurrences
        predposlog <- predict(gamNfin, samp, type="response", se.fit=T)         # predict eggs when present
        predpos   <- lnorm.mean(predposlog$fit, predposlog$se.fit)              # convert lognormal mean and SE to normal space
#        predposse <- lnorm.se(predposlog$fit, predposlog$se.fit)
#        co <- as.numeric(cor(predlogit$fit, predpos, method="pearson"))        # calculate covariance
#        predse <- (comb.var(predpos, predposse, predlogit$fit, predlogit$se.fit, co))^0.5    # calculate combined variance
        predind <-  predlogit$fit * predpos                                     # calculate combined index
    tempmat <- cbind(samp, predind)                                             # temporary save of prediction terms and predictions
    preddoy <- rbind(preddoy, tempmat)                }

nrow(tempmat) * length(seq(min(d$doy), max(d$doy), 6))
dim(preddoy)

par(mfrow=c(1,2))
plotSAmap(tempmat$predind/1000, tempmat$lon, tempmat$lat, cexnum=0.6, pchnum=15)
plotSAmap(predlogit$fit*10, tempmat$lon, tempmat$lat, cexnum=0.6, pchnum=15)

dat <- data.frame(cbind(as.numeric(preddoy$lon), as.numeric(preddoy$lat), 
                        as.numeric(preddoy$predind), as.numeric(preddoy$doy), as.numeric(preddoy$ocModDep)))
names(dat) <- c("lon", "lat", "N", "doy", "depest")

windows()
par(mfrow=c(3,5), mex=0.5)
for (i in unique(dat$doy)[1:15]) {
f <- dat[which(dat$doy==i),]
plotSAmap(f$N/1000, f$lon, f$lat, cexnum=0.6, pchnum=15)
mtext(side=3, line=0.5, paste("day of year", i)) }

##########################  assign polygon labels  #############################

co <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", header=F)                    # edited version to align with state boundaries

################################################################################

################  plot original recruitment habitat grid cells  ################
windows()
plot(1, xlim=c(-100,-75), ylim=c(24,36))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(mean(m[,1]), mean(m[,2]), m[1,3], cex=0.6)      }

points(dat$lon, dat$lat, pch=19, cex=0.25)             # release locations

################  get polygon numbers for release locations  ###################

pts = SpatialPoints(cbind(dat$lon, dat$lat))
m <- co[which(co[,3]==1), ]
m <- rbind(m, m[1,])
oldpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=1)
for (j in 2:max(co$V3))  {
  m <- co[which(co[,3]==j), ]
  m <- rbind(m, m[1,])
  newpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=j)
  oldpol <- c(oldpol, newpol)           }
pol = SpatialPolygons(oldpol)
polynams <- over(pts, pol)
rel <- cbind(polynams, dat)
rel                                          # file coded with polygon numbers

points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)
summary(is.na(rel$polynams))

##############  label release site with closest polygon number  ################
ctrslat <- tapply(co$V2, co$V3, mean)        # find center point of each polygon
ctrslon <- tapply(co$V1, co$V3, mean)        # acts as approximate method for assigning nearest polygon
points(ctrslon, ctrslat, pch=20, cex=2)

for (i in 1:nrow(rel)) {
  if (is.na(rel$polynams[i])) {
            rel$polynams[i] <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))  }}
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
summary(is.na(rel$polynams))

######################    plot new polygon assignments   #######################
plot(1, xlim=c(-85,-75), ylim=c(24,36))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, lwd=1)
  text(mean(m[,1]), mean(m[,2]), m[1,3], cex=0.6)      }
map('usa', add=T)
cols <- rainbow(117)
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])    #  check
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)          #  check

############## adjust depth 
head(rel)

points(rel$lon, rel$lat, pch=15, col=is.na(rel$depest))
rel$depest[is.na(rel$depest)] <- 0

table(rel$depest)
hist(rel$depest)
plot(rel$lon, rel$lat, pch=15, col=gray(1-(rel$depest/max(rel$depest+10)+0.04)))
plotSAmap(rel$depest, rel$lon, rel$lat, cexnum=0.6, pchnum=15)

length(which(rel$depest<45))
length(which(rel$depest<=0))
length(which(rel$depest>=45))

points(rel$lon[which(rel$depest<=0)], rel$lat[which(rel$depest<=0)], cex=1.5, pch=19, col=2)
points(rel$lon[which(rel$depest>45)], rel$lat[which(rel$depest>45)], cex=1.5, pch=19, col=3)

rel$spawndep <- rel$depest - 5                     # set spawning 5m above ocean floor
rel$spawndep[which(rel$spawndep >= 45)] <- 40       # for deeper than 45m, set at 40m

minspawndep <- 9                                  #  minimum spawning depth set at 15m based on literature

length(which(rel$depest<minspawndep))
dim(rel)
tapply(rel$N, rel$depest<minspawndep, sum)         # relative spawning biomass above 15m cutoff
rel <- rel[which(rel$depest>=minspawndep),]
dim(rel)

plot(rel$depest, rel$spawndep)                 # check assigned spawning depths
table(rel$spawndep < rel$depest)
table(rel$spawndep)
table(rel$depest, rel$spawndep)

plotSAmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)      # check on map

######################## convert temporal information ##########################

rel$mon <- as.numeric(format(strptime(rel$doy, format="%j"), format="%m"))
rel$day <- as.numeric(format(strptime(rel$doy, format="%j"), format="%d"))

lis <- 2008:2009               # loop over years in matrix below

##############  trim release file to cells with sufficient eggs   ##############

m <- rel[-c(5,6)]       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

x <- m$N /20000

m$N <-round(m$N / 2000)
mean(m$N); min(m$N); max(m$N)
table(m$N==0)
tapply(x, (m$N == 0), sum)/sum(x)          # what percentage of particles get lost by rounding?

m <- m[which(m$N>0),]

sum(m$N) * length(lis)
prod(length(lis), nrow(m))

###################### now, making the release file ############################

mat <- as.data.frame(matrix(data=NA, nrow=length(lis)*nrow(m), ncol=9))   # empty matrix to be filled
mat[,1] <- rep(m$polynams, length(lis))                                   # column 1: release polygon number
mat[,2] <- rep(m$lon, length(lis))                                        # column 2: release longitude
mat[,3] <- rep(m$lat, length(lis))                                        # column 3: release latitude
mat[,4] <- rep(m$spawndep, length(lis))                                   # column 4: release depth
mat[,5] <- rep(m$N, length(lis))                                          # column 5: number of particles per release
mat[,6] <- sort(rep(lis, nrow(m)))                                        # column 6: release year
mat[,7] <- rep(m$mon, length(lis))                                        # column 7: release month
mat[,8] <- rep(m$day, length(lis))                                        # column 8: release day
mat[,9] <- 0                                                              # column 9: release hour

sum(mat[,5])
head(mat)
table((mat$V5 > 0))

mean(mat$V5); min(mat$V5); max(mat$V5)                                  # distribution of particles

getwd()
dim(mat)
dim(mat)/8
dim(mat)/7
sum(mat$V5)

write.table(mat, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/RS_ATL_release.txt", sep="\t", col.names=F, row.names=F)          # WRITE FILE TO TXT

matN <- mat[which(mat$V3 > 34),]
matS <- mat[which(mat$V3 <= 34),]

write.table(matN, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/RS_ATL_releaseHatteras_HYCOM150.txt", sep="\t", col.names=F, row.names=F)          # WRITE FILE TO TXT
write.table(matS, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/RS_ATL_releaseMain_HYCOM150.txt", sep="\t", col.names=F, row.names=F)  

#  break down by year, then above and below 34N
#  >34N run thru global HYCOM
#  <34N run thru SABGOM

rm(list=ls())
d <- read.table("RS_ATL_releaseHatteras.txt", sep="\t")
d <- d[which(d$V6==2010 & d$V7 > 2 & d$V7 < 7),]
d$V6 <- 2017
d$V4 <- 7
write.table(d, file="releaseHatteras2017.txt", sep="\t", col.names=F, row.names=F) 

rm(list=ls())
d <- read.table("RS_ATL_releaseHatteras.txt", sep="\t")
d <- d[which(d$V6==2010 & d$V7 > 2 & d$V7 < 7),]
d$V4 <- 7
write.table(d, file="releaseHatteras2010.txt", sep="\t", col.names=F, row.names=F) 


d <- read.table("RS_ATL_releaseS.txt", sep="\t")
d <- d[which(d$V6==2010),]
write.table(d, file="RS_ATL_releaseS2010.txt", sep="\t", col.names=F, row.names=F) 

####################     end construction of release file    ###################

################    double check that matrix came out ok     ###################

matfin <- matS
table(matfin$V6, matfin$V7)       # numbers in columns should be same
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2008 & matfin$V7 == 4 & matfin$V8 == 22); length(f)
plotSAmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2008 & matfin$V7 == 07 & matfin$V8 == 21); length(f)
plotSAmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

##################################   END    ####################################


