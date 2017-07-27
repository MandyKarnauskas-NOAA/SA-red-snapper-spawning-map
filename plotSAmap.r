
###########################   MAPPING FUNCTION   ###############################
plotSAmap <- function(x, lon, lat, pchnum, cexnum, cutof=NA, startzero=F)  {
    pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 500, 1000)
  if (startzero==F)  {a <- floor(min(x)) } else { a <- 0 }
  b <- max((x)-a)*1.03
  pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
  cols <- rainbow(100, start=0.01, end=0.7)[100:1]; if(is.na(pchnum)) {pchnum=15}
#  cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]; if(is.na(pchnum)) {pchnum=15}

  map("state", interior = TRUE, xlim=c(-81.75, -75), ylim=c(26.8, 35.2))
  points(lon, lat, col=cols[pind], pch=pchnum, cex=cexnum)
  axis(1, at=c(-80, -78, -76), lab=c(expression(paste(80,degree,"W")), expression(paste(78,degree,"W")), expression(paste(76,degree,"W"))))
  axis(2, at=c(28, 30, 32, 34), lab=c(expression(paste(28,degree,"N")), expression(paste(30,degree,"N")), expression(paste(32,degree,"N")), expression(paste(34,degree,"N"))), las=1)
  box()
  yloc <- seq(28, 32, length.out=100)
 for (j in 1:100) {   polygon(c(-77, -76.5, -76.5, -77), c(yloc[j], yloc[j], yloc[j+1], yloc[j+1]), col=cols[j], border=NA) }
   w <- which.min(abs(((max(x)-min(x))/6) - pos))
  if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
  text(x=-76.6, y=yloc[round((xx-a)/b*100+1)], xx, pos=4)     }
################################################################################