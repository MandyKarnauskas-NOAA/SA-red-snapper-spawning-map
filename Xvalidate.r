#########################  X-fold X-Validation   ###############################
#########################  M. Karnauskas 7/10/2014  ############################
#########################  mandy.karnauskas@noaa.gov  ##########################
#                                                                              #
#   This code has not been thoroughly tested or peer-reviewed!                 #
#   Use with caution!!!                                                        #
#                                                                              #
################################################################################
# 
# Description: 10-fold cross-validation for binomial logistic regression models
#
# dat         original data
# model       model object to be X-validated
# kfold       number of K-fold cross validation groups
# sample ID   only needed if data are in long 1-column format (zeros and ones)
#               with sites having observations on multiple rows, and 
#               randomization by site is desired (in this case, ID is needed so 
#               that randomization is done by sites rather than data rows)
#
#  as of 7/20/2017, can take model outputs from GLM, GLMER, or GAM 
#
################################################################################

xvalid <- function(model, dat, sampleID = c(), kfold = 10)  {

if (!"pROC" %in% installed.packages()) install.packages("pROC", repos='http://cran.us.r-project.org')
library(pROC)
if (!"AICcmodavg" %in% installed.packages()) install.packages("AICcmodavg", repos='http://cran.us.r-project.org')
library(AICcmodavg)
if (!"formula.tools" %in% installed.packages()) install.packages("formula.tools", repos='http://cran.us.r-project.org')
library(formula.tools)

summaryTable <- data.frame(matrix(nrow=kfold, ncol=5, dimnames=list(c(), c("kfold", "Area Under Curve", "False Positive Rate", "False Negative Rate", "total FPR_FNR"))))

dat$rand <- runif(nrow(dat))   # random group assignment for kfold division
if (length(sampleID) > 0)   {    dat$ind <- dat[,which(names(dat)==sampleID)]
  for (i in unique(dat$ind))  {  dat$rand[which(dat$ind==i)] <- runif(1)  }}
dat$grp <- as.numeric(cut(dat$rand, quantile(unique(dat$rand), seq(0, 1, length.out = kfold + 1)), include.lowest = TRUE))

if (class(model)[1]=="glmerMod")  {  depen <- 1: (which(all.names(model@call)=="~") - 1)
                                  cols <- which(colnames(dat) %in% all.vars(model@call)[depen]) }
if (class(model)[1]=="glm")       {  depen <- 1 : (length(all.vars(model$formula)) - length(rhs.vars(model$formula))) 
                                  cols <- which(colnames(dat) %in% all.vars(model$formula)[depen]) }   # deduce # columns in dependent variable
if (class(model)[1]=="gam")       {  depen <- 1 : (length(all.vars(model$formula)) - length(rhs.vars(model$formula))) 
                                  cols <- which(colnames(dat) %in% all.vars(model$formula)[depen]) }                                  
                                  
# convert to single-column of data (0s and 1s) if necessary
   if (length(depen)==2)  {
      tpos <- dat[rep(1, dat[1,cols[1]]),]
      tneg <- dat[rep(1, dat[1,cols[2]]),]
      tpos$presence <- rep(1, nrow(tpos))
      tneg$presence <- rep(0, nrow(tneg))
      dat2 <- rbind(tpos, tneg)
    for (j in 2:nrow(dat))  {
      tpos <- dat[rep(j, dat[j,cols[1]]),]
      tneg <- dat[rep(j, dat[j,cols[2]]),]
      tpos$presence <- rep(1, nrow(tpos))
      tneg$presence <- rep(0, nrow(tneg))
      dat2 <- rbind(dat2, tpos, tneg)  }  }    
      
   if (length(depen)==1)   {  dat$presence <- dat[,cols];  dat2 <- dat   }

options(warn=2)

for (k in 1:kfold)  {
  summaryTable[k,1] <- k
  
  train <- dat2[which(dat2$grp!=k), ]
  test  <- dat2[which(dat2$grp==k), ]
  
modelFitting <- function()  {  
if (class(model)[1]=="glmerMod")   {   fmla <- as.formula(paste("presence ~ ", paste(unlist(strsplit(unlist(strsplit(as.character(model@call[2]), "~")[[1]])[2], "+", fixed=T)), collapse= "+"))) 
                                  res <- glmer(fmla, data=train, family=binomial(logit), control=glmerControl(optimizer="bobyqa"))    
                                  return(res)        }
if (class(model)[1]=="glm")     {   fmla <- as.formula(paste("presence ~ ", paste(rhs.vars(model$formula), collapse= "+")))
                                res <- glm(fmla, data=train, family="binomial")   
                                return(res)       } 
if (class(model)[1]=="gam")     {   fmla <- as.formula(paste("presence ~ ", paste(rhs.vars(model$formula), collapse= "+")))
                                res <- gam(fmla, data=train, family=binomial, method="REML")
                                return(res)       }                                                                 
                            }                 
res <- try(modelFitting(), TRUE) ; res

if (class(res)[1]=="try-error") {  plot(1,1, col=0, axes=F, xlab="", ylab=""); text(1,1, "NA", cex=2) }

if (class(res)[1]!="try-error") {

if (class(model)[1]=="glmerMod")   {   pred <- as.numeric(predictSE(res, train, type="response", se.fit=F))    }
if (class(model)[1]=="glm")        {   pred <- res$fitted.values     } 
if (class(model)[1]=="gam")        {   pred <- res$fitted.values     } 
 
#############################   ROC ANALYSIS   #################################
# true positive rate (Sensitivity) is plotted in function of the false positive rate (100-Specificity) for different cut-off points of a parameter
temproc <- roc(train$presence, pred, plot=TRUE, grid=TRUE)
summaryTable[k,2] <- temproc$auc       # CALCULATE AREA UNDER THE CURVE, CONSTRUCT MATRIX OF ROC INFORMATION FOR EACH CUTOFF ("thresholds")
roctable <- cbind(temproc$sensitivities, temproc$specificities, temproc$thresholds, temproc$sensitivities+temproc$specificities)
        # Sensitivity = proportion of actual positives which are correctly identified as such
        # Specificity = proportion of negatives which are correctly identified as such
cutoff <- roctable[roctable[,4] == max(roctable[,4]), 3]
cutoff <- min(cutoff)
a <- table(pred > cutoff, train$presence)
a[2,1] / sum(a[,1])   # FPR for training set
a[1,2] / sum(a[,2])   # FNR for training set

#######################  prediction on test data set  ##########################

if (class(model)[1]=="glmerMod") {   predlogit <- predictSE(res, test, type="response", se.fit=F)   } # predict occurrences for test data set
if (class(model)[1]=="glm")      {   predlogit <- predict.glm(res, test, type="response")  }          # predict occurrences for test data set
if (class(model)[1]=="gam")      {   predlogit <- predict.gam(res, test, type="response")  }          # predict occurrences for test data set

b <- table(predlogit > cutoff, test$presence)
   if(nrow(b)==1)  {   b <- rbind(b, c(0,0))  }
   if(ncol(b)==1)  {   b <- cbind(b, c(0,0))  }
summaryTable[k,3] <- b[2,1] / sum(b[,1])       # FPR
summaryTable[k,4] <- b[1,2] / sum(b[,2])       # FNR
summaryTable[k,5] <- sum(summaryTable[k,3:4])       # FNR
         }      # end if statement for non-error runs
         }      
         
options(warn=0)         
         
#########################   plot results   #####################################
barplot(cbind( colMeans(summaryTable, na.rm=T)[3:4], rbind(summaryTable[,3], summaryTable[,4])), args=c("top", horiz=T, bty="n"),
    beside=F, ylim=c(0,1.1), legend=chartr(".", " ", colnames(summaryTable)[3:4]), names.arg=c("MEAN", 1:kfold), axes=F)
    axis(2, at=seq(0,1,0.1), las=2)
return(summaryTable)

       }
       
##############################  end function  ##################################
################################################################################

                   

