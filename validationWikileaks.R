## For cluster Generating x BS samples, fit trees and predict and save it all
# Thomas Rusch, 2012 

###################################################
########### Sampling        #######################
###################################################

#naive RS 
naive_rs <- function(data,fraction)
       {
         ind <- 1:dim(data)[1]
         #set.seed(seed)
         rsamp <- sample(ind,size=round(fraction*length(ind)),replace=TRUE)
         list("resamp"=rsamp,"train"=data[rsamp,],"oob"=data[-rsamp,])
       }


stratified <- function(df, id, group, size, seed="NULL", ...) {
  ### Code by Ananda from http://news.mrdwab.com/2011/05/20/stratified-random-sampling-in-r-from-a-data-frame/
  #
  #  USE: * Specify your data frame, ID variable (as column number), and
  #         grouping variable (as column number) as the first three arguments.
  #       * Decide on your sample size. For a sample proportional to the
  #         population, enter "size" as a decimal. For an equal number of
  #         samples from each group, enter "size" as a whole number.
  #       * Decide on if you want to use a seed or not. If not, leave blank
  #         or type "NULL" (with quotes). 
  #
  #  Example 1: To sample 10% of each group from a data frame named "z", where
  #             the ID variable is the first variable, the grouping variable
  #             is the fourth variable, and the desired seed is "1", use:
  # 
  #                 > stratified(z, 1, 4, .1, 1)
  #
  #  Example 2: To run the same sample as above but without a seed, use:
  # 
  #                 > stratified(z, 1, 4, .1)
  #
  #  Example 3: To sample 5 from each group from a data frame named "z", where
  #             the ID variable is the first variable, the grouping variable
  #             is the third variable, and the desired seed is 2, use:
  #
  #                 > stratified(z, 1, 3, 5, 2)
  #
  k = unstack(data.frame(as.vector(df[id]), as.vector(df[group])))
  l = length(k)
  results = vector("list", l)

  if (seed == "NULL" & size < 1) {
      for (i in 1:length(k)) {
        N = k[[i]]
        n = round(length(N)*size)
        results[[i]] = list(sample(N, n,...))
      }
    } else if (seed == "NULL" & size >= 1) {
      for (i in 1:length(k)) {
        N = k[[i]]
        results[[i]] = list(sample(N, size, ...))
      }
    } else if (size < 1) {
      for (i in 1:length(k)) {
        set.seed(seed)
        N = k[[i]]
        n = round(length(N)*size)
        results[[i]] = list(sample(N, n, ...))
      }
    } else if (size >= 1) {
      for (i in 1:length(k)) {
        set.seed(seed)
        N = k[[i]]
        results[[i]] = list(sample(N, size, ...))
      }
    }
  z = data.frame(c(unlist(results)))
  names(z) = names(df[id])
  w = merge(df, z)
  w[order(w[group]), ]
}


pred_bs <- function(afg1,targetnumber,fraction=5/6,seed=4711,sampling=c("naive","stratified"), mcontrol=mob_control(alpha=1e-03,minsplit=200))
  #creates the bootstrap samples, fits the trees and predicts the training and oob observations
  #saves the trees and the training set as RRStreeSeed-run.rda (or SRS if stratified)
  #outputs a list of the partitions for all BS samples as elements
  #
  #afg1 .. data
  #targetnumber ... number of bs samples
  #fraction... relative size of the bs sample to the data set
  #seed ... seed
  #sampling.. choice of stratified or random sampling
  #mcontrol... controls parameters to pass to mob
  #
  #needs the original tree object in file fmAll.rda
  {
   res <- vector("list",1)
   i <- 1
   set.seed(seed)
   while(length(res[sapply(res, function(x) !inherits(x, "try-error"))])<targetnumber)
   {
   afgRRS <- naive_rs(afg1,fraction)
   if(sampling=="stratified")
        {
          load("fmAll.RData")
          afg1 <- cbind(afg1,unclass(fmAll)@where)
          afgRRS <- list()
          afgRRS[[1]] <- stratified(afg1,1,dim(afg1)[2],fraction,replace=TRUE)
          afgRRS[[2]] <- afg1[-unique(afgRRS[[1]][,"ids"]),]
          names(afgRRS)[1] <- "train"
          names(afgRRS)[2] <- "oob"
          rm(fmAll)
        }
    cat("###### RUN ",i," #######\n")
    cat("Size unique training: ",dim(unique(afgRRS[["train"]]))[1],"\n")
    cat("** Model fit **\n")
    rrs<- mob(kiaAll ~ 1 |attackOn+affiliation+dcolor+region+complexAttack+t1+t2+t3+t4+Topic_5+t6+t7+t8+t9+t10+t11+t12+t13+Topic_14+t15+t16+t17+Topic_18+Topic_19+t20+t21+t22+t23+t24+t25+t26+Topic_27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+Topic_61+t62+t63+t64+t65+t66+t67+t68+t69+t70+Topic_71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+Topic_85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100, data = afgRRS[["train"]], model = negbinModel, na.action=NULL, control=mcontrol)
   if(sampling=="naive")      save(afgRRS,rrs,file=paste("RRStree",sqrt(seed),"-",i,".rda",sep=""))
   if(sampling=="stratified") save(afgRRS,rrs,file=paste("SRStree",seed,"-",i,".rda",sep=""))                  
   cat("** Prediction **\n")
   res[[i]] <-try(predict(rrs,newdata=afg1,type="node"))
   i <- i+1
   rm(rrs)
   gc()
  }
  res
 }

###################################################
########### Calculating the Jaccard Index  ########
###################################################


segment_stab <- function(origpart, bspart)
  {
    #calculates overall and segment-wise stabilty statistics (Jaccard agreement and (adjusted) rand index)
    #
    #origpart ... original segmentation (all observations and their terminal node
    #bspart... bootstrap segmentation (all observations and their terminal node) 
    require("flexclust")
     #Contingency table
    tab <- table(origpart,bspart)
    #matrix of segment-wise jaccard agreements 
    Imat <- matrix(NA,ncol=dim(tab)[2],nrow=dim(tab)[1])
    for(i in 1:dim(tab)[1])
      {
      for(j in 1:dim(tab)[2])
        {
          Imat[i,j] <- tab[i,j]/(rowSums(tab)[i]+colSums(tab)[j]-tab[i,j]) #jaccard agreement 
        }                 
      }
      #rand index
      rando <- randIndex(tab,correct=FALSE)
      #adjusted rand index
      randa <- randIndex(tab)
      names(rando) <- "RandIndex"
      names(randa) <- "AdjRandIndex"
      #segments of the bootstrap segmentation with the highest Jaccard per original segment (the corresponding segments) 
      maxpos <- apply(Imat,1,function(x) colnames(tab)[which(x==max(x))]) 
     #the highest Jaccard per original segment of the bootstrap segmentation 
      maxval <- apply(Imat,1,max)
      corsis <- cbind(table(origpart),as.numeric(maxpos),maxval)
      #nice display
      colnames(corsis) <- c("original","corresponding","Index")
      out <- list("randIndex"=c(rando,randa),"contingency"=tab,"indices"=Imat,"maxindex"=maxval,"corresp"=corsis)
      out
    }

bs_segment_stab <- function(origpart,bsparts)
  #applies segment_stab() over all bootstrap samples 
  {
    out <- lapply(bsparts,function (x) segment_stab(origpart,x))
    out
  }

###################################################
########### Convenience functions #################
###################################################

get_corresptable <- function(bssegmentstab)
  #extract only the nice table
  {
    out <- lapply(bssegmentstab, function (x) x[["corresp"]])
    out 
  }

get_rands <- function(bssegmentstab)
  #extract the rand indices
  {
    out <- lapply(bssegmentstab, function (x) x[["randIndex"]])
    out 
  }

get_corresp <- function(bssegmentstab)
  #extract the corresponding segments
  {
    out <- lapply(bssegmentstab, function (x) as.data.frame(x[["corresp"]][,"corresponding"]))
    return(out)
  }

get_maxindex <- function(bssegmentstab)
  #extract the maxindex per given segment
  {
    out <- lapply(bssegmentstab, function (x) x[["maxindex"]])
    return(out)
  }


###########################################
#### Cluster Code    ######################
###########################################

library(party)
load("afgClean.RData")
source("fitfunctionsWikileaks.R")
source("stratify.R")

load("fmAll.RData")
afg1 <- afg[!is.na(afg$kiaAll),]
ids <- 1:dim(afg1)[1]
afg1 <- cbind(ids,afg1)


#stratified RS
#seed <- as.integer(Sys.getenv("SGE_TASK_ID"))
seed <- 4711
preds <- pred_bs(afg1,5,seed=seed,sampling="stratified",mcontrol=mob_control(alpha=1e-03,minsplit=250)) 
predRRS <- preds[sapply(preds, function(x) !inherits(x, "try-error"))]
save(preds,predRRS,file=paste("predictionsSRS-",seed,".rda"))


#naive resampling 
seed <- seed^2
preds <- pred_bs(afg1,5,sampling="naive",seed=seed) 
predRRS <- preds[sapply(preds, function(x) !inherits(x, "try-error"))]
save(preds,predRRS,file=paste("predictionsRRS-",seed,".rda"))

##########################################################
###### Execution for Stability Validation ################
###### WARNING: Messy interactive code ahead #############
##########################################################

    
load("fmAll.RData")
origpart <- unclass(fmAll)@where

allpred <- list()
for(i in 1:length(dir(".",pattern="predictionsRRS-*")))
  {
  load(paste("predictionsRRS-",i,".rda",sep=""))
  allpred[[i]] <- predRRS
  }


allpredRRS <- unlist(allpred,recursive=FALSE)

allpred <- list()
for(i in 1:length(dir(".",pattern="predictionsSRS-*")))
  {
  load(paste("predictionsSRS-",i,".rda",sep=""))
  allpred[[i]] <- predRRS
  }

allpredSRS <- unlist(allpred,recursive=FALSE)

stabRRS <- bs_segment_stab(origpart,allpredRRS)
stabSRS <- bs_segment_stab(origpart,allpredSRS)

#get_corresptable(stabRRS)
get_rands(stabRRS)
get_rands(stabSRS)


#average rand RRS
mean(unlist(lapply(get_rands(stabRRS), function(x) x[1])))
#average adj rand RRS
mean(unlist(lapply(get_rands(stabRRS), function(x) x[2])))

#average rand SRS
mean(unlist(lapply(get_rands(stabSRS), function(x) x[1])))
#average adj rand SRS
mean(unlist(lapply(get_rands(stabSRS), function(x) x[2])))


#average jaccard for segment
indmatRRS <- matrix(unlist(get_maxindex(stabRRS)),nrow=15)
meanRRS <- rowSums(indmatRRS)/dim(indmatRRS)[2]
meanRRS

#average jaccard for segment
indmatSRS <- matrix(unlist(get_maxindex(stabSRS)),nrow=15)
meanSRS <- rowSums(indmatSRS)/dim(indmatSRS)[2]
meanSRS

#median jaccard for segment
medRRS <- apply(indmatRRS,1,median)
meanRRS-medRRS


#######################################################################################
##################### MAX INDICES #####################################################
#######################################################################################

############### coinciding = 1
#coinciding RRS - all
coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(all.equal(i,1))))
#coinRRS
mean(colSums(coinRRS)/15)

rrs1 <- rowSums(coinRRS)/200

colSums(coinRRS)/15

#coinciding SRS - all 
coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(all.equal(i,1))))
#coinSRS

srs1 <- rowSums(coinSRS)/200

mean(colSums(coinSRS)/15)

#overall - all
0.5*(mean(colSums(coinSRS)/15)+mean(colSums(coinRRS)/15))


#coinciding RRS - describ
coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(all.equal(i,1))))
coinRRS <- coinRRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinRRS)/7)

coinRRS <- coinRRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinRRS)/5)
#mean(colSums(coinRRS)/7)

coinRRS <- coinRRS[(16-c(7,12)),]
mean(colSums(coinRRS)/2)


#coinciding SRS - describ
coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(all.equal(i,1))))
coinSRS <- coinSRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinSRS)/7)
coinSRS <- coinSRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinSRS)/5)

coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(all.equal(i,1))))
coinSRS <- coinSRS[(16-c(7,12)),]
mean(colSums(coinSRS)/2)

#overall - describ 
0.5*(mean(colSums(coinSRS)/7)+mean(colSums(coinRRS)/7))


############### Strong correspondence = 0.8 
#corresp RRS - all
corrRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>0.8)))
mean(colSums(corrRRS)/15)


rrs08 <- rowSums(corrRRS)/200
mean(rrs08)

#corresp SRS  - all
corrSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>0.8)))
mean(colSums(corrSRS)/15)

srs08 <- rowSums(corrSRS)/200


#overall  - all
0.5*(mean(colSums(corrRRS)/15)+mean(colSums(corrSRS)/15))

#coinciding RRS - describ
coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>.8)))
coinRRS <- coinRRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinRRS)/7)
coinRRS <- coinRRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinRRS)/5)
#mean(colSums(coinRRS)/7)

coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>.8)))
coinRRS <- coinRRS[(16-c(7,12)),]
mean(colSums(coinRRS)/2)

#coinciding SRS - describ
coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>.8)))
coinSRS <- coinSRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinSRS)/7)
coinSRS <- coinSRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinSRS)/5)

coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i>.8)))
coinSRS <- coinSRS[(16-c(7,12)),]
mean(colSums(coinSRS)/2)

#overall - describ 
0.5*(mean(colSums(coinSRS)/7)+mean(colSums(coinRRS)/7))


############### Weak correspondence < .25
#corresp RRS - all
corrRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.25)))
mean(colSums(corrRRS)/15)

rrs025 <- rowSums(corrRRS)/200


#corresp SRS  - all
corrSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.25)))
mean(colSums(corrSRS)/15)

srs025 <- rowSums(corrSRS)/200

#overall  - all
0.5*(mean(colSums(corrRRS)/15)+mean(colSums(corrSRS)/15))

#coinciding RRS - describ
coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.25)))
coinRRS <- coinRRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinRRS)/7)
coinRRS <- coinRRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinRRS)/5)
#mean(colSums(coinRRS)/7)


#coinciding SRS - describ
coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.25)))
coinSRS <- coinSRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinSRS)/7)
coinSRS <- coinSRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinSRS)/5)

############### Weak correspondence < .5
#corresp RRS - all
corrRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.5)))
mean(colSums(corrRRS)/15)

rrs05 <- rowSums(corrRRS)/200


#corresp SRS  - all
corrSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.5)))
mean(colSums(corrSRS)/15)

srs05 <- rowSums(corrSRS)/200


#overall  - all
0.5*(mean(colSums(corrRRS)/15)+mean(colSums(corrSRS)/15))

#coinciding RRS - describ
coinRRS <- sapply(get_corresptable(stabRRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.5)))
coinRRS <- coinRRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinRRS)/7)
coinRRS <- coinRRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinRRS)/5)
#mean(colSums(coinRRS)/7)


#coinciding SRS - describ
coinSRS <- sapply(get_corresptable(stabSRS), function(x) sapply(x[,"Index"],function(i) isTRUE(i<0.5)))
coinSRS <- coinSRS[(16-c(1,2,3,4,5,7,12)),]
mean(colSums(coinSRS)/7)

coinSRS <- coinSRS[(16-c(1,2,3,4,5)),]
mean(colSums(coinSRS)/5)


########################################################################################
############## Summary statistics tables ###############################################
########################################################################################

## Table for RRS
indmatRRS <- matrix(unlist(get_maxindex(stabRRS)),nrow=15)
summary(indmatRRS[1,])

summRRS <- apply(indmatRRS,1,summary)
St.Dev. <- apply(indmatRRS,1,sd)

tabRRS <- rbind(summRRS,rrs1,rrs08,rrs05,rrs025)
tabRRS <- t(tabRRS)
tabRRS <- tabRRS[order(15:1),]
tabRRS <- rbind(tabRRS,colSums(tabRRS)/15)

require(xtable)
xtable(tabRRS)


## Table for SRS
indmatSRS <- matrix(unlist(get_maxindex(stabSRS)),nrow=15)
summary(indmatSRS[1,])

summSRS <- apply(indmatSRS,1,summary)
St.Dev. <- apply(indmatSRS,1,sd)

tabSRS <- rbind(summSRS,srs1,srs08,srs05,srs025)
tabSRS <- t(tabSRS)
tabSRS <- tabSRS[order(15:1),]
tabSRS <- rbind(tabSRS,colSums(tabSRS)/15)
tabSRS

require(xtable)
xtable(tabSRS)


#jaccard boxplot for segments
boxplot.matrix(indmatSRS,use.cols=FALSE)

#####################################################
###### violinplots of the segment-wise jaccard#######
#####################################################
#dark are the ones we describe in the paper

pdf("jaccardviolins.pdf",height=20,width=2,pointsize=8)
require(vioplot)
par(mfrow=c(2,1),mai=c(0,0.6,0.5,0.25))
plot(seq(0,1,length.out=15),1:15,type="n",xlim=c(0,1),ylim=c(1,15),yaxt="n",xaxt="n",ylab="RRS")
axis(2,at=1:15,labels=paste("R",1:15,sep=""))
axis(3,at=seq(0,1,by=0.2))

for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(indmatRRS[i,], at=i, horizontal=TRUE ,col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(indmatRRS[i,], horizontal=TRUE, col="gray90",at=i,add=TRUE,rectCol="white",colMed="black")
 }
points(1,10,pch=19)

par(mai=c(0.5,0.6,0,0.25))
plot(seq(0,1,length.out=15), 1:15,  type="n",yaxt="n",xlim=c(0,1),ylim=c(1,15),ylab="SRS",xlab="")
axis(2,at=1:15,labels=paste("R",1:15,sep=""))
for(i in c(4,9,11,13,14,15))
  {
   vioplot(indmatSRS[i,], at=i, horizontal=TRUE, wex=0.5, col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(indmatSRS[i,], col="gray90",at=i, horizontal=TRUE, add=TRUE,rectCol="white",colMed="black")
 }
points(1,12,pch=19)
dev.off()

pdf("jaccardviolins2.pdf",width=12,height=7)
require(vioplot)
par(mfrow=c(2,1),mai=c(0.5,0.6,0.5,0.25),las=1)
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,1),xlim=c(1,15),yaxt="n",xaxt="n",main="RRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,1,by=0.2))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(indmatRRS[i,], at=i ,col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(indmatRRS[i,], col="gray90",at=i,add=TRUE,rectCol="white",colMed="black")
 }
points(10,1,pch=19)

par(mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15),  type="n",yaxt="n", xaxt="n",ylim=c(0,1),xlim=c(1,15),xlab="",main="SRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,1,by=0.2))
for(i in c(4,9,11,13,14,15))
  {
   vioplot(indmatSRS[i,], at=i, col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(indmatSRS[i,], col="gray90",at=i, add=TRUE,rectCol="white",colMed="black")
 }
points(12,1,pch=19)
dev.off()

###################################################
################# beanplots for the jaccard index
###################################################

require(beanplot)

pdf("beanplot_Jacc2.pdf",width=12,height=7)
par(mfrow=c(1,1),mai=c(0.5,0.6,0.5,0.25),las=1)
bw <- 0.075
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,1),xlim=c(1,15),yaxt="n",xaxt="n",)
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,1,by=0.2))
for(i in c(4,9,11,12,13,14## ,15
           ))
  {
   beanplot(indmatRRS[i,], bw=bw, cutmin=0, cutmax=1, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95)
 }
for(i in c(4,9,11,12,13,14,#15
           ))
  {
   beanplot(indmatSRS[i,], bw=bw, cutmin=0, cutmax=1, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(indmatRRS[i,], bw=bw, cutmin=0, cutmax=1,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0), col = list(c("gray95", "white","white","black")),grownage=50, maxwidth=.95)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(indmatSRS[i,], bw=bw, cutmin=0, what=c(0,1,1,0),cutmax=1,at=i,add=TRUE,side="second",beanlines="median",col = list(c("gray80", "white","white","black")),grownage=50, maxwidth=.95)
 }
    beanplot(indmatRRS[15,], bw=bw, cutmin=0.5, cutmax=1, at=15,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95)
    beanplot(indmatRRS[15,], bw=bw, cutmin=0.5, cutmax=1, at=15,add=TRUE,side="second",what=c(0,1,1,0), col = list(c("gray45", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95)
dev.off()


#This plot is optimized for bandwidth (as good as possible); screw the other one
#also tried to use beanlines="quantile" but sucks too

#pdf("beanplot_Jacc.pdf",width=12,height=7)
postscript("beanplot_Jacc.eps",width=12,height=7)
par(mfrow=c(1,1),mai=c(0.5,0.6,0.5,0.25),las=1)
bw <- 0.035
beanlines <- "median"
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,1),xlim=c(1,15),yaxt="n",xaxt="n",)
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,1,by=0.2))
for(i in c(4,9,11))
  {
   beanplot(indmatRRS[i,], cutmin=0, cutmax=1, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines=beanlines,grownage=50, maxwidth=.95)
 }
for(i in c(12,13,15))
  {
 beanplot(indmatRRS[i,],bw=bw,cutmin=0, cutmax=1, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines=beanlines,grownage=50, maxwidth=.95)
}
 beanplot(indmatRRS[14,],bw=0.02,cutmin=0, cutmax=1, at=14,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines=beanlines,grownage=50, maxwidth=.95) 
for(i in c(4,9,11,13,14)) 
  {
   beanplot(indmatSRS[i,], cutmin=0, cutmax=1, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),beanlines=beanlines,grownage=50, maxwidth=.95)
 }
for(i in c(12,15)) 
  {
   beanplot(indmatSRS[i,], bw=bw,cutmin=0, cutmax=1, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),beanlines=beanlines,grownage=50, maxwidth=.95)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(indmatRRS[i,], cutmin=0, cutmax=1,at=i,add=TRUE,side="first",beanlines=beanlines, what=c(0,1,1,0), col = list(c("gray95", "white","white","black")),grownage=50, maxwidth=.95)
 }
for(i in c(1,2,3,5,6,7,8))
  {
   beanplot(indmatSRS[i,], cutmin=0, what=c(0,1,1,0),cutmax=1,at=i,add=TRUE,side="second",beanlines=beanlines,col = list(c("gray80", "white","white","black")),grownage=50, maxwidth=.95)
 }
   beanplot(indmatSRS[10,], bw=0.01, cutmin=0, what=c(0,1,1,0),cutmax=1,at=10,add=TRUE,side="second",beanlines=beanlines,col = list(c("gray80", "white","white","black")),grownage=50, maxwidth=.95)
dev.off()

pdf("beanplot_JaccHorizontal.pdf",width=12,height=7)
par(mfrow=c(1,1),mai=c(0.5,0.6,0.5,0.25))
bw <- 0.075
plot(1:15, seq(0,1,length.out=15), type="n",xlim=c(0,1),ylim=c(1,15),yaxt="n",xaxt="n",)
axis(2,at=1:15,labels=paste("R",15:1,sep=""))
axis(1,at=seq(0,1,by=0.2))
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(indmatRRS[i,], bw=bw, cutmin=0, cutmax=1, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95, horizontal=TRUE)
 }
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(indmatSRS[i,], bw=bw, cutmin=0, cutmax=1, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),beanlines="median",grownage=50, maxwidth=.95, horizontal=TRUE)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(indmatRRS[i,], bw=bw, cutmin=0, cutmax=1,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0), col = list(c("gray95", "white","white","black")),grownage=50, maxwidth=.95, horizontal=TRUE)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(indmatSRS[i,], bw=bw, cutmin=0, what=c(0,1,1,0),cutmax=1,at=i,add=TRUE,side="second",beanlines="median",col = list(c("gray80", "white","white","black")),grownage=50, maxwidth=.95, horizontal=TRUE)
 }

dev.off()

   


###############################################################################################
########     Plots for segment-wise parameter varaibility     ################################
###############################################################################################
### Paul Hofmarcher, 2012


#### corresponding segments
get_corresp(stabRRS)
get_corresp(stabSRS)

#########################################################
# code snippet for getting the data and the coefficients
######################
seeds <- 1:40
runs <- 1:5

corRRS <- get_corresp(stabRRS)
pos <- paste(rep(seeds,each=5),runs, sep="-")
meanmat <- matrix(NA,nrow=15,ncol=200)
dispmat <- matrix(NA,nrow=15,ncol=200)
for (i in 1:200)#200
  {
    ## tmptree <- 
      load(paste("RRStree",pos[[i]],".rda",sep=""))
    tmptree <- rrs
    segs <- corRRS[[i]]
    meanmat[,i] <- sapply(nodes(tmptree,unlist(segs)),function(x) x$model$coefficients)
    dispmat[,i] <- sapply(nodes(tmptree,unlist(segs)),function(x) x$model$theta)
    rm(tmptree)
    gc()
  }

meanmatRRS <- meanmat
dispmatRRS <- dispmat

save(meanmatRRS, dispmatRRS, file="RRSmean_disp_matrix_validation.rda")

corSRS <- get_corresp(stabSRS)
pos <- paste(rep(seeds,each=5),runs, sep="-")
meanmat <- matrix(NA,nrow=15,ncol=200)
dispmat <- matrix(NA,nrow=15,ncol=200)
for (i in 1:200)#200
  {
    ## tmptree <- 
      load(paste("SRStree",pos[[i]],".rda",sep=""))
    tmptree <- rrs
    segs <- corSRS[[i]]
    meanmat[,i] <- sapply(nodes(tmptree,unlist(segs)),function(x) x$model$coefficients)
    dispmat[,i] <- sapply(nodes(tmptree,unlist(segs)),function(x) x$model$theta)
    rm(tmptree)
    gc()
  }

meanmatSRS <- meanmat
dispmatSRS <- dispmat

save(meanmatSRS, dispmatSRS, file="SRSmean_disp_matrix_validation.rda")

## original values
muk <- c(0.779, -0.399, 0.917,0.904, 0.215, 0.269, 0.114,-0.376,-1.804,-2.979, 0.269, 0.389, -0.329, -0.011, -1.282)

xval <- 15:1


## for exp mean
muke <- exp(muk)
meanmatSRSe <- exp(meanmatSRS)
meanmarRRSe <- exp(meanmatRRS)

pdf("mean.pdf",width=12,height=7)
require(vioplot)
par(mfrow=c(2,1),mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(-4,1.3),xlim=c(1,15),yaxt="n",xaxt="n",main="RRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-4,1,by=1))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(meanmatRRS[i,], at=i ,col="gray60",add=TRUE,rectCol="white",colMed="black")
   #segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], ltw=4, lty=3)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(meanmatRRS[i,], col="gray90",at=i,add=TRUE,rectCol="white",colMed="black")
  # segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], ltw=4, lty=3)
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], lwd=1.5, lty=3)
}

par(mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15),  type="n",yaxt="n", xaxt="n",ylim=c(-4,1.3),xlim=c(1,15),xlab="",main="SRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-4,1,by=1))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(meanmatSRS[i,], at=i, col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(meanmatSRS[i,], col="gray90",at=i, add=TRUE,rectCol="white",colMed="black")
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], lwd=1.5, lty=3)
}
dev.off()





## for exp mean
muke <- exp(muk)
meanmatSRSe <- exp(meanmatSRS)
meanmatRRSe <- exp(meanmatRRS)

pdf("meanexp.pdf",width=12,height=7)
require(vioplot)
par(mfrow=c(2,1),mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,4),xlim=c(1,15),yaxt="n",xaxt="n",main="RRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,4,by=1))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(meanmatRRSe[i,], at=i ,col="gray60",add=TRUE,rectCol="white",colMed="black")
   #segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], ltw=4, lty=3)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(meanmatRRSe[i,], col="gray90",at=i,add=TRUE,rectCol="white",colMed="black")
  # segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], ltw=4, lty=3)
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muke[i], y1=muke[i], lwd=1.5, lty=3)
}

par(mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15),  type="n",yaxt="n", xaxt="n",ylim=c(0,4),xlim=c(1,15),xlab="",main="SRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,4,by=1))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(meanmatSRSe[i,], at=i, col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(meanmatSRSe[i,], col="gray90",at=i, add=TRUE,rectCol="white",colMed="black")
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muke[i], y1=muke[i], lwd=1.5, lty=3)
}
dev.off()



thetak <- c(.089, .069, .096, .386, .468, .128, .275, .090, .043, .009, .205, .373, .117, .199, .047)
pdf("disp.pdf",width=12,height=7)
require(vioplot)
par(mfrow=c(2,1),mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,0.8),xlim=c(1,15),yaxt="n",xaxt="n",main="RRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-3,1,by=0.2))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(dispmatRRS[i,], at=i ,col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(dispmatRRS[i,], col="gray90",at=i,add=TRUE,rectCol="white",colMed="black")
 }


for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=thetak[i], y1=thetak[i], lwd=1.5, lty=3)
}

par(mai=c(0.5,0.6,0.5,0.25))
plot(1:15, seq(0,1,length.out=15),  type="n",yaxt="n", xaxt="n",ylim=c(0,0.8),xlim=c(1,15),xlab="",main="SRS")
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-3,1,by=0.2))
for(i in c(4,9,11,12,13,14,15))
  {
   vioplot(dispmatSRS[i,], at=i, col="gray60",add=TRUE,rectCol="white",colMed="black")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   vioplot(dispmatSRS[i,], col="gray90",at=i, add=TRUE,rectCol="white",colMed="black")
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=thetak[i], y1=thetak[i], lwd=1.5, lty=3)
}
dev.off()





################# beanplot

require(beanplot)

pdf("beanplotmean_ohnebw.pdf",width=12,height=7)
par(mfrow=c(1,1),mai=c(0.5,0.6,0.5,0.25))
bw <- 0.075
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(-4,1.5),xlim=c(1,15),yaxt="n",xaxt="n",)
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-4,1,by=1))
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(meanmatRRS[i,], cutmin=-4, cutmax=2, at=i,add=TRUE,side="first",what=c(0,1,1,0),grownage=50, maxwidth = 0.95,col = list(c("gray60", "white","white","black")),beanlines="median")
 
 }
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(meanmatSRS[i,], cutmin=-4, cutmax=2, at=i,add=TRUE,side="second",what=c(0,1,1,0),grownage=50,maxwidth = 0.95,col = list(c("gray45", "white","white","black")),beanlines="median")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(meanmatRRS[i,], cutmin=-4, cutmax=2,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0),grownage=50,maxwidth = 0.95, col = list(c("gray95", "white","white","black")))
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(meanmatSRS[i,], cutmin=-4, what=c(0,1,1,0),cutmax=2,at=i,add=TRUE,side="second",beanlines="median",grownage=50,maxwidth = 0.95,col = list(c("gray80", "white","white","black")))
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], lwd=1.5, lty=3)
}
##  legendbean <- seq(from=-4, to=-2.5, by=0.01)
## subplot(beanplot(legendbean, what=c(0,1,0,0), col=list(c("gray95", "white","white","black")),axis(1, at=seq(14, 15, by=1), labels = FALSE),axis(2, at=seq(-4, -3, by=1), labels = FALSE)), type="plt",x=15, y=-3.4)
#mtext(x=c(14, -3.4), "RRS")
## subplot(beanplot(legendbean, what=c(0,1,0,0), col=list(c("gray95", "white","white","black")), axes=FALSE, side="second"), type="plt",x=15, y=-3.4)

dev.off()


###########################################
pdf("beanplotdisp_ohnebw.pdf",width=12,height=4)
par(mfrow=c(1,1),mai=c(0.5,0.6,0.5,0.25))
bw <- 0.075
plot(1:15, seq(0,1,length.out=15), type="n",ylim=c(0,0.7),xlim=c(1,15),yaxt="n",xaxt="n",)
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(0,1,by=0.2))
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(dispmatRRS[i,], cutmin=0, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),grownage=50,maxwidth=0.95,beanlines="median", cutmax=0.7)
 }
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(dispmatSRS[i,], cutmin=0, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),grownage=50,maxwidth=0.95,beanlines="median", cutmax=0.7)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(dispmatRRS[i,], cutmin=0,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0), col = list(c("gray95", "white","white","black")),grownage=50,maxwidth=0.95, cutmax=0.7)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(dispmatSRS[i,], cutmin=0, what=c(0,1,1,0),at=i,add=TRUE,side="second",beanlines="median",col = list(c("gray80", "white","white","black")),grownage=50,maxwidth=0.95, cutmax=0.7)
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=thetak[i], y1=thetak[i], lwd=1.5, lty=3)
}
dev.off()

############################# two beanplots in one pdf#############################
load("RRSmean_disp_matrix_validation.rda")
load("SRSmean_disp_matrix_validation.rda")

muk <- c(0.779, -0.399, 0.917,0.904, 0.215, 0.269, 0.114,-0.376,-1.804,-2.979, 0.269, 0.389, -0.329, -0.011, -1.282)
thetak <- c(.089, .069, .096, .386, .468, .128, .275, .090, .043, .009, .205, .373, .117, .199, .047)
xval <- 15:1

require(beanplot)
#pdf("bean.pdf",width=12,height=7)
postscript("bean.eps",width=12,height=7)
 nf <- layout(matrix(c(1,1,2,2), 2, 2, byrow=TRUE), respect=TRUE, heights=c(1,0.5))
#split.screen(c(2,1))
par(mai=c(0,0.5,0.0,0.1), las=1)
#for pdf par(## mfrow=c(2,1),
    mai=c(0,0.1,0.0,0.1), las=1)
bw <- 0.075
#screen(1)
plot(1:15, ## seq(0,1,length.out=15)
      type="n",ylim=c(-4,1.5),xlim=c(1,15),yaxt="n",xaxt="n",ylab="" #, ylab=expression(hat(mu)[k]
     )
mtext(side=2, text=expression(log(hat(mu)[k])), line=3) 
#axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(2,at=seq(-4,1,by=1))

for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(meanmatRRS[i,], cutmin=-4, cutmax=2, at=i,add=TRUE,side="first",what=c(0,1,1,0),grownage=50, maxwidth = 0.95,col = list(c("gray60", "white","white","black")),beanlines="median")
 
 }
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(meanmatSRS[i,], cutmin=-4, cutmax=2, at=i,add=TRUE,side="second",what=c(0,1,1,0),grownage=50,maxwidth = 0.95,col = list(c("gray45", "white","white","black")),beanlines="median")
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(meanmatRRS[i,], cutmin=-4, cutmax=2,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0),grownage=50,maxwidth = 0.95, col = list(c("gray95", "white","white","black")))
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(meanmatSRS[i,], cutmin=-4, what=c(0,1,1,0),cutmax=2,at=i,add=TRUE,side="second",beanlines="median",grownage=50,maxwidth = 0.95,col = list(c("gray80", "white","white","black")))
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=muk[i], y1=muk[i], lwd=1.5, lty=3)
}

#title(ylab=expression("Estimates of"~mu))
#screen(2)

par(mai=c(0.35,0.5,0,0.1), las=1)
#for pdf par(mai=c(0.35,0.1,0,0.1), las=1)
plot(1:15, ## seq(0,1,length.out=15),
     type="n",ylim=c(0,0.7),xlim=c(1,15),yaxt="n",xaxt="n",ylab="" ## ylab=expression(hat(theta)[k])
     )
axis(1,at=1:15,labels=paste("R",15:1,sep=""))
axis(4,at=seq(0,1,by=0.2))
mtext(4, text=expression(hat(theta)[k]), line=3) 
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(dispmatRRS[i,], cutmin=0, at=i,add=TRUE,side="first",what=c(0,1,1,0), col = list(c("gray60", "white","white","black")),grownage=50,maxwidth=0.95,beanlines="median", cutmax=0.7)
 }
for(i in c(4,9,11,12,13,14,15))
  {
   beanplot(dispmatSRS[i,], cutmin=0, at=i,add=TRUE,side="second",what=c(0,1,1,0),col = list(c("gray45", "white","white","black")),grownage=50,maxwidth=0.95,beanlines="median", cutmax=0.7)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(dispmatRRS[i,], cutmin=0,at=i,add=TRUE,side="first",beanlines="median", what=c(0,1,1,0), col = list(c("gray95", "white","white","black")),grownage=50,maxwidth=0.95, cutmax=0.7)
 }
for(i in c(1,2,3,5,6,7,8,10))
  {
   beanplot(dispmatSRS[i,], cutmin=0, what=c(0,1,1,0),at=i,add=TRUE,side="second",beanlines="median",col = list(c("gray80", "white","white","black")),grownage=50,maxwidth=0.95, cutmax=0.7)
 }

for(i in 1:15){
segments(x0=xval[i]-0.5,x1=xval[i]+0.5, y0=thetak[i], y1=thetak[i], lwd=1.5, lty=3)
}
#close.screen(all = TRUE)
dev.off()
