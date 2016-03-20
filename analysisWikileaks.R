#!/usr/bin/env Rscript
#############################################################################################
###### Pre-processing Afghanistan war log with topic models
#################################################################
#P. Hofmarcher, 2011

#require("tm")
#require("topicmodels")
require("NLP",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("tm",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("topicmodels",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("SnowballC",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("slam",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("mvtnorm",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")
require("party",lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2/")


readAFG <- function(elem, language, id)
    PlainTextDocument(elem$content[7],
                      datetimestamp = elem$content[2],
                      id = id, language = language)

## prepare data
afg <- read.csv("afg.csv", stringsAsFactors = FALSE, header = FALSE)

## meta data == data without report summaries
meta_data <- afg
colnames(meta_data)<-c("reportkey", "date", "type", "category", "tracking_n", "title", "summary", "region", "attackon", "complex", "reportingu", "unitname", "typeofunit", "FriendlyWo", "FriendlyKi", "HostNatWou", "HostNatKil", "CivWounded", "CivKill", "EnemyWound", "EnemyKille", "EnemyDetai", "MGRS", "lat", "long", "Originator", "UpdatedByG", "CCIR", "sigact", "affiliatio", "dcolor", "classification")
meta_data$date <- as.Date(meta_data$date)

## Document Term Matrix
print(date())
print("Phase 1 - start")
afg_corp <- Corpus(DataframeSource(afg), readerControl = list(reader =readAFG))  #takes some time
## JLP
options(mc.cores=1)
##
afg_DTM <- DocumentTermMatrix(afg_corp, control = list(removePunctuation = TRUE, removeNumbers = TRUE, stemming = TRUE, stopwords = TRUE, minWordLength= 2))  #takes a LOT of time
index <- which(row_sums(afg_DTM)==0)      
#removing those without a report summary
afg_DTM <- afg_DTM[-index,]
meta_data<- meta_data[-index,]
print("Phase 1 - complete")
print(date())
print("Phase 2 - start")
## removing report out of DTM. Each Report summary contains the term "report"
afg_DTM <- afg_DTM[,-c(which(colnames(afg_DTM)=="report"))]

## LDA with k=100 topics (can be computationally very challenging)
afg_LDA <- LDA(afg_DTM, k=100, control=list(alpha=0.001, estimate.alpha=FALSE, verbose=4))      

#returns the topics
topic <- topics(afg_LDA)

#returns the 30 most frequent terms
terms_afg <- terms(afg_LDA, 30)

usedData <- cbind(meta_data,topic)
save.image(file="origData_cbind_topics.rda")
print("Phase 2 - complete")
print(date())

######################################################
###### Data cleaning
######################################################
#T. Rusch, 2011

load("origData_cbind_topics.rda")
newnam <- names(usedData)[-7] #remove the text summaries, as they are to big

#some new names 
newnam[c(27,28)] <- c("ccir","sigact")
newnam[22] <- "mgrs"
names(usedData) <- newnam

#change cat vars to factors
factorizers <- c(1,3,4,5,6,7,8,9,10,11,12,22,25,26,27,28,29,30,31,32)
for (i in factorizers) usedData[,i] <- factor(usedData[,i])

#create dummy matrix for the topics
m1 <- model.matrix(~-1+factor(usedData$topic))
summary(m1)

#name the dummy matrix
n1 <- rep(NA,100)
for(i in 1:100) n1[i] <- paste("t",i,sep="")
dimnames(m1)[[2]] <- n1
m1 <- as.data.frame(m1)
for(i in 1:100) m1[,i] <- factor(m1[,i])
summary(m1)

afg <- cbind(usedData,m1)
save(afg,file="afgClean.RData")
###################################################
####Analysis of the Wikileaks Afghanistan war logs
######################################################
#T. Rusch, 2011

#load data
load("afgClean.RData")

#source the function that are not in party and written by us
#contains a negative binomial S4 statsmodel  and a fitting function based on glm.nb()
source("fitfunctionsWikileaks.R") 

#setting up hyperparameters
control.me <- mob_control(alpha=1e-04,minsplit=300,verbose=TRUE)

#remove the documents with no entries in the fatality columns
afg1 <- afg[!is.na(afg$kiaAll),]
system.time(fmAll <- mob(kiaAll ~ 1 |attackOn+affiliation+dcolor+region+complexAttack+t1+t2+t3+t4+Topic_5+t6+t7+t8+t9+t10+t11+t12+t13+Topic_14+t15+t16+t17+Topic_18+Topic_19+t20+t21+t22+t23+t24+t25+t26+Topic_27+t28+t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+Topic_61+t62+t63+t64+t65+t66+t67+t68+t69+t70+Topic_71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+Topic_85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100, data = afg1, model = negbinModel, na.action=NULL, control=control.me))
summary(fmAll)

#source the plot functions that are not in party and written by us
#contains a binary tree with a mean-standard deviation terminal function plot
source("plotfunctionsWikileaks.R")

#plotting
plot.BinaryTree(fmAll,terminal_panel=node_meansddevplot3)

