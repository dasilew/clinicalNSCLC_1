#calculating follow-up time (http://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Advanced-statistics/SurvivalAnalysis.html)
library(prodlim)
library(survival)
library(Publish)
quantile(prodlim(Hist(survival,censoring)~1,data=NSCLC22,reverse=TRUE))


quantile(prodlim(Hist(survival,censoring)~1,data=NSCLC23,reverse=TRUE))
#to count number of death events
NSCLC22$censoring
count(NSCLC22$censoring)
table(NSCLC22$censoring)
sum <- 89 + 295

a = 89/384
a
b = 295/384
b
