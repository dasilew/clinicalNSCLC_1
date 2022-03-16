#07.09.2021
#updated 19.09.2021
#updated 02.10.2021
#updated 16.01.2022: only NSCLC
#updated 28.01.2022: only NSCLC
#updated 31.01.2022
#updated 03.02.2022: do it only for the two groups of interest
#updated 05.03.2022
# Laod required libraries
install.packages(c("survival", "survminer"))
library(gtsummary)
library(dplyr)
library(MatchIt)
library(Rcpp)
library(survival)
library(survminer)
library(gapminder)
library(dplyr)
library(tidyverse)
library(readxl)
library(swimplot)
glimpse(lung)
glimpse(NSCLC22)

NSCLC23 <- read_excel("~/Desktop/NSCLC22.xlsx")
View(NSCLC23)
#exploratory data analysis
head(NSCLC23)
tail(NSCLC23)
glimpse(NSCLC23)
length(levels(factor(NSCLC23$patient_ID)))

NSCLC22

my_table2 <- NSCLC22 %>% select(GPA2, no_of_BrMs, localization_of_BrMs, primary_tumor.Rx2, extracran.met, naive, adj.RTx3, adj.tx2,
                                censoring, survival)



my_table2$GPA2<- factor(my_table2$GPA2,       
                        levels = c("bad", "good"), 
                        labels = c("GPA below 2,5", "GPA at least 2,5"))


my_table2$no_of_BrMs <- factor(my_table2$no_of_BrMs, 
                               levels = c("0", "1", "2"), 
                               labels = c("one", "two", "more than two"))


my_table2$extracran.met <- factor(my_table2$extracran.met,       
                                  levels = c("0", "1"), 
                                  labels = c("No", "Yes"))


my_table2$localization_of_BrMs <- factor(my_table2$localization_of_BrMs, 
                                         levels = c("0", "1", "2"), 
                                         labels = c("supratentorial", "infratentorial", "both"))


my_table2$primary_tumor.Rx2 <- factor(my_table2$primary_tumor.Rx2, 
                                      levels = c("0", "1"), 
                                      labels = c("no primary tumor resection", 
                                                 "primary tumor resection"))

my_table2$naive <- factor(my_table2$naive, 
                          levels = c("0", "1"), 
                          labels = c("no pre-treatment", "pre-treatment"))


my_table2$adj.RTx3 <- factor(my_table2$adj.RTx3, 
                             levels = c("0", "1"), 
                             labels = c("other than SRS", "SRS"))


my_table2$adj.tx2 <- factor(my_table2$adj.tx2 , 
                            levels = c("0", "1"), 
                            labels = c("RTx + CTx", "RTx + ICI"))



# Load example dataset, modified version of survival::colon
colon_s3 <- my_table2[, c("survival", "censoring", "GPA2", "localization_of_BrMs", "no_of_BrMs", "extracran.met",
                          "primary_tumor.Rx2", "naive", "adj.RTx3", "adj.tx2")]

# Convert categorial variables into factor variables to avoid NAs in the output table
colon_s3$localization_of_BrMs <- as.factor(colon_s3$localization_of_BrMs)

colon_s3$GPA2 <- as.factor(colon_s3$GPA2)
colon_s3$no_of_BrMs <- as.factor(colon_s3$no_of_BrMs)
colon_s3$extracran.met <- as.factor(colon_s3$extracran.met)
#colon_s3$volume <- as.factor(colon_s3$volume)
colon_s3$primary_tumor.Rx2 <- as.factor(colon_s3$primary_tumor.Rx2)
colon_s3$naive <- as.factor(colon_s3$naive)
colon_s3$adj.RTx3 <- as.factor(colon_s3$adj.RTx3)
colon_s3$adj.tx2 <- as.factor(colon_s3$adj.tx2)


colon_s3 

#my_table2 <- NSCLC22 %>% select(no_of_BrMs, GPA2, localization_of_BrMs, primary_tumor.Rx2, extracran.met, adj.RTx3, adj.tx2,
 #                               censoring, survival)


res.cox <- coxph(Surv(survival, censoring) ~ GPA2 + no_of_BrMs + localization_of_BrMs + primary_tumor.Rx2 + extracran.met + adj.RTx3 + adj.tx2, data =  NSCLC22)
res.cox


#older script with different covariates:
#res.cox <- coxph(Surv(survival, censoring) ~ age + gender + no_of_BrMs + volume + extracran.met, data =  NSCLC23)
#res.cox

#adj.tx2, gender, age, KPS_group, extracran.met,
#no_of_BrMs, volume)


#
#res.cox <- coxph(Surv(survival, censoring) ~ age + gender + no_of_BrMs + no_of_BrM_Rxs2 + volume + primary_tumor.Rx2 + 
              #     extracran.met + UICC2 + GPA + KPS + localization_of_BrMs, NLR, data =  NSCLC22

#The proportional hazards (PH) assumption can be checked using statistical tests and graphical diagnostics based on the scaled Schoenfeld residuals.
#The function cox.zph() [in the survival package] provides a convenient solution to test the proportional hazards assumption for each covariate included in a Cox refression model fit.
#For each covariate, the function cox.zph() correlates the corresponding set of scaled Schoenfeld residuals with time, to test for independence between residuals and time. 
#Additionally, it performs a global test for the model as a whole.

#To test for the proportional-hazards (PH) assumption, type this:
  
  test.ph <- cox.zph(res.cox)
test.ph
p2 <- ggcoxzph(test.ph)

ggpar(p2, 
      font.title = c(6),
      font.subtitle = c(10,  "orange"),
      font.caption = c(10,  "orange"),
      font.x = c(8,  "blue"),
      font.y = c(8,  "#993333")
)








