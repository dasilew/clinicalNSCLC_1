#07.09.2021
#updated 19.09.2021
#updated 02.10.2021
#updated 16.01.2022: only NSCLC
#updated 28.01.2022: only NSCLC
#updated 31.01.2022
#updated 03.02.2022: do it only for the two groups of interest
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
res.cox <- coxph(Surv(survival, censoring) ~ age + gender + no_of_BrMs + volume + extracran.met, data =  NSCLC23)
res.cox

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











# Step 1: Select the data apply preprocessing steps 
#-------------------------------------------------------------------------------

# Only two possible values for the dependent grouping variable are possible 
# (otherwise the logistic regression would not work)
# Therefore, select only patients with adj.tx2 in (0, 1) for the analysis
# Another option would be to merge values 1 and 2 into one category
propensity.master <- subset(NSCLC22, adj.tx2 %in% c(0, 1))

# Recoding variables and converting into factor variables
propensity.master$gender <- factor(propensity.master$gender)
propensity.master$no_of_BrMs <- factor((propensity.master$no_of_BrMs > 1) * 1, labels = c("one", "more than one"))
propensity.master$no_of_BrM_Rxs2 <- factor(propensity.master$no_of_BrM_Rxs2, labels = c("one", "two"))
propensity.master$primary_tumor.Rx2 <- factor(propensity.master$primary_tumor.Rx2, labels = c("no primary tumor resection", "primary tumor resection"))
propensity.master$pre.tx <- factor((propensity.master$pre.tx >= 1) * 1, labels = c("no pre-treatment", "pre-treatment"))
propensity.master$adj.RTx2 <- factor(propensity.master$adj.RTx2, labels = c("not known", "WBRT", "local radiation of resection cavity", "SRS"))
propensity.master$dose_of_RTx <- factor((propensity.master$dose_of_RTx > 30) * 1, labels = c("< 30 Gy radiation dose", ">/= 30 Gy radiation dose"))
propensity.master$adj.tx2 <- factor(propensity.master$adj.tx2, labels = c("RTx + CTx", "RTx + ICI"))
propensity.master$KPS_group <- factor((propensity.master$KPS > 70) * 1, labels = c("< 70%", ">/= 70%"))
propensity.master$UICC2 <- factor((propensity.master$UICC2 >= 1) * 1, labels = c("other than IV", "IV"))
propensity.master$extracran.met <- factor(propensity.master$extracran.met, labels = c("No", "Yes"))
#propensity.master$localization_of_BrMs <- factor(propensity.master$localization_of_BrMs, labels = c("supratentorial", "infratentorial", "both"))
#propensity.master$histology2 <- factor(propensity.master$histology2, labels = c("ADC", "other than ADC"))

#not included into matching only into the summary table
propensity.master$PDL1_extracranial <- factor((propensity.master$PDL1_extracranial > 1) * 1, labels = c("< 1%", ">/= 1%"))
propensity.master$PDL1_intracranial <- factor((propensity.master$PDL1_intracranial > 1) * 1, labels = c("< 1%", ">/= 1%"))
propensity.master$Ki67_extracranial <- factor((propensity.master$Ki67_extracranial > 30) * 1, labels = c("< 1%", ">/= 1%"))
propensity.master$Ki67_intracranial <- factor((propensity.master$Ki67_intracranial > 30) * 1, labels = c("< 1%", ">/= 1%"))
propensity.master$NLR <- factor((propensity.master$NLR > 5) * 1, labels = c("< 5", ">/= 5"))
propensity.master$volume <- factor((propensity.master$volume > 5) * 1, labels = c("< 15", ">/= 15"))

#remove GPA from PSM

# Select only non-missing values since the matching algorithm does not allow missing values
m.out0 <- propensity.master  %>%
  select(gender, age, no_of_BrMs, no_of_BrM_Rxs2, primary_tumor.Rx2, pre.tx, adj.RTx2, 
         dose_of_RTx, adj.tx2, KPS_group, UICC2, extracran.met, survival, censoring, volume) %>%
  na.omit()
summary(m.out0)


# No matching; constructing a pre-match matchit object
m.out0 <- matchit(adj.tx2 ~ gender + age +  no_of_BrMs + primary_tumor.Rx2 + pre.tx + adj.RTx2 + 
                    dose_of_RTx + KPS_group + UICC2 + extracran.met + volume + censoring + survival, data = m.out0,
                  method = NULL, distance = "glm") %>% 
  na.omit()

# Checking balance prior to matching
summary(m.out0)


library(Hmisc)
library(xlsx)

write.xlsx(as.data.frame(m.out0$r), file = "results.xlsx")



# Step 2: Perform univariate statistics for each covariate
#-------------------------------------------------------------------------------

# Descriptive summary statistics
summary.tab <- propensity.master %>% select(gender, age, no_of_BrMs, no_of_BrM_Rxs2, primary_tumor.Rx2, 
                                            pre.tx, adj.RTx2, dose_of_RTx, adj.tx2, KPS_group, UICC2, 
                                            extracran.met, volume, PDL1_intracranial, Ki67_intracranial, 
                                            PDL1_extracranial, Ki67_extracranial, NLR, TTF1_status_BrM)

#add other variables here if needed...




summary.tab %>% tbl_summary(
  # Group table by adj.tx2
  by = adj.tx2) %>%
  bold_labels() %>% 
  # Add p values in last column (depending on the scale level of the data)
  add_p() %>%
  # Add the descriptive statistics for the total sample in the first data column
  add_overall()  %>%
  # Add label for the Grouping variable in the table
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment**")

# Remarks: 
# - The covariates PDL1i, Ki67i and NLR should not be taken into account in the propensity score matching because of high portion of missings
# - histology: Currently, there are many "sparse" categories. Either remove or merge categories together



# Step 3: Conduct propensity score matching with adj.tx2 as dependent variable
#-------------------------------------------------------------------------------

# Select only non-missing values since the matching algorithm does not allow missing values
match.data <- propensity.master  %>%
  select(gender, age, no_of_BrMs, volume, 
         adj.tx2, KPS_group,  
         extracran.met, volume, censoring, survival
  ) %>%
  na.omit()

#addendum 15.09.21:
#remove histology, or put together ADC vs. ohter; TTF1 status
#according to UICC IV vs. other (= non IV)
#match for pre-treatment (no pretx vs. pretx)

#addendum 16.11.2022: 
#all are NSCLC
#histology and entity are removed

#addendum january 2022, addenddum from 28.01.2022 remove NSCLC!
#31.01.2022: changed match.data

# Set seed value to avoid differing results for different runs
set.seed(1234)

# Propensity Score Matching algorithm
# Input parameters of the function matchit()
# - formula: is the model formula with grouping variable as dependent variable followed by all covariates connected by "+"
# - distance: the setting of "glm" ensures that a logistic regression is running in the background to derive the propensity score
#   The prop. score is the probability of each patient of getting selected into group "adjuvant immune therapy" based on the information of the covariates
# - method: "nearest" means that the matching is based on the nearest neighbor method. 
#   This means basically, that for each patient in group "adjuvant immune therapy" will matched with a patient in group "no adjuvant immune therapy" 
#   that is most similar. Statistically speaking, this means that for each patient in group "adjuvant immune therapy", 
#   the patient in group "no adjuvant immune therapy" is chosen where the propensity score is closest to the patient in group "adjuvant immune therapy".
prop.matching <- matchit(adj.tx2 ~ gender + age + no_of_BrMs + volume + KPS_group + extracran.met + survival + censoring
                    data = match.data, method = "nearest", distance = "glm")


#summary
summary(prop.matching)
plot(summary(prop.matching, interactions = FALSE),
     var.order = "unmatched")
#19.09.2021
#31.01.2022
#Plot balance creating a Love plot of the covariates.
#A Love plot is a clean way to visually summarize balance. Using plot on the output of a call to summary(),
##on a matchit object produces a Love plot of the standardized mean differences. 
#plot.summary.matchit() has several additional arguments that can be used to customize the plot.
#plot(summary(prop.matching))

# The output shows (amont others) two important statistics: 
# 1. An overview of means for both groups ("Summary of Balance for All Data")
#   - "Means Treated" is the mean value for each covariate for the smaller group ("adjuvant immune therapy")
#   - "Means Control" is the mean value for each covariate for the smaller group ("no adjuvant immune therapy") AFTER matching
#   - Example for age: mean of treated is 62.833 and mean of controll is now 63.45. Before matchting, the mean for controll was (slightly) higher
# 2. A general overview at the end of the output how man observatios are control and treated group and how many were matched

# Get dataset with only the matched observations 
# This dataset will contain than only the matched patients, all other unmatched patients in group "no adjuvant immune therapy" are removed
match.data.final <- match.data(prop.matching)
match.data.final
#test data set
test <- subset(match.data.final, subclass==1)

# Compare dimension of original dataset (match.data) and final matched dataset (match.data.final)
dim(match.data)
dim(match.data.final)

# The dataset match.data.final contains three new variables
# - distance: This is the propensity score, i.e. the outcome of the underlying logistic regression
#   Based on this distance, patients of both groups are matched
# - weights: This variable shows the weighting of patients, since no weighting is used, the value is 1 always
# - subclass: This variable shows the matched pairs, i.e. in each of the 60 subclasses, there are two patients. 
#   One patient is from group "adjuvant immune therapy", the other is from group "no adjuvant immune therapy"
#   When checking the distance between the two patients in one subclass, the values should be very close.
#   Example for subclass 1, see below: distance values are 0.282 and 0.285
match.data.final[match.data.final$subclass == 1, "distance"]

plot(prop.matching, type = "jitter", interactive = FALSE)
plot(prop.matching, type = "qq", interactive = FALSE,
     which.xs = c("KPS_group"))
#take different covariates
plot(prop.matching, type = "qq", interactive = FALSE,
     which.xs = c("KPS_group", "age"))
plot(prop.matching, type = "qq", interactive = FALSE,
     which.xs = c("volume"))

# Re-run descriptive statistics
summary.tab.prop <- match.data.final %>% select(adj.tx2, gender, age, KPS_group, extracran.met,
                                                no_of_BrMs, volume)

summary.tab.prop %>% tbl_summary(
  # Group table by adj.tx2
  by = adj.tx2) %>%
  bold_labels() %>% 
  # Add p values in last column (depending on the scale level of the data)
  add_p() %>%
  # Add the descriptive statistics for the total sample in the first data column
  add_overall()  %>%
  # Add label for the Grouping variable in the table
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment**")

# The means and percentage values for both groups are now much closer than before the matching. 
# All p-values are now > 0.05, i.e. there is no significant differences between both groups compared to all covariates. 
# Regarding the covariates, both groups are now "equal"



# Step 4: Re-perform survivial analysis with matched dataset
#-------------------------------------------------------------------------------

surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match <- survfit(surv_object ~ adj.tx2, data = match.data.final)
ggsurvplot(fit.match,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("RT + CT", "RT + ICI"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match

#trying to get confidence intervals for survival for different treatment cohorts
fit.match

#addition from 19.09.2021; KMS stratified according to covariates within the matched groups





#I stopped here





#Fig2a extracranial mets
surv_object <- Surv(time = match.data.final$survival_months, event = match.data.final$censoring_OS)
fit.match2 <- survfit(surv_object ~ extracran.met, data = match.data.final)
ggsurvplot(fit.match2,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("without exracranial metastases", "extracranial metastases"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match2

#Fig2b localization
surv_object <- Surv(time = match.data.final$survival_months, event = match.data.final$censoring_OS)
fit.match3 <- survfit(surv_object ~ localization_of_BrMs, data = match.data.final)
ggsurvplot(fit.match3,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("supratentorial", "infratentorial", "both"),
           pval = TRUE, pval.coord = c(20, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,48),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match3

#Fig2c GPA_group
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match4 <- survfit(surv_object ~ GPA_group, data = match.data.final)
ggsurvplot(fit.match4,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("GPA >/=2", "GPA < 2"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match4

#Fig2d primary tumor resection
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match5 <- survfit(surv_object ~ primary_tumor.Rx_2, data = match.data.final)
ggsurvplot(fit.match5,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("no PTRx","PTRx before or after BrMRx"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match5

#Fig2f nr.brm resections (no_of_BrM_Rxs2)
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match6 <- survfit(surv_object ~ no_of_BrM_Rxs2, data = match.data.final)
ggsurvplot(fit.match6,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("1 brain met resect", ">/= 2 brain met resect"),
           pval = TRUE, pval.coord = c(50, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match6

#Fig2g no_of_BrMs
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match7 <- survfit(surv_object ~ no_of_BrMs, data = match.data.final)
ggsurvplot(fit.match7,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("1 brain metastasis", "2 or more brain metastases"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match7


#Fig2h UICC
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match8 <- survfit(surv_object ~ UICC2, data = match.data.final)
ggsurvplot(fit.match8,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("other than IV", "IV"),
           pval = TRUE, pval.coord = c(54, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match8

#Fig2i 
surv_object <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
fit.match8 <- survfit(surv_object ~ UICC2, data = match.data.final)
ggsurvplot(fit.match8,
           #title    = "Survival Analysis with matched dataset",
           data = match.data.final,
           legend.labs = c("IV", "other than IV"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)
#confidence interval
fit.match8









match.data.final
typeof(match.data.final)
#You can view each component’s structure using the str() function.
str(match.data.final)
###############################
#using base R: 
#see also: https://dabblingwithdata.wordpress.com/2018/01/02/my-favourite-r-package-for-summarising-data/
#summary(match.data.final)
#new try
by(match.data.final, match.data.final$adj.tx3, summary)


match.data.final_selected <- match.data.final %>%
  select(adj.tx3, survival, censoring)
match.data.final_selected

#new try
match.data.final
surv_objectxxx <- Surv(time = match.data.final$survival, event = match.data.final$censoring)
surv_objectxxx

fitxxx <- survfit(surv_objectxxx ~ adj.tx3, data = match.data.final)
summary(fitxxx)
fitxxx
#now for cox

fitcoxph <- coxph(surv_objectxxx ~ adj.tx3, data = match.data.final)
summary(fitcoxph)
ggforest(fitcoxph, data = match.data.final)

library(dplyr)
match.data.final %>%
  group_by(adj.tx3) %>%
  summarize(median = median(survival))

library(tidyverse)

install.packages("srvyr")
# or for development version
# devtools::install_github("gergness/srvyr")
library(srvyr, warn.conflicts = FALSE)
data(api, package = "survey")

t.test(Bwt ~ Sex, data=cats, conf.level=.99)

######################################################
#NOT MATCHED DATA:

surv_object <- Surv(time = match.data$survival, event = match.data$censoring)
fit.NOTmatch <- survfit(surv_object ~ adj.tx3, data = match.data)
ggsurvplot(fit.NOTmatch,
           #title    = "Survival Analysis with unmatched dataset",
           data = match.data,
           legend.labs = c("RT + CT", "RT + ICI"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",
)

fit.NOTmatch
#https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R7_LogisticRegression-Survival/R7_LogisticRegression-Survival4.html

match.data

# Useful links for propensity score matching in R: 
# https://sejdemyr.github.io/r-tutorials/statistics/tutorial8.html
# https://datascienceplus.com/how-to-use-r-for-matching-samples-propensity-score/
# https://stackoverflow.com/questions/29672088/speedup-matchit


install.packages("xlsx")
library("xlsx")
# Write the first data set in a new workbook
write.xlsx(match.data, file = "myworkbook.xlsx",
           sheetName = "NOT MATCHED DATA", append = FALSE)
# Add a second data set in a new worksheet
write.xlsx(match.data.final, file = "myworkbook.xlsx", 
           sheetName="MATCHED DATA", append=TRUE)



#02.10.2021
#perform uni- and multivariable analysis
#1. Summarise variables/factors by a categorical variable
library(finalfit)
library(dplyr)
match.data.final
#do not use the unmatched dataset (summary.tab) but use as for the KMS analysis the matched dataset: data = match.data.final
# Table 1 - Patient demographics by variable of interest ----
colon_s3 <- match.data.final[c("gender", "age", "KPS_group", "GPA_group", "extracran.met", "no_of_BrMs", "no_of_BrM_Rxs2", "localization_of_BrMs", 
                        "UICC2", "primary_tumor.Rx_2", "pre.tx2", "adj.tx3", "survival", "censoring")]
colon_s3


explanatory = c("extracran.met", "no_of_BrMs", "no_of_BrMs",
                 "UICC2", "pre.tx2", "adj.tx3", "primary_tumor.Rx_2")

dependent = "Surv(survival, censoring)"
colon_s3 %>%
  finalfit(dependent, explanatory) -> t3
t3

