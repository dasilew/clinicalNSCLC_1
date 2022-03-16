#1. Summarise variables/factors by a categorical variable
#updated 05.03.2022
library(finalfit)
library(dplyr)

# Table 1 - Patient demographics by variable of interest ----
#NSCLC15 <- read.csv2(file.choose())

#TRY AS IN THE MANUSCRIPT
#localization_of_BrMs
#no_of_BrMs
#extracran.met
#primary_tumor.Rx_2
#no_of_BrM_Rxs2
#pre.tx
#adj.tx3
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

my_table2 <- NSCLC22 %>% select(no_of_BrMs, GPA2, localization_of_BrMs, primary_tumor.Rx2, extracran.met, adj.RTx3, adj.tx2,
                                censoring, survival)


# Issue with results (see warning message): For adj.tx2 and adj.tx3, values are NA for some categories
# Reason is multicolinearity: Both variables highly correlate with each other
# See the following cross-table: category 2 for adj.tx3 has in adj.tx2 only values in category 3 and so on
# Recommendation: use only one of the in the multivariate analysis
#table(colon_s3$adj.tx2, colon_s3$adj.tx3)

# Fix is shown with table using only adj.tx2 and omitting adj.tx3
explanatory = c("localization_of_BrMs", "GPA2", "no_of_BrMs", "extracran.met", "primary_tumor.Rx2", 
                 "adj.RTx3", "adj.tx2")
dependent = "Surv(survival, censoring)"

colon_s3 %>%
  finalfit(dependent, explanatory) -> t3.mod
knitr::kable(t3.mod, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"))
t3.mod

# Cox Proportional Hazards multivariable analysis. 
#https://finalfit.org/reference/ff_merge.html

explanatory = c("localization_of_BrMs", "GPA2", "no_of_BrMs", "extracran.met", "primary_tumor.Rx2", 
                "adj.RTx3", "adj.tx2")
dependent = "Surv(survival, censoring)"


# Create separate tables
colon_s3 %>%
  summary_factorlist(dependent, explanatory, fit_id=TRUE)
#> Note: dependent includes missing data. These are dropped.

a <- colon_s3 %>%
  coxphuni(dependent, explanatory) %>%
  fit2df(estimate_suffix=" (CPH univar. analysis)")
a
b <- colon_s3 %>%
  coxphmulti(dependent, explanatory) %>%
  fit2df(estimate_suffix=" (CPH multivar.analysis)")
b
class(b)

#merge a and b

c <- right_join(a, b, by = 'explanatory')
c
knitr::kable(c, row.names=FALSE, align=c("l", "l", "r", "r", "r", "r"))
c



# Create separate tables
axample.summary <- colon_s3 %>%
  summary_factorlist(dependent, explanatory, fit_id=TRUE)
#> Note: dependent includes missing data. These are dropped.

example.univariable <- colon_s %>%
  glmuni(dependent, explanatory) %>%
  fit2df(estimate_suffix=" (univariable)")
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...

colon_s %>%
  glmmulti(dependent, explanatory) %>%
  fit2df(estimate_suffix=" (multivariable)") -> example.multivariable
#> Waiting for profiling to be done...

colon_s %>%
  glmmixed(dependent, explanatory, random_effect) %>%
  fit2df(estimate_suffix=" (multilevel)") -> example.multilevel

# Pipe together

