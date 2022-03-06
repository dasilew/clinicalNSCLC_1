#NSCLC brain metastasis table Sunday 30.01.2022, updated 13.02.2022
#updated 15.02.2022
#update 24.02.2022
#update 01.03.2022
#update 05.03.2022
#Figure 1
install.packages("gtsummary")
# install.packages("gtsummary")
library(gtsummary)
library(boot) 
# Load required packages
library(survival)
library(survminer)
library(dplyr)
NSCLC22
glimpse(NSCLC22)
colnames(NSCLC22) #Names of Columns which are the names of predictors and outcome variables
str(NSCLC22) # Structure of the dataset
#tabel1a
my_table1 <- NSCLC22 %>% select(gender, age, anatomical_localization, hydrocephalus, hemorrhage, LMM, 
                                no_of_BrM_Rxs, no_of_BrMs, localization_of_BrMs, volume, extracran.met, UICC2, KPS, GPA)
my_table1 %>% mutate(age_group = ifelse(age >=70, "70 years or older", "less than 70 years"))
my_table1 %>% mutate(no_of_BrM_Rxs = ifelse(no_of_BrM_Rxs >1, "more than one", "one"))

my_table1$no_of_BrM_Rxs <- factor(my_table1$no_of_BrM_Rxs)
my_table1 %>% mutate(no_of_BrMs= ifelse(no_of_BrMs >1, "more than one", "one"))
my_table1$no_of_BrMs <- factor(my_table1$no_of_BrMs)
my_table1$localization_of_BrMs <- factor(my_table1$localization_of_BrMs, 
                                         levels = c("0", "1", "2"), 
                                         labels = c("supratentorial", "infratentorial", "both"))

my_table1$no_of_BrMs <- factor(my_table1$no_of_BrMs, 
                               levels = c("0", "1", "2"), 
                               labels = c("one", "two", "more than two"))

my_table1 %>% mutate(volume = ifelse(volume >15, "> 15 mL", "< 15 mL"))



my_table1$extracran.met <- factor(my_table1$extracran.met, 
                                  levels = c("0", "1"), 
                                  labels = c("No", "Yes"))

my_table1$extracran.met <- factor(my_table1$UICC2, 
                                  levels = c("0", "1"), 
                                  labels = c("UICC stage IV", "UICC stage I-III"))

my_table1$anatomical_localization <- factor(my_table1$anatomical_localization, 
                                            levels = c("1", "2", "3", "4", "5", "6"), 
                                            labels = c("frontal", "parietal", "temporal", "occipital", "cerebellar", "other"))

#my_table1$leading_symptom <- factor(my_table1$leading_symptom, 
                                            #levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
                                            #labels = c("headache", "sensory-motor symptoms or hemiparesis", "seizures", 
                                                     #  "alterations behaviour or desorientation", "aphasia", 
                                                     #  "vertigo and dyscoordination", "headache, nausea, vomiting",
                                                      # "cerebellar symptoms", "visual impairment", "incidental"))



my_table1 <- my_table1 %>% mutate(KPS_group = ifelse(KPS >70, ">/= 70%", "< 70%"))
my_table1$KPS_group <- factor(my_table1$KPS_group)
my_table1 <- my_table1 %>% mutate(GPA_group = ifelse(GPA >=2, " >/=2", "<2"))
my_table1$GPA_group  <- factor(my_table1$GPA_group )

theme_gtsummary_journal(journal = "jama")
#> Setting theme `JAMA`
theme_gtsummary_compact()
#> Setting theme `Compact`
#SUpplementary table 1
#Table 1. Patient Characteristics
my_table1
my_table1a <- my_table1 %>% select(gender, age, anatomical_localization, 
                                   hydrocephalus, hemorrhage, LMM, no_of_BrMs, localization_of_BrMs, volume, 
                                   extracran.met, UICC2, KPS_group, GPA_group)
my_table1a
my_table1a %>% tbl_summary(
  label = list(age ~ "Age", 
               gender ~ "Gender", 
               anatomical_localization  ~ "anatomical localization",
               hydrocephalus ~ "hydrocephalus at baseline", 
               hemorrhage ~ "hemorrhage at baseline", 
               LMM ~ "LMM at baseline by means of MRI", 
               #no_of_BrM_Rxs ~ "Number of brain metastasis resections",
               no_of_BrMs  ~ "Number of brain metastases at baseline",
               localization_of_BrMs ~ "Brain region affected of brain metastasis",
               volume ~ "volume of dominant (resected) brain metastasis (mL)",
               UICC2 ~ "UICC stage",
               extracran.met ~ "Extracranial metastasis at baseline", 
               KPS_group ~ "KPS after operation", 
               GPA_group ~ "GPA after operation "
               
  )
) %>%
  modify_caption("**Table 1: Patient characteristics**") %>%
  bold_labels()

#Table 2 Patient Characteristics - treatment-related characteristics
my_table2 <- NSCLC22 %>% select(no_of_BrM_Rxs, rx_of_independent_BrM, rx_due_to_local_recurrence, time_to_BrM_recurrence, 
                                primary_tumor.Rx, primary_tumor.Rx2, pre.tx, naive, adj.RTx2, adj.RTx, dose_of_RTx, adj.tx.modalities,
                                adj.tx2, adj.tx3)


my_table2$primary_tumor.Rx <- factor(my_table2$primary_tumor.Rx , 
                                     levels = c("0", "1", "2"), 
                                     labels = c("no primary tumor resection", 
                                                "primary tumor resection before brain metastasis resection",
                                                "primary tumor resection after brain metastasis resection"))

my_table2$primary_tumor.Rx2 <- factor(my_table2$primary_tumor.Rx2 , 
                                       levels = c("0", "1"), 
                                       labels = c("no primary tumor resection", 
                                                  "primary tumor resection"))

my_table2$pre.tx <- factor(my_table2$pre.tx, 
                          levels = c("0", "2", "1", "3"), 
                          labels = c("no neoadjuvant treatment", "RTx + CTx", "RTx + ICI", 
                                      "RTx + SMI"))
my_table2$naive <- factor(my_table2$naive, 
                                    levels = c("0", "1"), 
                                    labels = c("no pre-treatment", "pre-treatment"))


my_table2$adj.RTx <- factor(my_table2$adj.RTx, 
                            levels = c("wbr", "rtx", "SRS", "NA"), 
                            labels = c("WBRT", "local irradation of resection cavity", "SRS", "not known"))


my_table2 %>% mutate(dose_of_RTx = ifelse(dose_of_RTx >29, "30 Gy or more ", "less than 30 Gy"))


my_table2$adj.tx2 <- factor(my_table2$adj.tx2 , 
                            levels = c("0", "1"), 
                            labels = c("RTx + CTx", "RTx + ICI"))

my_table2$adj.tx3 <- factor(my_table2$adj.tx3 , 
                            levels = c("0", "1", "2"), 
                            labels = c("RTx + CTx", "RTx + ICI", "RTx + SMIs"))

#my_table2$censoring <- factor(my_table2$censoring, 
                     #         levels = c("0", "1"), 
                      #        labels = c("censored", "death occured"))


#my_table1b <- my_table2 %>% select(primary_tumor.Rx, primary_tumor.Rx_2, neoadj.Tx, naive_at_BrM.Rx, adj.tx, adj.tx2, time.to.BrM.Rx, censoring, survival)

my_table2 <- my_table2 %>% select(no_of_BrM_Rxs, rx_of_independent_BrM, rx_due_to_local_recurrence, time_to_BrM_recurrence, 
primary_tumor.Rx, primary_tumor.Rx2, pre.tx, adj.RTx, adj.RTx2, dose_of_RTx, adj.tx2, adj.tx3, adj.tx.modalities)



my_table2 %>% tbl_summary(
  missing = "no", ###NIIIIIICE
  label = list(no_of_BrM_Rxs ~ "number of brain metastasis resections", 
               rx_of_independent_BrM ~ "resection of independent lesions",
               no_of_BrM_Rxs ~ "primary tumor resection in context of brain metastasis resection",
               rx_due_to_local_recurrence ~ "resection due to local recurrence",
               time_to_BrM_recurrence ~ "time to local re-resection",
               primary_tumor.Rx ~ "primary tumor resection in context of brain metastasis resection", 
               primary_tumor.Rx2 ~ "primary tumor resection",
               no_of_BrM_Rxs ~ "Number of brain metastasis resections",
               pre.tx ~ "pre-treatment at baseline", 
               adj.RTx ~ "mode of cranial radiation after first brain metastasis resection", 
               adj.RTx2 ~ "fraction of stereotactically irradiated patients", 
               dose_of_RTx ~ "dose of cranial radiation after first brain metastasis resection", 
               adj.tx2 ~ "adjuvant treatment groups of interest: RTx + CTx vs. RTx + ICI",
               adj.tx3 ~ "main adjuvant treatment groups",
               adj.tx.modalities ~ "adjuvant treatment modalities"
            
               
               
  )
) %>%
  modify_caption("**Table 2: Treatment-related patient characteristics**") %>%
  bold_labels() ###supplementary table 2

##Table 3: Histopathological and biomarker-related patient characteristics
NSCLC22
my_table3 <- NSCLC22 %>% select(NLR, entity, histology, PDL1_intracranial, Ki67_intracranial, 
                                PDL1_extracranial, Ki67_extracranial, 
                                EGFR_status_BrM, ALK_status_BrM, TTF1_status_BrM)
my_table3

my_table3 <-my_table3 %>% mutate(NL_ratio = ifelse(NLR >=5, ">/= 5", "< 5"))
my_table_1c <- my_table3
my_table_1c
my_table_1c <-my_table3 %>% mutate(PDL1_intracranial = ifelse(PDL1_intracranial >=1, ">/= 1", "< 1"))
my_table_1c
my_table_1cc <-my_table_1c %>% mutate(PDL1_extracranial = ifelse(PDL1_extracranial >=1, ">/= 1", "< 1"))
my_table_1cc


my_table_1ccc <-my_table_1cc %>% mutate(Ki67_intracranial = ifelse(Ki67_intracranial >=30, ">/= 30", "< 30"))
my_table_1ccc

my_table_1cccc <-my_table_1ccc %>% mutate(Ki67_extracranial = ifelse(Ki67_extracranial >=30, ">/= 30", "< 30"))
my_table_1cccc

my_table_1cccc$EGFR_status_BrM <- factor(my_table_1cccc$EGFR_status_BrM , 
                                         levels = c("0", "1"), 
                                         labels = c("wildtype", "mutated"))
my_table_1cccc$ALK_status_BrM <- factor(my_table_1cccc$ALK_status_BrM , 
                                        levels = c("0", "1"), 
                                        labels = c("wildtype", "mutated"))
my_table_1cccc

my_table_1ccccc <- my_table_1cccc %>% select(NL_ratio, 
                                             entity, histology, 
                                             PDL1_intracranial, Ki67_intracranial, 
                                             PDL1_extracranial, Ki67_extracranial,
                                             EGFR_status_BrM, 
                                             ALK_status_BrM, TTF1_status_BrM)
my_table_1ccccc
my_table_1ccccc %>% tbl_summary(
  #missing = "no", ###NIIIIIICE
  label = list(NL_ratio ~ "neutrophil-to-lymphocyte ratio", 
               entity ~ "tumor subtype", 
               histology ~ "histologic subtype", 
               PDL1_intracranial ~ "PDL1 status brain metastasis tissue", 
               Ki67_intracranial ~ "Ki67 brain metastasis tissue", 
               PDL1_extracranial ~ "PDL1 extracranial tissue", 
               Ki67_extracranial ~ "Ki67 extracranial tissue", 
               EGFR_status_BrM ~ "EGFR status of brain metastasis tissue", 
               ALK_status_BrM ~ "ALK status of brain metastasis tissue", 
               TTF1_status_BrM ~ "TTF1 status of brain metastasis tissue"
               
  )
) %>%
  modify_caption("**Table 3: Histopathological and biomarker-related patient characteristics**") %>%
  bold_labels()

#Fit survival curves
#Fig1
require("survival")
fit <- survfit(Surv(survival, censoring) ~ 1, data = NSCLC22)
# Drawing curves
ggsurvplot(fit, color = "#2E9FDF")
# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, 
           #title = "Fig 1",
           data = NSCLC22,
           #pval = TRUE, 
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.y.text.col = TRUE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit
NSCLC22

##############################################################

################################################################
#other KM curves

pending

sumtable <- NSCLC22 %>% select(gender, age, anatomical_localization, 
                                   hydrocephalus, hemorrhage, LMM, no_of_BrM_Rxs, no_of_BrMs, localization_of_BrMs, volume, 
                                   extracran.met, UICC2, KPS, GPA, no_of_BrM_Rxs, rx_of_independent_BrM, rx_due_to_local_recurrence, 
                               time_to_BrM_recurrence, 
                                  primary_tumor.Rx, primary_tumor.Rx2, pre.tx, adj.RTx2, dose_of_RTx, adj.tx2, adj.tx3, adj.tx.modalities, 
                               NLR, entity, histology, PDL1_intracranial, Ki67_intracranial, 
                                PDL1_extracranial, Ki67_extracranial, 
                                EGFR_status_BrM, ALK_status_BrM, TTF1_status_BrM, survival)




sumtable <- NSCLC22 %>% select(gender, age, anatomical_localization, 
                               hydrocephalus, hemorrhage, LMM, no_of_BrM_Rxs, no_of_BrMs, localization_of_BrMs, volume, 
                               extracran.met, UICC2, KPS, GPA, no_of_BrM_Rxs, rx_of_independent_BrM, rx_due_to_local_recurrence, 
                               time_to_BrM_recurrence, 
                               primary_tumor.Rx, primary_tumor.Rx2, pre.tx, adj.RTx2, dose_of_RTx, adj.tx2, adj.tx3, adj.tx.modalities, 
                               NLR, entity, histology, PDL1_intracranial, Ki67_intracranial, 
                               PDL1_extracranial, Ki67_extracranial, 
                               EGFR_status_BrM, ALK_status_BrM, TTF1_status_BrM, survival)




############################################


#Fig2a extracranial mets
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit2 <- survfit(Surv(survival, censoring) ~ 1, data = NSCLC22)
summary(fit2)
print(fit2)
ggsurvplot(fit2,
           #title    = "Fig 2a",
           data = summarytable,
           legend.labs = c("without exracranial metastases", "extracranial metastases"),
           pval = TRUE, pval.coord = c(15, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit2














############################



summarytable$no_of_BrM_Rxs2 <- factor(summarytable$no_of_BrM_Rxs2, 
                                      levels = c("1", "2"), 
                                      labels = c("one", "two or more"))

summarytable$no_of_BrMs <- factor(summarytable$no_of_BrMs, 
                                  levels = c("0", "1", "2"), 
                                  labels = c("one", "two", "more than two"))

summarytable$primary_tumor.Rx2 <- factor(summarytable$primary_tumor.Rx2, 
                                         levels = c("0", "1", "2"), 
                                         labels = c("no primary tumor removal", "primary tumor removal before brain metastasis resection", 
                                                    "primary tumor removal after brain metastasis resection"))

summarytable$localization_of_BrMs <- factor(summarytable$localization_of_BrMs, 
                                            levels = c("0", "1", "2"), 
                                            labels = c("supratentorial", "infratentorial", "both"))
summarytable
summarytable$extracran.met <- factor(summarytable$extracran.met, 
                                     levels = c("0", "1"), 
                                     labels = c("No", "Yes"))
summarytable<- summarytable %>% mutate(KPS_group = ifelse(KPS >70, ">/= 70%", "< 70%"))
summarytable$KPS_group <- factor(summarytable$KPS_group)
summarytable <- summarytable %>% mutate(NLR_group = ifelse(NLR >=5, " >/=5", "<5"))
summarytable$NLR_group  <- factor(summarytable$NLR_group)
summarytable <- summarytable %>% mutate(PDL1_intracranial_group = ifelse(PDL1_intracranial >=1, " >/=1", "<1"))
summarytable$PDL1_intracranial_group  <- factor(summarytable$PDL1_intracranial_group)

#############################

#Fig2a extracranial mets
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit2 <- survfit(surv_object ~ extracran.met, data = summarytable)
summary(fit2)
print(fit2)
ggsurvplot(fit2,
           #title    = "Fig 2a",
           data = summarytable,
           legend.labs = c("without exracranial metastases", "extracranial metastases"),
           pval = TRUE, pval.coord = c(15, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit2

#Fig2b localization
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit3 <- survfit(surv_object ~ localization_of_BrMs, data = summarytable)
summary(fit3)
print(fit3)
ggsurvplot(fit3,
           #title    = "Fig 2b",
           data = summarytable,
           legend.labs = c("supratentorial", "infratentorial", "both"),
           pval = TRUE, pval.coord = c(15, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit3

#Fig2c GPA_group
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit4 <- survfit(surv_object ~ KPS_group, data = summarytable)
summary(fit4)
print(fit4)
ggsurvplot(fit4,
           #title    = "Fig 2c",
           data = summarytable,
           legend.labs = c("KPS < 70", "KPS >/=70"),
           pval = TRUE, pval.coord = c(15, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit4
#Fig2d primary tumor resection
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit5 <- survfit(surv_object ~ adj.RTx, data = summarytable)
summary(fit5)
print(fit5)
ggsurvplot(fit5,
           #title    = "Fig 2d",
           data = summarytable,
           legend.labs = c("no PTRx", "PTRx before or after BrMRx"),
           pval = TRUE, pval.coord = c(35, 0.8),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit5

#display CI and median OS
#fit.list2 <- surv_fit(Surv(time, status) ~ sex, data = colon,
                      #group.by = "rx")
#surv_median(fit.list2)
library(survival)
surv_median(fit5)



#Fig2e primary tumor resection
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit6 <- survfit(surv_object ~ primary_tumor.Rx2, data = summarytable)
summary(fit6)
print(fit6)
ggsurvplot(fit6,
           #title    = "Fig 2e",
           data = summarytable,
           legend.labs = c("no PTRx", "PTRx", ),
           pval = TRUE, pval.coord = c(33, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit6

#Fig2f nr.brm resections
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
summarytable
fit7 <- survfit(surv_object ~ no_of_BrM_Rxs2, data = summarytable)
summary(fit7)
print(fit7)
ggsurvplot(fit7,
           #title    = "Fig 2f",
           data = summarytable,
           legend.labs = c("1 brain metastasis resection", ">/= 2 brain metastasis resections"),
           pval = TRUE, pval.coord = c(35, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit7

#Fig2g nr.brms
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
fit8 <- survfit(surv_object ~ no_of_BrMs, data = summarytable)
summary(fit8)
print(fit8)
ggsurvplot(fit8,
           #title    = "Fig 2g",
           data = summarytable,
           legend.labs = c("1 brain metastasis", "2 brain metastases", "> 2 brain metastases"),
           pval = TRUE, pval.coord = c(26, 0.6),
           palette = "Dark2",
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit8

#Fig2h adj.tx2 (updated 30.08.2021)
surv_object <- Surv(time = summarytable$survival, event = summarytable$censoring)
surv_object 
summarytable
fit8 <- survfit(surv_object ~ adj.tx3, data = summarytable)
summary(fit8)
ggsurvplot(fit8,
           #title    = "Fig 3",
           data = summarytable,
           legend.labs = c("CT + RT", 
                           "ICI + RT"
                           ),
           pval = TRUE, pval.coord = c(40, 0.8),
           palette = "Dark2", 
           conf.int = TRUE,
           risk.table = FALSE, 
           risk.table.y.text.col = FALSE,
           xlim = c(0,84),
           break.time.by = 12,
           xlab = "Time in months", 
           surv.median.line = "hv",  # add the median survival pointer.
)
fit8

----------------------------------------------------

summarytable
colnames(summarytable) 
str(summarytable)


devtools::install_github("vandomed/tab")
library("tab")
glm_v(
  censoring ~ poly(age, 2, raw = TRUE) + gender * KPS_group, 
  data = summarytable, 
  family = binomial
)





#NOT WORKING





#univariate and multivariate analysis
install.packages("finalfit")
install.packages("rstan")
install.packages("boot")
library(finalfit)
library(dplyr)
# Load example dataset, modified version of survival::colon
colon_s <- summarytable

#survival::coxph(dependent ~ explanatory)
explanatory = c("localization_of_BrMs", "no_of_BrMs", "extracran.met", "neoadj.tx", "adj.tx3", "PDL1_intracranial", "no_of_BrM_Rxs2", "primary_tumor.Rx_2")
dependent = "Surv(survival, censoring)"
colon_s %>%
  finalfit(dependent, explanatory) -> t3
t3
# Save objects for knitr/markdown
save(t3, dependent, explanatory, file = "out.rda")

options(kableExtra.auto_format = FALSE)
library(kableExtra)
t3 %>%
  kbl() %>%
  kable_paper("hover", full_width = F, font_size = 10)
t3 <-na.omit(t3)
t3 %>%
  kbl() %>%
  kable_paper("hover", full_width = F, font_size = 10)


dt %>%
  kbl() %>%
  kable_paper("hover", full_width = F)

