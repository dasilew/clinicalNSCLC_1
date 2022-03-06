#NSCLC brain metastasis table Sunday 30.01.2022, updated 13.02.2022
#Correlation plots
#summarytable merged covariates from the above tables
#update 05.03.2022
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

#Correlation plots
#summarytable merged covariates from the above tables

sumtable2 <- sumtable %>% select(
  no_of_BrMs, localization_of_BrMs, volume, extracran.met, 
  NLR, PDL1_intracranial, Ki67_intracranial, PDL1_extracranial, Ki67_extracranial, TTF1_status_BrM)

library(readxl)
library(tidyverse)
library(flextable)
install.packages("gtsummary")
# install.packages("gtsummary")
library(gtsummary)
library(boot) 
# Load required packages
library(survival)
library(survminer)
library(dplyr)
library(readxl)
volume

sumtable3 <- sumtable2 %>% 
  
  dplyr::mutate(PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                extracran.met = as.numeric(factor(extracran.met, exclude = "NA")),
                
                localization_of_BrMs = as.numeric(factor(localization_of_BrMs, exclude = "NA")),
                
                volume = as.numeric(factor(volume, exclude = "NA")),
                
                NLR = as.numeric(factor(NLR, exclude = "NA")),
                
                PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                Ki67_intracranial = as.numeric(factor(Ki67_intracranial, exclude = "NA")),
                
                PDL1_extracranial = as.numeric(factor(PDL1_extracranial, exclude = "NA")),
                
                Ki67_extracranial = as.numeric(factor(Ki67_extracranial, exclude = "NA")),
                
                
                TTF1_status_BrM = as.numeric(factor(TTF1_status_BrM, exclude = "NA"))
                
  )

glimpse(sumtable3)

sumtable4 <- sumtable3 %>% 
  
  dplyr::mutate(PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                localization_of_BrMs = as.numeric(factor(localization_of_BrMs, exclude = "NA")),
                
                volume = as.numeric(factor(volume, exclude = "NA")),
                
                NLR = as.numeric(factor(NLR, exclude = "NA")),
                
                PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                Ki67_intracranial = as.numeric(factor(Ki67_intracranial, exclude = "NA")),
                
                PDL1_extracranial = as.numeric(factor(PDL1_extracranial, exclude = "NA")),
                
                Ki67_extracranial = as.numeric(factor(Ki67_extracranial, exclude = "NA")),
                
                
                TTF1_status_BrM = as.numeric(factor(TTF1_status_BrM, exclude = "NA"))
                
  )

glimpse(sumtable4)





library(GGally)

ggcorr(data = sumtable4, method = c("pairwise","kendall"),
       
       label = T, label_round = 2, label_size = 3, size = 3, hjust = 1)

sumtable3

sumtable4<-sumtable3 %>% tidyr::gather(key = "variable", value = "value",
                                       
                                       c(1:7, 9:12))

glimpse(sumtable4)


library(corrplot)

## generate correlation matrix

cor_matrix <- cor(sumtable4, method = "kendall",
                  
                  use = "pairwise")

cor_matrix
corrplot(cor_matrix, method = 'number', tl.cex = 0.75, number.cex = 0.75)

corrplot(cor_matrix, tl.cex = 0.75, 
         
         addCoef.col = "black", number.cex = 0.75)



test_res <- cor.mtest(sumtable4, conf.level = 0.95, method = "kendall")

corrplot(cor_matrix, p.mat = test_res$p)

corrplot(cor_matrix, p.mat = test_res$p, type = "lower", insig= "blank",
         
         method = "color", addCoef.col = "white", number.cex = 0.75, 
         
         diag = F,
         
         addgrid.col = "black")



corrplot(cor_matrix, p.mat = test_res$p, type = "lower", insig= "label_sig",
         
         method = "color", diag = F, addgrid.col = "black", pch.cex = 0.9)




####


#summarytable merged covariates from the above tables

sumtable22 <- sumtable %>% select(
  no_of_BrMs, localization_of_BrMs, volume, extracran.met, 
  NLR, PDL1_intracranial, Ki67_intracranial, PDL1_extracranial, Ki67_extracranial, survival)

library(readxl)
library(tidyverse)
library(flextable)
install.packages("gtsummary")
# install.packages("gtsummary")
library(gtsummary)
library(boot) 
# Load required packages
library(survival)
library(survminer)
library(dplyr)
library(readxl)

sumtable23 <- sumtable22 %>% 
  
  dplyr::mutate(PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                extracran.met = as.numeric(factor(extracran.met, exclude = "NA")),
                
                localization_of_BrMs = as.numeric(factor(localization_of_BrMs, exclude = "NA")),
                
                volume = as.numeric(factor(volume, exclude = "NA")),
                
                NLR = as.numeric(factor(NLR, exclude = "NA")),
                
                PDL1_intracranial = as.numeric(factor(PDL1_intracranial, exclude = "NA")),
                
                Ki67_intracranial = as.numeric(factor(Ki67_intracranial, exclude = "NA")),
                
                PDL1_extracranial = as.numeric(factor(PDL1_extracranial, exclude = "NA")),
                
                Ki67_extracranial = as.numeric(factor(Ki67_extracranial, exclude = "NA")),
                
                
                survival = as.numeric(factor(survival, exclude = "NA"))
                
  )

glimpse(sumtable23)






library(GGally)

ggcorr(data = sumtable23, method = c("pairwise","kendall"),
       
       label = T, label_round = 2, label_size = 3, size = 3, hjust = 1)

sumtable23


library(corrplot)

## generate correlation matrix

cor_matrix <- cor(sumtable23, method = "spearmann",
                  
                  use = "pairwise")

cor_matrix
corrplot(cor_matrix, method = 'number', tl.cex = 0.75, number.cex = 0.75)

corrplot(cor_matrix, tl.cex = 0.75, 
         
         addCoef.col = "black", number.cex = 0.75)



test_res <- cor.mtest(sumtable23, conf.level = 0.95, method = "spearmann")

corrplot(cor_matrix, p.mat = test_res$p)

corrplot1 <- corrplot(cor_matrix, p.mat = test_res$p, type = "lower", insig= "blank",
         
         method = "color", addCoef.col = "white", number.cex = 0.75, 
         
         diag = F,
         
         addgrid.col = "black")

corrplot1



corrplot2 <- corrplot(cor_matrix, p.mat = test_res$p, type = "lower", insig= "label_sig",
         
         method = "color", diag = F, addgrid.col = "black", pch.cex = 0.9)

corrplot2

