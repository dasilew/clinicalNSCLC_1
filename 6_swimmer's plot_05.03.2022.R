#swimmer's plot
#data called NSCLC 22.2

#updated 06.03.2022: do it only for the two groups of interest
# Laod required libraries
#install.packages(c("survival", "survminer"))
library(gtsummary)
#library(MatchIt)
library(Rcpp)
library(survival)
library(survminer)
library(gapminder)
library(dplyr)
library(swimplot)
library(tidyverse)
library(readxl)
library(swimplot)
library(flextable)
library(janitor)

# original column names
dat<-NSCLC22_2
names(dat)
# data structure
glimpse(dat)
#View(dat)
#exploratory data analysis
head(dat)
tail(dat)

### Remove duplicated columns

dat2<-dat[!duplicated(as.list(dat))]

glimpse(dat2)

### Remove empty columns and rows

dat3<-dat2 %>% remove_empty()

glimpse(dat3)

#We select the important columns
#patient_ID, patient_ID2,
#date_of_first_operative_brain_metastasis_resection_OP1...35,
#date_of_second_operative_brain_metastasis_resection_OP... 36,
#date_of_third_operative_brain_metastasis_resection_OP3... 37,
#date_of_fourth_operative_brain_metastasis_resection_OP4...38,
#date_of death_according_to_Gießener_Tumordokumentationssystem...39,
#last_contact....40,
#DIA_DAT  ... 41,
#CTx_start 50,
#CTx_end ...51,
#ICI_start.. 56,
#ICI_end -...57,
#adj.tx2 ...60
dat4<-dat3 %>% select(5,6,35:41,50,51,55,56,59)

names(dat4)

# rename the columns

names(dat4)<-c("patient_ID", "patient_ID2", 
               "OP1_date","OP2_date","OP3_date",
               "OP4_date","death_date","last_contact","diagnosis_date",
               
               "CTx_start", "CTx_end", "ICI_start", "ICI_end", "adj_tx2")


#glimpse(dat4)

#names(dat4)


dat5<-dat4 %>% 
  
  dplyr::mutate(death = as.numeric(difftime(death_date,diagnosis_date,
                                            
                                            units = "days")),
                OP1 = as.numeric(difftime(OP1_date, diagnosis_date, 
                                          units = "days")),
                
                OP2 =as.numeric(difftime(OP2_date, diagnosis_date, 
                                         units = "days")),
                
                OP3 =as.numeric(difftime(OP3_date, diagnosis_date, 
                                         units = "days")),
                
                OP4 =as.numeric(difftime(OP4_date, diagnosis_date, 
                                         units = "days")),
                
                follow_up = as.numeric(difftime(last_contact, diagnosis_date, 
                                                units = "days")),
                
                Chemo_start = as.numeric(difftime(CTx_start, diagnosis_date,
                                                  
                                                  units = "days")),
                
                Chemo_end = as.numeric(difftime(CTx_end, diagnosis_date,
                                                
                                                units = "days")),
                
                Immuno_start = as.numeric(difftime(ICI_start, diagnosis_date,
                                                   
                                                   units = "days")),
                
                Immuno_end = as.numeric(difftime(ICI_end, diagnosis_date,
                                                 
                                                 units = "days"))) %>%
  
  group_by(patient_ID) %>% 
  
  mutate(max_event = max(death, OP1,OP2, OP3,OP4,Chemo_start,Chemo_end,
                         Immuno_start, Immuno_end,
                         na.rm = T)) %>% ungroup() %>%
  
  dplyr::mutate(follow_up = ifelse(follow_up < max_event,
                                   
                                   max_event, follow_up)) %>%
  
  dplyr::mutate(death = ifelse(death < follow_up,
                               
                               follow_up, death))

glimpse(dat5)
summary(dat5$Chemo_start)

summary(dat5$Chemo_end)

summary(dat5$Immuno_start)

summary(dat5$Immuno_end)

dat5 <-dat5 %>% dplyr::mutate(Therapy = ifelse(is.na(Chemo_start),
                                               "Immunotherapy",
                                               
                                               "Chemotherapy"))

table(dat5$Therapy)
#We have 108 chemotherapy patients and 63 immunotherapy patients, so a total of 172 patients.

#Our data is now composed of 174 rows and 26 columns (glimpse function)

library(rstatix)
library(flextable)


swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", fill= "lightblue",
             
             width=0.9, col= "black")+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))
#data of the plot

dat5 %>% arrange(desc(follow_up)) %>% 
  
  dplyr::select(patient_ID2,follow_up) %>%
  
  flextable() %>% align(align = "center", part = "all")


#stratify by therapy
swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))


dat5 %>% arrange(desc(follow_up)) %>% 
  
  dplyr::select(patient_ID2,follow_up, Therapy) %>%
  
  pivot_wider(names_from = "Therapy", values_from = "follow_up") %>%
  
  flextable() %>% align(align = "center", part = "all")

# some summary statistics

dat5 %>% group_by(Therapy) %>% 
  
  get_summary_stats(follow_up, type = "common") %>%
  
  flextable() %>% align(align = "center", part = "all")



dat6<- dat5 %>% tidyr::gather(key = "Event", value = "value", death:OP4)

# looking at the new columns

dat6 %>% dplyr::select(patient_ID2, Event, value)

# data dimensions

dim(dat6)

# frequency table of events

table(dat6$Event)



swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))+
  
  swimmer_points(df_points = dat6 %>% data.frame(stringsAsFactors = F),
                 
                 id = "patient_ID2", name_shape = "Event", time = "value",
                 
                 size = 3)+
  
  guides(fill = guide_legend(override.aes = list(shape = NA)))









swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))+
  
  swimmer_points(df_points = dat6 %>% data.frame(stringsAsFactors = F),
                 
                 id = "patient_ID2", time = "value",
                 
                 size = 3, name_col = "Event")+
  
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  
  scale_color_manual(name = "Event", values = c("#404040", "#9480fc",
                                                
                                                "#c296c0","#ebdf2e","#377b6d"))








swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))+
  
  swimmer_points(df_points = dat6 %>% data.frame(stringsAsFactors = F),
                 
                 id = "patient_ID2", time = "value",
                 
                 size = 3, name_col = "Event", name_shape = "Event")+
  
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  
  scale_color_manual(name = "Event", values = c("#404040", "#9480fc",
                                                
                                                "#c296c0","#ebdf2e","#377b6d"))










# gather start columns and remove NA from start value column

dat7<-dat5 %>% tidyr::gather(key = "start", value = "Therapy_start", 
                             
                             c(Chemo_start,Immuno_start)) %>%
  
  drop_na("Therapy_start")

head(dat7$start)

# gather end columns and remove NA from end value column

# add start columns

# add another column for the therapy period

dat8<-dat5 %>% tidyr::gather(key = "end", value = "Therapy_end", 
                             
                             c(Chemo_end,Immuno_end)) %>%
  
  drop_na("Therapy_end") %>% full_join(dat7) %>%
  
  dplyr::mutate(Period = ifelse(start=="Chemo_start","Chemotherapy period",
                                
                                "Immunotherapy period"))

head(dat8$end)

table(dat8$Therapy)

table(dat8$Period)



swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 7, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))+
  
  swimmer_points(df_points = dat6 %>% data.frame(stringsAsFactors = F),
                 
                 id = "patient_ID2", name_shape = "Event", time = "value",
                 
                 size = 3)+
  
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  
  
  swimmer_lines(df_lines = dat8 %>% data.frame(stringsAsFactors = F),
                
                id="patient_ID2",name_linetype = "Period",
                
                start = "Therapy_start", end="Therapy_end",
                
                size=1)









sp1 <- swimmer_plot(df=dat5 %>% data.frame(stringsAsFactors = F), id= "patient_ID2",
             
             end= "follow_up", name_fill= "Therapy", col = "black",
             
             width=0.9, stratify = c("Therapy"))+
  
  theme(axis.text.y = element_text(size = 6, face = "bold"))+
  
  scale_fill_manual(values = c("#66CC33","#FF6666"))+
  
  swimmer_points(df_points = dat6 %>% data.frame(stringsAsFactors = F),
                 
                 id = "patient_ID2", time = "value",
                 
                 size = 2, name_col = "Event")+
  
  
  swimmer_lines(df_lines = dat8 %>% data.frame(stringsAsFactors = F),
                
                id="patient_ID2",name_linetype = "Period",
                
                start = "Therapy_start", end="Therapy_end",
                
                size=1)+
  scale_color_manual(name = "Event", values = c("#404040", "#9480fc",
                                                
                                                "#c296c0","#ebdf2e","#377b6d"))+
  
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         
         color = guide_legend(override.aes = list(linetype = NA)))
sp1

