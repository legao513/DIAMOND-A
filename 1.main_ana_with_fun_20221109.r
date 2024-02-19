library(lubridate)
library(stringr)
library(dplyr)
library('readxl')
library("survival")
library("survminer")
library(dplyr)
library(tidyr)
library(haven)

options(digits = 12)
777/365.25


#step1: read 2001to2018 mother_child linkage data
mc_01to18 <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mc_01to18.rds")


#step2: select records with child's ref_id, dob between 2001 to 2018
mc_01to15 <- mc_01to18 %>% filter(!is.na(Baby.s.Reference.Key.)) %>% 
  mutate(Reference.Key.=as.numeric(Reference.Key.),
         Maternity.Episode..Delivery.Date..yyyy.mm.dd..=as.Date(Maternity.Episode..Delivery.Date..yyyy.mm.dd..),
         Maternity.Episode..Gestation.weeks.=as.numeric(Maternity.Episode..Gestation.weeks.),
         Date.of.Birth..yyyy.mm.dd..=as.Date(Date.of.Birth..yyyy.mm.dd..)) %>% 
  filter(year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)<=2018,
         year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)>=2001,
         !is.na(Date.of.Birth..yyyy.mm.dd..))#535924
# saveRDS(mc_01to15,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15.rds")
mc_01to15 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15.rds")
length(unique(mc_01to15$Reference.Key.))  #395457


#step3: read mother alldx
mc_01to18_mm_alldx <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mc_01to18_mm_alldx.rds")
mc_01to18_mm_alldx1 <- mc_01to18_mm_alldx %>% 
  mutate(Reference.Key.=as.numeric(Reference.Key.),Reference.Date.=ymd(Reference.Date.)) %>% 
  filter(Reference.Key.%in%mc_01to15$Reference.Key.)
#step3.1: select mother with PGDM (ICD9:250)
mm_pgdm_dx <- mc_01to18_mm_alldx1 %>% 
  filter(grepl("^250",All.Diagnosis.Code..ICD9..)==T)#8410
length(unique(mm_pgdm_dx$Reference.Key.)) #2461
# saveRDS(mm_pgdm_dx,"D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx.rds")
mm_pgdm_dx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx.rds")

#step3.2: select mother with GDM (ICD9:648.83)
mm_gdm_dx <- mc_01to18_mm_alldx1 %>% 
  filter(grepl("^648.83",All.Diagnosis.Code..ICD9..)==T)#71549
length(unique(mm_gdm_dx$Reference.Key.)) #35974
# saveRDS(mm_gdm_dx,"D:/data_no_share_20200506/6.DIAMOND/mm_gdm_dx.rds")
mm_gdm_dx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_gdm_dx.rds")




#step4: read mother allrx
mc_01to18_mm_allrx <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mc_01to18_mm_allrx.rds")
mc_01to18_mm_allrx1 <- mc_01to18_mm_allrx %>% 
  mutate(Reference.Key.=as.numeric(Reference.Key.),
         Dispensing.Date..yyyy.mm.dd..=ymd(Dispensing.Date..yyyy.mm.dd..),
         Prescription.Start.Date.=ymd(Prescription.Start.Date.),
         Prescription.End.Date.=ymd(Prescription.End.Date.),
         Drug.Name.=toupper(Drug.Name.))%>% 
  filter(Reference.Key.%in%mc_01to15$Reference.Key.)
#step4.1: select mother with PGDM (BNF:6.1)
mm_anti_drug <- mc_01to18_mm_allrx1 %>% 
  filter(grepl("^6.1",Therapeutic.Classification..BNF..Principal..)==T) %>% 
  mutate(Rxst=ifelse(is.na(Prescription.Start.Date.),Dispensing.Date..yyyy.mm.dd..,Prescription.Start.Date.),
         Rxed=ifelse(is.na(Prescription.End.Date.),NA,Prescription.End.Date.),
         Rxst=as.Date(Rxst,origin="1970-01-01"),Rxed=as.Date(Rxed,origin="1970-01-01"),
         Dosage.=gsub(' ',"",Dosage.))#122461
length(unique(mm_anti_drug$Reference.Key.)) #5993
#step4.2: missing value imputation
mm_anti_drug_1 <- mm_anti_drug %>% filter(!is.na(Prescription.End.Date.))
mm_anti_drug_2 <- mm_anti_drug %>% filter(is.na(Prescription.End.Date.))
table(factor(mm_anti_drug_2$Drug.Frequency.)) %>% as.data.frame()
table(factor(mm_anti_drug_2$Dosage.)) %>% as.data.frame()

mm_anti_drug_2$Drug_Frequency_Fix <- NA
mm_anti_drug_2$Drug_Frequency_Fix [which(grepl("^AT|^BEFORE|^DAILY|^ONCE|^WITH|AS DIRECTED|^BEFORE",mm_anti_drug_2$Drug.Frequency.)==T)] <- 1
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE AFTERNOON","IN THE EVENING","IN THE MORNING","WITH BREAKFAST","WITH EVENING MEAL"))] <- 1
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("EVERY SIX HOURS"))] <- 4
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("HOURLY"))] <- 24
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [0.5] AT NOON", "IN THE MORNING AND [0.5] IN THE AFTERNOON"))] <- 1.5
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [1.5] AT NOON"))] <- 2.5
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [1] AT NOON", "IN THE MORNING AND [1] IN THE AFTERNOON"))] <- 2
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [2] AT NIGHT","IN THE MORNING AND [2] AT NOON","IN THE MORNING AND [2] IN THE AFTERNOON"))] <- 3
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [3] AT NIGHT"))] <- 4
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("IN THE MORNING AND [4] AT NIGHT"))] <- 5
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("TWICE DAILY"))] <- 2
mm_anti_drug_2$Drug_Frequency_Fix [which(mm_anti_drug_2$Drug.Frequency. %in% c("THREE TIMES DAILY"))] <- 3
mm_anti_drug_2$Drug_Frequency_Fix [which(grepl("WHEN NECESSARY",mm_anti_drug_2$Drug.Frequency.)==T)] <- NA



mm_anti_drug_2_1 <- mm_anti_drug_2 %>% mutate(dosagefix=substring(Dosage.,1,3)) %>% 
  mutate(dosagefix=str_extract(dosagefix, "\\d+\\.*\\d*")) %>% 
  mutate(duration1=round(as.numeric(Quantity..Named.Patient..)/(as.numeric(dosagefix) * as.numeric(Drug_Frequency_Fix)))) %>% 
  mutate(Rxed=ymd(Rxst)+days(duration1)-1) %>% 
  select(-c(Drug_Frequency_Fix,dosagefix,duration1))
mm_anti_drug_new <- rbind(mm_anti_drug_1,mm_anti_drug_2_1)
x <- mm_anti_drug_new %>% filter(is.na(Rxed))
table(x$Drug.Name.)
x1 <- mm_anti_drug_new %>% filter(!is.na(Rxed)) %>% group_by(Drug.Name.) %>% 
  mutate(du=median(Rxed-Rxst+1)) %>% filter(Drug.Name.%in%x$Drug.Name.) %>% 
  filter(duplicated(Drug.Name.)==F) %>% select(Drug.Name.,du)
mm_anti_drug_new1 <- merge(mm_anti_drug_new,x1,by="Drug.Name.",all.x = T) %>% 
  mutate(Rxed1=ifelse(is.na(Rxed),Rxst+du-1,Rxed),Rxed1=as.Date(Rxed1,origin="1970-01-01")) %>% 
  select(-c(Rxed,du))
x <- mm_anti_drug_new1 %>% filter(is.na(Rxed1))

x2 <- mm_anti_drug_new1 %>% filter(!is.na(Rxed1)) %>% group_by(Drug.Item.Code.) %>% 
  mutate(du=median(Rxed1-Rxst+1)) %>% filter(Drug.Item.Code.%in%x$Drug.Item.Code.) %>% 
  filter(duplicated(Drug.Item.Code.)==F) %>% select(Drug.Item.Code.,du)
mm_anti_drug_new2 <- merge(mm_anti_drug_new1,x2,by="Drug.Item.Code.",all.x = T) %>% 
  mutate(Rxed2=ifelse(is.na(Rxed1),Rxst+du-1,Rxed1),Rxed2=as.Date(Rxed2,origin="1970-01-01")) %>% 
  select(-c(Rxed1,du))
# saveRDS(mm_anti_drug_new2,"D:/data_no_share_20200506/6.DIAMOND/mm_anti_drug_new2.rds")
mm_anti_drug_new2 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_anti_drug_new2.rds")



#step5: read lab test data
#remove space in LIS.Test.Description. and keep only number in LIS.Result.
mc_01to18_mm_lab  <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mc_01to18_mm_lab.rds")
mc_01to18_mm_lab1 <- mc_01to18_mm_lab %>% 
  mutate(Reference.Key.=as.numeric(Reference.Key.)) %>% 
  mutate(LIS.Reference.Datetime.=as.Date(LIS.Reference.Datetime.)) %>% 
  mutate(LIS.Result..Numeric.Result.=as.numeric(LIS.Result..Numeric.Result.)) %>% 
  mutate(LIS.Test.Description.1=toupper(gsub(' ',"",LIS.Test.Description.))) %>% 
  filter(Reference.Key.%in%mc_01to15$Reference.Key.)

#step5.1: select lab test for PGDM
#HbA1c >= 48mmol/mol (>=6.5 %)
mm_hba1cDCCT_pgdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("A1C_HBA1C(HPLC)_%","HAEMOGLOBINA1C","HBA1C","HEMOGLOBINA1C",
                                        "HEMOGLOBINA1C(IM)","HEMOGLOBINA1CRATIO","WHOLEBLOODHBA1C")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=6.5)
mm_hba1cIFCC_pgdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("HAEMOGLOBINA1C(IFCC)","HBA1C(IFCC)","HEMOGLOBINA1C(IFCC)",
                                        "WHOLEBLOODHBA1C(IFCC)")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=48)
#fasting plasma glucose concentration >= 7.0 mmol/L
mm_fpg_pgdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("FASTINGGLUCOSE","FASTINGPLASMAGLUCOSE","GLUCOSE(FASTING)",
                                        "GLUCOSE,FASTING","GLUCOSEFASTING","GLUCOSEFASTING,PLASMA",
                                        "MBI_GLUCOSE,FASTING_MMOL/L","MOG_GLUCOSE,FASTING_MMOL/L")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=7.0)
#2hplasma glucose >= 11.1mmol/L
mm_2hpg_pgdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("2-HOURPOSTPRANDIALPLASMAGLUCOSE","GLUCOSE(2HRS)",
                                        "GLUCOSE,2HOURSPOSTPRANDIAL","GLUCOSE,2HR",
                                        "GLUCOSE,2HRPOSTPRANDIAL","GLUCOSE,2HRPP",
                                        "GLUCOSE2HRPOSTPRANDIAL","GLUCOSETOLERANCETEST,2HR,PLASMA",
                                        "GLUCOSETOLERANCETEST2HOUR","GTT2HR",
                                        "MOG_GLUCOSE,2HOURPOSTLOAD_MMOL/L","OGTT2HR,PLASMA")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=11.1)

#combine all lab test meet the criteria for PGDM indentification
mm_lab_pgdm <- rbind(mm_hba1cDCCT_pgdm,mm_hba1cIFCC_pgdm,mm_fpg_pgdm,mm_2hpg_pgdm)#18516
length(unique(mm_lab_pgdm$Reference.Key.)) #3706
# saveRDS(mm_lab_pgdm,"D:/data_no_share_20200506/6.DIAMOND/mm_lab_pgdm.rds")
mm_lab_pgdm <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_lab_pgdm.rds")


#step5.2: select lab test for GDM
#fasting plasma glucose concentration >= 5.1 mmol/L
mm_fpg_gdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("FASTINGGLUCOSE","FASTINGPLASMAGLUCOSE","GLUCOSE(FASTING)",
                                        "GLUCOSE,FASTING","GLUCOSEFASTING","GLUCOSEFASTING,PLASMA",
                                        "MBI_GLUCOSE,FASTING_MMOL/L","MOG_GLUCOSE,FASTING_MMOL/L")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=5.1)
#2hplasma glucose >= 8.5mmol/L
mm_2hpg_gdm <- mc_01to18_mm_lab1 %>% 
  filter((LIS.Test.Description.1 %in% c("2-HOURPOSTPRANDIALPLASMAGLUCOSE","GLUCOSE(2HRS)",
                                        "GLUCOSE,2HOURSPOSTPRANDIAL","GLUCOSE,2HR",
                                        "GLUCOSE,2HRPOSTPRANDIAL","GLUCOSE,2HRPP",
                                        "GLUCOSE2HRPOSTPRANDIAL","GLUCOSETOLERANCETEST,2HR,PLASMA",
                                        "GLUCOSETOLERANCETEST2HOUR","GTT2HR",
                                        "MOG_GLUCOSE,2HOURPOSTLOAD_MMOL/L","OGTT2HR,PLASMA")==T)) %>% 
  filter(LIS.Result..Numeric.Result.>=8.5)

#combine all lab test meet the criteria for PGDM indentification
mm_lab_gdm <- rbind(mm_fpg_gdm,mm_2hpg_gdm)#31875
length(unique(mm_lab_gdm$Reference.Key.)) #18235
# saveRDS(mm_lab_gdm,"D:/data_no_share_20200506/6.DIAMOND/mm_lab_gdm.rds")
mm_lab_gdm <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_lab_gdm.rds")





#step6: find the exposure of mother
#step6.1 merge gdm dx rx and lab result, PGDM dx rx and lab
# mm_gdm_dx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_gdm_dx.rds")
# mm_lab_gdm <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_lab_gdm.rds")
# mm_anti_drug_new2 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_anti_drug_new2.rds")
# mm_pgdm_dx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx.rds")
# mm_lab_pgdm <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_lab_pgdm.rds")

mm_gdm_dx_col2 <- mm_gdm_dx %>% 
  select(Reference.Key.,Reference.Date.) %>% 
  `colnames<-` (c("Reference.Key.","gdm_date"))
mm_gdm_lab_col2 <- mm_lab_gdm %>% 
  select(Reference.Key.,LIS.Reference.Datetime.)%>% 
  `colnames<-` (c("Reference.Key.","gdm_date"))
mm_anti_rx_gdm_col2 <- mm_anti_drug_new2 %>% 
  select(Reference.Key.,Rxst) %>% 
  `colnames<-` (c("Reference.Key.","gdm_date"))
mm_pgdm_dx_col2 <- mm_pgdm_dx %>% 
  select(Reference.Key.,Reference.Date.) %>% 
  `colnames<-` (c("Reference.Key.","gdm_date"))
mm_pgdm_lab_col2 <- mm_lab_pgdm %>% 
  select(Reference.Key.,LIS.Reference.Datetime.)%>% 
  `colnames<-` (c("Reference.Key.","gdm_date"))

mm_gdm_dx_lab_rx <- rbind(mm_gdm_dx_col2,mm_gdm_lab_col2,mm_anti_rx_gdm_col2,mm_pgdm_dx_col2,mm_pgdm_lab_col2)
# saveRDS(mm_gdm_dx_lab_rx,"D:/data_no_share_20200506/6.DIAMOND/mm_gdm_dx_lab_rx.rds")
mm_gdm_dx_lab_rx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_gdm_dx_lab_rx.rds")

#combine with gdm dx, rx and lab data, find out gdm cases
median(mc_01to15$Maternity.Episode..Gestation.weeks.[which(!is.na(mc_01to15$Maternity.Episode..Gestation.weeks.))])#39
mc_01to15_1 <- mc_01to15 %>% 
  mutate(Maternity.Episode..Gestation.weeks.1=ifelse(is.na(Maternity.Episode..Gestation.weeks.),39,Maternity.Episode..Gestation.weeks.),
         pregancy_date=ymd(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)-weeks(as.numeric(Maternity.Episode..Gestation.weeks.1))) %>% #calculate the start date of the pregnancy
  #mutate(m_pgdm_dx_lab=NA,m_gdm_dx_lab=NA,m_pgdm_rx=NA,m_gdm_rx=NA) %>% 
  mutate(id=1:length(Reference.Key.))#give each record a new id

mc_01to15_2 <- merge(mc_01to15_1,mm_gdm_dx_lab_rx,by="Reference.Key.",all.x = T) %>% 
  mutate(gdm=ifelse(!is.na(gdm_date)&gdm_date>=pregancy_date&gdm_date<=Maternity.Episode..Delivery.Date..yyyy.mm.dd..,1,0)) %>% 
  #select(-c(Reference.Date.)) %>% 
  arrange(id,-gdm,gdm_date) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup()
table(mc_01to15_2$gdm) #46404
t1 <-mc_01to15_2 %>% filter(Reference.Key.==407865)
length(unique(t1$Baby.s.Reference.Key.))
write.csv(mc_01to15_2 %>% filter(gdm==1) %>% select(Baby.s.Reference.Key.),"M:/Personal/Gao Le/DIAMOND_crosscheck/GDM_baby_46404.csv")


#step6.1merge pgdm dx rx and lab result
# mm_pgdm_dx <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx.rds")
# mm_lab_pgdm <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_lab_pgdm.rds")
# mm_anti_drug_new2 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_anti_drug_new2.rds")
mm_pgdm_dx_col2 <- mm_pgdm_dx %>% 
  select(Reference.Key.,Reference.Date.) %>% 
  `colnames<-` (c("Reference.Key.","pgdm_date"))
mm_pgdm_lab_col2 <- mm_lab_pgdm %>% 
  select(Reference.Key.,LIS.Reference.Datetime.)%>% 
  `colnames<-` (c("Reference.Key.","pgdm_date"))
mm_anti_rx_pgdm_col2 <- mm_anti_drug_new2 %>% 
  select(Reference.Key.,Rxst) %>% 
  `colnames<-` (c("Reference.Key.","pgdm_date"))
mm_pgdm_dx_lab_rx <- rbind(mm_pgdm_dx_col2,mm_pgdm_lab_col2,mm_anti_rx_pgdm_col2)
# saveRDS(mm_pgdm_dx_lab_rx,"D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx_lab_rx.rds")
mm_pgdm_dx_lab_rx <-readRDS("D:/data_no_share_20200506/6.DIAMOND/mm_pgdm_dx_lab_rx.rds")


# find the delivery date and pregnancy date of each ppl
pd_date <- mc_01to15_2 %>% 
  select(Reference.Key.,pregancy_date,Maternity.Episode..Delivery.Date..yyyy.mm.dd..) %>% 
  arrange(Reference.Key.,pregancy_date) %>% 
  group_by(Reference.Key.) %>% mutate(id_new=1:length(Reference.Key.)) %>% 
  `colnames<-` (c("Reference.Key.","pd","dd","id_new"))
table(pd_date$id_new)
pddate_wide <- pd_date %>% pivot_wider(names_from = id_new, values_from = c("pd","dd"))


mc_01to15_3 <- merge(mc_01to15_2,pddate_wide,by="Reference.Key.") %>% 
  arrange(Reference.Key.,pregancy_date) %>% 
  group_by(Reference.Key.) %>% mutate(id_new=1:length(Reference.Key.)) %>% ungroup() %>% 
  mutate(pd_2=ifelse(pregancy_date>=pd_2,pd_2,NA),dd_2=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_2,dd_2,NA)) %>% 
  mutate(pd_3=ifelse(pregancy_date>=pd_3,pd_3,NA),dd_3=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_3,dd_3,NA)) %>% 
  mutate(pd_4=ifelse(pregancy_date>=pd_4,pd_4,NA),dd_4=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_4,dd_4,NA)) %>% 
  mutate(pd_5=ifelse(pregancy_date>=pd_5,pd_5,NA),dd_5=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_5,dd_5,NA)) %>% 
  mutate(pd_6=ifelse(pregancy_date>=pd_6,pd_6,NA),dd_6=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_6,dd_6,NA)) %>% 
  mutate(pd_7=ifelse(pregancy_date>=pd_7,pd_7,NA),dd_7=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_7,dd_7,NA)) %>% 
  mutate(pd_8=ifelse(pregancy_date>=pd_8,pd_8,NA),dd_8=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_8,dd_8,NA)) %>% 
  mutate(pd_9=ifelse(pregancy_date>=pd_9,pd_9,NA),dd_9=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_9,dd_9,NA)) %>% 
  mutate(pd_2=as.Date(pd_2,origin="1970-01-01"),dd_2=as.Date(dd_2,origin="1970-01-01")) %>% 
  mutate(pd_3=as.Date(pd_3,origin="1970-01-01"),dd_3=as.Date(dd_3,origin="1970-01-01")) %>% 
  mutate(pd_4=as.Date(pd_4,origin="1970-01-01"),dd_4=as.Date(dd_4,origin="1970-01-01")) %>% 
  mutate(pd_5=as.Date(pd_5,origin="1970-01-01"),dd_5=as.Date(dd_5,origin="1970-01-01")) %>% 
  mutate(pd_6=as.Date(pd_6,origin="1970-01-01"),dd_6=as.Date(dd_6,origin="1970-01-01")) %>% 
  mutate(pd_7=as.Date(pd_7,origin="1970-01-01"),dd_7=as.Date(dd_7,origin="1970-01-01")) %>% 
  mutate(pd_8=as.Date(pd_8,origin="1970-01-01"),dd_8=as.Date(dd_8,origin="1970-01-01")) %>% 
  mutate(pd_9=as.Date(pd_9,origin="1970-01-01"),dd_9=as.Date(dd_9,origin="1970-01-01"))


#combine with pgdm dx rx and lab data, find out pgdm
mc_01to15_4 <- merge(mc_01to15_3,mm_pgdm_dx_lab_rx,by="Reference.Key.",all.x = T)

mc_01to15_4$pgdm <- 0

mc_01to15_4$pgdm[which(mc_01to15_4$id_new==1&
                         mc_01to15_4$pgdm_date<mc_01to15_4$pd_1)] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==2&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==3&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==4&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==5&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_4&mc_01to15_4$pgdm_date<mc_01to15_4$pd_5)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==6&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_4&mc_01to15_4$pgdm_date<mc_01to15_4$pd_5)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_5&mc_01to15_4$pgdm_date<mc_01to15_4$pd_6)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==7&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_4&mc_01to15_4$pgdm_date<mc_01to15_4$pd_5)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_5&mc_01to15_4$pgdm_date<mc_01to15_4$pd_6)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_6&mc_01to15_4$pgdm_date<mc_01to15_4$pd_7)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==8&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_4&mc_01to15_4$pgdm_date<mc_01to15_4$pd_5)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_5&mc_01to15_4$pgdm_date<mc_01to15_4$pd_6)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_6&mc_01to15_4$pgdm_date<mc_01to15_4$pd_7)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_7&mc_01to15_4$pgdm_date<mc_01to15_4$pd_8)))] <- 1
mc_01to15_4$pgdm[which(mc_01to15_4$id_new==9&
                         (mc_01to15_4$pgdm_date<mc_01to15_4$pd_1|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_1&mc_01to15_4$pgdm_date<mc_01to15_4$pd_2)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_2&mc_01to15_4$pgdm_date<mc_01to15_4$pd_3)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_3&mc_01to15_4$pgdm_date<mc_01to15_4$pd_4)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_4&mc_01to15_4$pgdm_date<mc_01to15_4$pd_5)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_5&mc_01to15_4$pgdm_date<mc_01to15_4$pd_6)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_6&mc_01to15_4$pgdm_date<mc_01to15_4$pd_7)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_7&mc_01to15_4$pgdm_date<mc_01to15_4$pd_8)|
                            (mc_01to15_4$pgdm_date>mc_01to15_4$dd_8&mc_01to15_4$pgdm_date<mc_01to15_4$pd_9)))] <- 1




#keep once pgdm
mc_01to15_5 <- mc_01to15_4 %>% 
  arrange(id,-pgdm,pgdm_date) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() %>% 
  select(-pd_1,-pd_2,-pd_3,-pd_4,-pd_5,-pd_6,-pd_7,-pd_8,-pd_9,-dd_1,-dd_2,-dd_3,-dd_4,-dd_5,-dd_6,-dd_7,-dd_8,-dd_9)
table(mc_01to15_5$pgdm) #2573
write.csv(mc_01to15_5 %>% filter(pgdm==1) %>% select(Reference.Key.,Baby.s.Reference.Key.),"M:/Personal/Gao Le/DIAMOND_crosscheck/pgdm_2573.csv")

#combine pgdm and gdm as mdm
mc_01to15_6 <- mc_01to15_5 %>% mutate(dm=pgdm+gdm)
mc_01to15_6$dm1 <- NA
mc_01to15_6$dm1[which(mc_01to15_6$pgdm==0&mc_01to15_6$gdm==0)] <- " Non DM"
mc_01to15_6$dm1[which(mc_01to15_6$pgdm==1)] <- "PGDM"
mc_01to15_6$dm1[which(mc_01to15_6$pgdm==0&mc_01to15_6$gdm==1)] <- "GDM"
mc_01to15_6$dm1 <- as.factor(mc_01to15_6$dm1)
table(mc_01to15_6$dm1)


mc_01to15_6$dm2 <- NA
mc_01to15_6$dm2 <- ifelse(mc_01to15_6$dm1==" Non DM"," Non DM","MDM")
table(mc_01to15_6$dm2)
# saveRDS(mc_01to15_6,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_6.rds")


#summarise the date of dm, categarise to BFP<=-91 -90<=BFP<=- 0<=AFP<=90 91<AFP<=180 and AFP>180, reference non MDM(0)
mc_01to15_7 <- mc_01to15_6
mc_01to15_7$mdm_period <- NA
mc_01to15_7$mdm_period[which(mc_01to15_7$dm1==" Non DM")] <- "  Non DM"
mc_01to15_7$mdm_period[which((mc_01to15_7$dm1=="PGDM")&(mc_01to15_7$pregancy_date-mc_01to15_7$pgdm_date+1>90))] <- "0_90BFP"
mc_01to15_7$mdm_period[which((mc_01to15_7$dm1=="PGDM")&(mc_01to15_7$pregancy_date-mc_01to15_7$pgdm_date+1>0)&(mc_01to15_7$pregancy_date-mc_01to15_7$pgdm_date+1<=90))] <- "1_0_90BFP"
mc_01to15_7$mdm_period[which((mc_01to15_7$dm1=="GDM")&(mc_01to15_7$gdm_date-mc_01to15_7$pregancy_date+1<=90))] <- "2_90AFP"
mc_01to15_7$mdm_period[which((mc_01to15_7$dm1=="GDM")&(mc_01to15_7$gdm_date-mc_01to15_7$pregancy_date+1>90)&(mc_01to15_7$gdm_date-mc_01to15_7$pregancy_date+1<=180))] <- "3_180AFP"
mc_01to15_7$mdm_period[which((mc_01to15_7$dm1=="GDM")&(mc_01to15_7$gdm_date-mc_01to15_7$pregancy_date+1>180))] <- "4_up180AFP"
table(mc_01to15_7$mdm_period)
# saveRDS(mc_01to15_7,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_7.rds")



#find T1DM AND T2DM
mm_pgdm_dx_t1 <- mm_pgdm_dx %>% 
  filter(grepl('^250.01|^250.03|^250.11|^250.13|^250.21|^250.23|^250.31|^250.33|^250.41|^250.43|^250.51|^250.53|^250.61|^250.63|^250.71|^250.73|^250.81|^250.83|^250.91|^250.93',All.Diagnosis.Code..ICD9..)==T)
mm_pgdm_rx_insulin <- mm_anti_drug %>% 
  arrange(Reference.Key.,ymd(Rxst)) %>% 
  filter(duplicated(Reference.Key.)==F,
         grepl("INSULIN",Drug.Name.)==T)

mc_01to15_7$pgdm_type <- NA
mc_01to15_7$pgdm_type[which(mc_01to15_7$dm1=="PGDM"&(mc_01to15_7$Reference.Key.%in%mm_pgdm_rx_insulin$Reference.Key.|mc_01to15_7$Reference.Key.%in%mm_pgdm_dx_t1$Reference.Key.))] <- "T1DM"
mc_01to15_7$pgdm_type[which(mc_01to15_7$dm1=="PGDM"&is.na(mc_01to15_7$pgdm_type))] <- "T2DM"
mc_01to15_7$pgdm_type[which(mc_01to15_7$dm1=="GDM")] <- "GDM"
mc_01to15_7$pgdm_type[which(mc_01to15_7$dm1==" Non DM")] <- " Non DM"

table(mc_01to15_7$pgdm_type)

#findout gdm treated 
mc_01to15_7_1 <- merge(mc_01to15_7,mm_anti_drug_new2 %>% select(Reference.Key.,Rxst),by="Reference.Key.",all.x = T) %>% 
  mutate(gdm_treat=ifelse(!is.na(Rxst)&Rxst>=pregancy_date&Rxst<=Maternity.Episode..Delivery.Date..yyyy.mm.dd..,1,0)) %>% 
  arrange(id,-gdm_treat) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() %>% 
  select(-Rxst)


#findout t2dm treated
mc_01to15_7_2 <- merge(mc_01to15_7_1,mm_anti_drug_new2 %>% select(Reference.Key.,Rxst),by="Reference.Key.",all.x = T) %>% 
  mutate(t2dm_dp=ifelse(!is.na(Rxst)&Rxst>=pregancy_date&Rxst<=Maternity.Episode..Delivery.Date..yyyy.mm.dd..,1,0)) %>% 
  arrange(id,-t2dm_dp) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() %>% 
  select(-Rxst)

mc_01to15_7_3 <- merge(mc_01to15_7_2,mm_anti_drug_new2 %>% select(Reference.Key.,Rxst),by="Reference.Key.",all.x = T) %>% 
  mutate(t2dm_bp=ifelse(!is.na(Rxst)&Rxst<pregancy_date,1,0)) %>% 
  arrange(id,-t2dm_bp) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() %>% 
  select(-Rxst)

mc_01to15_7_3$treat_t2dm <- NA
mc_01to15_7_3$treat_t2dm[which(mc_01to15_7_3$pgdm_type=="T2DM"&mc_01to15_7_3$t2dm_dp==1&mc_01to15_7_3$t2dm_bp==0)] <- "tod"
mc_01to15_7_3$treat_t2dm[which(mc_01to15_7_3$pgdm_type=="T2DM"&mc_01to15_7_3$t2dm_dp==0&mc_01to15_7_3$t2dm_bp==1)] <- "tob"
mc_01to15_7_3$treat_t2dm[which(mc_01to15_7_3$pgdm_type=="T2DM"&mc_01to15_7_3$t2dm_dp==1&mc_01to15_7_3$t2dm_bp==1)] <- "tb"
mc_01to15_7_3$treat_t2dm[which(mc_01to15_7_3$pgdm_type=="T2DM"&mc_01to15_7_3$t2dm_dp==0&mc_01to15_7_3$t2dm_bp==0)] <- "non_treat"
table(mc_01to15_7_3$treat_t2dm)




#step7: find outcome in children
#step7.1 read raw data_MPH
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/1.MPH"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})

All1 <- lapply(All, function(x) x[c('Reference Key\n','Dispensing Date (yyyy-mm-dd)\n')])
mph <- do.call(rbind.data.frame, All1)
colnames(mph) <- make.names(colnames(mph)) #875008
mph1 <- mph %>% `colnames<-` (c("Reference.Key.","Reference.Date.")) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))
length(unique(mph1$Reference.Key.))#43841

#step7.2 read raw data_ADHD
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/2.ADHD"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})
adhd <- do.call(rbind.data.frame, All)
colnames(adhd) <- make.names(colnames(adhd)) #52928
adhd1 <- adhd %>% select(Reference.Key.,Reference.Date.) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))
length(unique(adhd1$Reference.Key.))#39403

#step7.3 combine MPH and ADHD data and select first episode
mph_adhd <- rbind(mph1,adhd1) %>% 
  arrange(Reference.Key.,Reference.Date.) %>% 
  filter(duplicated(Reference.Key.)==F) %>% 
  `colnames<-` (c('Reference.Key.','adhd_date'))#51599
# saveRDS(mph_adhd,"D:/data_no_share_20200506/6.DIAMOND/mph_adhd.rds")
mph_adhd <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mph_adhd.rds")
write.csv(mph_adhd,"M:/Personal/Gao Le/DIAMOND_crosscheck/mph_adhd_combine.csv")


#step7.5 merge mother and baby data
mc_01to15_9 <- merge(mc_01to15_7_3,mph_adhd,by.x = "Baby.s.Reference.Key.",by.y = "Reference.Key.",all.x = T)
# saveRDS(mc_01to15_9,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_9.rds")

#step7.6merge baby death data
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/3.death"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})

bb_death <- do.call(rbind.data.frame, All)
colnames(bb_death) <- make.names(colnames(bb_death)) #535925
length(unique(bb_death$Reference.Key.)) #535925
# saveRDS(bb_death,"M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/bb_death.rds")
bb_death <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/bb_death.rds")

bb_death1 <- bb_death %>% 
  select(Reference.Key.,Date.of.Registered.Death.) %>% 
  mutate(Date.of.Registered.Death.=ymd(Date.of.Registered.Death.)) %>% 
  `colnames<-` (c("Baby.s.Reference.Key.","bb_dod"))
mc_01to15_9_1 <- merge(mc_01to15_9,bb_death1,by="Baby.s.Reference.Key.",all.x = T) %>% 
  mutate(obs_end=pmin(adhd_date, ymd('2020-12-31'), bb_dod, na.rm= T))

mc_01to15_10 <- mc_01to15_9_1 %>% 
  mutate(child.adhd=ifelse(is.na(adhd_date),0,1)) %>% 
  mutate(time=obs_end-Maternity.Episode..Delivery.Date..yyyy.mm.dd..+1) %>% 
  mutate(y=as.numeric(year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)))
# saveRDS(mc_01to15_10,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_10.rds")
mc_01to15_10 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_10.rds")

# #calculate the ADHD incidence information (rerun following code)
# #total cohort
# table(year(mc_01to15_10$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))
# y <- mc_01to15_10 %>% filter(child.adhd==1) %>% select(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)
# table(year(y$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))
# #non dm
# y <- mc_01to15_10 %>% filter(child.adhd==1,dm1==" Non DM") %>% select(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)
# table(year(y$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))
# #any dm
# y <- mc_01to15_10 %>% filter(child.adhd==1,dm1!=" Non DM") %>% select(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)
# table(year(y$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))
# #gdm
# y <- mc_01to15_10 %>% filter(child.adhd==1,dm1=="GDM") %>% select(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)
# table(year(y$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))
# #pgdm
# y <- mc_01to15_10 %>% filter(child.adhd==1,dm1=="PGDM") %>% select(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)
# table(year(y$Maternity.Episode..Delivery.Date..yyyy.mm.dd..))




# #8. calculate the baseline character
# #8.1 maternal age
mc_01to15_11 <- mc_01to15_10 %>% mutate(m_age=as.numeric((Maternity.Episode..Delivery.Date..yyyy.mm.dd..-Date.of.Birth..yyyy.mm.dd..+1)/365.25))
# #non_dm group
# mean(mc_01to15_11 %>% filter(dm1==" Non DM") %>% filter(!is.na(m_age))%>% pull(m_age))#31.26693
# sd(mc_01to15_11 %>% filter(dm1==" Non DM") %>% filter(!is.na(m_age))%>% pull(m_age))#5.004657
# #any mdm group
# mean(mc_01to15_11 %>% filter(dm1!=" Non DM") %>% filter(!is.na(m_age))%>% pull(m_age))#33.39326
# sd(mc_01to15_11 %>% filter(dm1!=" Non DM") %>% filter(!is.na(m_age))%>% pull(m_age))#4.892947
# #gdm group
# mean(mc_01to15_11 %>% filter(dm1=="GDM") %>% filter(!is.na(m_age))%>% pull(m_age))#33.36918
# sd(mc_01to15_11 %>% filter(dm1=="GDM") %>% filter(!is.na(m_age))%>% pull(m_age))#4.899491
# #pgdm group
# mean(mc_01to15_11 %>% filter(dm1=="PGDM") %>% filter(!is.na(m_age))%>% pull(m_age))#34.27699
# sd(mc_01to15_11 %>% filter(dm1=="PGDM") %>% filter(!is.na(m_age))%>% pull(m_age))#4.560879
# 
# #8.2 delivery year
# #non_dm group
# y <- mc_01to15_11 %>% filter(dm1==" Non DM") %>% mutate(y=year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)) %>% select(y)
# table(y$y)
# #any mdm group
# y <- mc_01to15_11 %>% filter(dm1!=" Non DM") %>% mutate(y=year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)) %>% select(y)
# table(y$y)
# #gdm group
# y <- mc_01to15_11 %>% filter(dm1=="GDM") %>% mutate(y=year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)) %>% select(y)
# table(y$y)
# #pgdm group
# y <- mc_01to15_11 %>% filter(dm1=="PGDM") %>% mutate(y=year(Maternity.Episode..Delivery.Date..yyyy.mm.dd..)) %>% select(y)
# table(y$y)


#8.3 SES status
table(mc_01to15_11$District.of.Residence.on.Latest.Selected.Encounter..district..)
mc_01to15_11$ses <- NA
mc_01to15_11$ses[which(mc_01to15_11$District.of.Residence.on.Latest.Selected.Encounter..district.. %in% c("CENTRAL & WESTE","WANCHAI","SAI KUNG excl.","TSEUNG KWAN O","EASTERN"))] <- "HIGH"
mc_01to15_11$ses[which(mc_01to15_11$District.of.Residence.on.Latest.Selected.Encounter..district.. %in% c("SOUTHERN","KOWLOON CITY","MONGKOK","YAU TSIM","TSUEN WAN","TAI PO"))] <- "Midium1"
mc_01to15_11$ses[which(mc_01to15_11$District.of.Residence.on.Latest.Selected.Encounter..district.. %in% c("SHATIN","OTHERS","ISLANDS excl. N","NORTH LANTAU","YUEN LONG","NORTH"))] <- "Midium2"
mc_01to15_11$ses[which(mc_01to15_11$District.of.Residence.on.Latest.Selected.Encounter..district.. %in% c("WONG TAI SIN","TUEN MUN","SHAM SHUI PO","KWAI TSING","KWUN TONG"))] <- " LOW"
# #LOW group
# y <- mc_01to15_11 %>% filter(ses==" LOW") %>% select(dm1)
# table(y$dm1)
# #medium2 group
# y <- mc_01to15_11 %>% filter(ses=="Midium2") %>% select(dm1)
# table(y$dm1)
# #medium1 group
# y <- mc_01to15_11 %>% filter(ses=="Midium1") %>% select(dm1)
# table(y$dm1)
# #HIGH group
# y <- mc_01to15_11 %>% filter(ses=="HIGH") %>% select(dm1)
# table(y$dm1)
saveRDS(mc_01to15_11,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11.rds")
mc_01to15_11 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11.rds")



#8.4 medication use function
rx_code <- list(c("antihypertensive","^2.[24]|^2.5.[1245]|^2.6.2"),
                c("antipsychotics","^4.2"),
                c("sedatives","^4.1"),
                c("antidepressants","^4.3"),
                c("antiepileptics","^4.8"),
                c("antiparkinson","^4.9"),
                c("stimulants","^4.4"),
                c("opioids","^4.7.2"))

func_rx <- function(z){
  x <- z[1]
  y <- z[2]
  
  #find the medication using BNF code
  mm_drug <- mc_01to18_mm_allrx1 %>% 
    filter(grepl(y,Therapeutic.Classification..BNF..Principal..)==T) %>% 
    mutate(Dispensing.Date..yyyy.mm.dd..==as.Date(Dispensing.Date..yyyy.mm.dd..))
  #merge rx data with main database, select mother with x rx use before pregnancy date
  mc_01to15_rx <- merge(mc_01to15_11,mm_drug %>% select(Reference.Key.,Dispensing.Date..yyyy.mm.dd..),by="Reference.Key.",all.x = T) %>% 
    mutate(x1=ifelse((!is.na(Dispensing.Date..yyyy.mm.dd..)&Dispensing.Date..yyyy.mm.dd..<=pregancy_date),1,0)) %>% 
    select(-c(Dispensing.Date..yyyy.mm.dd..)) %>% 
    arrange(id,-x1) %>% 
    group_by(id) %>% 
    filter(duplicated(id)==F) %>% ungroup() %>% 
    arrange(id) %>% 
    select(id,x1) %>%
    `colnames<-` (c("id",x))
    
  
  return(mc_01to15_rx)
}
all_mc_01to15_rx <- bind_cols(lapply(rx_code,func_rx))
all_mc_01to15_rx1 <- all_mc_01to15_rx %>% 
  select(-id...3,-id...5,-id...7,-id...9,-id...11,-id...13,-id...15)
mc_01to15_11_1 <- merge(mc_01to15_11,all_mc_01to15_rx1,by.x="id",by.y = "id...1")

triptan_drug <- mc_01to18_mm_allrx1 %>% 
  filter(grepl("TRIPTAN",Drug.Name.)==T) %>% 
  mutate(Dispensing.Date..yyyy.mm.dd..==as.Date(Dispensing.Date..yyyy.mm.dd..))%>% 
  select(Reference.Key.,Dispensing.Date..yyyy.mm.dd..)
mc_01to15_11_1_triptan <- merge(mc_01to15_11_1,triptan_drug,by="Reference.Key.",all.x = T) %>% 
  mutate(triptan=ifelse((!is.na(Dispensing.Date..yyyy.mm.dd..)&Dispensing.Date..yyyy.mm.dd..<=pregancy_date),1,0)) %>% 
  select(-c(Dispensing.Date..yyyy.mm.dd..)) %>% 
  arrange(id,-triptan) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup()
  
# test <- teratogenic_drug <- mc_01to18_mm_allrx1 %>%
#   filter(grepl("^8.1|^2.5.5|^6.2.2",Therapeutic.Classification..BNF..Principal..)==T|grepl('WARFARIN|ORFARIN|METHOTREXAT|EMTHEXATE|ISOTRETINOIN|ACNOTIN|NIMEGEN|ORATANE|REDUCAR|ROACCUTANE|MISOPRISTOL|THALIDOMIDE|MYCOPHENOLIC ACID|MYFORTIC|LEFLUNOMIDE|ARAVA|ARAY|LENALIDOMIDE|REVLIMID|POMALIDOMIDE|POMALYST|ERGOT ALKALOIDS',Drug.Name.)) %>%
#   mutate(Dispensing.Date..yyyy.mm.dd..==as.Date(Dispensing.Date..yyyy.mm.dd..)) %>%
#   group_by(Drug.Name.) %>%
#   mutate(ct=length(unique(Reference.Key.)))%>%
#   filter(duplicated(Drug.Name.)==F)
# 


teratogenic_drug <- mc_01to18_mm_allrx1 %>% 
  filter(grepl("^8.1|^6.2.2",Therapeutic.Classification..BNF..Principal..)==T|grepl('WARFARIN|ORFARIN|METHOTREXAT|EMTHEXATE|ISOTRETINOIN|ACNOTIN|NIMEGEN|ORATANE|REDUCAR|ROACCUTANE|MISOPRISTOL|THALIDOMIDE|MYCOPHENOLIC ACID|MYFORTIC|LEFLUNOMIDE|ARAVA|ARAY|LENALIDOMIDE|REVLIMID|POMALIDOMIDE|POMALYST|ERGOT ALKALOIDS',Drug.Name.)) %>% 
  mutate(Dispensing.Date..yyyy.mm.dd..==as.Date(Dispensing.Date..yyyy.mm.dd..))%>% 
  select(Reference.Key.,Dispensing.Date..yyyy.mm.dd..)
mc_01to15_11_1_teratogenic <- merge(mc_01to15_11_1_triptan,teratogenic_drug,by="Reference.Key.",all.x = T) %>% 
  mutate(teratogenic=ifelse((!is.na(Dispensing.Date..yyyy.mm.dd..)&Dispensing.Date..yyyy.mm.dd..<=pregancy_date),1,0)) %>% 
  select(-c(Dispensing.Date..yyyy.mm.dd..)) %>% 
  arrange(id,-teratogenic) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup()
table(mc_01to15_11_1_teratogenic$teratogenic)

#8.5 dx function
dx_code <- list(c("pcos","^256.4"),
                c("anxiety","^293.84|^300"),
                c("depression","^296.[23]|^311|^309.[01]"),
                c("epilepsy","^345"),
                c("migraine","^346|^784.0"),
                c("hypertension","^40[1-5]"),
                c("renal","^58[0-9]|^59[0-3]"),
                c("crohn","^55[56]"),
                c("sleep","^780.5|^347|^327|^307.4"),
                c("scz","^295"),
                c("bp","296.[0145678]"),
                c("personality","^301"),
                c("intellectual","^31[789]"),
                c("development","^315"),
                c("illicit","^304|^305.[2-9]|^648.3|^760.7[235]|^292|^965.0|^967.[0689]|^969.[67]|^970.1|^970.81|^E850.[012]|^E852|^E854.[12]|^E935.[0-2]|^E937|^E939.[67]|^E940.1|^E950.1|^E980.2"),
                c("somke_dx","^305.1|^V15.82|^649.0"),
                c("alcohol_dx","^291|^303|^305.0|^357.5|^425.5|^535.3|^571.0|^571.1|^571.2|^571.3|^980|^V11.3"),
                c("adhd_asd","^314|^299"),
                c("thyroid","^24[0-6]"),
                c("Rheumatoid","^714"),
                c("headache","^339.0")
                #c("obe_dx","^278.00|^278.01|^278.03|^V85.3|^V85.4")
                )

func_dx <- function(z){
  x <- z[1]
  y <- z[2]
  
  #find the dx using icd code
  mm_dx <- mc_01to18_mm_alldx1 %>% 
    filter(grepl(y,All.Diagnosis.Code..ICD9..)==T) %>% 
    mutate(Reference.Date.==as.Date(Reference.Date.))
  #merge dx data with main database, select mother with x dx before pregnancy date
  mc_01to15_dx <- merge(mc_01to15_11,mm_dx %>% select(Reference.Key.,Reference.Date.),by="Reference.Key.",all.x = T) %>% 
    mutate(x2=ifelse((!is.na(Reference.Date.)&Reference.Date.<=pregancy_date),1,0)) %>% 
    select(-c(Reference.Date.)) %>% 
    arrange(id,-x2) %>% 
    group_by(id) %>% 
    filter(duplicated(id)==F) %>% ungroup() %>% 
    arrange(id) %>% 
    select(id,x2) %>%
    `colnames<-` (c("id",x)) %>% 
  
  return(mc_01to15_dx)
}

all_mc_01to15_dx <- bind_cols(lapply(dx_code,func_dx))
all_mc_01to15_dx1 <- all_mc_01to15_dx %>% 
  select(-id...3,-id...5,-id...7,-id...9,-id...11,-id...13,-id...15,-id...17,
         -id...19,-id...21,-id...23,-id...25,-id...27,-id...29,-id...31,
         -id...33,-id...35,-id...37,-id...39,-id...41)
mc_01to15_11_2 <- merge(mc_01to15_11_1_teratogenic,all_mc_01to15_dx1,by="id",by.y = "id...1",all.x = T)
saveRDS(mc_01to15_11_2,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_2.rds")





#8.6 preterm
mc_01to15_11_3 <- mc_01to15_11_2 %>% 
  mutate(preterm=ifelse(Maternity.Episode..Gestation.weeks.1<37,1,0))
saveRDS(mc_01to15_11_3,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_3.rds")
table(mc_01to15_11_3$preterm)#45212
# y <- mc_01to15_11_3 %>% filter(preterm==1) %>% select(dm1)
# table(y$dm1)




#8.22 smoke data and dx
#find smoke information from social nursing patient form
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/2.mother/social_nursing"
filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})
mm_social_nursing <- do.call(rbind.data.frame, All)
colnames(mm_social_nursing) <- make.names(colnames(mm_social_nursing))
mm_social_nursing$PAS..Created.Date..yyyy.mm.dd.. <- as.Date(mm_social_nursing$PAS..Created.Date..yyyy.mm.dd..)
# saveRDS(mm_social_nursing,"M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_social_nursing.rds")
mm_social_nursing <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_social_nursing.rds")

#find out somkers and ex-smokers from social nursing form
mm_social_nursing_s <- mm_social_nursing %>% 
  filter(PAS..Smoking.status.=="Ex-smoker"|PAS..Smoking.status.=="Smoker") %>% 
  select(Reference.Key.,PAS..Created.Date..yyyy.mm.dd..)

mc_01to15_11_4 <- merge(mc_01to15_11_3,mm_social_nursing_s,by="Reference.Key.",all.x = T) %>% 
  mutate(smoke_sn=ifelse((!is.na(PAS..Created.Date..yyyy.mm.dd..)&PAS..Created.Date..yyyy.mm.dd..<=pregancy_date),1,0)) %>% 
  select(-c(PAS..Created.Date..yyyy.mm.dd..)) %>%
  arrange(id,-smoke_sn) %>%
  group_by(id) %>%
  filter(duplicated(id)==F) %>% ungroup() %>% 
  mutate(smoke=ifelse(smoke_sn+somke_dx==0,0,1)) %>% 
  select(-somke_dx,-smoke_sn)
table(mc_01to15_11_4$smoke)#1907
saveRDS(mc_01to15_11_4,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_4.rds")
y <- mc_01to15_11_4 %>% filter(smoke==1) %>% select(dm1)
table(y$dm1)




#8.7 alcohol data and dx
#find alcohol information from social nursing patient form
mm_social_nursing <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_social_nursing.rds")

#find out drinker and ex-drinker from social nursing form
mm_social_nursing_d <- mm_social_nursing %>% 
  filter(PAS..Alcohol.use.status.=="Ex-drinker"|PAS..Alcohol.use.status.=="Drinker") %>% 
  select(Reference.Key.,PAS..Created.Date..yyyy.mm.dd..)

mc_01to15_11_5 <- merge(mc_01to15_11_4,mm_social_nursing_d,by="Reference.Key.",all.x = T) %>% 
  mutate(drink_sn=ifelse((!is.na(PAS..Created.Date..yyyy.mm.dd..)&PAS..Created.Date..yyyy.mm.dd..<=pregancy_date),1,0)) %>% 
  select(-c(PAS..Created.Date..yyyy.mm.dd..)) %>%
  arrange(id,-drink_sn) %>%
  group_by(id) %>%
  filter(duplicated(id)==F) %>% ungroup() %>% 
  mutate(drink=ifelse(drink_sn+alcohol_dx==0,0,1)) %>% 
  select(-drink_sn,-alcohol_dx)
table(mc_01to15_11_5$drink)#3107
saveRDS(mc_01to15_11_5,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_5.rds")
mc_01to15_11_5 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_5.rds")
y <- mc_01to15_11_5 %>% filter(drink==1) %>% select(dm1)
table(y$dm1)



#8.8 BMI data and dx

#find out BMI information from social nursing form
mm_social_nursing_BMI <- mm_social_nursing %>% 
  select(Reference.Key.,PAS..Body.Mass.Index.,PAS..Created.Date..yyyy.mm.dd..) %>% 
  filter(!is.na(PAS..Body.Mass.Index.)) %>% 
  `colnames<-` (c("Reference.Key.","FM_result","FM_date"))
#find BMI information from FM
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/2.mother/BMI"
filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})
mm_FM <- do.call(rbind.data.frame, All)
colnames(mm_FM) <- make.names(colnames(mm_FM))
mm_FM_BMI <- mm_FM %>% 
  filter(FM.Module..Health.Status...Laboratory.Test.Name.=="BMI") %>% 
  mutate(FM.Module..Health.Status...Examination...Investigation.Result.=as.numeric(FM.Module..Health.Status...Examination...Investigation.Result.)) %>% 
  mutate(FM.Module..Health.Status...Examination...Investigation.Result.Date.=as.Date(FM.Module..Health.Status...Examination...Investigation.Result.Date.)) %>% 
  select(Reference.Key.,FM.Module..Health.Status...Examination...Investigation.Result.,FM.Module..Health.Status...Examination...Investigation.Result.Date.) %>% 
  `colnames<-` (c("Reference.Key.","FM_result","FM_date"))
# saveRDS(mm_FM_BMI,"M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_FM_BMI.rds")
mm_FM_BMI <- readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_FM_BMI.rds")

#combine BMI data from social nursing and FM module
mm_bmi_sc_fm <- rbind(mm_social_nursing_BMI,mm_FM_BMI) %>% mutate(FM_result=as.numeric(FM_result))
# saveRDS(mm_bmi_sc_fm,"M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_bmi_sc_fm.rds")
# mm_bmi_sc_fm <-readRDS("M:/Cohort Raw Data (do not edit)/DIAMOND/3.combined/mm_bmi_sc_fm.rds")



#find out underweight from FM module (BMI <19)
mm_FM_underw <- mm_bmi_sc_fm %>% 
  filter(FM_result<19) %>% 
  select(Reference.Key.,FM_date) 

#merge data with long data of pregnancy records
mc_01to15_11_6 <- merge(mc_01to15_11_5,pddate_wide,by="Reference.Key.") %>% 
  arrange(Reference.Key.,pregancy_date) %>% 
  group_by(Reference.Key.) %>% mutate(id_new=1:length(Reference.Key.)) %>% ungroup() %>% 
  mutate(pd_2=ifelse(pregancy_date>=pd_2,pd_2,NA),dd_2=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_2,dd_2,NA)) %>% 
  mutate(pd_3=ifelse(pregancy_date>=pd_3,pd_3,NA),dd_3=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_3,dd_3,NA)) %>% 
  mutate(pd_4=ifelse(pregancy_date>=pd_4,pd_4,NA),dd_4=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_4,dd_4,NA)) %>% 
  mutate(pd_5=ifelse(pregancy_date>=pd_5,pd_5,NA),dd_5=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_5,dd_5,NA)) %>% 
  mutate(pd_6=ifelse(pregancy_date>=pd_6,pd_6,NA),dd_6=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_6,dd_6,NA)) %>% 
  mutate(pd_7=ifelse(pregancy_date>=pd_7,pd_7,NA),dd_7=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_7,dd_7,NA)) %>% 
  mutate(pd_8=ifelse(pregancy_date>=pd_8,pd_8,NA),dd_8=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_8,dd_8,NA)) %>% 
  mutate(pd_9=ifelse(pregancy_date>=pd_9,pd_9,NA),dd_9=ifelse(Maternity.Episode..Delivery.Date..yyyy.mm.dd..>=dd_9,dd_9,NA)) %>% 
  mutate(pd_2=as.Date(pd_2,origin="1970-01-01"),dd_2=as.Date(dd_2,origin="1970-01-01")) %>% 
  mutate(pd_3=as.Date(pd_3,origin="1970-01-01"),dd_3=as.Date(dd_3,origin="1970-01-01")) %>% 
  mutate(pd_4=as.Date(pd_4,origin="1970-01-01"),dd_4=as.Date(dd_4,origin="1970-01-01")) %>% 
  mutate(pd_5=as.Date(pd_5,origin="1970-01-01"),dd_5=as.Date(dd_5,origin="1970-01-01")) %>% 
  mutate(pd_6=as.Date(pd_6,origin="1970-01-01"),dd_6=as.Date(dd_6,origin="1970-01-01")) %>% 
  mutate(pd_7=as.Date(pd_7,origin="1970-01-01"),dd_7=as.Date(dd_7,origin="1970-01-01")) %>% 
  mutate(pd_8=as.Date(pd_8,origin="1970-01-01"),dd_8=as.Date(dd_8,origin="1970-01-01")) %>% 
  mutate(pd_9=as.Date(pd_9,origin="1970-01-01"),dd_9=as.Date(dd_9,origin="1970-01-01"))
  

#combine with underweight data, find out the BMI date within one year before the pregnancy
mc_01to15_11_7 <- merge(mc_01to15_11_6,mm_FM_underw,by="Reference.Key.",all.x = T) %>% 
  mutate(underw_fm=0)

mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==1)&(mc_01to15_11_7$FM_date>=(mc_01to15_11_7$pd_1-365)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_1))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==2)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_2-365,mc_01to15_11_7$dd_1+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_2))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==3)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_3-365,mc_01to15_11_7$dd_2+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_3))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==4)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_4-365,mc_01to15_11_7$dd_3+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_4))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==5)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_5-365,mc_01to15_11_7$dd_4+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_5))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==6)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_6-365,mc_01to15_11_7$dd_5+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_6))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==7)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_7-365,mc_01to15_11_7$dd_6+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_7))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==8)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_8-365,mc_01to15_11_7$dd_7+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_8))] <- 1
mc_01to15_11_7$underw_fm[which((mc_01to15_11_7$id_new==9)&(mc_01to15_11_7$FM_date>=pmax(mc_01to15_11_7$pd_9-365,mc_01to15_11_7$dd_8+1)&mc_01to15_11_7$FM_date<=mc_01to15_11_7$pd_9))] <- 1

#keep once underweight
mc_01to15_11_8 <- mc_01to15_11_7 %>% 
  select(-c(FM_date)) %>%
  arrange(id,-underw_fm) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() 
table(mc_01to15_11_8$underw_fm) #2038

saveRDS(mc_01to15_11_8,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_8.rds")
# mc_01to15_11_8 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_8.rds")



#find out overweight from FM module (25<=BMI <30)
mm_FM_overw <- mm_bmi_sc_fm %>% 
  filter(FM_result<30,FM_result>=25) %>% 
  select(Reference.Key.,FM_date)

mc_01to15_11_9 <- merge(mc_01to15_11_8,mm_FM_overw,by="Reference.Key.",all.x = T) %>% 
  mutate(overw_fm=0)

mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==1)&(mc_01to15_11_9$FM_date>=(mc_01to15_11_9$pd_1-365)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_1))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==2)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_2-365,mc_01to15_11_9$dd_1+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_2))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==3)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_3-365,mc_01to15_11_9$dd_2+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_3))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==4)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_4-365,mc_01to15_11_9$dd_3+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_4))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==5)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_5-365,mc_01to15_11_9$dd_4+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_5))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==6)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_6-365,mc_01to15_11_9$dd_5+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_6))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==7)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_7-365,mc_01to15_11_9$dd_6+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_7))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==8)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_8-365,mc_01to15_11_9$dd_7+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_8))] <- 1
mc_01to15_11_9$overw_fm[which((mc_01to15_11_9$id_new==9)&(mc_01to15_11_9$FM_date>=pmax(mc_01to15_11_9$pd_9-365,mc_01to15_11_9$dd_8+1)&mc_01to15_11_9$FM_date<=mc_01to15_11_9$pd_9))] <- 1

#keep once overweight
mc_01to15_11_10 <- mc_01to15_11_9 %>% 
  select(-c(FM_date)) %>%
  arrange(id,-overw_fm) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() 
table(mc_01to15_11_10$overw_fm) #1698

saveRDS(mc_01to15_11_10,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_10.rds")



#find out obesity using icd code
mm_obe_dx <- mc_01to18_mm_alldx1 %>% 
  filter(grepl("^278.00|^278.01|^278.03|^V85.3|^V85.4",All.Diagnosis.Code..ICD9..)==T) %>% 
  mutate(Reference.Date.==as.Date(Reference.Date.)) %>% 
  select(Reference.Key.,Reference.Date.)

#find out obesity from FM module (BMI >= 30)
mm_FM_obes <- mm_bmi_sc_fm %>% 
  filter(FM_result>=30) %>% 
  select(Reference.Key.,FM_date)
colnames(mm_obe_dx) <- colnames(mm_FM_obes)
mm_obe_fm_dx <- rbind(mm_FM_obes,mm_obe_dx)

mc_01to15_11_11 <- merge(mc_01to15_11_10,mm_obe_fm_dx,by="Reference.Key.",all.x = T) %>% 
  mutate(obe=0)

mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==1)&(mc_01to15_11_11$FM_date>=(mc_01to15_11_11$pd_1-365)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_1))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==2)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_2-365,mc_01to15_11_11$dd_1+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_2))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==3)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_3-365,mc_01to15_11_11$dd_2+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_3))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==4)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_4-365,mc_01to15_11_11$dd_3+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_4))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==5)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_5-365,mc_01to15_11_11$dd_4+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_5))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==6)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_6-365,mc_01to15_11_11$dd_5+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_6))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==7)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_7-365,mc_01to15_11_11$dd_6+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_7))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==8)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_8-365,mc_01to15_11_11$dd_7+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_8))] <- 1
mc_01to15_11_11$obe[which((mc_01to15_11_11$id_new==9)&(mc_01to15_11_11$FM_date>=pmax(mc_01to15_11_11$pd_9-365,mc_01to15_11_11$dd_8+1)&mc_01to15_11_11$FM_date<=mc_01to15_11_11$pd_9))] <- 1

#keep once underweight
mc_01to15_11_12 <- mc_01to15_11_11 %>% 
  select(-c(FM_date)) %>%
  arrange(id,-obe) %>% 
  group_by(id) %>% 
  filter(duplicated(id)==F) %>% ungroup() 
table(mc_01to15_11_12$obe) #662

saveRDS(mc_01to15_11_12,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_12.rds")
# mc_01to15_11_12 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_12.rds")


mc_01to15_11_12$BMI_fac <- NA
mc_01to15_11_12$BMI_fac[which(mc_01to15_11_12$obe==1)] <- "obesity"
mc_01to15_11_12$BMI_fac[which(mc_01to15_11_12$underw_fm==1&mc_01to15_11_12$obe==0)] <- "underweight"
mc_01to15_11_12$BMI_fac[which(mc_01to15_11_12$overw_fm==1&mc_01to15_11_12$obe==0&mc_01to15_11_12$underw_fm==0)] <- "overweight"
# mc_01to15_11_12$BMI_fac[which(mc_01to15_11_12$overw_fm==1&mc_01to15_11_12$underw_fm==0)] <- "overweight"
mc_01to15_11_12$BMI_fac[which(mc_01to15_11_12$overw_fm==0&mc_01to15_11_12$obe==0&mc_01to15_11_12$underw_fm==0)] <- " normal"
table(mc_01to15_11_12$BMI_fac)
#change some covariates to factor for cox regression
median(mc_01to15_11_12$m_age[which(!is.na(mc_01to15_11_12$m_age))]) #31.90691
mc_01to15_11_13 <- mc_01to15_11_12 %>%
  arrange(Reference.Key.,id_new) %>% 
  mutate(preg_diff=ifelse(id_new==1,1000,pregancy_date-lag(pregancy_date))) %>% 
  mutate(BMI_fac=ifelse(preg_diff==0,lag(BMI_fac),BMI_fac)) %>% 
  mutate(ses =as.factor(ses),y = as.factor(y)) %>% 
  mutate(m_age1 = ifelse(!is.na(m_age),m_age,31.90691)) %>% 
  select(-pd_1,-pd_2,-pd_3,-pd_4,-pd_5,-pd_6,-pd_7,-pd_8,-pd_9,-dd_1,-dd_2,-dd_3,-dd_4,-dd_5,-dd_6,-dd_7,-dd_8,-dd_9)
table(mc_01to15_11_13$BMI_fac,mc_01to15_11_13$dm2)

# saveRDS(mc_01to15_11_13,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_13.rds")
mc_01to15_11_13 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_13.rds")




# #baseline of comorbidities
# d_code <- c("antihypertensive","antipsychotics","sedatives","antidepressants","antiepileptics","antiparkinson",
#             "stimulants","opioids","triptan","teratogenic","pcos","anxiety","depression","epilepsy","migraine",
#             "hypertension","renal","crohn","sleep","scz","bp","personality","intellectual","development",
#             "illicit","adhd_asd","thyroid","Rheumatoid","headache","preterm","smoke","drink")
# 
# fun_com <- function(y){
#   n <- data.frame(mc_01to15_11_13[,c("dm1",y)])
#   n1 <- n[n[ ,2] == 1, ]
#   n2 <- data.frame(table(n1$dm1)) %>% 
#     mutate(disease=y)
#   return(n2)
# }
# sum_com <- bind_rows(lapply(d_code, fun_com))


#save data for pssw non_MDM vs MDM
names(mc_01to15_11_13)[11]<-paste("dob")
names(mc_01to15_11_13)[25]<-paste("sex")


#step7: find outcome in children
#step7.1 read raw data_MPH
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/1.MPH"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})

All1 <- lapply(All, function(x) x[c('Reference Key\n','Dispensing Date (yyyy-mm-dd)\n')])
mph <- do.call(rbind.data.frame, All1)
colnames(mph) <- make.names(colnames(mph)) #875008
mph1 <- mph %>% `colnames<-` (c("Reference.Key.","Reference.Date.")) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))
length(unique(mph1$Reference.Key.))#43841

#step7.2 read raw data_ADHD
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/2.ADHD"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})
adhd <- do.call(rbind.data.frame, All)
colnames(adhd) <- make.names(colnames(adhd)) #52928
adhd1 <- adhd %>% select(Reference.Key.,Reference.Date.) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))
length(unique(adhd1$Reference.Key.))#39403

#step7.3 combine MPH and ADHD data and select first episode
mph_adhd <- rbind(mph1,adhd1) %>% 
  arrange(Reference.Key.,Reference.Date.) %>% 
  `colnames<-` (c('Reference.Key.','adhd_date'))#51599
mph_adhd_3 <- merge(mph_adhd,mc_01to15_11_13 %>% select(Baby.s.Reference.Key.,dob),by.x = "Reference.Key.",by.y = "Baby.s.Reference.Key.") %>% 
  filter(adhd_date-dob>365.25*3) %>% 
  arrange(Reference.Key.,adhd_date) %>% 
  filter(duplicated(Reference.Key.)==F)



mc_01to15_11_14 <- merge(mc_01to15_11_13 %>% select(-adhd_date,-time),mph_adhd_3 %>% select(-dob),by.x = "Baby.s.Reference.Key.",by.y = "Reference.Key.",all.x = T) %>% 
  mutate(child.adhd=ifelse(is.na(adhd_date),0,1),
         obs_end=pmin(adhd_date, ymd('2020-12-31'), bb_dod, na.rm= T)) %>% 
  mutate(time=obs_end-dob+1) %>% 
  mutate(y=as.numeric(year(dob)))
saveRDS(mc_01to15_11_14,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_14.rds")
mc_01to15_11_14 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_14.rds")



#step8.1 cox assumption
res.cox <- coxph(Surv(time, child.adhd) ~ dm1+y, data =  mc_01to15_11_14)
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
ggcoxzph(test.ph)
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

a <- ggsurvplot(
  fit = survfit(Surv(as.numeric(time)/365.25, child.adhd) ~ dm1, data =  mc_01to15_11_14), 
  fun = "cumhaz",
  xlab = "Years after birth", 
  ylab = "Cumulative hazard",
  ggtheme = theme_classic2(base_size=20)
)

a$plot + geom_vline(xintercept = 6)+
  annotate("text", size=6,x=6, y=0.08, label="age=6")
  



 # options(digits = 12)
mc_01to15_sas <- mc_01to15_11_14 %>% 
  mutate(m_age1=as.numeric(dob-Date.of.Birth..yyyy.mm.dd..)) %>% 
  filter(!is.na(Baby.s.Reference.Key.),
         No..of.Abortion.==0,
         No..of.Perinatal.Deaths.==0,
         sex!="U",
         !is.na(Maternity.Episode..Gestation.weeks.),
         !is.na(Date.of.Birth..yyyy.mm.dd..)) %>% 
  mutate(dm_sas=ifelse(dm2=="MDM",1,0),
         # m_age1=m_age2,
         child_adhd=child.adhd,
         reference_key=Reference.Key.,
         multi_pre=Maternity.Episode..OBS.Hx...Cx..Multiple.Pregnancy..Y.N..,
         parity=Maternity.Episode..OBS.Hx...Cx..Parity.,
         insti=Institution..IPAS..,
         b_id=as.numeric(Baby.s.Reference.Key.)) %>% 
  select(reference_key,dm_sas,child_adhd,time,m_age1,ses,y,antihypertensive,antipsychotics,sedatives,
         antidepressants,antiepileptics,antiparkinson,stimulants,opioids,triptan,teratogenic,pcos,anxiety,
         depression,epilepsy,migraine,hypertension,renal,crohn,sleep,scz,bp,personality,intellectual,
         development,illicit,adhd_asd,thyroid,Rheumatoid,headache,smoke,drink,BMI_fac,dm1,mdm_period,
         treat_t2dm,gdm_treat,pgdm_type,sex,dob,preterm,multi_pre,parity,insti,b_id,adhd_date,bb_dod) %>% 
  mutate(dm1=ifelse(dm1==" Non DM","non_dm",as.character(dm1)),
         mdm_period=ifelse(mdm_period==" Non DM","non_dm",as.character(mdm_period)),
         pgdm_type=ifelse(pgdm_type==" Non DM","non_dm",as.character(pgdm_type)))

mc_01to15_sas$parity[which(mc_01to15_sas$parity==">= 5")] <- 5
# ad_ter <- read.csv("M:/Personal/Gao Le/DIAMOND_crosscheck/teratogenic_ac.csv")
# gl_ter <- mc_01to15_sas %>% 
#   filter(dm_sas==1,teratogenic==1)
# 
# ter_dif <- ad_ter %>% filter(!b_id%in%gl_ter$Baby.s.Reference.Key.)
# ter_dif <- gl_ter %>% filter(!Baby.s.Reference.Key.%in%ad_ter$b_id) %>% 
#   select(Baby.s.Reference.Key.)
# 
# ad_obe <- read.csv("M:/Personal/Gao Le/DIAMOND_crosscheck/cohort_5.csv")
# table(ad_obe$MDM,ad_obe$Obesity)
# table(ad_obe$MDM,ad_obe$Normal)
# table(ad_obe$MDM,ad_obe$Overweight)
# table(ad_obe$MDM,ad_obe$Underweight)
# ad_obe1 <- ad_obe %>%
#   filter(GDM_treat==1)
# gl_obe <- mc_01to15_11_13 %>%
#   filter(dm1=="GDM",gdm_treat==1)
# ad_gl <- merge(ad_obe1,gl_obe,by.x = 'b_id',by.y = 'Baby.s.Reference.Key.') %>% 
#   mutate(dif=as.numeric(time)-as.numeric(studydur))


# obe_dif <- ad_obe1 %>% filter(!b_id%in%gl_obe$Baby.s.Reference.Key.)
# ter_dif <- gl_obe %>% filter(!Baby.s.Reference.Key.%in%ad_obe1$b_id) %>%
#   select(Baby.s.Reference.Key.,Reference.Key.)
# 

mc_01to15_sas$mdm_period_1 <- NA
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="0_90BFP")] <- 1
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="1_0_90BFP")] <- 2
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="2_90AFP")] <- 3
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="3_180AFP")] <- 4
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="4_up180AFP")] <- 5
mc_01to15_sas$mdm_period_1[which(is.na(mc_01to15_sas$mdm_period_1))] <- 0

mc_01to15_sas$dm1_cat <- NA
mc_01to15_sas$dm1_cat[which(mc_01to15_sas$dm1=="PGDM")] <- 1
mc_01to15_sas$dm1_cat[which(mc_01to15_sas$dm1=="GDM")] <- 2
mc_01to15_sas$dm1_cat[which(is.na(mc_01to15_sas$dm1_cat))] <- 0


mc_01to15_sas$mdm_period_sen <- NA
mc_01to15_sas$mdm_period_sen[which(mc_01to15_sas$mdm_period=="0_90BFP")] <- 0
mc_01to15_sas$mdm_period_sen[which(is.na(mc_01to15_sas$mdm_period_sen))] <- 1

mc_01to15_sas$dm1_catsen <- NA
mc_01to15_sas$dm1_catsen[which(mc_01to15_sas$dm1=="PGDM")] <- 0
mc_01to15_sas$dm1_catsen[which(is.na(mc_01to15_sas$dm1_catsen))] <- 1

mc_01to15_sas$treat_t2dm1 <- NA
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tob")] <- 1
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tod")] <- 2
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tb")] <- 3
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="non_treat")] <- 0


mc_01to15_sas1 <- mc_01to15_sas %>% 
  mutate(nor=ifelse(BMI_fac==" normal",1,0),
         obe=ifelse(BMI_fac=="obesity",1,0),
         over=ifelse(BMI_fac=="overweight",1,0),
         under=ifelse(BMI_fac=="underweight",1,0),
         b_id_new=1:535016)

test <- mc_01to15_sas1 %>% 
  filter(as.numeric(as.character(y))<=2014) %>% 
  filter(is.na(adhd_date)|adhd_date>=(dob+years(6))) %>% 
  mutate(obs_end=pmin(adhd_date, ymd('2020-12-31'), bb_dod, na.rm= T),
         time1=obs_end-(dob+years(6))+1)
  
# PH assumption
res.cox <- coxph(Surv(time1, child_adhd) ~ dm_sas+m_age1+y+ses+sex+multi_pre+parity+insti+antihypertensive+antipsychotics+sedatives+antidepressants+antiepileptics+antiparkinson+stimulants+opioids+triptan +teratogenic+pcos+anxiety+depression+epilepsy+migraine +hypertension+renal+crohn+sleep+scz+bp+personality +intellectual+development+illicit+adhd_asd+thyroid +Rheumatoid+headache+smoke+drink+nor+obe+over+under, data =  test)
res.cox <- coxph(Surv(time1, child_adhd) ~ dm_sas+m_age1+y+ses+sex+multi_pre+parity+insti, data =  mc_01to15_sas1%>% filter(as.numeric(as.character(y))<=2014) )
res.cox <- coxph(Surv(time1, child_adhd) ~ dm_sas,data=test)
summary(res.cox)
test.ph <- cox.zph(res.cox)
test.ph
ggcoxzph(test.ph)
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())


nor obe over under
write.csv(mc_01to15_sas1,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas.csv")
write.csv(mc_01to15_sas1 %>% filter(as.numeric(as.character(y))<=2014) %>% mutate(time1=time-2190),
          "D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_2014.csv")
write.csv(mc_01to15_sas1 %>% filter(as.numeric(as.character(y))<=2011) %>% mutate(time1=time-2190),
          "D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_2011.csv")
write_sas(mc_01to15_sas1 %>% filter(as.numeric(as.character(y))<=2014) %>% mutate(time1=time-2190),
          "D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_2014.sas7bdat")
write.csv(mc_01to15_sas1 %>% filter(sex=="F"),"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_F.csv")
write.csv(mc_01to15_sas1 %>% filter(sex=="M"),"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_M.csv")






#check no. of different drug users
#step7: find outcome in children
mc_01to15_11_13 <- readRDS("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_11_13.rds")
names(mc_01to15_11_13)[11]<-paste("dob")
names(mc_01to15_11_13)[25]<-paste("sex")

#step7.1 read raw data_MPH
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/1.MPH"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})

All1 <- lapply(All, function(x) x[c('Reference Key\n','Dispensing Date (yyyy-mm-dd)\n','Drug Name\n',"Drug Item Code\n")])
mph <- do.call(rbind.data.frame, All1)
colnames(mph) <- make.names(colnames(mph)) #875008


#find adhd count
table(mph$Drug.Name.)
table(mph$Drug.Item.Code.)
mph1 <- mph %>% 
  filter(grepl("METH|ATOM",Drug.Name.)==T|grepl("METH|ATOM|DEXM",Drug.Item.Code.)==T) %>% 
  select(Reference.Key.,Dispensing.Date..yyyy.mm.dd..) %>% 
  `colnames<-` (c("Reference.Key.","Reference.Date.")) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))


#step7.2 read raw data_ADHD
path <- "M:/Cohort Raw Data (do not edit)/DIAMOND/2.convert/3.child/2.ADHD"

filenames_list <- list.files(path=path, pattern="*.xlsx",full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read_xlsx(filename,sheet = "Data")
})
adhd <- do.call(rbind.data.frame, All)
colnames(adhd) <- make.names(colnames(adhd)) #52928
adhd1 <- adhd %>% select(Reference.Key.,Reference.Date.) %>% 
  mutate(Reference.Date.=as.Date(Reference.Date.))
length(unique(adhd1$Reference.Key.))#39403

#step7.3 combine MPH and ADHD data and select first episode
mph_adhd <- rbind(mph1,adhd1) %>% 
  arrange(Reference.Key.,Reference.Date.) %>% 
  `colnames<-` (c('Reference.Key.','adhd_date'))#51599
mph_adhd_3 <- merge(mph_adhd,mc_01to15_11_13 %>% select(Baby.s.Reference.Key.,dob),by.x = "Reference.Key.",by.y = "Baby.s.Reference.Key.") %>% 
  filter(adhd_date-dob>365.25*3) %>% 
  arrange(Reference.Key.,adhd_date) %>% 
  filter(duplicated(Reference.Key.)==F)



mc_01to15_11_14 <- merge(mc_01to15_11_13 %>% select(-adhd_date,-time),mph_adhd_3 %>% select(-dob),by.x = "Baby.s.Reference.Key.",by.y = "Reference.Key.",all.x = T) %>% 
  mutate(child.adhd=ifelse(is.na(adhd_date),0,1),
         obs_end=pmin(adhd_date, ymd('2020-12-31'), bb_dod, na.rm= T)) %>% 
  mutate(time=obs_end-dob+1) %>% 
  mutate(y=as.numeric(year(dob)))


# options(digits = 12)
mc_01to15_sas <- mc_01to15_11_14 %>% 
  mutate(m_age1=as.numeric(dob-Date.of.Birth..yyyy.mm.dd..)) %>% 
  filter(!is.na(Baby.s.Reference.Key.),
         No..of.Abortion.==0,
         No..of.Perinatal.Deaths.==0,
         sex!="U",
         !is.na(Maternity.Episode..Gestation.weeks.),
         !is.na(Date.of.Birth..yyyy.mm.dd..)) %>% 
  mutate(dm_sas=ifelse(dm2=="MDM",1,0),
         # m_age1=m_age2,
         child_adhd=child.adhd,
         reference_key=Reference.Key.,
         multi_pre=Maternity.Episode..OBS.Hx...Cx..Multiple.Pregnancy..Y.N..,
         parity=Maternity.Episode..OBS.Hx...Cx..Parity.,
         insti=Institution..IPAS..,
         b_id=as.numeric(Baby.s.Reference.Key.)) %>% 
  select(reference_key,dm_sas,child_adhd,time,m_age1,ses,y,antihypertensive,antipsychotics,sedatives,
         antidepressants,antiepileptics,antiparkinson,stimulants,opioids,triptan,teratogenic,pcos,anxiety,
         depression,epilepsy,migraine,hypertension,renal,crohn,sleep,scz,bp,personality,intellectual,
         development,illicit,adhd_asd,thyroid,Rheumatoid,headache,smoke,drink,BMI_fac,dm1,mdm_period,
         treat_t2dm,gdm_treat,pgdm_type,sex,dob,preterm,multi_pre,parity,insti,b_id) %>% 
  mutate(dm1=ifelse(dm1==" Non DM","non_dm",as.character(dm1)),
         mdm_period=ifelse(mdm_period==" Non DM","non_dm",as.character(mdm_period)),
         pgdm_type=ifelse(pgdm_type==" Non DM","non_dm",as.character(pgdm_type)))
# ad_ter <- read.csv("M:/Personal/Gao Le/DIAMOND_crosscheck/teratogenic_ac.csv")
# gl_ter <- mc_01to15_sas %>% 
#   filter(dm_sas==1,teratogenic==1)
# 
# ter_dif <- ad_ter %>% filter(!b_id%in%gl_ter$Baby.s.Reference.Key.)
# ter_dif <- gl_ter %>% filter(!Baby.s.Reference.Key.%in%ad_ter$b_id) %>% 
#   select(Baby.s.Reference.Key.)
# 
# ad_obe <- read.csv("M:/Personal/Gao Le/DIAMOND_crosscheck/cohort_5.csv")
# table(ad_obe$MDM,ad_obe$Obesity)
# table(ad_obe$MDM,ad_obe$Normal)
# table(ad_obe$MDM,ad_obe$Overweight)
# table(ad_obe$MDM,ad_obe$Underweight)
# ad_obe1 <- ad_obe %>%
#   filter(GDM_treat==1)
# gl_obe <- mc_01to15_11_13 %>%
#   filter(dm1=="GDM",gdm_treat==1)
# ad_gl <- merge(ad_obe1,gl_obe,by.x = 'b_id',by.y = 'Baby.s.Reference.Key.') %>% 
#   mutate(dif=as.numeric(time)-as.numeric(studydur))


# obe_dif <- ad_obe1 %>% filter(!b_id%in%gl_obe$Baby.s.Reference.Key.)
# ter_dif <- gl_obe %>% filter(!Baby.s.Reference.Key.%in%ad_obe1$b_id) %>%
#   select(Baby.s.Reference.Key.,Reference.Key.)
# 

mc_01to15_sas$mdm_period_1 <- NA
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="0_90BFP")] <- 1
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="1_0_90BFP")] <- 2
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="2_90AFP")] <- 3
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="3_180AFP")] <- 4
mc_01to15_sas$mdm_period_1[which(mc_01to15_sas$mdm_period=="4_up180AFP")] <- 5
mc_01to15_sas$mdm_period_1[which(is.na(mc_01to15_sas$mdm_period_1))] <- 0

mc_01to15_sas$dm1_cat <- NA
mc_01to15_sas$dm1_cat[which(mc_01to15_sas$dm1=="PGDM")] <- 1
mc_01to15_sas$dm1_cat[which(mc_01to15_sas$dm1=="GDM")] <- 2
mc_01to15_sas$dm1_cat[which(is.na(mc_01to15_sas$dm1_cat))] <- 0


mc_01to15_sas$mdm_period_sen <- NA
mc_01to15_sas$mdm_period_sen[which(mc_01to15_sas$mdm_period=="0_90BFP")] <- 0
mc_01to15_sas$mdm_period_sen[which(is.na(mc_01to15_sas$mdm_period_sen))] <- 1

mc_01to15_sas$dm1_catsen <- NA
mc_01to15_sas$dm1_catsen[which(mc_01to15_sas$dm1=="PGDM")] <- 0
mc_01to15_sas$dm1_catsen[which(is.na(mc_01to15_sas$dm1_catsen))] <- 1

mc_01to15_sas$treat_t2dm1 <- NA
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tob")] <- 1
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tod")] <- 2
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="tb")] <- 3
mc_01to15_sas$treat_t2dm1[which(mc_01to15_sas$treat_t2dm=="non_treat")] <- 0


mc_01to15_sas1 <- mc_01to15_sas %>% 
  mutate(nor=ifelse(BMI_fac==" normal",1,0),
         obe=ifelse(BMI_fac=="obesity",1,0),
         over=ifelse(BMI_fac=="overweight",1,0),
         under=ifelse(BMI_fac=="underweight",1,0),
         b_id_new=1:535016)

nor obe over under
write.csv(mc_01to15_sas1 %>% filter(as.numeric(as.character(y))<=2014) %>% mutate(time1=time-2190),
          "D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_2014_license.csv")




#step8.2 cox regression
data_cox <- mc_01to15_11_13


#calculate person year
py <- mc_01to15_11_13 %>% 
  group_by(dm1) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(dm1)==F) %>% 
  select(dm1,adhd_count,p_year,pop)

py1 <- mc_01to15_11_13 %>% 
  group_by(dm2) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(dm2)==F) %>% 
  select(dm2,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13 %>% 
  group_by(mdm_period) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(mdm_period)==F) %>% 
  select(mdm_period,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13 %>% 
  group_by(pgdm_type) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(pgdm_type)==F) %>% 
  select(pgdm_type,adhd_count,p_year,pop)






#sen for at 6 y follow up
mc_01to15_11_13_6y <- mc_01to15_sas1%>% 
  filter(as.numeric(as.character(y))<=2014)

data_cox <- mc_01to15_11_13_6y



#filter mothers with GMD
mc_01to15_12_GDM_ref <- mc_01to15_11_13_6y %>% filter(dm1=="GDM") %>% 
  select(reference_key)
#find out all delivery records of those mother with at least one GDM
mc_01to15_12_1 <- mc_01to15_11_13_6y %>% 
  filter(reference_key%in%mc_01to15_12_GDM_ref$reference_key)
#find out mothers with non-mdm from the last step
mc_01to15_12_1_nonmdm_ref <- mc_01to15_12_1 %>% 
  filter(dm_sas==0) %>% select(reference_key)
#keep those mother with both MDM and NON-MDM records
mc_01to15_12_sibling <- mc_01to15_12_1 %>% 
  filter(reference_key%in%mc_01to15_12_1_nonmdm_ref$reference_key) %>% 
  filter(dm1_cat!=1)
table(mc_01to15_12_sibling$dm1)

mc_01to15_sas_sibling <- mc_01to15_12_sibling %>% 
  # mutate(dm_sas=ifelse(dm2=="MDM",1,0),
  #        dm_sas=ifelse(dm2=="MDM",1,0),
  #        # m_age1=m_age2,
  #        child_adhd=child.adhd,
  #        reference_key=Reference.Key.,
  #        multi_pre=Maternity.Episode..OBS.Hx...Cx..Multiple.Pregnancy..Y.N..,
  #        parity=Maternity.Episode..OBS.Hx...Cx..Parity.,
  #        insti=Institution..IPAS..,
  #        b_id=as.numeric(Baby.s.Reference.Key.)
  #        # m_age1=ifelse(!is.na(m_age),m_age,31.76),
  # ) %>% 
  select(reference_key,dm_sas,child_adhd,time,m_age1,ses,y,antihypertensive,antipsychotics,sedatives,
         antidepressants,antiepileptics,antiparkinson,stimulants,opioids,triptan,teratogenic,pcos,anxiety,
         depression,epilepsy,migraine,hypertension,renal,crohn,sleep,scz,bp,personality,intellectual,
         development,illicit,adhd_asd,thyroid,Rheumatoid,headache,smoke,drink,BMI_fac,dm1,mdm_period,
         treat_t2dm,gdm_treat,pgdm_type,sex,dob,preterm,multi_pre,parity,insti,b_id) %>% 
  mutate(nor=ifelse(BMI_fac==" normal",1,0),
         obe=ifelse(BMI_fac=="obesity",1,0),
         over=ifelse(BMI_fac=="overweight",1,0),
         under=ifelse(BMI_fac=="underweight",1,0))

test1 <- mc_01to15_sas_sibling %>% 
  group_by(reference_key,dm_sas) %>% 
  mutate(no_b=length(reference_key),no_adhd=sum(child_adhd),rate=no_adhd/no_b) %>% 
  select(reference_key,dm_sas,child_adhd,no_b,no_adhd,rate) %>% 
  ungroup()
test1_dm <- test1 %>% filter(dm_sas==1) %>% 
  filter(duplicated(reference_key)==F)
test1_no <- test1 %>% filter(dm_sas==0) %>% 
  filter(duplicated(reference_key)==F)

test2 <- merge(test1_dm,test1_no,by="reference_key") %>% 
  filter(rate.x!=rate.y)

test3 <- mc_01to15_sas_sibling %>% 
  filter(reference_key%in%test2$reference_key)

write.csv(mc_01to15_sas_sibling,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_sibling_6y.csv")

