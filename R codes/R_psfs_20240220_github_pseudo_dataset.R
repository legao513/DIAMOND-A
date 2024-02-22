#generate random dataset
# set.seed(5530)
# data1 <- data.frame(baby_id=1:20000,
#                     mother_id=round(runif(20000, min=20001, max=35000)),
#                    mdm=rbinom(20000, 1, 0.1),
#                    child_adhd=rbinom(20000, 1, 0.1),
#                    followup_time=round(runif(20000, min=1, max=7305)),
#                    mother_age=round(runif(20000, min=4825, max=21223)),
#                    ses=sample(c("LOW","HIGH","MEDIUMLOW","MEDIUMHIGH"), 20000, replace=TRUE),
#                    birth_year=sample(2001:2014, 20000, replace=TRUE),
#                    antihypertensive=rbinom(20000, 1, 0.01),
#                    antipsychotics=rbinom(20000, 1, 0.01),
#                    sedatives=rbinom(20000, 1, 0.01),
#                    antidepressants=rbinom(20000, 1, 0.01),
#                    antiepileptics=rbinom(20000, 1, 0.01),
#                    antiparkinson=rbinom(20000, 1, 0.01),
#                    stimulants=rbinom(20000, 1, 0.01),
#                    opioids=rbinom(20000, 1, 0.01),
#                    triptan=rbinom(20000, 1, 0.005),
#                    teratogenic=rbinom(20000, 1, 0.01),
#                    pcos=rbinom(20000, 1, 0.005),
#                    anxiety=rbinom(20000, 1, 0.01),
#                    depression=rbinom(20000, 1, 0.01),
#                    epilepsy=rbinom(20000, 1, 0.01),
#                    migraine=rbinom(20000, 1, 0.01),
#                    hypertension=rbinom(20000, 1, 0.01),
#                    renal_disease=rbinom(20000, 1, 0.01),
#                    crohn_disease=rbinom(20000, 1, 0.01),
#                    sleep_disorder=rbinom(20000, 1, 0.01),
#                    scz=rbinom(20000, 1, 0.01),
#                    bipolar=rbinom(20000, 1, 0.01),
#                    personality_disorder=rbinom(20000, 1, 0.005),
#                    intellectual_disability=rbinom(20000, 1, 0.01),
#                    psy_development_disorder=rbinom(20000, 1, 0.01),
#                    illicit_drug=rbinom(20000, 1, 0.01),
#                    adhd_asd=rbinom(20000, 1, 0.001),
#                    thyroid_disorder=rbinom(20000, 1, 0.01),
#                    rheumatoid_arthristis=rbinom(20000, 1, 0.01),
#                    headache=rbinom(20000, 1, 0.01),
#                    smoke=rbinom(20000, 1, 0.01),
#                    drink=rbinom(20000, 1, 0.01),
#                    baby_sex=sample(c("F","M"), 20000, replace=TRUE),
#                    preterm_birth=rbinom(20000, 1, 0.01),
#                    multi_pregnancy=sample(c("N","Y"), 20000, replace=TRUE),
#                    parity=sample(c(">= 5","0","1","2","3","4"), 20000, replace=TRUE),
#                    birth_institution=sample(c("AHN","KWH","PMH","PYN","QMH","UCH"), 20000, replace=TRUE),
#                    normal_bmi=rbinom(20000, 1, 0.999)) %>% 
#   mutate(bmi_fac=ifelse(normal_bmi==0,sample(c("obe","over","under"),length(normal_bmi),replace = T),0),
#          obesity=ifelse(bmi_fac=="obe",1,0),
#          overweight=ifelse(bmi_fac=="over",1,0),
#          underweight=ifelse(bmi_fac=="under",1,0),
#          gdm=ifelse(mdm==1,rbinom(sum(data1$mdm), 1, 0.93),0),
#          pgdm=ifelse(mdm==1&gdm==0,1,0),
#          t1dm=ifelse(pgdm==1,rbinom(sum(data1$pgdm), 1, 0.45),0),
#          t2dm=ifelse(pgdm==1& t1dm==0,1,0))
# 
# write.csv(data1,"D:/OneDrive - connect.hku.hk/Other-share/5.Phd project/10.mother-child linkage/Github DIAMOND-A/DIAMOND-A/R codes/data1.csv")



library(rlang)
library(WeightIt)
library(haven)
library(dplyr)
library(tableone)
library(data.table)
library(smd)
library(ggplot2)
library(survival)
library(survminer)
library(EValue)



#1.step 1 Calculate the propensity score for t-PA mdm for the study population using logistic regression
## This is the analysis for mdm vs non-mdm, for the analysis of gdm, pgdm, t1dm and t2dm, please change the variable 'mdm' to 'gdm', 'pgdm', 't1dm' and 't2dm' respectively
## For the sibling analysis, please change the variable 'mdm' to 'gdm', and add a strata function in the cox regression model
#read the cleaned database
data1 <- read.csv("D:/OneDrive - connect.hku.hk/Other-share/5.Phd project/10.mother-child linkage/Github DIAMOND-A/DIAMOND-A/R codes/data1.csv")

### calculate the person year, number of total population and number of evenets
data_count <- data1 %>% 
  group_by(mdm) %>% 
  mutate(py_total=round(sum(followup_time)/365),
         ppl_total=length(baby_id),
         event_total=sum(child_adhd)) %>% 
  filter(duplicated(mdm)==F) %>% 
  dplyr::select(mdm,child_adhd,py_total,ppl_total,event_total) %>% 
  mutate(rate=round(event_total/py_total*100,2))

data1_1 <- data1
#generate the propensity score using logistic regression
#pls note, except for maternal age, all other covariates were considered as caterogy vaiable

ps.fit <- glm(mdm ~ mother_age + as.factor(birth_year)+as.factor(ses)+baby_sex+multi_pregnancy+
                birth_institution+parity+antihypertensive+antipsychotics+sedatives+antidepressants+
                antiepileptics+antiparkinson+stimulants+opioids+triptan+teratogenic+pcos+
                anxiety+depression+epilepsy+migraine+hypertension+renal_disease+crohn_disease+
                sleep_disorder+scz+bipolar+personality_disorder+intellectual_disability+
                psy_development_disorder+illicit_drug+adhd_asd+thyroid_disorder+rheumatoid_arthristis+
                headache+smoke+drink+normal_bmi+obesity+overweight+underweight,data = data1_1,family = binomial)

data1_1$pscore <- predict(ps.fit, type='response')


#some data test --- just igonore
# data1_1$pscore <- NA
# data1_1$pscore <- as.numeric(substring(data1_1$pscore1,1,11))
# data1_new <- data1_1 %>%
#   arrange(baby_id) %>%
#   select(baby_id,pscore)
# 
# data2 <- read_sas("D:/data_no_share_20200506/6.DIAMOND/ps_fs.sas7bdat") %>%
#   arrange(baby_id) %>%
#   select(baby_id,psweight,strata,ps,mdm)
# 



#2.step 2- trim the non-overlapping regions of the PS
#find out the min and max ps score of the exposure group
options(scipen=999)

test_exp <- data1_1 %>% filter(mdm==1)
min_exp <- min(test_exp$pscore)
min_exp <- round(min_exp,6) #currently I found if keep 6 decimals, the result is closest to sas

max_exp <- max(test_exp$pscore)
max_exp <- round(max_exp,6)
#find out the min and max ps score of the non exposure group
test_noexp <- data1_1 %>% filter(mdm==0)
min_non <- min(test_noexp$pscore)
min_non <- round(min_non,6)

max_non <- max(test_noexp$pscore)
max_non <- round(max_non,6)

#select the max of the (exposure min score and non exposure min score) and the min of (exposure max score and non exposure max score)
#keep records between the min and max as the trim 
test_trim <- data1_1 %>% filter(pscore>=max(min_exp,min_non)&pscore<=min(max_exp,max_non))
# y <- test_trim %>% filter(mdm==1)
# x <- test_trim %>% filter(mdm==0)



#3. Stratified analysis: categorize the propensity score into quatiles and perform a stratified analysis. (For pss method=cohort )
test_stra <- test_trim %>% 
  arrange(pscore) %>% 
  mutate(rank_raw=1:length(baby_id)) %>% 
  mutate(strata=floor(rank_raw*50/(nrow(test_trim)+1))+1)





#step 4 Finally, assign weights to unexposed in each strata based on the distribution of the exposed
#as we use ATT method, all treated group has a weight of 1
strata_n <- test_stra

strata_n1 <- strata_n %>% filter(mdm==1) %>% 
  group_by(strata) %>% 
  mutate(ct=length(baby_id)) %>%
  arrange(strata,-ct) %>% 
  filter(duplicated(strata)==F) %>% 
  ungroup() %>% 
  select(strata,mdm,ct) %>% 
  `colnames<-` (c("strata","exp","ct_exp"))

strata_n2 <- strata_n %>% filter(mdm==0) %>% 
  group_by(strata) %>% 
  mutate(ct=length(baby_id)) %>%
  arrange(strata,-ct) %>% 
  filter(duplicated(strata)==F) %>% 
  ungroup() %>% 
  select(strata,mdm,ct) %>% 
  `colnames<-` (c("strata","nonexp","ct_nonexp"))


strata_nall <- merge(strata_n1,strata_n2,by="strata") %>% 
  filter(!is.na(ct_exp),!is.na(ct_nonexp),ct_exp!=0,ct_nonexp!=0) %>% 
  mutate(w_unexp=(ct_exp/sum(ct_exp))/(ct_nonexp/sum(ct_nonexp)),w_exp=1)
w_exp <- strata_nall %>% select(strata,w_exp,exp) %>% 
  `colnames<-` (c("strata","w","mdm"))

w_unexp <- strata_nall %>% select(strata,w_unexp,nonexp) %>% 
  `colnames<-` (c("strata","w","mdm"))
w_all <- rbind(w_exp,w_unexp)




#output
# unweighted database -- SMD (Appendix 5 and 6)
#pls change all category varibale to binary, i did not list all of them here
data_unweighted <- data1_1 %>% 
  mutate(y2001=ifelse(birth_year==2001,1,0),
         y2002=ifelse(birth_year==2002,1,0),
         y2003=ifelse(birth_year==2003,1,0),
         y2004=ifelse(birth_year==2004,1,0),
         y2005=ifelse(birth_year==2005,1,0),
         y2006=ifelse(birth_year==2006,1,0),
         y2007=ifelse(birth_year==2007,1,0),
         y2008=ifelse(birth_year==2008,1,0),
         y2009=ifelse(birth_year==2009,1,0),
         y2010=ifelse(birth_year==2010,1,0),
         y2011=ifelse(birth_year==2011,1,0),
         y2012=ifelse(birth_year==2012,1,0),
         y2013=ifelse(birth_year==2013,1,0),
         y2014=ifelse(birth_year==2014,1,0),
         seslow=ifelse(ses==" LOW",1,0),
         seshigh=ifelse(ses=="HIGH",1,0),
         sesm1=ifelse(ses=="Midium1",1,0),
         sesm2=ifelse(ses=="Midium2",1,0))


varsp <- c("mother_age","y2001","y2002","y2003","y2004","y2005","y2006","y2007","y2008","y2009",
           "y2010","y2011","y2012","y2013","y2014",
           "seslow","seshigh","sesm1","sesm2",
           "baby_sex","multi_pregnancy","birth_institution","parity","antihypertensive","antipsychotics",
           "sedatives","antidepressants","antiepileptics","antiparkinson","stimulants",
           "opioids","triptan","teratogenic","pcos","anxiety","depression","epilepsy","migraine",
           "hypertension","renal_disease","crohn_disease","sleep_disorder","scz","bipolar","personality_disorder","intellectual_disability",
           "psy_development_disorder","illicit_drug","adhd_asd","thyroid_disorder","rheumatoid_arthristis","headache","smoke","drink",
           "normal_bmi","obesity","overweight","underweight")
catvarsp<- c("y2001","y2002","y2003","y2004","y2005","y2006","y2007","y2008","y2009",
             "y2010","y2011","y2012","y2013","y2014",
             "seslow","seshigh","sesm1","sesm2",
             "baby_sex","multi_pregnancy","birth_institution","parity","antihypertensive","antipsychotics",
             "sedatives","antidepressants","antiepileptics","antiparkinson","stimulants",
             "opioids","triptan","teratogenic","pcos","anxiety","depression","epilepsy","migraine",
             "hypertension","renal_disease","crohn_disease","sleep_disorder","scz","bipolar","personality_disorder","intellectual_disability",
             "psy_development_disorder","illicit_drug","adhd_asd","thyroid_disorder","rheumatoid_arthristis","headache","smoke","drink",
             "normal_bmi","obesity","overweight","underweight")

basTable<-CreateTableOne( data = data_unweighted,vars = varsp, factorVars = catvarsp, 
                          strata = "mdm",includeNA=TRUE)

baseline<-as.data.frame(as.data.table(print(basTable,test = F,noSpaces = T,smd = TRUE),
                                      keep.rownames = T)) #SMD*100 and then ill in appendix 6



#code for weighted table after trimming (Appendix 6, SMD*100 and then fill in appendix 6)
#count * weight for each covariate, and round the sum data
test_trim1 <- merge(test_stra,w_all,by=c("strata","mdm")) %>% 
  mutate(y2001=ifelse(birth_year==2001,1,0),
         y2002=ifelse(birth_year==2002,1,0),
         y2003=ifelse(birth_year==2003,1,0),
         y2004=ifelse(birth_year==2004,1,0),
         y2005=ifelse(birth_year==2005,1,0),
         y2006=ifelse(birth_year==2006,1,0),
         y2007=ifelse(birth_year==2007,1,0),
         y2008=ifelse(birth_year==2008,1,0),
         y2009=ifelse(birth_year==2009,1,0),
         y2010=ifelse(birth_year==2010,1,0),
         y2011=ifelse(birth_year==2011,1,0),
         y2012=ifelse(birth_year==2012,1,0),
         y2013=ifelse(birth_year==2013,1,0),
         y2014=ifelse(birth_year==2014,1,0),
         seslow=ifelse(ses==" LOW",1,0),
         seshigh=ifelse(ses=="HIGH",1,0),
         sesm1=ifelse(ses=="Midium1",1,0),
         sesm2=ifelse(ses=="Midium2",1,0))

unloadNamespace("tableone")
library(survey)

trim_weight <- svydesign(ids = ~ baby_id, strata = ~ mdm, weights = ~ w,
                         data = test_trim1,nest=TRUE)
library(tableone)
tab1 <- svyCreateTableOne(vars = varsp,
                          strata = "mdm", data = trim_weight,
                          factorVars = catvarsp)
baseline<-as.data.frame(as.data.table(print(tab1,test = F,noSpaces = T,smd = TRUE),
                                      keep.rownames = T))









#cox regression
#calculate the total number of events, no of patients and no of person year in each exposure group (unweighted dataset)
#No of events (tables)
table(data1_1$mdm,data1_1$child_adhd)
#No of individuals (figures)
table(data1_1$mdm)
#person years (tables)
py <- data1_1 %>% 
  group_by(mdm) %>% 
  mutate(pyt=sum(followup_time)/365) %>% 
  filter(duplicated(mdm)==F) %>% 
  select(mdm,pyt)

#cox regression, pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the cox model
cox_data <- merge(test_stra,w_all,by=c("strata","mdm")) %>% 
  mutate(followup_time=ifelse(followup_time==0,followup_time+1,followup_time))
res.cox1 <- coxph(Surv(followup_time, child_adhd) ~ mdm, 
                  weights=w,data =  cox_data)
summary(res.cox1) #(keep three decimals pls, it will be easier for us to meta analyse)

# for the sibling analysis, change the variable 'mdm' to 'gdm', and add a strata of mother_id in the cox regression
res.cox1 <- coxph(Surv(followup_time, child_adhd) ~ gdm+strata(mother_id), 
                  weights=w, data =  cox_data)
summary(res.cox1) #(keep three decimals pls, it will be easier for us to meta analyse)




km_sz<-survfit(Surv(followup_time, child_adhd) ~ mdm, 
               weights = w,data=cox_data)

a <- ggsurvplot(km_sz,
           fun = "event")
#could you please save the data after trimming, with columns of followup_time, child_adhd, mdm, and weight for us to generate the cumulative plot
#for MDM vs Non-MDMD, GDM vs Non-MDMD and PGDM vs Non-MDM only of the main analysis

a$data.survplot

output0<- a$data.survplot %>% 
  filter(mdm==0) %>% 
  mutate(incidence_0=n.event/n.risk) %>% 
  select(incidence_0,time)

output1<- a$data.survplot %>% 
  filter(mdm==1) %>% 
  mutate(incidence_1=n.event/n.risk) %>% 
  select(incidence_1,time)

output <- merge(output0,output1,by="time") %>% 
  mutate(rd=incidence_1-incidence_0)





# poisson regression (pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the fomular)
#point estimate
exp(summary(glm(child_adhd ~ mdm,weights = w,offset = log(followup_time),family="poisson",data=cox_data))$coefficients[2,1])
#confindence interval
exp(confint(glm(child_adhd ~ mdm,weights = w,offset = log(followup_time),family="poisson",data=cox_data)))


# negative binomial regression (pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the fomular)
exp(MASS::glm.nb(child_adhd ~ mdm+offset(log(followup_time)) ,weights = w,data=cox_data %>% 
                   mutate(followup_time=as.numeric(followup_time)))$coefficients[2])
#confindence interval
exp(confint(MASS::glm.nb(child_adhd ~ mdm+offset(log(followup_time)) ,weights = w,data=cox_data %>% 
                           mutate(followup_time=as.numeric(followup_time)))))





options(digits = 8)
# calculate the risk difference using coxph
res.cox1 <- coxph(Surv(followup_time, child_adhd) ~ mdm, 
                  weights=w,data =  cox_data)
#exposure==1
summary(survfit(res.cox1,newdata = data.frame(mdm=1)), times = seq(6,18,2)*360)
#exposure==0
summary(survfit(res.cox1,newdata = data.frame(mdm=0)), times = seq(6,18,2)*360)




#calculate e value
evalues.RR(est = 1.16, lo = 1.08, hi = 1.24)






