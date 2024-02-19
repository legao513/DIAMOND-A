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


#1.step 1 Calculate the propensity score for t-PA dm_sas for the study population using logistic regression
#read the cleaned databse
data1 <- read.csv("D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_2014.csv")

### calculate the person year, number of total population and number of evenets
data_count <- data1 %>% 
  group_by(dm_sas) %>% 
  mutate(py_total=round(sum(time)/365),
         ppl_total=length(b_id_new),
         event_total=sum(child_adhd)) %>% 
  filter(duplicated(dm_sas)==F) %>% 
  dplyr::select(dm_sas,child_adhd,py_total,ppl_total,event_total) %>% 
  mutate(rate=round(event_total/py_total*100,2))

data1_1 <- data1
#generate the propensity score using logistic regression
#pls note, except for maternal age, all other covariates were considered as caterogy vaiable

ps.fit <- glm(dm_sas ~ m_age1 + as.factor(y)+as.factor(ses)+sex+multi_pre+
                insti+parity+antihypertensive+antipsychotics+sedatives+antidepressants+
                antiepileptics+antiparkinson+stimulants+opioids+triptan+teratogenic+pcos+
                anxiety+depression+epilepsy+migraine+hypertension+renal+crohn+sleep+scz+bp+personality+
                intellectual+development+illicit+adhd_asd+thyroid+Rheumatoid+headache+smoke+
                drink+nor+obe+over+under,data = data1_1,family = binomial)

data1_1$pscore <- predict(ps.fit, type='response')


#some data test --- just igonore
# data1_1$pscore <- NA
# data1_1$pscore <- as.numeric(substring(data1_1$pscore1,1,11))
# data1_new <- data1_1 %>%
#   arrange(b_id_new) %>%
#   select(b_id_new,pscore)
# 
# data2 <- read_sas("D:/data_no_share_20200506/6.DIAMOND/ps_fs.sas7bdat") %>%
#   arrange(b_id_new) %>%
#   select(b_id_new,psweight,strata,ps,dm_sas)
# 



#2.step 2- trim the non-overlapping regions of the PS
#find out the min and max ps score of the exposure group
options(scipen=999)

test_exp <- data1_1 %>% filter(dm_sas==1)
min_exp <- min(test_exp$pscore)
min_exp <- round(min_exp,6) #currently I found if keep 6 decimals, the result is closest to sas

max_exp <- max(test_exp$pscore)
max_exp <- round(max_exp,6)
#find out the min and max ps score of the non exposure group
test_noexp <- data1_1 %>% filter(dm_sas==0)
min_non <- min(test_noexp$pscore)
min_non <- round(min_non,6)

max_non <- max(test_noexp$pscore)
max_non <- round(max_non,6)

#select the max of the (exposure min score and non exposure min score) and the min of (exposure max score and non exposure max score)
#keep records between the min and max as the trim 
test_trim <- data1_1 %>% filter(pscore>=max(min_exp,min_non)&pscore<=min(max_exp,max_non))
# y <- test_trim %>% filter(dm_sas==1)
# x <- test_trim %>% filter(dm_sas==0)



#3. Stratified analysis: categorize the propensity score into quatiles and perform a stratified analysis. (For pss method=cohort )
test_stra <- test_trim %>% 
  arrange(pscore) %>% 
  mutate(rank_raw=1:length(b_id_new)) %>% 
  mutate(strata=floor(rank_raw*50/(nrow(test_trim)+1))+1)
  




#step 4 Finally, assign weights to unexposed in each strata based on the distribution of the exposed
#as we use ATT method, all treated group has a weight of 1
strata_n <- test_stra

strata_n1 <- strata_n %>% filter(dm_sas==1) %>% 
  group_by(strata) %>% 
  mutate(ct=length(b_id_new)) %>%
  arrange(strata,-ct) %>% 
  filter(duplicated(strata)==F) %>% 
  ungroup() %>% 
  select(strata,dm_sas,ct) %>% 
  `colnames<-` (c("strata","exp","ct_exp"))

strata_n2 <- strata_n %>% filter(dm_sas==0) %>% 
  group_by(strata) %>% 
  mutate(ct=length(b_id_new)) %>%
  arrange(strata,-ct) %>% 
  filter(duplicated(strata)==F) %>% 
  ungroup() %>% 
  select(strata,dm_sas,ct) %>% 
  `colnames<-` (c("strata","nonexp","ct_nonexp"))


strata_nall <- merge(strata_n1,strata_n2,by="strata") %>% 
  filter(!is.na(ct_exp),!is.na(ct_nonexp),ct_exp!=0,ct_nonexp!=0) %>% 
  mutate(w_unexp=(ct_exp/sum(ct_exp))/(ct_nonexp/sum(ct_nonexp)),w_exp=1)
w_exp <- strata_nall %>% select(strata,w_exp,exp) %>% 
  `colnames<-` (c("strata","w","dm_sas"))

w_unexp <- strata_nall %>% select(strata,w_unexp,nonexp) %>% 
  `colnames<-` (c("strata","w","dm_sas"))
w_all <- rbind(w_exp,w_unexp)




#output
# unweighted database -- SMD (Appendix 5 and 6)
#pls change all category varibale to binary, i did not list all of them here
data_unweighted <- data1_1 %>% 
  mutate(y2001=ifelse(y==2001,1,0),
         y2002=ifelse(y==2002,1,0),
         y2003=ifelse(y==2003,1,0),
         y2004=ifelse(y==2004,1,0),
         y2005=ifelse(y==2005,1,0),
         y2006=ifelse(y==2006,1,0),
         y2007=ifelse(y==2007,1,0),
         y2008=ifelse(y==2008,1,0),
         y2009=ifelse(y==2009,1,0),
         y2010=ifelse(y==2010,1,0),
         y2011=ifelse(y==2011,1,0),
         y2012=ifelse(y==2012,1,0),
         y2013=ifelse(y==2013,1,0),
         y2014=ifelse(y==2014,1,0),
         y2015=ifelse(y==2015,1,0),
         y2016=ifelse(y==2016,1,0),
         y2017=ifelse(y==2017,1,0),
         y2018=ifelse(y==2018,1,0),
         seslow=ifelse(ses==" LOW",1,0),
         seshigh=ifelse(ses=="HIGH",1,0),
         sesm1=ifelse(ses=="Midium1",1,0),
         sesm2=ifelse(ses=="Midium2",1,0))


varsp <- c("m_age1","y2001","y2002","y2003","y2004","y2005","y2006","y2007","y2008","y2009",
           "y2010","y2011","y2012","y2013","y2014","y2015","y2016","y2017","y2018",
           "seslow","seshigh","sesm1","sesm2",
           "sex","multi_pre","insti","parity","antihypertensive","antipsychotics",
           "sedatives","antidepressants","antiepileptics","antiparkinson","stimulants",
           "opioids","triptan","teratogenic","pcos","anxiety","depression","epilepsy","migraine",
           "hypertension","renal","crohn","sleep","scz","bp","personality","intellectual",
           "development","illicit","adhd_asd","thyroid","Rheumatoid","headache","smoke","drink",
           "nor","obe","over","under")
catvarsp<- c("y2001","y2002","y2003","y2004","y2005","y2006","y2007","y2008","y2009",
             "y2010","y2011","y2012","y2013","y2014","y2015","y2016","y2017","y2018",
             "seslow","seshigh","sesm1","sesm2",
             "sex","multi_pre","insti","parity","antihypertensive","antipsychotics",
             "sedatives","antidepressants","antiepileptics","antiparkinson","stimulants",
             "opioids","triptan","teratogenic","pcos","anxiety","depression","epilepsy","migraine",
             "hypertension","renal","crohn","sleep","scz","bp","personality","intellectual",
             "development","illicit","adhd_asd","thyroid","Rheumatoid","headache","smoke","drink",
             "nor","obe","over","under")

basTable<-CreateTableOne( data = data_unweighted,vars = varsp, factorVars = catvarsp, 
                          strata = "dm_sas",includeNA=TRUE)

baseline<-as.data.frame(as.data.table(print(basTable,test = F,noSpaces = T,smd = TRUE),
                                      keep.rownames = T)) #SMD*100 and then ill in appendix 6



#code for weighted table after trimming (Appendix 6, SMD*100 and then fill in appendix 6)
#count * weight for each covariate, and round the sum data
test_trim1 <- merge(test_stra,w_all,by=c("strata","dm_sas")) %>% 
  mutate(y2001=ifelse(y==2001,1,0),
         y2002=ifelse(y==2002,1,0),
         y2003=ifelse(y==2003,1,0),
         y2004=ifelse(y==2004,1,0),
         y2005=ifelse(y==2005,1,0),
         y2006=ifelse(y==2006,1,0),
         y2007=ifelse(y==2007,1,0),
         y2008=ifelse(y==2008,1,0),
         y2009=ifelse(y==2009,1,0),
         y2010=ifelse(y==2010,1,0),
         y2011=ifelse(y==2011,1,0),
         y2012=ifelse(y==2012,1,0),
         y2013=ifelse(y==2013,1,0),
         y2014=ifelse(y==2014,1,0),
         y2015=ifelse(y==2015,1,0),
         y2016=ifelse(y==2016,1,0),
         y2017=ifelse(y==2017,1,0),
         y2018=ifelse(y==2018,1,0),
         seslow=ifelse(ses==" LOW",1,0),
         seshigh=ifelse(ses=="HIGH",1,0),
         sesm1=ifelse(ses=="Midium1",1,0),
         sesm2=ifelse(ses=="Midium2",1,0))

unloadNamespace("tableone")
library(survey)

trim_weight <- svydesign(ids = ~ b_id_new, strata = ~ dm_sas, weights = ~ w,
                       data = test_trim1)
library(tableone)
tab1 <- svyCreateTableOne(vars = varsp,
                          strata = "dm_sas", data = trim_weight,
                          factorVars = catvarsp)
baseline<-as.data.frame(as.data.table(print(tab1,test = F,noSpaces = T,smd = TRUE),
                                      keep.rownames = T))









#cox regression
#calculate the total number of events, no of patients and no of person year in each exposure group (unweighted dataset)
#No of events (tables)
table(data1_1$dm_sas,data1_1$child_adhd)
#No of individuals (figures)
table(data1_1$dm_sas)
#person years (tables)
py <- data1_1 %>% 
  group_by(dm_sas) %>% 
  mutate(pyt=sum(time)/365) %>% 
  filter(duplicated(dm_sas)==F) %>% 
  select(dm_sas,pyt)

#cox regression, pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the cox model
cox_data <- merge(test_stra,w_all,by=c("strata","dm_sas")) %>% 
  mutate(time=ifelse(time==0,time+1,time))
res.cox1 <- coxph(Surv(time, child_adhd) ~ dm_sas, 
                  weights=w,data =  cox_data)
summary(res.cox1) #(keep three decimals pls, it will be easier for us to meta analyse)



km_sz<-survfit(Surv(time, child_adhd) ~ dm_sas, 
               weights = w,data=cox_data)

ggsurvplot(km_sz,
           fun = "event")
#could you please save the data after trimming, with columns of time, child_adhd, dm_sas, and weight for us to generate the cumulative plot
#for MDM vs Non-MDMD, GDM vs Non-MDMD and PGDM vs Non-MDM only of the main analysis

a$data.survplot

output0<- a$data.survplot %>% 
  filter(dm_sas==0) %>% 
  mutate(incidence_0=n.event/n.risk) %>% 
  select(incidence_0,time)

output1<- a$data.survplot %>% 
  filter(dm_sas==1) %>% 
  mutate(incidence_1=n.event/n.risk) %>% 
  select(incidence_1,time)

output <- merge(output0,output1,by="time") %>% 
  mutate(rd=incidence_1-incidence_0)





# poisson regression (pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the fomular)
#point estimate
exp(summary(glm(child_adhd ~ dm_sas,weights = w,offset = log(time),family="poisson",data=cox_data))$coefficients[2,1])
#confindence interval
exp(confint(glm(child_adhd ~ dm_sas,weights = w,offset = log(time),family="poisson",data=cox_data)))


# negative binomial regression (pls check if the smd (%) of the weighted table is more than 10, if so, pls add the covariates in the fomular)
exp(MASS::glm.nb(child_adhd ~ dm_sas+offset(log(time)) ,weights = w,data=cox_data %>% 
                   mutate(time=as.numeric(time)))$coefficients[2])
#confindence interval
exp(confint(MASS::glm.nb(child_adhd ~ dm_sas+offset(log(time)) ,weights = w,data=cox_data %>% 
                           mutate(time=as.numeric(time)))))





options(digits = 8)
# calculate the risk difference using coxph
res.cox1 <- coxph(Surv(time, child_adhd) ~ dm_sas, 
                  weights=w,data =  cox_data)
#exposure==1
summary(survfit(res.cox1,newdata = data.frame(dm_sas=1)), times = seq(6,18,2)*360)
#exposure==0
summary(survfit(res.cox1,newdata = data.frame(dm_sas=0)), times = seq(6,18,2)*360)
