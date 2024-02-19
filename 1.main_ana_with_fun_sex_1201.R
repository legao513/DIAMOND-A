#sen for at 6 y follow up
mc_01to15_11_13_6y <- mc_01to15_11_13 %>% 
  filter(ymd('2020-12-31')-ymd(dob)>=365.25*6)

data_cox <- mc_01to15_11_13_6y
options(scipen=999)


#calculate person year
py <- mc_01to15_11_13_6y %>% 
  group_by(dm1) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(dm1)==F) %>% 
  select(dm1,adhd_count,p_year,pop)

py1 <- mc_01to15_11_13_6y %>% 
  group_by(dm2) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(dm2)==F) %>% 
  select(dm2,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13_6y %>% 
  group_by(mdm_period) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(mdm_period)==F) %>% 
  select(mdm_period,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13_6y %>% 
  group_by(pgdm_type) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(pgdm_type)==F) %>% 
  select(pgdm_type,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13_6y %>% 
  filter(dm1=="GDM") %>% 
  group_by(gdm_treat) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(pgdm_type)==F) %>% 
  select(gdm_treat,adhd_count,p_year,pop)

py2 <- mc_01to15_11_13_6y %>% 
  filter(dm1=="PGDM") %>% 
  group_by(treat_t2dm) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(pgdm_type)==F) %>% 
  select(treat_t2dm,adhd_count,p_year,pop)


library(finalfit)
library(dplyr)

#crude HR
options(digits = 10)

res.cox1 <- coxph(Surv(time, child.adhd) ~ dm2, data =data_cox)
summary(res.cox1)

res.cox2 <- coxph(Surv(time, child.adhd) ~ dm1, data = data_cox)
summary(res.cox2)

res.cox3 <- coxph(Surv(time, child.adhd) ~ dm1, data = data_cox, 
                  subset = dm1 %in% c("PGDM","GDM"))
summary(res.cox3)

res.cox4 <- coxph(Surv(time, child.adhd) ~ mdm_period, data = data_cox)
summary(res.cox4)

res.cox5 <- coxph(Surv(time, child.adhd) ~ mdm_period, data = data_cox,
                  subset = mdm_period %in% c("0_90BFP","1_0_90BFP","2_90AFP","3_180AFP","4_up180AFP"))
summary(res.cox5)

res.cox6 <- coxph(Surv(time, child.adhd) ~ pgdm_type, data = data_cox,
                  subset = pgdm_type %in% c(" Non DM","T1DM","T2DM"))
summary(res.cox6)

res.cox7 <- coxph(Surv(time, child.adhd) ~ pgdm_type, data = data_cox,
                  subset = pgdm_type %in% c("T1DM","T2DM"))
summary(res.cox7)

res.cox8 <- coxph(Surv(time, child.adhd) ~ gdm_treat, data = data_cox %>% filter(dm1=="GDM"))
summary(res.cox8)

res.cox9 <- coxph(Surv(time, child.adhd) ~ treat_t2dm, data = data_cox %>% filter(pgdm_type=="T2DM"))
summary(res.cox9)




#filter mothers with GMD
mc_01to15_12_GDM_ref <- mc_01to15_11_13_6y %>% filter(dm1=="GDM") %>% 
  select(Reference.Key.)
#find out all delivery records of those mother with at least one GDM
mc_01to15_12_1 <- mc_01to15_11_13_6y %>% 
  filter(Reference.Key.%in%mc_01to15_12_GDM_ref$Reference.Key.)
#find out mothers with non-mdm from the last step
mc_01to15_12_1_nonmdm_ref <- mc_01to15_12_1 %>% 
  filter(dm1==" Non DM") %>% select(Reference.Key.)
#keep those mother with both MDM and NON-MDM records
mc_01to15_12_sibling <- mc_01to15_12_1 %>% 
  filter(Reference.Key.%in%mc_01to15_12_1_nonmdm_ref$Reference.Key.) %>% 
  filter(dm1!="PGDM") %>% 
  arrange(Reference.Key.,pregancy_date)
table(mc_01to15_12_sibling$dm1)

mc_01to15_sas_sibling <- mc_01to15_12_sibling %>% 
  mutate(dm_sas=ifelse(dm2=="MDM",1,0),
         # m_age1=ifelse(!is.na(m_age),m_age,31.76),
         child_adhd=child.adhd,reference_key=Reference.Key.) %>% 
  select(reference_key,dm_sas,child_adhd,time,m_age1,ses,y,antihypertensive,antipsychotics,sedatives,
         antidepressants,antiepileptics,antiparkinson,stimulants,opioids,triptan,teratogenic,
         pcos,anxiety,depression,epilepsy,migraine,hypertension,crohn,renal,sleep,scz,bp,
         personality,intellectual,development,illicit,adhd_asd,thyroid,Rheumatoid,
         headache,smoke,drink,BMI_fac,mdm_period,preterm)

mc_01to15_sas_sibling$mdm_period_1 <- NA
mc_01to15_sas_sibling$mdm_period_1[which(mc_01to15_sas_sibling$mdm_period=="0_90BFP")] <- 1
mc_01to15_sas_sibling$mdm_period_1[which(mc_01to15_sas_sibling$mdm_period=="1_0_90BFP")] <- 2
mc_01to15_sas_sibling$mdm_period_1[which(mc_01to15_sas_sibling$mdm_period=="2_90AFP")] <- 3
mc_01to15_sas_sibling$mdm_period_1[which(mc_01to15_sas_sibling$mdm_period=="3_180AFP")] <- 4
mc_01to15_sas_sibling$mdm_period_1[which(mc_01to15_sas_sibling$mdm_period=="4_up180AFP")] <- 5
mc_01to15_sas_sibling$mdm_period_1[which(is.na(mc_01to15_sas_sibling$mdm_period_1))] <- 0


mc_01to15_sas_sibling$mdm_period_sen <- NA
mc_01to15_sas_sibling$mdm_period_sen[which(mc_01to15_sas_sibling$mdm_period=="2_90AFP")] <- 0
mc_01to15_sas_sibling$mdm_period_sen[which(is.na(mc_01to15_sas_sibling$mdm_period_sen))] <- 1



write.csv(mc_01to15_sas_sibling,"D:/data_no_share_20200506/6.DIAMOND/mc_01to15_sas_sibling_6y.csv")



#calculate person year
py <- mc_01to15_12_sibling %>% 
  group_by(dm1) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(dm1)==F) %>% 
  select(dm1,adhd_count,p_year,pop)

py1 <- mc_01to15_12_sibling %>% 
  group_by(mdm_period) %>% 
  mutate(p_year=as.numeric(sum(time)/365.25),adhd_count=sum(child.adhd),pop=length(Baby.s.Reference.Key.)) %>% 
  filter(duplicated(mdm_period)==F) %>% 
  select(mdm_period,adhd_count,p_year,pop)



#sibling crude
res.cox1 <- coxph(Surv(time, child.adhd) ~ dm1, data =mc_01to15_12_sibling)
summary(res.cox1)

res.cox1 <- coxph(Surv(time, child.adhd) ~ mdm_period, data =mc_01to15_12_sibling)
summary(res.cox1)
