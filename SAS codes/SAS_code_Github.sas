*Propensity score (PS) fine-stratification weighting was used to address the differences in baseline covariates.; 
*In SAS, this was achieved by applying the SAS macro developed by Desai, Rishi, 2020, 
	"Propensity score fine stratification SAS macro", https://doi.org/10.7910/DVN/U8JLCW, Harvard Dataverse, V5;

*Step 1: Assign directory containing analytic file (containing information on patient identifiers,
exposure, confounding varaibles, and binary outcomes/censoring time) and PS fine-stratification macro respectively;
libname a "D:\SAS\MDM_ADHD\analysis";
libname pss "D:\SAS\MDM_ADHD\ps";

*Step 2: Run PSFS;
%let path= D:\SAS\MDM_ADHD\ps; 		/*Please put the directory here*/ 
%let data= a.data1;		/*The dataset name*/
%let exp= MDM; 			/*The exposure variable*/
%let class= birth_institution parity multi_pregnancy baby_sex
birth_year SES antihypertensive antipsychotics sedatives antidepressants antiepileptics
antiparkinson stimulants opioids teratogenic triptan
PCOS Smoke drink Anxiety Depression
Epilepsy Migraine sleep_disorder hypertension renal_disease crohn_disease rheumatoid_arthristis
scz Bipolar intellectual_disability psy_development_disorder personality_disorder
thyroid_disorder headache ADHD_ASD
Overweight Obesity normal_bmi
Underweight illicit_drug; 		/*The list of class covariates*/
%let con= mother_age ; 		/*The list of continuous variable covariates*/
%let patientid=baby_id; 	/*Variable for patient id*/
%let outcome=child_adhd; 		/*Outcome variable name*/
%let followup=followup_time;		/*Follow-up time variable name*/
%let outexcel="D:\SAS\MDM_ADHD\results";		/*Output excel folder*/

%include "&path\PSS weighted analysis2.sas"; ** the stratification and analysis macro- %fine_stratification; 
%include "&path\Weighted Table 1s.sas"; ** macro for creating patient characteristic tables- %table1
										NOTE- this macro is invoked by calling the %fine_stratification macro ; 
options mprint mlogic;

%fine_stratification (in_data= &data, 
exposure= &exp, 
PS_provided= no, 
ps_var= ps, 
ps_class_var_list=  &class, 
ps_cont_var_list=  &con , 
PSS_method=cohort, 
time_unit= days,
estimand= att,
n_of_strata= 50, 
out_data= a.cohort_ps, 
id_var= b_id, effect_estimate= HR,
outcome= &outcome, survival_time= &followup, 
out_excel= D:\SAS\MDM_ADHD\results\MDMresults);

*Step 3: Run Cox PHReg with assigned weights;
*Note: Factor(s) with standardised differences greater than 10% was further adjusted in the Cox models.;
proc phreg data=a.cohort_ps COVS(AGGREGATE);
weight psweight;
model &followup*&outcome(0)=&exp /*add variables if SMD>0.1*/ /rl ties=efron;
id &patientid; /*if sibling analysis then replace with "strata m_id;" */
title "Cox model with robust standard error";
ods output ParameterEstimates=results1;
run;
