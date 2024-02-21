
%macro fine_stratification (in_data= , exposure= , PS_provided= no , ps_var= , ps_class_var_list= , ps_cont_var_list= ,
							interactions= , PSS_method=exposure,
						    n_of_strata= , out_data= PS_FS, id_var= , estimand= , effect_estimate= ,
							outcome= , survival_time= , time_unit= ,
							out_excel= , work_lib=clean);

*************************************************************************************************************

** step 0- Create the crude Table 1 **;

*************************************************************************************************************;

%table1 (in_for_table1= &in_data, treatment_var= &exposure, 
categorical_var_list= &ps_class_var_list, continuous_var_list= &ps_cont_var_list, out_table1= Crude_T1);

*************************************************************************************************************

** step 1- build and calculate the PS- if not already done **;

*************************************************************************************************************;
title 'PS model'; 

%let ps_provided1 = %upcase(&PS_provided); %if "&PS_provided1" = "NO" %then %do;

ods output Association=save_c; 

proc logistic data=&in_data desc ;
class &ps_class_var_list;
model &exposure= &ps_cont_var_list &ps_class_var_list &interactions; 
output out=sample_ps p=&ps_var;
run;

ods output clear; 
data save_c; set save_c; where label2='c'; keep label2 nvalue2; run;

%end;

%let ps_provided1 = %upcase(&PS_provided); %if "&PS_provided1" = "YES" %then %do;

data sample_ps; set &in_Data; run;

%end;

title; 

*************************************************************************************************************

** step 2- trim the non-overlapping regions of the PS **;

*************************************************************************************************************;

	*figure out the region of overlap;
	proc sql noprint;
	create table bounds_pscore as
	select max(&ps_var) as max_p, min(&ps_var) as min_p, &exposure
	from sample_ps
	group by &exposure;
	quit;

	*Save upper and lower bound of overlapping PS regions;

	proc sql noprint;
	select max(min_p) into: overlap_lower_bound from bounds_pscore; quit;
	proc sql noprint;
	select min(max_p) into: overlap_upper_bound from bounds_pscore; quit;

	*delete non-overlapping observations;

	data sample_ps_trimmed; set sample_ps;
	where &overlap_lower_bound<=&ps_var<=&overlap_upper_bound;
	run;
	
   data trim_summary; merge sample_ps (in=r) sample_ps_trimmed (in=s); 
   by &id_var; 
   if r=1 and s^=1 then Trimmed= 'Yes'; else trimmed='No'; 
   run; 

*************************************************************************************************************

** step 3- Conduct propensity score stratification 

*************************************************************************************************************;

%let PSS_method1 = %upcase(&PSS_method); %if "&PSS_method1" = "EXPOSURE" %then %do;

		** separate exposed only **;

		data only_exp; set sample_ps_trimmed; where &exposure=1; run;

		** Create strata based on the PS distribution of the exposed **;

		proc rank data=only_exp groups=&n_of_strata out=fs_exp;
		ranks strata_sas;
		var &ps_var;
		run;

		data fs_exp; set fs_exp; strata=1+strata_sas; * number of strata indicator, since proc rank assigns strata beginning from 0, we add 1 here **;
		drop strata_sas;
		run;

		*Get bounds of the PS for strata created;

				proc means data=fs_exp noprint;  
				 class strata; 
				 var &ps_var; 
				 output out=bounds_for_fs n=n 
				 min(&ps_var)=PS_min max(&ps_var)=PS_max; 
				run; 

				data bounds_for_fs; set bounds_for_fs;
				where strata ne .; *delete overall mean; 
				run;

		** Now, using the bounds of PS from the exposed distribution, assign strata in the unexposed patients;

		%let n=&n_of_strata; 

				%macro strata_in_unexp;

				** lets save the bounds for each separate strata into SAS datasets, which will be merged to the unexposed only 
				   cohort for the ease of strata creation in the next step **;

				%do strata= 1 %to &n;

					proc sql ; create table PSvalue&strata as
					select (ps_min) as ps_low&strata
					from bounds_for_fs
				    where strata=&strata; 
					quit;

				%end;

				data only_unexp; set sample_ps_trimmed; where &exposure=0; run;

				data only_unexp1; 
		              merge only_unexp psvalue1-psvalue&n;
				run;

				** the above merge creates a dataset with extra columns with each column containing 
				value of the lower bound of the stratum boundary calculated from the ps distribution of the exposed (eg 5 strata= 5 extra columns). 

				This information is used below to determine stratum membership of each unexposed subject. however, before
				moving on, we fill all the rows with the values of this variable by using missing=mean because the merge only fills one row since
				the datasets psvalue1 through n only has one row; 

				proc stdize out=only_unexp2 reponly missing=mean; var ps_low:; run; 

				** create strata **;

				data unexp_strata; set only_unexp2;
				array abc(*) ps_low: ;
				do i=3 to dim(abc) while (strata=.);
				if 0<=&ps_var<ps_low2 then strata=1; ** the first strata between 0 and the minimum PS for the second stratum (stratum 1) **;
				if abc(i-1)<=&ps_var<abc(i) then strata=(i-1); ** rest of the strata based on the PS between minimum PS for two neighbouring 
				                                             strata. Note that since the last strata will be defined based on PS values
				                                             greater than the min PS in the penultimate strata, this step will code that
				                                             value as . In the next step, we give this strata the value equal to the 
															 highest strata**;
				end;
				run;

				data fs_unexp; set unexp_strata; if strata=. then strata=&n; 
				drop i ps_low:;
				run;

				/*
				proc means; var &ps_var; class strata; run;
				*/
				%mend;

				 %strata_in_unexp;

		** now concatenate the exposed and unexposed patients **;

		proc sort data=fs_exp; by &id_var; run;
		proc sort data=fs_unexp; by &id_var; run;

		data stratified;
		set fs_exp fs_unexp;
		by &id_var;
		run;

%end;

%let PSS_method1 = %upcase(&PSS_method); %if "&PSS_method1" = "COHORT" %then %do;

		proc rank data=sample_ps_trimmed groups=&n_of_strata out=fs;
		ranks strata_sas;
		var &ps_var;
		run;

		data stratified; set fs; 
		strata=1+strata_sas; 
		* number of strata indicator, since proc rank assigns strata beginning from 0, we add 1 here **;
		drop strata_sas;
		run;
			*Get bounds of the PS for strata created;

				proc means data=stratified noprint;  
				 class strata; 
				 var &ps_var; 
				 output out=bounds_for_fs n=n 
				 min(&ps_var)=PS_min max(&ps_var)=PS_max; 
				run; 

				data bounds_for_fs; set bounds_for_fs;
				where strata ne .;
				run;

%end;

*************************************************************************************************************

** step 4 Finally, assign weights in each strata to account for stratum membership

*************************************************************************************************************;

		proc sort data=stratified out=matched_ps1; by strata; run;

		*total exposed; 

		data total_exposed; set matched_ps1; where &exposure=1; run;

		data total_exposed1; set total_exposed; by strata; 
		if first.strata then total_Exp=0; total_Exp+1;
		if first.strata then events_exp=0; events_exp+&outcome;
		if last.strata then output;
		keep strata total_exp events_exp;
		
		run;

		*total unexposed; 

		data total_unexposed; set matched_ps1; where &exposure=0; run;

		data total_unexposed1; set total_unexposed; by strata; 
		if first.strata then total_unExp=0; total_unExp+1;
		if first.strata then events_unexp=0; events_unexp+&outcome;
		if last.strata then output;
		keep strata total_unexp events_unexp;
		run;

		*merge the two files;

		data  strata_n; 
		merge  total_exposed1  total_unexposed1 ;
		by strata;
		if total_exp=. then total_exp=0; if total_unexp=. then total_unexp=0;
		run;

		*Drop the strata with either 0 exposed or 0 unexposed;

		data  strata_n; set  strata_n; if total_exp=0 or total_unexp=0 then delete; run;

		proc sql noprint; select sum(total_exp) into: sum_exp from  strata_n; quit; 
		
		proc sql noprint; select sum(total_unexp) into: sum_unexp from  strata_n; quit; 

		proc sql noprint; select sum(total_unexp)+sum(total_exp) into: sum_total from  strata_n; quit; 

		**Create weights **;

		data  strata_n; set  strata_n; 
		att_unexp_weight=(total_exp/&sum_exp)/(total_unexp/&sum_unexp); *** ATT Weighting ;
		ate_exp_weight= ((total_exp+total_unexp)/&sum_total)/(total_exp/&sum_exp); **ate weights for the exposed group;
		ate_unexp_weight= ((total_exp+total_unexp)/&sum_total)/(total_unexp/&sum_unexp);**ate weights for the unexposed group;
		run;

		data strata_n; set strata_n; 
		att_weighted_unexp= total_unexp*att_unexp_weight;
		att_weighted_unexp_events= events_unexp*att_unexp_weight; 

		ate_weighted_unexp= total_unexp*ate_unexp_weight; 
		ate_weighted_unexp_events= events_unexp*ate_unexp_weight; 
		ate_weighted_exp= total_exp*ate_exp_weight; 
		ate_weighted_exp_events= events_exp*ate_exp_weight; 
		run;

	   data strata_n1;
  		merge strata_n (in=r) bounds_for_fs ( in=s keep = strata ps_min ps_max);
 		 by strata;
		if r; 
 		run;


		proc sql; create table strata_n2 as select 'Whole sample' as strata,  sum(total_Exp) as n_total_exp, sum(total_unexp) as n_total_unexp, 
		sum(events_Exp) as n_total_exp_events, sum(events_unexp) as n_total_unexp_events,
		round (sum (att_weighted_unexp)) as total_att_weighted_unexp, round (sum (att_weighted_unexp_events)) as total_att_weighted_unexp_events,
		round (sum (ate_weighted_unexp)) as total_ate_weighted_unexp, round (sum (ate_weighted_unexp_events)) as total_ate_weighted_unexp_events,
		round (sum (ate_weighted_exp)) as total_ate_weighted_exp, round (sum (ate_weighted_exp_events)) as total_ate_weighted_exp_events
		from strata_n1;
		quit; 

		*/

*** merge these weights with original cohort***;
		
proc sql;
create table with_weights as
select stratified.*, strata_n.att_unexp_weight, strata_n.ate_unexp_weight, strata_n.ate_exp_weight
from stratified
left join strata_n
on stratified.strata=strata_n.strata; 
quit;

%let estimand1 = %upcase(&estimand); 

%if "&estimand1" = "ATT" %then %do;

		*** these weights only apply to unexposed, every exposed person get the weight of 1 ***;

		data &out_data; set with_weights;
		psweight=att_unexp_weight; 
		if &exposure=1 and att_unexp_weight ne . then psweight=1; 
		drop att_unexp_weight ate_unexp_weight ate_exp_weight;
		run;

		title; 

%end;

%if "&estimand1" = "ATE" %then %do;

		*** these weights apply to both groups ***;

		data &out_data; set with_weights;
		psweight=ate_unexp_weight; 
		if &exposure=1 then psweight=ate_exp_weight; 
		drop att_unexp_weight ate_unexp_weight ate_exp_weight;
		run;

		title; 

%end;

*************************************************************************************************************

** step 5- Create the propensity score weighted Table 1 & calculate post-weighting c statistic
		   as a measure of improvement in balance between the two groups in the covariate vector**;

*************************************************************************************************************;

%table1 (in_for_table1= &out_data, treatment_var= &exposure, 
categorical_var_list= &ps_class_var_list, continuous_var_list= &ps_cont_var_list, weight=psweight, out_table1= Weighted_T1);

*title 'Post-weighting PS model'; 

ods graphics off;           
ods exclude all;             
ods output Association=save_c1; 

proc logistic data=&out_data desc rocoptions(weighted);
class &ps_class_var_list;
model &exposure= &ps_cont_var_list &ps_class_var_list; 
weight psweight;
run;

ods output clear; 

ods exclude none; 

data save_c1; set save_c1; where label2='c'; keep label2 nvalue2; run;

%if "&PS_provided1" = "NO" %then %do;
data balance; length measure $30; set save_c (in=r) save_c1 (in=s); 
if r=1 then measure='PS model c-statistic'; else measure='Post weighting c-statistic';
run;
%end; 

%if "&PS_provided1" = "YES" %then %do;
data balance; length measure $30; set save_c1; 
measure='Post weighting c-statistic';
run;
%end; 

*************************************************************************************************************

** step 6- Prepare datasets for Diagnostic plots **;

*************************************************************************************************************;

proc sort data=&out_data out=for_boxplot; by &exposure; run;

data for_balance_plot; length sample $12.;
set crude_t1 (in=r) weighted_t1 (in=s);
if r=1 then sample='Unweighted';
if s=1 then sample='PS weighted';
run;

proc sql noprint; select count (variable) into: total_rows from for_balance_plot; quit; 

title;

*************************************************************************************************************

** step 7- Estimate the effect of treatment using weighted genmod or weighted PH model

*************************************************************************************************************;

** generalized linear model with identity link for risk differences; 

%let effect_estimate1 = %upcase(&effect_estimate); %if "&effect_estimate1" = "RD" %then %do;

* crude and weighted event counts; 

proc sql; create table crude_counts as 
select &exposure, sum (&outcome) as n_crude_events, count (&id_var) as n_crude_total, (calculated n_crude_events/calculated n_crude_total)*100 as crude_risk_per100,
(cinv(0.025,2*calculated n_crude_events)/(2*calculated n_crude_total))*100 as crude_risk_lcl,
(cinv(0.975,2*(1+calculated n_crude_events))/(2*calculated n_crude_total))*100 as crude_risk_ucl
from &in_data
group by &exposure;
quit;

proc sql; create table weighted_counts as 
select &exposure, round(sum (&outcome*psweight)) as n_weighted_events, sum(psweight) as n_weighted_total, (calculated n_weighted_events/calculated n_weighted_total)*100 as weighted_risk_per100,
(cinv(0.025,2*calculated n_weighted_events)/(2*calculated n_weighted_total))*100 as weighted_risk_lcl,
(cinv(0.975,2*(1+calculated n_weighted_events))/(2*calculated n_weighted_total))*100 as weighted_risk_ucl
from &out_data
group by &exposure;
quit;

data event_counts; merge crude_counts weighted_counts; by &exposure; run;

*crude analysis; 

title 'Unweighted outcome model'; 

ods output ParameterEstimates=crude_beta;
proc genmod data=&in_data desc;
class &exposure/param=ref ref=first;
model &outcome= &exposure/dist=bin link=identity;
run;
ods output clear; 

data result1; set crude_beta; where parameter not in ('Intercept','Scale'); 
RD=(Estimate); LCL=(Estimate-1.96*StdErr); UCL=(Estimate+1.96*StdErr);
keep parameter estimate stderr RD LCL UCL;
run;


*ps-adjusted analysis; 

title 'PS-weighted outcome model'; 

ods output GEEEmpPEst=weighted_beta ;
proc genmod data=&out_data desc;
class &exposure &id_var/param=ref ref=first;
model &outcome= &exposure/dist=bin link=identity;
weight psweight;
repeated subject= &id_var/ type=ind;
run;
ods output clear; 

data result2; set weighted_beta; where parm not in ('Intercept','Scale'); 
RD=(Estimate); LCL=(Estimate-1.96*StdErr); UCL=(Estimate+1.96*StdErr);
keep parameter estimate stderr RD LCL UCL;
run;

data results; length method $20;  set result1 (in=r) result2; 
if r=1 then Method= 'Unweighted'; else Method='PS-weighted'; 
run;

title; 

%end;

** generalized linear model with log link for risk ratios; 

%let effect_estimate1 = %upcase(&effect_estimate); %if "&effect_estimate1" = "RR" %then %do;

* Unweighted and weighted event counts; 

proc sql; create table crude_counts as 
select &exposure, sum (&outcome) as n_crude_events, count (&id_var) as n_crude_total, (calculated n_crude_events/calculated n_crude_total)*100 as crude_risk_per100,
(cinv(0.025,2*calculated n_crude_events)/(2*calculated n_crude_total))*100 as crude_risk_lcl,
(cinv(0.975,2*(1+calculated n_crude_events))/(2*calculated n_crude_total))*100 as crude_risk_ucl
from &in_data
group by &exposure;
quit;

proc sql; create table weighted_counts as 
select &exposure, round(sum (&outcome*psweight)) as n_weighted_events, sum(psweight) as n_weighted_total, (calculated n_weighted_events/calculated n_weighted_total)*100 as weighted_risk_per100,
(cinv(0.025,2*calculated n_weighted_events)/(2*calculated n_weighted_total))*100 as weighted_risk_lcl,
(cinv(0.975,2*(1+calculated n_weighted_events))/(2*calculated n_weighted_total))*100 as weighted_risk_ucl
from &out_data
group by &exposure;
quit;

data event_counts; merge crude_counts weighted_counts; by &exposure; run;

** crude analysis; 

title 'Unweighted outcome model'; 

ods output ParameterEstimates=crude_beta;
proc genmod data=&in_data desc;
class &exposure/param=ref ref=first;
model &outcome= &exposure/dist=bin link=log;
run;
ods output clear; 

data result1; set crude_beta; where parameter not in ('Intercept','Scale'); 
RR=exp(Estimate); LCL=exp(Estimate-1.96*StdErr); UCL=exp(Estimate+1.96*StdErr);
keep parm estimate stderr RR LCL UCL;
run;

** ps adjusted; 

title 'PS-weighted outcome model'; 

ods output GEEEmpPEst=weighted_beta ;
proc genmod data=&out_data desc;
class &exposure &id_var/param=ref ref=first;
model &outcome= &exposure/dist=bin link=log;
weight psweight;
repeated subject= &id_var/ type=ind;
run;
ods output clear; 

data result2; set weighted_beta; where parm not in ('Intercept','Scale'); 
RR=exp(Estimate); LCL=exp(Estimate-1.96*StdErr); UCL=exp(Estimate+1.96*StdErr);
keep parm estimate stderr RR LCL UCL;
run;

data results; length method $20;  set result1 (in=r) result2; 
if r=1 then Method= 'Unweighted'; else Method='PS-weighted'; 
run;

title;

%end;

** generalized linear model with logit link for odds ratios; 

%let effect_estimate1 = %upcase(&effect_estimate); %if "&effect_estimate1" = "OR" %then %do;

* crude and weighted event counts; 

proc sql; create table crude_counts as 
select &exposure, sum (&outcome) as n_crude_events, count (&id_var) as n_crude_total, (calculated n_crude_events/calculated n_crude_total)*100 as crude_risk_per100,
(cinv(0.025,2*calculated n_crude_events)/(2*calculated n_crude_total))*100 as crude_risk_lcl,
(cinv(0.975,2*(1+calculated n_crude_events))/(2*calculated n_crude_total))*100 as crude_risk_ucl
from &in_data
group by &exposure;
quit;

proc sql; create table weighted_counts as 
select &exposure, round(sum (&outcome*psweight)) as n_weighted_events, sum(psweight) as n_weighted_total, (calculated n_weighted_events/calculated n_weighted_total)*100 as weighted_risk_per100,
(cinv(0.025,2*calculated n_weighted_events)/(2*calculated n_weighted_total))*100 as weighted_risk_lcl,
(cinv(0.975,2*(1+calculated n_weighted_events))/(2*calculated n_weighted_total))*100 as weighted_risk_ucl
from &out_data
group by &exposure;
quit;

data event_counts; merge crude_counts weighted_counts; by &exposure; run;

**Crude analysis;

title 'Unweighted outcome model'; 

ods output ParameterEstimates=crude_beta;
proc genmod data=&in_data desc;
class &exposure/param=ref ref=first;
model &outcome= &exposure/dist=bin link=logit;
run;
ods output clear; 

data result1; set crude_beta; where parameter not in ('Intercept','Scale'); 
OR=exp(Estimate); LCL=exp(Estimate-1.96*StdErr); UCL=exp(Estimate+1.96*StdErr);
keep parameter estimate stderr OR LCL UCL;
run;

** ps adjusted; 

title 'PS-weighted outcome model'; 
ods output GEEEmpPEst =weighted_beta ;
proc genmod data=&out_data desc;
class &exposure &id_var/param=ref ref=first;
model &outcome= &exposure/dist=bin link=logit;
weight psweight;
repeated subject= &id_var/ type=ind;
run;
ods output clear; 

data result2; set weighted_beta; where parm not in ('Intercept','Scale'); 
OR=exp(Estimate); LCL=exp(Estimate-1.96*StdErr); UCL=exp(Estimate+1.96*StdErr);
keep parm estimate stderr OR LCL UCL;
run;

data results; length method $20;  set result1 (in=r) result2; 
if r=1 then Method= 'Unweighted'; else Method='PS-weighted'; 
run;

title; 

%end;

** cox proportional hazard regression model with time to event data for hazard ratios; 

%let effect_estimate1 = %upcase(&effect_estimate); %if "&effect_estimate1" = "HR" %then %do;

* crude and weighted event counts; 

%let time_unit1 = %upcase(&time_unit); %if "&time_unit1" = "DAYS" %then %do;
 
proc sql; create table crude_counts as 
select &exposure,
sum (&outcome) as n_crude_events, round(sum(&survival_time/365)) as cumulative_pyears_crude, (calculated n_crude_events/calculated cumulative_pyears_crude)*100 as crude_IR_100py,
(cinv(0.025,2*calculated n_crude_events)/(2*calculated cumulative_pyears_crude))*100 as crude_IR_lcl,
(cinv(0.975,2*(1+calculated n_crude_events))/(2*calculated cumulative_pyears_crude))*100 as crude_IR_ucl
from &in_data
group by &exposure;
quit;

proc sql; create table weighted_counts as 
select &exposure, 
round(sum (&outcome*psweight)) as n_weighted_events, round(sum(&survival_time*psweight)/365) as cumulative_pyears_weighted, 
(calculated n_weighted_events/calculated cumulative_pyears_weighted)*100 as weighted_IR_100py,
(cinv(0.025,2*calculated n_weighted_events)/(2*calculated cumulative_pyears_weighted))*100 as weighted_IR_lcl,
(cinv(0.975,2*(1+calculated n_weighted_events))/(2*calculated cumulative_pyears_weighted))*100 as weighted_IR_ucl
from &out_data
group by &exposure;
quit;

%end;

%let time_unit1 = %upcase(&time_unit); %if "&time_unit1" = "YEARS" %then %do;
 
proc sql; create table crude_counts as 
select &exposure,
sum (&outcome) as n_crude_events, round(sum(&survival_time)) as cumulative_pyears_crude, (calculated n_crude_events/calculated cumulative_pyears_crude)*100 as crude_IR_100py,
(cinv(0.025,2*calculated n_crude_events)/(2*calculated cumulative_pyears_crude))*100 as crude_IR_lcl,
(cinv(0.975,2*(1+calculated n_crude_events))/(2*calculated cumulative_pyears_crude))*100 as crude_IR_ucl
from &in_data
group by &exposure;
quit;

proc sql; create table weighted_counts as 
select &exposure, 
round(sum (&outcome*psweight)) as n_weighted_events, round(sum(&survival_time*psweight)) as cumulative_pyears_weighted, 
(calculated n_weighted_events/calculated cumulative_pyears_weighted)*100 as weighted_IR_100py,
(cinv(0.025,2*calculated n_weighted_events)/(2*calculated cumulative_pyears_weighted))*100 as weighted_IR_lcl,
(cinv(0.975,2*(1+calculated n_weighted_events))/(2*calculated cumulative_pyears_weighted))*100 as weighted_IR_ucl
from &out_data
group by &exposure;
quit;

%end;

data event_counts; merge crude_counts weighted_counts; by &exposure; run;

** crude analysis**;

title 'Unweighted outcome model'; 

ods output ParameterEstimates=crude_beta; 
  proc phreg data=&in_data;
  class &exposure/param=ref ref=first;
  model &survival_time*&outcome(0) = &exposure/rl;
  run;
ods output clear; 

** adjusted analysis **;

title 'PS-weighted outcome model'; 

ods output ParameterEstimates=weighted_beta; 

  proc phreg data=&out_data covs;
  class &exposure/param=ref ref=first;
     model &survival_time*&outcome(0) = &exposure/rl;
     weight psweight;
     run;
ods output clear; 

data results; length method $20;  set crude_beta (in=r) weighted_beta; 
if r=1 then Method= 'Unweighted'; else Method='PS-weighted'; 
drop chisq probchisq classval0 df;
run;

title; 

%end;

*************************************************************************************************************

** step 8- Save all the result outputs in an excel file

*************************************************************************************************************;

ods listing close;
ods excel file="&out_excel..xlsx" style=htmlblue;

title 'Non-overlap trimming summary';

ods excel options(SHEET_NAME='Non-overlap trimming summary');

   proc freq data=trim_summary ; table trimmed*&exposure/nopct norow; run; 
   title; 

title 'Stratification weighting summary';

ods excel options(SHEET_NAME='Stratification weighting summary');

%if "&estimand1" = "ATT" %then %do;

		proc print data=strata_n1 noobs label;
		var strata ps_min ps_max total_exp events_exp total_unexp events_unexp att_unexp_weight att_weighted_unexp att_weighted_unexp_events;
		label strata= "PS Stratum"
		ps_min= "PS Stratum boundary- low"
		ps_max= "PS Stratum boundary- high"
		total_exp= "Total n exposed"
		events_exp= "Event count in exposed"
		total_unexp= "Total n reference"
		events_unexp= "Event count in reference"
		att_unexp_weight = "Stratum-specific ATT weight in reference"
		att_weighted_unexp = "Total n reference after weighting"
		att_weighted_unexp_events= "Event count in reference after weighting";
		run; 

%end;

%if "&estimand1" = "ATE" %then %do;

		proc print data=strata_n1 noobs label;
		var strata ps_min ps_max total_exp events_exp total_unexp events_unexp
		ate_exp_weight ate_unexp_weight ate_weighted_exp ate_weighted_exp_events ate_weighted_unexp ate_weighted_unexp_events;
		label strata= "PS stratum"
		ps_min= "PS Stratum boundary- low"
		ps_max= "PS Stratum boundary- high"
		total_exp= "Total n exposed"
		events_exp= "Event count in exposed"
		total_unexp= "Total n reference"
		events_unexp= "Event count in reference"

		ate_exp_weight = "Stratum-specific ATE weight in exposed"
		ate_weighted_exp = "Total n exposed after weighting"
		ate_weighted_exp_events= "Event count in exposed after weighting"

		ate_unexp_weight = "Stratum-specific ATE weight in reference"
		ate_weighted_unexp = "Total n reference after weighting"
		ate_weighted_unexp_events= "Event count in reference after weighting";
		run; 

%end;
title 'Unweighted and weighted counts';

ods excel options(SHEET_NAME='Unweighted and weighted counts');

%if "&estimand1" = "ATT" %then %do;

		proc print data=strata_n2 noobs label;
		var strata n_total_exp n_total_exp_events n_total_unexp n_total_unexp_events total_att_weighted_unexp total_att_weighted_unexp_events;
		label strata= "Population"
		n_total_exp= "Total n exposed"
		n_total_exp_events= "Total exposed events"
		n_total_unexp= "Total n reference"
		n_total_unexp_events= "Total reference events"
		total_att_weighted_unexp= "Total weighted reference"
		total_att_weighted_unexp_events= "Total weighted reference events";
		run; 

%end;

%if "&estimand1" = "ATE" %then %do;

		proc print data=strata_n2 noobs label;
		var strata n_total_exp n_total_exp_events n_total_unexp n_total_unexp_events 
		total_ate_weighted_exp total_ate_weighted_exp_events total_ate_weighted_unexp total_ate_weighted_unexp_events;
		label strata= "Population"
		n_total_exp= "Total n exposed"
		n_total_exp_events= "Total exposed events"
		n_total_unexp= "Total n reference"
		n_total_unexp_events= "Total reference events"
		total_ate_weighted_exp= "Total weighted exposed"
		total_ate_weighted_exp_events= "Total weighted exposed events"
		total_ate_weighted_unexp= "Total weighted reference"
		total_ate_weighted_unexp_events= "Total weighted reference events";
		run; 
%end;


title 'Unweighted PS Distribution';

ods excel options(SHEET_NAME='Unweighted PS Distribution');

proc sgplot data=sample_ps noborder;
  histogram &ps_var / group=&exposure name='a' transparency=0.5;
  density &ps_var / group=&exposure;
  keylegend 'a' / location=inside position=topright across=1 noborder  title="Exposure";
  yaxis offsetmin=0 display=(noline noticks) grid label= "% patients";
  xaxis min=0 max=1 ;
run;

title 'Weighted PS Distribution';

ods excel options(SHEET_NAME='Weighted PS Distribution');

proc sgplot data=&out_data noborder;
  histogram &ps_var / group=&exposure weight=psweight name='a' transparency=0.5;
  density &ps_var / group=&exposure weight=psweight;
  keylegend 'a' / location=inside position=topright across=1 noborder  title="Exposure";
  yaxis offsetmin=0 display=(noline noticks) grid label= "% patients" ;
  xaxis min=0 max=1 ;
run;

title 'PS weight distribution';

ods excel options(SHEET_NAME= 'PS weight distribution');

proc sgplot data=for_boxplot;
  vbox psweight/ group=&exposure groupdisplay=cluster
    lineattrs=(pattern=solid) whiskerattrs=(pattern=solid); 
  keylegend / location=inside position=topright across=1;
run;


ods excel options(SHEET_NAME='Balance plot');

%if &total_rows<=60 %then %do;  ** for up to 30 variables; 

	proc sgplot data=for_balance_plot ;
	title 'Balance plot';
	dot variable/response=std_diff group=sample ;
	REFLINE -10 10 /axis=x lineattrs=(pattern=2);
	refline 0/axis=x lineattrs=(pattern=1);
    xaxis label='Standardized differences'   valueattrs=(size=8pt) min=-50 max=50;
    yaxis label='Variable' ;  
	run;
	quit;

%end; 

%if &total_rows>60 %then %do; 

ods graphics on / width=10in;
ods graphics on / height=20in;

	proc sgplot data=for_balance_plot ;
	title 'Balance plot';
	dot variable/response=std_diff group=sample ;
	REFLINE -10 10 /axis=x lineattrs=(pattern=2);
	refline 0/axis=x lineattrs=(pattern=1);
    xaxis label='Standardized differences'   valueattrs=(size=8pt) min=-50 max=50;
    yaxis label='Variable' valueattrs=(size=7pt); 
	run;
	quit;

ods graphics off; 
%end; 

title 'Unweighted Table 1';

ods excel options(SHEET_NAME='Unweighted Table 1');

proc print data=Crude_T1 noobs;
var variable variable_type Treatment_column Ref_column crude_diff std_diff;
run;

title 'Weighted Table 1';

ods excel options(SHEET_NAME='Weighted Table 1');

proc print data=Weighted_T1 noobs;
var variable variable_type Treatment_column Ref_column crude_diff std_diff;
run;

title 'Overall balance';

ods excel options(SHEET_NAME='Overall Balance');

proc print data=balance noobs;
run;

title 'Incidence rates';

ods excel options(SHEET_NAME='Incidence rates');

proc print data=event_counts noobs;
run;

title 'Effect estimates';

ods excel options(SHEET_NAME='Effect estimates');

proc print data=results noobs;
run;

ods excel close;

title; 


%let work_lib1 = %upcase(&work_lib); %if "&work_lib1" = "CLEAN" %then %do;

proc datasets lib=work nolist kill;
quit;
run;

%end;

%mend; 
