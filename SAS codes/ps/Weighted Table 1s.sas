
%macro table1 (in_for_table1= , treatment_var= , categorical_var_list= , continuous_var_list= , weight=dummy_weight, out_table1= ); 

* we dont need to see each 2 by 2 table or proc means output, so suppressing it. this information is saved as sas datasets 
and provided in excel format later; 

ods graphics off;            /* or use the %ODSOff macro */
ods exclude all;             /* suspend all open destinations */

data s1; set &in_for_table1; dummy_weight=1; *create dummy weights ;
trt=&treatment_var; * rename treatment variable to a simpler form for use within this macro; 
run;

************************** Crude table- categorical covariates *****************************************;


data fortable; retain &categorical_var_list; set s1;
keep &categorical_var_list; 
run;

proc contents data = fortable 
out = vars(keep = name varnum)
noprint;
run; 

proc sort data=vars; by varnum; run; 

proc sql noprint; select max(varnum) into: n from vars; quit;

%do i= 1 %to &n;

proc sql noprint; select NAME into: variable from vars where varnum=&i; quit;

ods output CrossTabFreqs=tabs_&i; * see in reults explorer-> results tab of interest-> properties-> name;

proc freq data=s1;
tables (&variable)*trt/nopct norow chisq;
weight &weight; 
run;

ods output clear;
run;

proc sort data=tabs_&i;
by trt;
run;

data tabs_&i;
set tabs_&i; where trt ne . ; *delete totals;
run;

data tabs_&i; set tabs_&i; rename &variable=cat; run;
data tabs_&i; set tabs_&i; category=put(cat, 10.); drop cat _LABEL_; run;

	proc append base= freqs data = tabs_&i; 
	run;

%end;

data cat_freqs; set freqs; where ColPercent ne .; 
trt_cat=catx('_', 'trt_cat', trt);
run;

*The following statements delete the base version and rename the youngest historical version to the base version. This is done because otherwise proc append keeps 
adding output rows to your output dataset across different runs, which is often times not desirable;

proc datasets NOLIST;
 delete freqs (gennum=0);
quit;

data temp (index = (trt_cat));
  set cat_freqs;
  obs=_n_;
run;

data _null_ ;
  dcl hash hh   (             ) ;
  hh.definekey  ('k'          ) ;
  hh.definedata ('table', 'trt', 'category','frequency', 'colpercent', 'obs') ;
  hh.definedone () ;
  do k = 1 by 1 until ( last.trt_cat ) ;
    set temp;
    by trt_cat ;
    hh.add () ;
  end ;
  hh.output (dataset: trt_cat) ;
run ;

proc sort data= trt_cat_0; by table category; run; 
proc sort data= trt_cat_1; by table category; run; 

data t1a; merge trt_cat_0 (rename= (frequency= Ref_freq colpercent=Ref_percent)) 
trt_cat_1 (rename= (frequency= Treated_freq colpercent=Treated_percent)); 
by table category; 
drop trt;
run;

proc datasets NOLIST;
   modify t1a; 
     attrib _all_ label=' '; *removing labels to preserve ordering; 
run;

proc sort data=t1a; by obs ; run;

data t1a; set t1a; drop obs; run;

************************** Crude table- continuous covariates *****************************************;

data formeans; retain &continuous_var_list; set s1;
keep &continuous_var_list; 
run;

proc contents data = formeans
out = contvars(keep = name varnum)
noprint;
run; 

proc sql noprint; select max(varnum) into: n from contvars; quit;

%do i= 1 %to &n;

proc sql noprint; select NAME into: variable from contvars where varnum=&i; quit;

proc means data=s1 noprint;
var &variable;
class trt;
weight &weight; 
OUTPUT OUT=c1;
run;

data c1; set c1; where trt ne . and _stat_ in ('MEAN', 'STD'); run;

proc transpose data=c1 out=means_&i; 
by trt; *subject id;
		id _stat_; *prescription id: Goes to the heading of the column;
		var &variable; *variable that goes in the column;
		run;

data means_&i; set means_&i; category=put(_name_, 10.); drop _name_ _LABEL_; run;

	proc append base= means data = means_&i; 
	run;

%end;

data cont_means; set means;  
n_percent= catx('/', MEAN, STD);
trt_cat=catx('_', 'trt_cat', trt);
rename mean=frequency std=colpercent ;
run;

*The following statements delete the base version and rename the youngest historical version to the base version. This is done because otherwise proc append keeps 
adding output rows to your output dataset across different runs, which is often times not desirable;

proc datasets NOLIST;
 delete means (gennum=0);
quit;

data temp1 (index = (trt_cat));
  set cont_means;
run;

data _null_ ;
  dcl hash hh   (             ) ;
  hh.definekey  ('k'          ) ;
  hh.definedata ('trt', 'category','frequency', 'colpercent', 'n_percent') ;
  hh.definedone () ;
  do k = 1 by 1 until ( last.trt_cat ) ;
    set temp1;
    by trt_cat ;
    hh.add () ;
  end ;
  hh.output (dataset: trt_cat) ;
run ;

data t1b; merge trt_cat_0 (rename= (frequency= Ref_freq colpercent=Ref_percent n_percent= Reference_n_percent)) 
trt_cat_1 (rename= (frequency= Treated_freq colpercent=Treated_percent n_percent= Treated_n_percent)); 
drop trt;
run;

proc datasets NOLIST;
   modify t1b; 
     attrib _all_ label=' '; 
run;

data t1_combined; set t1b (in=r) t1a (in=s); length variable_type $25.;
if r=1 then variable_type='Continuous [mean(sd)]'; 
if s=1 then variable_type='Categorical [n(%)]';
run; 

** concatenate categorical and continuous covariate tables; 

data t1_combined1; set t1_combined;
pre=' (';
post=')';
rounded_ref_freq= round(ref_freq, 1); if table= '' then rounded_ref_freq=round(ref_freq, 0.1);
rounded_treated_freq= round(treated_freq, 1); if table= '' then rounded_treated_freq=round(treated_freq, 0.1);
rounded_ref_percent= round(ref_percent, 0.1);
rounded_treated_percent= round(treated_percent, 0.1);
** quantities needed for std diff calculation; 
pt=(treated_percent/100); pc= (ref_percent/100); ** proportion of treated and reference for categorial variables;
xt=(treated_freq); xc= (Ref_freq); *means in treated and reference for continuous variables; 
st2= treated_percent*treated_percent; sc2=ref_percent*ref_percent; 
*variances in treated and reference for continuous variables; 
run;

** calculate the crude and standardized difference and prepare columns for output; 

data t1_combined2; length table $50;
set t1_combined1;
length ref_column $25 Treatment_column $25;
crude_diff= round((treated_percent-ref_percent), 0.1); ** percent differences for categorical covariates; 
if table= '' then crude_diff= round((treated_freq-ref_freq), 0.1); ** mean differences for continuous covariates; 
std_diff= round((100*(pt-pc))/ (sqrt((pt*(1-pt)+ pc*(1-pc))/2)), 0.1); ** categorical covariates; 
if table= '' then std_diff= round((100*(xt-xc))/sqrt((st2+sc2)/2), 0.1); ** continuous covariates; 
Ref_column= cat(rounded_ref_freq, pre, rounded_ref_percent, post);
Treatment_column= cat(rounded_treated_freq, pre, rounded_treated_percent, post);
run;

data t1_combined2; set t1_combined2; 
cat=strip(category); 
v=tranwrd(table, "Table",'');
v=tranwrd(v,"* trt",'');
run; 

data t1_combined3; set t1_combined2; where cat ne "0"; run;

data t1_combined3; set t1_combined3; length variable $50.;
variable= cats (v, "=", cat);
if v= '' then variable=cat; 
run;

** total row; 

proc freq data=s1; table &treatment_var/out=total_row1; 
weight &weight; 
run; 

proc sql noprint; create table total_treat as select COUNT as Treatment_column1 from total_row1 where &treatment_var=1 ; quit; 

proc sql noprint; create table total_ref as select COUNT as ref_column1 from total_row1 where &treatment_var=0 ; quit; 

data total_row; merge total_treat total_ref; run; 

data total_row; set total_row; 
treatment_column  = put(treatment_column1, 25. -L);
ref_column  = put(ref_column1, 25. -L);
run;

data t1_combined4; set total_row t1_combined3;
drop treatment_column1 ref_column1;
run; 

** final table;

data &out_table1; set t1_combined4; 
keep variable variable_type Treated_freq Treated_percent Ref_freq Ref_percent Treatment_column Ref_column crude_diff std_diff;
if variable= '' then variable='Total';
run; 

ods exclude none; 

proc datasets lib=work nolist;
 delete tabs: means:;
quit;
run;

%mend;
