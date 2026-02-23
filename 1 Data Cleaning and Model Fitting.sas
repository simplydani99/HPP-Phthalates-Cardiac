/*****************************************************************************************************************************************
  Programmer: Danielle Stevens
  Purpose of syntax: Data curation of cardiac measures in HPP-3D Study
  Last updated on: 11/25/2025
******************************************************************************************************************************************/

libname origdat "J:\StevensLab\The HPP3D Study\Original Data";
libname hpp "J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Derived_data";

options nofmterr;

/*Read in dataset*/
data cardiac_1; set origdat.hpp_hpp3d_deid_20250813; run;
proc contents data=cardiac_1; run;
proc sort data=cardiac_1; by record ga_us; run;


/*Data cleaning*/
data cardiac_2; set cardiac_1;
/*Recalculate CRT*/
CARDIOTHORACIC_RATIO = (HEART_AREA/CHEST_AREA) * 100;

/*Re-calculate EFW using Hadlock's Formula*/
/*EFW equations have been slightly altered to accomodate conversion from mm to cm (data is provided in mm and original equations are in cm)*/
if 10 LE gaweeks then do;
bio_efw = 10**(1.5662 - (0.0108*(bio_headcircumf/10)) + (0.0468*(bio_abcircumf/10)) + 
                                 (0.171*(bio_femurlength/10)) + (0.00034*((bio_headcircumf/10)**2)) - 
                                 (0.003685*((bio_abcircumf/10)*(bio_femurlength/10)))); end;


/*Count missing cardiac US*/
miss_cout = nmiss(CARDIOTHORACIC_RATIO, MMODE_RV_TAPSE, MMODE_LV_MAPSE, MMODE_RV_FS, MMODE_LV_FS, STD_RV_TEI, STD_LVENT_TEI, TRICUSPID_EARATIO, MITRAL_EARATIO);

/*Flag poor quality cardiac ultrasound*/
flag_US = 0;
if (ultrasound_unverified = 1 OR image_unverified = 1) AND miss_cout ne 9 then flag_US = 1;

/*Flag missing gestational age at cardiac ultrasound*/
miss_ga = 0;
if miss_cout ne 9 AND ga_us = . then miss_ga = 1;

where visit in (1,4,7); /*Limit to when cardiac ultrasounds were performed*/
run;


/*Examine exclusions*/
proc freq data=cardiac_2;
tables visit * (miss_cout) / missing; /*8 missing at visit 1, 44 missing at visit 4, 53 missing at visit 7*/
tables miss_cout flag_US miss_ga / missing;
tables miss_cout*miss_ga / missing; /*None missing gestational age*/
tables miss_cout*flag_us / missing; /*10 poor-quality US (separate from those missing cardiac US)*/
run;
proc freq data=cardiac_2 noprint;
tables record / out=out;
run; *303 participants, 909 visits;


/*Examine covariates Enrolled Sample*/
proc sort data=cardiac_2;
by record visit;
run;
data first; set cardiac_2;
by record visit;
if first.record then output;
run;
proc freq data=first;
tables sex_nb drugs_new smoking drugs_new alcohol_new  insurance_r 
     clinicsite highest_edu_r race_eth any_prior_preg race other_race ethnicity employment income marital_status_r / missing;
run;
proc freq data=first;
tables race*other_race ethnicity / missing;
where race_eth = 4;
run;
proc means data=first;
var age bmi ;
run;


/*Limit to the analytic dataset*/
data cardiac_analytic; set cardiac_2;
where miss_cout ne 9 AND flag_us = 0 AND miss_ga = 0;
run;


/*Count number of participants*/
proc freq data=cardiac_analytic noprint;
tables record / out=out;
run;
proc freq data=out;
tables count;
run;
/*n=295
  21 participants (10%) have 1 visit
  49 participants (16%) have 2 visits
  225 participants (74%) have 3 visits
*/

proc freq data=cardiac_analytic;
tables visit ;
run;
/*289 (36.4%) participants have data at visit 1
  256 (32.2%) participants have data at visit 4
  249 (31.4%) participants have data at visit 7
Total of 794 long'l assessments
*/
 

/*Examine covariates: Analytic Sample*/
proc sort data=cardiac_analytic;
by record visit;
run;
data first; set cardiac_analytic;
by record visit;
if first.record then output;
run;
proc freq data=first;
tables sex_nb drugs_new smoking drugs_new alcohol_new  insurance_r 
     clinicsite highest_edu_r race_eth any_prior_preg race other_race ethnicity employment income marital_status_r / missing;
run;
proc freq data=first;
tables race*other_race ethnicity / missing;
where race_eth = 4;
run;
proc means data=first;
var age bmi ;
run;


/*Get z-scores for annulus displacement, chest area, heart area*/
ods csv file="J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Derived_data\Referent for zscores.csv";
proc means data=cardiac_2;
class gaweeks;
var mmode_rv_tapse mmode_lv_mapse HEART_AREA CHEST_AREA rv_inlet_length lv_inlet_length mitral_e mitral_a tricuspid_e tricuspid_a; *non-normally distributed cardiac measures;
where cohort = 1; *low risk cohort;
output out=meandat n(mmode_rv_tapse)=n_rv_tapse n(mmode_lv_mapse)=n_lv_mapse n(heart_area)=n_heart_area n(chest_area)=n_chest_area n(rv_inlet_length)=n_rv_inlet_length n(lv_inlet_length)=n_lv_inlet_length n(mitral_e)=n_mitral_e n(mitral_a)=n_mitral_a n(tricuspid_e)=n_tricuspid_e n(tricuspid_a)=n_tricuspid_a n(aortic_bfvelocity)=n_aortic_bfvelocity n(pulm_bfvelocity)=n_pulm_bfvelocity  
                   mean(mmode_rv_tapse)=mean_rv_tapse mean(mmode_lv_mapse)=mean_lv_mapse mean(heart_area)=mean_heart_area mean(chest_area)=mean_chest_area mean(rv_inlet_length)=mean_rv_inlet_length mean(lv_inlet_length)=mean_lv_inlet_length mean(mitral_e)=mean_mitral_e mean(mitral_a)=mean_mitral_a mean(tricuspid_e)=mean_tricuspid_e mean(tricuspid_a)=mean_tricuspid_a mean(aortic_bfvelocity)=mean_aortic_bfvelocity mean(pulm_bfvelocity)=mean_pulm_bfvelocity  
                   std(mmode_rv_tapse)=sd_rv_tapse std(mmode_lv_mapse)=sd_lv_mapse std(heart_area)=sd_heart_area std(chest_area)=sd_chest_area std(rv_inlet_length)=sd_rv_inlet_length std(lv_inlet_length)=sd_lv_inlet_length std(mitral_e)=sd_mitral_e std(mitral_a)=sd_mitral_a std(tricuspid_e)=sd_tricuspid_e std(tricuspid_a)=sd_tricuspid_a std(aortic_bfvelocity)=sd_aortic_bfvelocity std(pulm_bfvelocity)=sd_pulm_bfvelocity;
run;
ods csv close;



/*Ok so I'm going to exclude any with <5 observations*/
proc sort data=cardiac_analytic; by gaweeks; run;
data meandat; set meandat; if gaweeks=. then DELETE; run;
data cardiac_analytic; merge cardiac_analytic meandat; by gaweeks; run;
proc print data=cardiac_analytic (obs=200);
var record visit gaweeks ga_us mmode_rv_tapse n_rv_tapse mean_rv_tapse sd_rv_tapse;
run;

data cardiac_analytic; set cardiac_analytic;

/*z-score annulus displacement, heart, chest, area - base on z-scores from low risk group*/
z_rv_tapse = .;
  if n_rv_tapse GE 5 then z_rv_tapse = (mmode_rv_tapse - mean_rv_tapse) / sd_rv_tapse;
z_lv_mapse = .;
  if n_lv_mapse GE 5 then z_lv_mapse = (mmode_lv_mapse - mean_lv_mapse) / sd_lv_mapse;
z_heart_area = .;
  if n_heart_area GE 5 then z_heart_area = (heart_area - mean_heart_area) / sd_heart_area;
z_chest_area = .;
  if n_chest_area GE 5 then z_chest_area = (chest_area - mean_chest_area) / sd_chest_area;
z_rv_inlet_length = .;
  if n_rv_inlet_length GE 5 then z_rv_inlet_length = (rv_inlet_length - mean_rv_inlet_length) / sd_rv_inlet_length;
z_lv_inlet_length = .;
  if n_lv_inlet_length GE 5 then z_lv_inlet_length = (lv_inlet_length - mean_lv_inlet_length) / sd_lv_inlet_length;
z_mitral_e = .;
  if n_mitral_e GE 5 then z_mitral_e = (mitral_e - mean_mitral_e) / sd_mitral_e;
z_mitral_a = .;
  if n_mitral_a GE 5 then z_mitral_a = (mitral_a - mean_mitral_a) / sd_mitral_a;
z_tricuspid_e = .;
  if n_tricuspid_e GE 5 then z_tricuspid_e = (tricuspid_e - mean_tricuspid_e) / sd_tricuspid_e;
z_tricuspid_a = .;
  if n_tricuspid_a GE 5 then z_tricuspid_a = (tricuspid_a - mean_tricuspid_a) / sd_tricuspid_a;

/*quadratic*/
ga_sq = ga_us*ga_us;

/*Cubic*/
ga_cub = ga_us*ga_us*ga_us;
run;

%macro plotmac(input);
proc sgplot data=cardiac_analytic;
scatter x=ga_us y=&input / markerattrs=(color=blue);
run;
%mend;
%plotmac(chest_area);
%plotmac(heart_area);
%plotmac(z_heart_area);
%plotmac(z_chest_area);
%plotmac(CARDIOTHORACIC_RATIO);

%plotmac(MMODE_RV_TAPSE);
%plotmac(MMODE_LV_MAPSE);
%plotmac(z_RV_TAPSE);
%plotmac(z_LV_MAPSE);

%plotmac(MMODE_RV_FS);
%plotmac(MMODE_LV_FS);

%plotmac(STD_RV_TEI);
%plotmac(STD_LVENT_TEI);

%plotmac(tricuspid_e);
%plotmac(tricuspid_a);
%plotmac(z_tricuspid_e);
%plotmac(z_tricuspid_a);
%plotmac(TRICUSPID_EARATIO);

%plotmac(mitral_e);
%plotmac(mitral_a);
%plotmac(z_mitral_e);
%plotmac(z_mitral_a);
%plotmac(MITRAL_EARATIO);

%plotmac(rv_inlet_length);
%plotmac(lv_inlet_length);
%plotmac(z_rv_inlet_length);
%plotmac(z_lv_inlet_length);
/*NOTE that, upon checking the distributions, the non-z-scored e and a values for the mitral and tricuspid valves are fairly normal with limited evidence of heteroskedasticity so I will stick with the non-z-scored values moving forward*/

/*Check out outliers*/
proc sort data=cardiac_analytic; by record ga_us; run;
proc print data=cardiac_analytic; 
where z_heart_area > 6 OR z_chest_area > 6 OR z_rv_tapse > 6 OR z_lv_mapse >6 OR z_rv_inlet_length > 6 OR z_lv_inlet_length>6;
run;

/*Flags 2 participants.
  Record 102 has elevated heart and chest areas 
  Record 177 has elevated inlet lengths

  No indication of congenital anomaly or other complication  and the US was not flagged as unreliable. 
  I am going to leave as is, given that it is usually multiple values that are abnormal.
  We will include a sensitivity analysis using the non-z-scored data */


/*Add a category variable for plots in R*/
%macro catvar(input);
proc freq data=cardiac_analytic noprint;
tables visit / missing out=out;
where &input ne .;
run;
proc sort data=cardiac_analytic;
by visit record;
run;
data out; set out;
cat_&input = "Visit "||TRIM(LEFT(PUT(visit, 9.0)))||" (n="||TRIM(LEFT(PUT(COUNT, 9.0)))||")";
run;
data cardiac_analytic; merge cardiac_analytic out;
by visit;
run;
%mend;
%catvar(chest_area);
%catvar(heart_area);
%catvar(z_heart_area);
%catvar(z_chest_area);
%catvar(CARDIOTHORACIC_RATIO);

%catvar(MMODE_RV_TAPSE);
%catvar(MMODE_LV_MAPSE);
%catvar(z_RV_TAPSE);
%catvar(z_LV_MAPSE);

%catvar(MMODE_RV_FS);
%catvar(MMODE_LV_FS);

%catvar(STD_RV_TEI);
%catvar(STD_LVENT_TEI);

%catvar(tricuspid_e);
%catvar(tricuspid_a);
%catvar(TRICUSPID_EARATIO);

%catvar(mitral_e);
%catvar(mitral_a);
%catvar(MITRAL_EARATIO);

%catvar(rv_inlet_length);
%catvar(lv_inlet_length);
%catvar(z_rv_inlet_length);
%catvar(z_lv_inlet_length);


/**********************/
/* View distributions */
/**********************/
%macro meantabs(data,input);
proc means data=&data noprint;
class visit;
var &input;
output out=out n=n mean=mean std=std;
run;
data out; set out;
mean_std = TRIM(LEFT(PUT(mean,9.2)))||" ± "||TRIM(LEFT(PUT(std,9.2)))||" (N="||TRIM(LEFT(PUT(n, 9.0)))||")";
where visit in(1,4,7,9);
run;
proc transpose data=out out=mean_&input prefix=visit_;
id visit;
var mean_std;
run;
data mean_&input; set mean_&input;
exposure = "&input";
run;

proc means data=&data noprint;
class visit;
var &input;
output out=out median=median q1=q1 q3=q3;
run;
data out; set out;
med_iqr = TRIM(LEFT(PUT(median,9.2)))||" ("||TRIM(LEFT(PUT(q1,9.2)))||", "||TRIM(LEFT(PUT(q3, 9.2)))||")";
where visit in(1,4,7,9);
run;
proc transpose data=out out=med_&input prefix=visit_;
id visit;
var med_iqr;
run;
data med_&input; set med_&input;
exposure = "&input";
run;

proc means data=&data noprint;
class visit;
var &input;
output out=out min=min max=max;
run;
data out; set out;
min_max = TRIM(LEFT(PUT(min,9.2)))||", "||TRIM(LEFT(PUT(max,9.2)));
where visit in(1,4,7,9);
run;
proc transpose data=out out=min_max_&input prefix=visit_;
id visit;
var min_max;
run;
data min_max_&input; set min_max_&input;
exposure = "&input";
run;
%mend;
%meantabs(cardiac_analytic, gaweeks);
%meantabs(cardiac_analytic, chest_area);
%meantabs(cardiac_analytic, heart_area);
%meantabs(cardiac_analytic, z_heart_area);
%meantabs(cardiac_analytic, z_chest_area);
%meantabs(cardiac_analytic, CARDIOTHORACIC_RATIO);

%meantabs(cardiac_analytic, MMODE_RV_TAPSE);
%meantabs(cardiac_analytic, MMODE_LV_MAPSE);
%meantabs(cardiac_analytic, z_RV_TAPSE);
%meantabs(cardiac_analytic, z_LV_MAPSE);

%meantabs(cardiac_analytic, MMODE_RV_FS);
%meantabs(cardiac_analytic, MMODE_LV_FS);

%meantabs(cardiac_analytic, STD_RV_TEI);
%meantabs(cardiac_analytic, STD_LVENT_TEI);

%meantabs(cardiac_analytic, tricuspid_e);
%meantabs(cardiac_analytic, tricuspid_a);
%meantabs(cardiac_analytic, TRICUSPID_EARATIO);

%meantabs(cardiac_analytic, mitral_e);
%meantabs(cardiac_analytic, mitral_a);
%meantabs(cardiac_analytic, MITRAL_EARATIO);

%meantabs(cardiac_analytic, rv_inlet_length);
%meantabs(cardiac_analytic, lv_inlet_length);
%meantabs(cardiac_analytic, z_rv_inlet_length);
%meantabs(cardiac_analytic, z_lv_inlet_length);


data cardiac_means; length exposure $ 40.; 
merge mean_gaweeks 
     mean_chest_area
mean_heart_area
mean_z_heart_area
mean_z_chest_area
mean_CARDIOTHORACIC_RATIO

mean_MMODE_RV_TAPSE
mean_MMODE_LV_MAPSE
mean_z_RV_TAPSE
mean_z_LV_MAPSE

mean_MMODE_RV_FS
mean_MMODE_LV_FS

mean_STD_RV_TEI
mean_STD_LVENT_TEI

mean_tricuspid_e
mean_tricuspid_a
mean_TRICUSPID_EARATIO

mean_mitral_e
mean_mitral_a
mean_MITRAL_EARATIO

mean_rv_inlet_length
mean_lv_inlet_length
mean_z_rv_inlet_length
mean_z_lv_inlet_length;
by exposure;
run;
data cardiac_meds; length exposure $ 40.; 
merge med_gaweeks  
     med_chest_area
med_heart_area
med_z_heart_area
med_z_chest_area
med_CARDIOTHORACIC_RATIO

med_MMODE_RV_TAPSE
med_MMODE_LV_MAPSE
med_z_RV_TAPSE
med_z_LV_MAPSE

med_MMODE_RV_FS
med_MMODE_LV_FS

med_STD_RV_TEI
med_STD_LVENT_TEI

med_tricuspid_e
med_tricuspid_a
med_TRICUSPID_EARATIO

med_mitral_e
med_mitral_a
med_MITRAL_EARATIO

med_rv_inlet_length
med_lv_inlet_length
med_z_rv_inlet_length
med_z_lv_inlet_length;
by exposure;
run;
data cardiac_min_max; length exposure $ 40.; 
merge min_max_gaweeks  
     min_max_chest_area
min_max_heart_area
min_max_z_heart_area
min_max_z_chest_area
min_max_CARDIOTHORACIC_RATIO

min_max_MMODE_RV_TAPSE
min_max_MMODE_LV_MAPSE
min_max_z_RV_TAPSE
min_max_z_LV_MAPSE

min_max_MMODE_RV_FS
min_max_MMODE_LV_FS

min_max_STD_RV_TEI
min_max_STD_LVENT_TEI

min_max_tricuspid_e
min_max_tricuspid_a
min_max_TRICUSPID_EARATIO

min_max_mitral_e
min_max_mitral_a
min_max_MITRAL_EARATIO

min_max_rv_inlet_length
min_max_lv_inlet_length
min_max_z_rv_inlet_length
min_max_z_lv_inlet_length;
by exposure;
run;

ods excel file = "J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Results\Distribution of Cardiac Outcomes in analytic dataset.xlsx";
proc print data=cardiac_means; run;
proc print data=cardiac_meds; run;
proc print data=cardiac_min_max; run;
ods excel close;




/***************/
/*Model fitting*/
/***************/
proc sort data=cardiac_analytic; by record ga_us; run;
proc print data=cardiac_analytic (obs=10);
var ga_us;
run;


/* Random intercept and slope models */
%let outcome_list1 = chest_area heart_area z_heart_area z_chest_area CARDIOTHORACIC_RATIO
MMODE_RV_TAPSE MMODE_LV_MAPSE z_RV_TAPSE z_LV_MAPSE
MMODE_RV_FS MMODE_LV_FS
STD_RV_TEI STD_LVENT_TEI
tricuspid_e tricuspid_a TRICUSPID_EARATIO
mitral_e mitral_a MITRAL_EARATIO
rv_inlet_length lv_inlet_length z_rv_inlet_length z_lv_inlet_length;

%macro modeltesting;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

ods output parameterestimates=parameterestimates;
proc glimmix data=cardiac_analytic;
class record ;
model &next_outcome = ga_us / solution ddfm=bw;
random intercept ga_us / subject=record type=UN ; 
run; 
quit;

ods output parameterestimates=parameterestimates;
proc glimmix data=cardiac_analytic;
class record ;
model &next_outcome = ga_us ga_sq/ solution ddfm=bw;
random intercept ga_us ga_sq / subject=record type=UN ; 
run; 
quit;
data quad; set parameterestimates;
if effect = "ga_sq" AND Probt < 0.05 then model_select = 1;
   else model_select = 0;
run;

ods output parameterestimates=parameterestimates;
proc glimmix data=cardiac_analytic;
class record ;
model &next_outcome = ga_us ga_sq ga_cub/ solution ddfm=bw;
random intercept ga_us ga_sq ga_cub / subject=record type=UN ; 
run; 
data cubic; set parameterestimates;
if effect = "ga_cub" AND Probt < 0.05 then model_select = 2;
   else model_select = 0;
run;

data &next_outcome; set quad cubic; outcome = "&next_outcome"; run;
%let i = %eval(&i + 1);
%end;
%mend;
%modeltesting;
/*While not all models converged, they did converge for those with significant gestational ages, so I think this these models are the best option*/

/*Next, let's test how well these models fit*/
/*This section of code creates the dataset for us to output our predicted values into*/
proc freq data=cardiac_analytic;
tables record / out=out;
run;
data out2; set out;
by record;
ga_us = 91; output; *week 13;
ga_us = 231; output; *week 33;
keep record ga_us;
run;
data out3;
merge out2 cardiac_analytic;
by record ga_us;
run;

%macro linear(outcome);
proc glimmix data=out3;
class record ;
model &outcome = ga_us / solution ddfm=bw;
random intercept ga_us / subject=record type=UN solution; 
output out=pred_&outcome pred(BLUP ILINK)=pred; /* output predicted values taking into account the random effects */
run; 
quit;
proc print data=pred_&outcome (obs=100);
var record visit ga_us &outcome pred;
run;
proc sgplot data=pred_&outcome;
scatter x=ga_us y=&outcome / markerattrs=(color=blue);
scatter x=ga_us y=pred;
run;
proc sgplot data=pred_&outcome;
series x=ga_us y=&outcome / group=record;
run;
proc sgplot data=pred_&outcome;
series x=ga_us y=pred / group=record ;
run;

proc transpose data=pred_&outcome out=&outcome;
by record;
var pred;
run;

data &outcome; set &outcome;
diff_&outcome = col3-col1;
keep record diff:;
run;

%mend;
%linear(z_rv_tapse); /*Doesn't look like the best fit, tbh. It really restricts the variability. Let's look at one more*/
%linear(MMODE_RV_FS); /*Doesn't look like the best fit either*/
/*Just these ones because they're the ones with a linear fit selected above*/


/*Next, let's test a simple bspline using the effect statement* - split at 21 weeks (147 days)*/
%macro modeltesting;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

proc glimmix data=out3;
class record ;
  effect spl = spline(ga_us / details basis=bspline knotmethod=list(147));
model &next_outcome = spl / solution ddfm=bw;
random intercept spl / subject=record type=UN solution; 
output out=pred_&next_outcome pred(BLUP ILINK)=pred; /* output predicted values taking into account the random effects */
run; 
quit;

proc sgplot data=pred_&next_outcome;
scatter x=ga_us y=&next_outcome / markerattrs=(color=blue);
scatter x=ga_us y=pred;
run;
proc sgplot data=pred_&next_outcome;
series x=ga_us y=&next_outcome / group=record;
run;
proc sgplot data=pred_&next_outcome;
series x=ga_us y=pred / group=record ;
run;
%let i = %eval(&i + 1);
%end;
%mend;
%modeltesting;
/* Could not fit tricuspid_ea and MMODE_RV_FS */



/*Manual spline at individual mid-point?*/
data vis4; set out3;
where visit=4;
vis4ga=ga_us;
keep record vis4ga;
run;
data out4; merge out3 vis4; by record; run;

data out4; set out4;
if vis4ga=. then vis4ga = 147;
spl1=ga_us;
spl2=(ga_us > vis4ga)*(ga_us-vis4ga);
run;

%macro modeltesting;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

proc glimmix data=out4;
class record ;
model &next_outcome = spl1 spl2 / solution ddfm=bw;
random intercept spl1 spl2  / subject=record type=UN; 
output out=pred_&next_outcome pred(BLUP ILINK)=pred_&next_outcome; /* output predicted values taking into account the random effects */
run; 
quit;

proc sgplot data=pred_&next_outcome;
scatter x=ga_us y=&next_outcome / markerattrs=(color=blue);
scatter x=ga_us y=pred_&next_outcome;
run;
proc sgplot data=pred_&next_outcome;
series x=ga_us y=&next_outcome / group=record;
run;
proc sgplot data=pred_&next_outcome;
series x=ga_us y=pred_&next_outcome / group=record ;
run;

data pred_&next_outcome; set pred_&next_outcome; 
where ga_us in (91,231);
keep record visit ga_us pred:;
run;
%let i = %eval(&i + 1);
%end;
%mend;
%modeltesting;
/* Though not always positive definite, all of these ran and look to be ok fits to the original data*/
/* Let's move forward with this spline for analyses */

/*Calculate weekly velocity*/
data mergepred; merge 
     pred_chest_area
pred_heart_area
pred_z_heart_area
pred_z_chest_area
pred_CARDIOTHORACIC_RATIO

pred_MMODE_RV_TAPSE
pred_MMODE_LV_MAPSE
pred_z_RV_TAPSE
pred_z_LV_MAPSE

pred_MMODE_RV_FS
pred_MMODE_LV_FS

pred_STD_RV_TEI
pred_STD_LVENT_TEI

pred_tricuspid_e
pred_tricuspid_a
pred_TRICUSPID_EARATIO

pred_mitral_e
pred_mitral_a
pred_MITRAL_EARATIO

pred_rv_inlet_length
pred_lv_inlet_length
pred_z_rv_inlet_length
pred_z_lv_inlet_length;
by record ga_us; run;

%macro modeltesting;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

data mergepred; set mergepred;
by record;
v_&next_outcome = pred_&next_outcome - lag(pred_&next_outcome);
if first.record then do; 
   v_&next_outcome = .; end;
run;

%let i = %eval(&i + 1);
%end;
%mend;
%modeltesting;

proc sort data=mergepred;
by record visit;
run;

/*Save first measurement as separate variable*/
data first_pred; set mergepred;
keep record pred:;
where ga_us = 91;
run;
%macro modeltesting;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

data first_pred; set first_pred;
first_&next_outcome = pred_&next_outcome;
run;

%let i = %eval(&i + 1);
%end;
%mend;
%modeltesting;

data mergepred; set mergepred;
by record visit;
if last.record then output;
keep record v_:;
run;

/*Analytic dataset to use in R*/
data hpp_cardiac_analytic; merge out4 mergepred first_pred(keep=record first:);
by record;
run;
data hpp_cardiac_analytic; set hpp_cardiac_analytic;
where visit in (1,4,7); /*Limit to when cardiac ultrasounds were performed*/
drop MEHHTP MBP MBzP MCNP MCOCH MCOP MCPP MECPP MECPTP MEHHP MEHP MEOHP MEP MHBP MHiBP MHiNCH MONP MiBP; /*Drop original variables - overwrote below with SGADJ and averaged variables from the gm_exp_long dataset*/
run;
proc print data=hpp_cardiac_analytic (obs=50);
var record visit ga_us gaweeks preg_mep tri1_mep tri2_mep tri3_mep cardiothoracic_ratio  v_cardiothoracic_ratio mehhtp;
run;
proc print data=origdat.gm_exp_long (obs=50); run;

/*Finally, merge in longform exposure data for secondary analyses*/
proc sort data=hpp_cardiac_analytic;
by record visit;
run;
proc sort data=origdat.gm_exp_long;
by record visit;
run;
proc contents data=hpp_cardiac_analytic; run;
proc contents data=origdat.gm_exp_long; run;
data hpp.hpp_cardiac_analytic; merge hpp_cardiac_analytic(in=aaa) origdat.gm_exp_long;
by record visit;
if aaa;
run;
proc print data=hpp.hpp_cardiac_analytic (obs=50);
var record visit ga_us gaweeks preg_mep tri1_mep tri2_mep tri3_mep mep cardiothoracic_ratio v_cardiothoracic_ratio ;
run;

/*Distribution of pregnancy outcomes*/
proc sort data=hpp.hpp_cardiac_analytic;
by record visit;
run;
data first; set hpp.hpp_cardiac_analytic;
by record visit;
if first.record then output;
where visit ne .;
run;
proc freq data=first;
tables sga lga preterm death: fgr_nb preeclampsia_mom gdm_mom ghypertension_mom/ missing;
run; 






/****************************************************************************************************/
/* Look at plots of changes over time                                                               */
/* Any associations with bw-for-ga_us?                                                              */
/* This analysis was conducted to understand how these measures change with pregnancy complications */
/****************************************************************************************************/



/*This section of code creates the dataset for us to output our predicted values into*/
proc freq data=hpp.hpp_cardiac_analytic;
tables record / out=out;
run;
data out2; set out;
by record;
do i=91 to 231; 
ga_us = i;
output; 
end;
keep record ga_us;
run;
data out3;
merge out2 hpp.hpp_cardiac_analytic(keep=chest_area heart_area z_heart_area z_chest_area CARDIOTHORACIC_RATIO
MMODE_RV_TAPSE MMODE_LV_MAPSE z_RV_TAPSE z_LV_MAPSE
MMODE_RV_FS MMODE_LV_FS
STD_RV_TEI STD_LVENT_TEI
tricuspid_e tricuspid_a TRICUSPID_EARATIO
mitral_e mitral_a MITRAL_EARATIO
rv_inlet_length lv_inlet_length z_rv_inlet_length z_lv_inlet_length ga_us visit record) ;
by record ga_us;
run;


/*Splines*/
data vis4; set out3;
where visit=4;
vis4ga=ga_us;
keep record vis4ga;
run;
data out4; merge out3 vis4; by record; run;
data out5; set out4;
if vis4ga=. then vis4ga = 147;
spl1=ga_us;
spl2=(ga_us > 147)*(ga_us-vis4ga);
run;
proc print data=out5(obs=400);
var record vis4: ga_us spl1 spl2;
run;


data cardiac_bw; merge out5 hpp.hpp_cardiac_analytic(keep=record SGA LGA bwga_zscore);
by record; run;
data cardiac_bw; set cardiac_bw;
bwga_zcat = "AGA";
  if SGA = 1 then bwga_zcat = "SGA";
  if LGA = 1 then bwga_zcat = "LGA";
  if bwga_zscore = . then bwga_zcat = " ";
run;


%macro bwmodels;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

proc sort data=cardiac_bw;
by record spl1;
run;

ods output estimates=est_&next_outcome;
proc glimmix data=cardiac_bw;
class record bwga_zcat(ref="AGA");
model &next_outcome = spl1 spl2 bwga_zcat / solution ddfm=bw;
random intercept spl1 spl2  / subject=record type=UN; 
lsmeans bwga_zcat / cl;
estimate 'SGA' bwga_zcat 0 1 -1 / cl; 
estimate 'LGA' bwga_zcat 1 0 -1 / cl;  
  output out=&next_outcome pred(BLUP ILINK)=fit_&next_outcome lcl(blup ilink)=lcl_&next_outcome ucl(blup ilink)=ucl_&next_outcome;  /* output predicted values taking into account the random effects */
run; 
quit;

data est_&next_outcome; set est_&next_outcome;
outcome = "&next_outcome";
run;

proc sort data=&next_outcome; by record ga_us; run;
%let i = %eval(&i + 1);
%end;
%mend;
%bwmodels;

proc print data=cardiothoracic_ratio (obs=50);
var record ga_us spl1 spl2 bwga_zcat cardiothoracic_ratio fit: ;
run;
proc print data=STD_RV_TEI (obs=500);
var record ga_us spl1 spl2 bwga_zcat STD_RV_TEI fit: lcl: ucl:;
run;

data hpp.hpp_cardiac_pred;  merge  chest_area heart_area z_heart_area z_chest_area CARDIOTHORACIC_RATIO
MMODE_RV_TAPSE MMODE_LV_MAPSE z_RV_TAPSE z_LV_MAPSE
MMODE_RV_FS MMODE_LV_FS
STD_RV_TEI STD_LVENT_TEI
tricuspid_e tricuspid_a TRICUSPID_EARATIO
mitral_e mitral_a MITRAL_EARATIO
rv_inlet_length lv_inlet_length z_rv_inlet_length z_lv_inlet_length; by record ga_us; run;
proc print data=hpp.hpp_cardiac_pred (obs=500);
var record ga_us spl1 spl2 bwga_zcat STD_RV_TEI fit: lcl: ucl:;
run;


data bwga_zcat_diffs;  set  est_CARDIOTHORACIC_RATIO
est_chest_area est_heart_area est_z_heart_area est_z_chest_area 
est_MMODE_RV_TAPSE est_MMODE_LV_MAPSE est_z_RV_TAPSE est_z_LV_MAPSE
est_MMODE_RV_FS est_MMODE_LV_FS
est_STD_RV_TEI est_STD_LVENT_TEI
est_tricuspid_e est_tricuspid_a est_TRICUSPID_EARATIO
est_mitral_e est_mitral_a est_MITRAL_EARATIO
est_rv_inlet_length est_lv_inlet_length est_z_rv_inlet_length est_z_lv_inlet_length; run;

ods csv file="J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Results\Differences in Cardiac Outcomes by Birthweight Zscores.csv";
proc print data=bwga_zcat_diffs; run;
ods csv close;
