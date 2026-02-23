/*****************************************************************************************************************************************
  Programmer: Danielle Stevens
  Purpose of syntax: Single chemical analysis for cardiac measures in HPP-3D Study
  Last updated on: 12/08/2025
******************************************************************************************************************************************/


/*********************************/
/* Read in imputed dataset       */
/* This dataset was created in R */
/*********************************/
libname hpp "J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Derived_data";

proc import datafile="J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Derived_data\Imputed HPP Cardiac Dataset.csv"
            dbms=csv
            out=hpp.hpp_cardiac_imp
            replace;
     getnames=yes;
	 guessingrows=ALL;
run;



/****************************/
/* Single-chemical analysis */
/* Adjusted analysis        */
/****************************/
data cardiac_analytic; set hpp.hpp_cardiac_imp;
where _imp ne 0; /*Limit to the imputed dataset for unadjusted analysis*/
run; 
proc freq data=cardiac_analytic noprint; table record/ out=out; run; /*n=295 - correct dataset*/
proc sort data=cardiac_analytic;
by _imp record ga_us;
run;
proc glimmix data=cardiac_analytic;
class record gender_nb any_prior_preg race_new highest_edu smoking alcohol_new drugs_new clinicsite;
model CARDIOTHORACIC_RATIO = iqr_preg_dehp spl1 spl2
                  gender_nb age bmi any_prior_preg race_new highest_edu smoking alcohol_new drugs_new clinicsite
                   / solution ddfm=bw;
random intercept  spl1 spl2 / subject=record type=UN ; 
estimate 'iqr_preg_dehp' iqr_preg_dehp 1  /  cl;
ods output estimates=estimates; 
by _imp;
run; 
quit;
proc sort data=estimates; by label _imp; run;
ods output parameterestimates=estimates2;
proc mianalyze data=estimates;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data iqr_preg_dehp; set estimates2;
out1=TRIM(LEFT(PUT(estimate, 9.1)))||" ("||TRIM(LEFT(PUT(lclmean, 9.1)))||", "||TRIM(LEFT(PUT(uclmean, 9.1)))||")";
out2=TRIM(LEFT(PUT(estimate, 9.2)))||" ("||TRIM(LEFT(PUT(lclmean, 9.2)))||", "||TRIM(LEFT(PUT(uclmean, 9.2)))||")";
keep label out1 out2 Probt estimate lclmean uclmean; 
run;



/*Macro*/
%let exposure_list1 = iqr_preg_MEP iqr_preg_DNBP iqr_preg_DIBP iqr_preg_MBZP iqr_preg_MCPP iqr_preg_DEHP iqr_preg_DINP iqr_preg_MCNP iqr_preg_DEHTP iqr_preg_DINCH;
%macro loopexp(outcome);
%local i next_exposure;
%let i=1;
%do %while (%scan(&exposure_list1,&i) ne );
%let next_exposure = %scan(&exposure_list1,&i);

proc glimmix data=cardiac_analytic;
class record gender_nb any_prior_preg race_new highest_edu smoking alcohol_new drugs_new clinicsite;
model &outcome = &next_exposure spl1 spl2 
                  gender_nb age bmi any_prior_preg race_new highest_edu smoking alcohol_new drugs_new clinicsite
                   / solution ddfm=bw;
random intercept spl1 spl2 / subject=record type=UN ; 
estimate "&next_exposure" &next_exposure 1  /  cl;
ods output estimates=estimates; 
by _imp;
run; 
quit;
proc sort data=estimates; by label _imp; run;
ods output parameterestimates=estimates2;
proc mianalyze data=estimates;
	by label;
	modeleffects Estimate;
    stderr StdErr;
    run;
data &next_exposure; set estimates2;
out1=TRIM(LEFT(PUT(estimate, 9.1)))||" ("||TRIM(LEFT(PUT(lclmean, 9.1)))||", "||TRIM(LEFT(PUT(uclmean, 9.1)))||")";
out2=TRIM(LEFT(PUT(estimate, 9.2)))||" ("||TRIM(LEFT(PUT(lclmean, 9.2)))||", "||TRIM(LEFT(PUT(uclmean, 9.2)))||")";
out3=TRIM(LEFT(PUT(estimate, 9.3)))||" ("||TRIM(LEFT(PUT(lclmean, 9.3)))||", "||TRIM(LEFT(PUT(uclmean, 9.3)))||")";
keep label out1 out2 out3 Probt estimate lclmean uclmean; 
run;
%let i = %eval(&i + 1);
%end;
%mend;

%let outcome_list1 = chest_area heart_area z_heart_area z_chest_area CARDIOTHORACIC_RATIO
MMODE_RV_TAPSE MMODE_LV_MAPSE z_RV_TAPSE z_LV_MAPSE
MMODE_RV_FS MMODE_LV_FS
STD_RV_TEI STD_LVENT_TEI
tricuspid_e tricuspid_a TRICUSPID_EARATIO
mitral_e mitral_a MITRAL_EARATIO
rv_inlet_length lv_inlet_length z_rv_inlet_length z_lv_inlet_length;
%macro loopout;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

%loopexp(&next_outcome);
data &next_outcome; set iqr_preg_DINCH iqr_preg_MEP iqr_preg_DNBP iqr_preg_DIBP iqr_preg_MBZP iqr_preg_MCPP iqr_preg_DEHP iqr_preg_DINP iqr_preg_MCNP iqr_preg_DEHTP ; run;

ods csv file="J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Results\Single Chemical Models - &next_outcome..csv";
proc print data=&next_outcome noobs; run;
ods csv close;

%let i = %eval(&i + 1);
%end;
%mend;
%loopout;




/*Determine sample sizes*/
proc freq data=cardiac_analytic noprint;
tables record / out=out;
where _imp = 1 AND iqr_preg_dehp ne . AND cardiothoracic_ratio ne .;
run;
proc summary data=out;
var count;
output out=n_visit sum=;
run;
data n_visit; set n_visit;
n_visit = TRIM(LEFT(PUT(_FREQ_, 9.0)))||", "||TRIM(LEFT(PUT(COUNT, 9.0)));
exposure = "iqr_preg_dehp";
outcome = "cardiothoracic_ratio";
keep exposure outcome n_visit;
run;


%macro loopexp(outcome);
%local i next_exposure;
%let i=1;
%do %while (%scan(&exposure_list1,&i) ne );
%let next_exposure = %scan(&exposure_list1,&i);

proc freq data=cardiac_analytic noprint;
tables record / out=out;
where _imp = 1 AND &next_exposure ne . AND &outcome ne .;
run;
proc summary data=out;
var count;
output out=n_visit sum=;
run;
data &next_exposure; set n_visit;
n_visit = TRIM(LEFT(PUT(_FREQ_, 9.0)))||", "||TRIM(LEFT(PUT(COUNT, 9.0)));
exposure = "&next_exposure";
outcome = "&outcome";
keep exposure outcome n_visit;
run;

%let i = %eval(&i + 1);
%end;
%mend;
%macro loopout;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

%loopexp(&next_outcome);
data &next_outcome; set iqr_preg_DINCH iqr_preg_MEP iqr_preg_DNBP iqr_preg_DIBP iqr_preg_MBZP iqr_preg_MCPP iqr_preg_DEHP iqr_preg_DINP iqr_preg_MCNP iqr_preg_DEHTP ; run;

%let i = %eval(&i + 1);
%end;
%mend;
%loopout;


ods excel file="J:\StevensLab\The HPP3D Study\HPP3D Phthalates and Cardiac_DRS\Analysis\Results\Single Chemical Models N and Num Visits for Preg Average Models for all Outcomes.xlsx";
%macro loopout;
%local i next_outcome;
%let i=1;
%do %while (%scan(&outcome_list1,&i) ne );
%let next_outcome = %scan(&outcome_list1,&i);

proc print data=&next_outcome; run;

%let i = %eval(&i + 1);
%end;
%mend;
%loopout;
ods excel close;


