//
//  code_check.cpp
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-06-17.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#include "code_check.h"


void CODECHECK_mandatory()
{
	/// Mandatory code checks
	/// and also used for visualization of input parameters
	
	// Population
	
	Population P;
	unsigned long founder_size = 500;
	double founder_femprop = 0.5;
	double founder_cswprop = 0.01;
	string folder_inputs = "../inputs/";
	P.setup_for_simulation(founder_size,
						   founder_femprop,
						   founder_cswprop,
						   folder_inputs,
						    "in_STI.csv",
						    "in_STI_SFincrease.csv",
						    "in_HIVrebound.csv",
						    "in_STItreatment.csv",
						    "in_STI_vaccine.csv",
						   false);
	
	// Partnerships
	CODECHECK_partnerships(P);
	CODECHECK_sexActivity(P);
	// Infectivity curves plots (used in documentation)
	CODECHECK_STI_infectivity_curves();
	
	
	
	string cmdstring = "cd " + _DIR_VIZINPUT + " ; Rscript visualize_ALL.R; cd ..";
	system(cmdstring.c_str());
	
	system("Rscript ./CODECHECK/codecheck_STI_IC.R");
}


void CODECHECK_population_growth()
{
	/// CHECK IF POPULATION GROWTH IS AS EXPECTED
	
	
	// Initialize empty Population object
	
	Population P(0);
	
	// Initial population STI free (and will remain so, just interested in population growth)
	
	P.initFromFile("./CODECHECK/startPopulation_noSTI.csv",   //./CODECHECK/startPopulation_noSTI.csv
				   "in_STI.csv",
				   "in_STI_SFincrease.csv","in_HIVrebound.csv");
	
	// Set all parameters
	P.setAllParameters("../inputs");
	

	
	// Re-define demographic parameters
	// Only births, no deaths
	
	double target_growth = 0.05;
	P.set_birthRate(target_growth);
	
	P.set_infantMortality(0.0);
	P.set_childMortality(0.0);
	
	P.set_deathParam_Weibull_scale(0.0);
	P.set_deathParam_Weibull_shape(0.0);
	P.set_deathParam_Weibull_scale_hiv(0.0);
	P.set_deathParam_Weibull_shape_hiv(0.0);
	
	P.displayInfo(false);
	
	// initial size of the population
	double P_init = P.get_size();
	
	// Horizon of the simulation
	double horizon = 20.0;
	double timestep = 1.0/12.0;
	
	
	// Simulation
	Simulation S(horizon, timestep, P, 0);
	
	bool logIndivInfo = false;
	bool dosex = false;
	bool traceNetwork = false;
	int displayProgress = 0;
	
	// Form partnerships first (no sex acts)
	
	S.set_horizon(horizon);
	S.runAllEvents_horizon(dosex,logIndivInfo,
						   traceNetwork,displayProgress,1);
	
	
	double P_final = S.get_population().census_alive();
	
	double growth_rate = (P_final/P_init);
	
	
	cout << " Target growth = "<<exp(target_growth*horizon)<<endl;
	cout << " Simulated growth = "<<growth_rate<<endl;
	
	double diff = (exp(target_growth*horizon)-growth_rate)/growth_rate;
	
	string msg = (diff<0.02)?"OK":"Problem!";
	cout << endl << "CODECHECK_population_growth: " << msg <<endl;
	
}





void CODECHECK_STI_infectivity_curves()
{
	
	/// CHECK ALL INFECTIVITY CURVES
	/// will be plotted in R
	
	ofstream f("./CODECHECK/codecheck_STI_IC.out");
	
	double timestep = 1.0 / 365.0;
	
	Gender defaultgender = female;
	
	// --- Syphilis ---
	
	cout << "CODECHECK: IC Tp (syphilis)"<<endl;
	
	STI s(Tp,"./in_STI.csv");
	
	s.set_is_recurrent(true);
	s.set_is_recurrent_2(true);
	s.set_is_symptomatic(false);

	for (double t=0; t<2; t+=timestep)
		f << "Tp," << t<<","<<s.infectivityCurve(t,defaultgender)<<endl;
	
	s.set_is_symptomatic(true);
	
	for (double t=0; t<2; t+=timestep)
		f << "Tp_condy," << t<<","<<s.infectivityCurve(t,defaultgender)<<endl;

	
	// --- Ct ---
	
	cout << "CODECHECK: IC Ct"<<endl;
	
	STI ct(Ct,"./in_STI.csv");
	
	ct.set_is_symptomatic(false);
	
	for (double t=0; t<2; t+=timestep)
		f << "Ct asymptom," << t<<","<<ct.infectivityCurve(t,defaultgender)<<endl;
	
	ct.set_is_symptomatic(true);
	
	for (double t=0; t<2; t+=timestep)
		f << "Ct symptom," << t<<","<<ct.infectivityCurve(t,defaultgender)<<endl;
	
	
	// --- Ng ---
	
	cout << "CODECHECK: IC Ng"<<endl;
	
	STI ng(Ng,"./in_STI.csv");
	
	//ng.set_is_symptomatic(true);
	
	for (double t=0; t<2; t+=timestep)
		f << "Ng," << t<<","<<ng.infectivityCurve(t,defaultgender)<<endl;
	
	
	// --- Hd ---
	
	cout << "CODECHECK: IC Hd"<<endl;
	
	STI hd(Hd,"./in_STI.csv");
	
	hd.set_is_symptomatic(false);
	
	for (double t=0; t<1; t+=timestep)
		f << "Hd asymptom," << t<<","<<hd.infectivityCurve(t,defaultgender)<<endl;
	
	hd.set_is_symptomatic(true);
	
	for (double t=0; t<1; t+=timestep)
		f << "Hd symptom," << t<<","<<hd.infectivityCurve(t,defaultgender)<<endl;

	
	
	// --- HPV ---
	
	cout << "CODECHECK: IC HPV"<<endl;
	
	STI hpv(HPV,"./in_STI.csv");
	
	hpv.set_is_recurrent(false);
	
	for (double t=0; t<1; t+=timestep)
		f << "HPV not recurrent," << t<<","<<hpv.infectivityCurve(t,defaultgender)<<endl;

	hpv.set_is_recurrent(true);
	
	for (double t=0; t<3; t+=timestep)
		f << "HPV recurrent," << t<<","<<hpv.infectivityCurve(t,defaultgender)<<endl;

	
	// --- HSV2 ---
	
	cout << "CODECHECK: IC HSV2"<<endl;
	
	STI hsv2(HSV2,"./in_STI.csv");
	
	hsv2.set_is_symptomatic(false);
	
	for (double t=0; t<5; t+=timestep)
		f << "HSV2," << t<<","<<hsv2.infectivityCurve(t,defaultgender)<<endl;
	
	hsv2.set_is_symptomatic(true);
	
	for (double t=0; t<5; t+=timestep)
		f << "HSV2 symptomatic," << t<<","<<hsv2.infectivityCurve(t,defaultgender)<<endl;
	
	
	// --- HIV ---
	
	cout << "CODECHECK: IC HIV"<<endl;
	
	STI hiv(HIV,"./in_STI.csv");
	
	for (double t=0; t<12; t+=timestep)
		f << "HIV," << t<<","<<hiv.infectivityCurve(t,defaultgender)<<endl;
	
	
	// --- Tv ---
	
	cout << "CODECHECK: IC Tv"<<endl;
	
	STI tv(Tv,"./in_STI.csv");
	Gender gmale = male;
	
	tv.set_is_symptomatic(false);
	
	for (double t=0; t<2; t+=timestep)
		f << "Tv female asympt," << t<<","<<tv.infectivityCurve(t,defaultgender)<<endl;

	for (double t=0; t<2; t+=timestep)
		f << "Tv male asympt," << t<<","<<tv.infectivityCurve(t,gmale)<<endl;
	
	
	tv.set_is_symptomatic(true);
	
	for (double t=0; t<2; t+=timestep)
		f << "Tv female symptom," << t<<","<<tv.infectivityCurve(t,defaultgender)<<endl;
	
	for (double t=0; t<2; t+=timestep)
		f << "Tv male symptom," << t<<","<<tv.infectivityCurve(t,gmale)<<endl;

}


void CODECHECK_partnerships(Population &P)
{
	
	// --- Formation: age
	
	vector<double> age = vector_seq(9.0, 80.0, 200);
	vector<double> fage(age.size());
	
	for (int i=0; i<age.size(); i++) {
		fage[i] = P.form_age(age[i]);
	}
	
	dcMatrix Mfage(age);
	Mfage.addColVector(fage);
	Mfage.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_form_age.out");
	
	
	// --- Formation: age gap
	
	vector<double> gap = vector_seq(-30, 50, 200);
	vector<double> fagegap(gap.size());
	
	for (int i=0; i<gap.size(); i++) {
		fagegap[i] = P.form_ageGap(gap[i]);
	}
	
	dcMatrix Mgap(gap);
	Mgap.addColVector(fagegap);
	Mgap.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_form_agegap.out");

	
	// --- Formation: risk group
	
	int rmx = P.get_maxRiskGroup();
	vector<double> r1;
	vector<double> r2;
	vector<double> frisk;
	
	for (int i=0; i<=rmx; i++)
	{
		for (int j=0; j<=rmx; j++)
		{
			r1.push_back(i);
			r2.push_back(j);
			frisk.push_back(P.form_riskGroup(i, j));
		}
	}
	dcMatrix Mfrisk(r1);
	Mfrisk.addColVector(r2);
	Mfrisk.addColVector(frisk);
	Mfrisk.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_form_risk.out");

	
	// --- Formation: deficit
	
	vector<double> x1 = vector_seq(0.0, 1.0, 10);
	vector<double> x2 = vector_seq(0.0, 1.0, 10);
	
	vector<double> fdeficit;
	vector<double> xx1;
	vector<double> xx2;
	
	
	for (int i=0; i<x1.size(); i++)
	{
		for (int j=0; j<x2.size(); j++)
		{
			xx1.push_back(x1[i]);
			xx2.push_back(x2[j]);
			fdeficit.push_back(P.form_deficit(x1[i], x2[j]));
		}
	}
	dcMatrix Mfdeficit(xx1);
	Mfdeficit.addColVector(xx2);
	Mfdeficit.addColVector(fdeficit);
	Mfdeficit.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_form_deficit.out");
	
	
	
	// --- Spousal progress: age and age gap
	
	vector<double> agef = vector_seq(10.0, 80.0, 50);
	vector<double> agegap = vector_seq(-30, 50, 50);
	
	vector<double> sp_age;
	vector<double> agef_;
	vector<double> agegap_;
	
	
	for (int i=0; i<agef.size(); i++)
	{
		for (int j=0; j<agegap.size(); j++)
		{
			agef_.push_back(agef[i]);
			agegap_.push_back(agegap[j]);
			
			double s = expGauss(agef[i],P.get_spousalProgress_meanAge_f(), P.get_spousalProgress_varAge_f())
			*expGauss(agegap_[j],P.get_spousalProgress_meanGap(), P.get_spousalProgress_varGap());
			
			sp_age.push_back(s);
		}
	}
	dcMatrix Msp_age(agef_);
	Msp_age.addColVector(agegap_);
	Msp_age.addColVector(sp_age);
	Msp_age.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_sp_age.out");
	
	
	// --- Dissolution: age
	
	vector<double> dage;
	vector<double> age_1;
	vector<double> age_2;
	
	for (int i=0; i<age.size(); i++) {
		for (int j=0; j<age.size(); j++)
		{
			age_1.push_back(age[i]);
			age_2.push_back(age[j]);
			dage.push_back(P.dissolve_age_fct(age[i],age[j]));
		}
		
	}
	
	dcMatrix Dage(age_1);
	Dage.addColVector(age_2);
	Dage.addColVector(dage);
	Dage.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_dissol_age.out");
	
	
	// --- Dissolution: risk group
	
	r1.clear();
	r2.clear();
	vector<double> drisk;
	
	for (int i=0; i<=rmx; i++)
	{
		for (int j=0; j<=rmx; j++)
		{
			r1.push_back(i);
			r2.push_back(j);
			drisk.push_back(P.dissolve_riskGroup_fct(i, j));
		}
	}
	dcMatrix Dfrisk(r1);
	Dfrisk.addColVector(r2);
	Dfrisk.addColVector(drisk);
	Dfrisk.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_dissol_risk.out");
	
	
	// --- Dissolution: partnership duration
	
	vector<double> pdur = vector_seq(0.0,70.0, 200);
	vector<double> ddur(age.size());
	
	for (int i=0; i<pdur.size(); i++) {
		ddur[i] = P.dissolve_duration_fct(pdur[i]);
	}
	
	dcMatrix Dpdur(pdur);
	Dpdur.addColVector(ddur);
	Dpdur.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_dissol_pduration.out");
	
	
	// --- Dissolution: deficit
	
	
	vector<double> ddeficit;
	xx1.clear();
	xx2.clear();
	
	
	for (int i=0; i<x1.size(); i++)
	{
		for (int j=0; j<x2.size(); j++)
		{
			xx1.push_back(x1[i]);
			xx2.push_back(x2[j]);
			ddeficit.push_back(P.dissolve_partnerDeficit_fct(x1[i], x2[j]));
		}
	}
	dcMatrix Ddeficit(xx1);
	Ddeficit.addColVector(xx2);
	Ddeficit.addColVector(ddeficit);
	Ddeficit.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_prtnr_dissol_deficit.out");
	
}



void CODECHECK_sexActivity(Population &P)
{
	vector<double> age = vector_seq(10.0, 80.0, 200);
	vector<double> rskgrp = vector_seq(0, 2, 3);
	vector<double> nprtn = vector_seq(0, 9, 10);
	
	// --- Age

	vector<double> fage(age.size());
	
	for (int i=0; i<age.size(); i++) {
		fage[i] = P.sexAct_reduce_age(age[i]);
	}
	
	dcMatrix Mage(age);
	Mage.addColVector(fage);
	Mage.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_age.out");

	
	// --- Risk
	
	vector<double> frisk(rskgrp.size());
	
	for (int i=0; i<rskgrp.size(); i++) {
		frisk[i] = P.sexAct_reduce_riskGroup((int)(rskgrp[i]));
	}
	
	dcMatrix Mrisk(rskgrp);
	Mrisk.addColVector(frisk);
	Mrisk.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_risk.out");

	
	// --- Number of partners
	
	vector<double> fnp(nprtn.size());
	
	for (int i=0; i<nprtn.size(); i++) {
		fnp[i] = P.sexAct_reduce_nPartners((int)(nprtn[i]));
	}
	
	dcMatrix Mnp(nprtn);
	Mnp.addColVector(fnp);
	Mnp.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_npartner.out");
	
	
	// --- Maximum number of concurrent partners
	
	vector<double> maxp(rskgrp.size());
	
	for (int i=0; i<rskgrp.size(); i++) {
		maxp[i] = P.proba_nMaxCurrSexPartner(male, (int)(rskgrp[i]));
	}
	
	dcMatrix Mmaxp(rskgrp);
	Mmaxp.addColVector(maxp);
	Mmaxp.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_maxpartner_male.out");
	
	
	maxp.clear();
	maxp.resize(rskgrp.size());
	
	for (int i=0; i<rskgrp.size(); i++) {
		maxp[i] = P.proba_nMaxCurrSexPartner(female, (int)(rskgrp[i]));
	}
	
	dcMatrix Mmaxp2(rskgrp);
	Mmaxp2.addColVector(maxp);
	Mmaxp2.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_maxpartner_female.out");

	
	// --- Probability sex with CSW
	
	vector<double> pcsw(rskgrp.size());
	
	for (int i=0; i<rskgrp.size(); i++) {
		pcsw[i] = P.probaSex_sexworker((int)(rskgrp[i]));
	}
	
	dcMatrix Mpcsw(rskgrp);
	Mpcsw.addColVector(pcsw);
	Mpcsw.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_probacsw.out");

	
	// --- Probability sex type 0
	
	vector<double> psex0;
	vector<double> r1;
	vector<double> r2;
	
	for (int i=0; i<=rskgrp.size(); i++) {
		for (int j=0; j<=rskgrp.size(); j++) {
			r1.push_back(i);
			r2.push_back(j);
			psex0.push_back(P.probaSex_type0(i, j));
		}
	}
	
	dcMatrix Msex0(r1);
	Msex0.addColVector(r2);
	Msex0.addColVector(psex0);
	Msex0.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_probatype0.out");
	
	// --- Probability sex type 1
	
	vector<double> psex1;
	r1.clear();
	r2.clear();
	
	for (int i=0; i<=rskgrp.size(); i++) {
		for (int j=0; j<=rskgrp.size(); j++) {
			r1.push_back(i);
			r2.push_back(j);
			psex1.push_back(P.probaSex_type1(i, j));
		}
	}
	
	dcMatrix Msex1(r1);
	Msex1.addColVector(r2);
	Msex1.addColVector(psex1);
	Msex1.WriteToFileCSV(_DIR_VIZINPUT + "viz_input_sexact_probatype1.out");

	
}




















