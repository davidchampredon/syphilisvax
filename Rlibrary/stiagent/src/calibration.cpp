/*
 *  calibration.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 2013-10-06.
 *  Copyright 2013 David Champredon. All rights reserved.
 *
 */

#include "calibration.h"


void calibration_set_all_target(Simulation& S,
								string target_filename_ageDistrib,
								string target_filename_ageGapDistrib,
								string target_filename_ageFirstSexDistrib_f,
								string target_filename_ageFirstSexDistrib_m,
								string target_filename_ageGapFirstSexSpouseDistrib_f,
								string target_filename_ageGapFirstSexSpouseDistrib_m,
								string target_filename_singleRatio,
								string target_filename_malesVisitCSW,
								string target_filename_ageMalesVisitCSWDistrib,
								string target_filename_nLifeSexPrtnrDistrib_f,
								string target_filename_nLifeSexPrtnrDistrib_m,
								string target_filename_STIprevalences_names,
								string target_filename_STIprevalences_global,
								string target_filename_STIprevalences_matrix,
								string target_filename_HIV_prev_age_f,
								string target_filename_HIV_prev_age_m)
{
	// ==========================================
	/// ==== RETRIEVE AND SET THE ALL TARGETS ====
	// ==========================================
	
	
	
	// age distribution target
	dcMatrix TAD = calib_getTarget_Distribution(target_filename_ageDistrib);
	S.set_target_ageDistribution(TAD);
	
	// age gap distribution target
	dcMatrix TAGD = calib_getTarget_AgeDistribution(target_filename_ageGapDistrib);
	S.set_target_ageGapDistribution(TAGD);
	
	// First sex (female) age distribution target
	dcMatrix TAFS_f = calib_getTarget_Distribution(target_filename_ageFirstSexDistrib_f);
	S.set_target_ageFirstSexDistribution_f(TAFS_f);
	
	// First sex (male) age distribution target
	dcMatrix TAFS_m = calib_getTarget_Distribution(target_filename_ageFirstSexDistrib_m);
	S.set_target_ageFirstSexDistribution_m(TAFS_m);
	
	// Gap First sex and first Spouse (female) age distribution target
	dcMatrix TGFSS_f = calib_getTarget_Distribution(target_filename_ageGapFirstSexSpouseDistrib_f);
	S.set_target_ageGapFirstSexSpouseDistribution_f(TGFSS_f);
	
	// Gap First sex and first Spouse (male) age distribution target
	dcMatrix TGFSS_m = calib_getTarget_Distribution(target_filename_ageGapFirstSexSpouseDistrib_m);
	S.set_target_ageGapFirstSexSpouseDistribution_m(TGFSS_m);
	
	
	// Single ratio target
	
	double TSR_f = getParameterFromFile("female", target_filename_singleRatio);
	double TSR_m = getParameterFromFile("male", target_filename_singleRatio);
	S.set_target_singleRatio_f(TSR_f);
	S.set_target_singleRatio_m(TSR_m);
	
	
	double TMVC;
	dcMatrix TAMVC;
	
	if (!target_filename_malesVisitCSW.empty())
	{	// Males visiting CSW
		TMVC = getParameterFromFile("prop.paidsex", target_filename_malesVisitCSW);
		S.set_target_malesVisitCSW(TMVC);
		
		// Age distribution of Males visiting CSW
		TAMVC = calib_getTarget_Distribution(target_filename_ageMalesVisitCSWDistrib);
		S.set_target_ageMalesVisitCSWDistribution(TAMVC);
	}
	
	
	dcMatrix TNLS_f;
	dcMatrix TNLS_m;
	
	if (!target_filename_nLifeSexPrtnrDistrib_f.empty())
	{
		// Number of lifetime sex partners distribution (females & males)
		TNLS_f = calib_getTarget_Distribution(target_filename_nLifeSexPrtnrDistrib_f);
		S.set_target_nLftSexPrtnrDistribution_f(TNLS_f);
		
		TNLS_m = calib_getTarget_Distribution(target_filename_nLifeSexPrtnrDistrib_m);
		S.set_target_nLftSexPrtnrDistribution_m(TNLS_m);
	}
	
	
	// =-=-=-=-=-=-=-=-=-=
	//   STI prevalences
	// =-=-=-=-=-=-=-=-=-=
	
	
	// Names of the STIs calibrated
	
	vector<string> stinames = calib_getTarget_STIprevalences_names(target_filename_STIprevalences_names);
	S.set_target_STIprevalence_names(StringToSTIname_vector(stinames));
	
	
	// STI Prevalence global (across risk groups, genders,etc.)
	vector<double> stiprev = calib_getTarget_STIprevalences_values(target_filename_STIprevalences_global);
	S.set_target_STIprevalence(stiprev);
	
	
	// STI Prevalence by risk group
	dcMatrix STIprevTargetRiskGroup(stinames.size(),
								  S.get_population().get_maxRiskGroup()+2);
	
	MatrixFromCSVfile(STIprevTargetRiskGroup,
					  target_filename_STIprevalences_matrix,
					  S.get_population().get_maxRiskGroup()+2);
	
	S.set_target_STIprevalence_by_riskGroup(STIprevTargetRiskGroup);
	
	// HIV prevalence by age and gender
	
	dcMatrix THIVAGE_f = calib_getTarget_AgeDistribution(target_filename_HIV_prev_age_f);
	dcMatrix THIVAGE_m = calib_getTarget_AgeDistribution(target_filename_HIV_prev_age_m);
	S.set_target_HIV_prev_age_f(THIVAGE_f);
	S.set_target_HIV_prev_age_m(THIVAGE_m);
	
	// -------------------------------------------------------------------------
	// -------------------------------------------------------------------------
	// For use outside C++ (optional):
	// Save target into files (double checks)
	
	TAD.WriteToFileCSV(_DIR_CALIB + "calib_AgeDistribution_target.out");
	TAGD.WriteToFileCSV(_DIR_CALIB + "calib_AgeGapDistribution_target.out");
	TAFS_f.WriteToFileCSV(_DIR_CALIB + "calib_AgeFirstSexDistribution_f_target.out");
	TAFS_m.WriteToFileCSV(_DIR_CALIB + "calib_AgeFirstSexDistribution_m_target.out");
	TGFSS_f.WriteToFileCSV(_DIR_CALIB + "calib_AgeGapFirstSexSpouseDistribution_f_target.out");
	TGFSS_m.WriteToFileCSV(_DIR_CALIB + "calib_AgeGapFirstSexSpouseDistribution_m_target.out");
	
	if (!target_filename_malesVisitCSW.empty())
		TAMVC.WriteToFileCSV(_DIR_CALIB + "calib_AgeMalesVisitCSWDistribution_target.out");
	
	
	if (!target_filename_nLifeSexPrtnrDistrib_f.empty())
	{
		TNLS_f.WriteToFileCSV(_DIR_CALIB + "calib_nLifeSexPrtnrDistribution_f_target.out");
		TNLS_m.WriteToFileCSV(_DIR_CALIB + "calib_nLifeSexPrtnrDistribution_m_target.out");
	}
	
	ofstream fSR_f(_DIR_CALIB + "calib_singleRatio_f_target.out");
	fSR_f<<TSR_f<<endl;
	
	ofstream fSR_m(_DIR_CALIB + "calib_singleRatio_m_target.out");
	fSR_m<<TSR_m<<endl;
	
	
	if (!target_filename_malesVisitCSW.empty())
	{
		ofstream fMVC(_DIR_CALIB + "calib_malesVisitCSW_target.out");
		fMVC<<TMVC<<endl;
	}
	
	vectorToFile(stinames, _DIR_CALIB + "calib_STI_names_target.out");
	vectorToFile(stiprev, _DIR_CALIB + "calib_STI_prevGlobal_target.out");
	STIprevTargetRiskGroup.WriteToFileCSV(_DIR_CALIB + "calib_STIprevalenceSRiskGroup_target.out");
	// -------------------------------------------------------------------------
	// -------------------------------------------------------------------------
	
	
}



void calibration_set_all_target_wrap(Simulation &S,
									 string target_filename_wrapper,
									 int DHSphase)
{
	
	/// WRAPS THE CALL OF 'calibration_set_all_target()' FOR CONVENIENCE
	
	
	// Retrieve file names of calibratino targets
	
	string target_filename_ageDistrib  = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageDistrib", target_filename_wrapper);
	string target_filename_ageGapDistrib = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageGapDistrib", target_filename_wrapper);
	string target_filename_ageFirstSexDistrib_f = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageFirstSexDistrib_f", target_filename_wrapper);
	string target_filename_ageFirstSexDistrib_m = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageFirstSexDistrib_m", target_filename_wrapper);
	string target_filename_ageGapFirstSexSpouseDistrib_f = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageGapFirstSexSpouseDistrib_f", target_filename_wrapper);
	string target_filename_ageGapFirstSexSpouseDistrib_m = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageGapFirstSexSpouseDistrib_m", target_filename_wrapper);
	string target_filename_singleRatio = _DIR_CALIB +  getParameterFromFile_string("target_filename_singleRatio", target_filename_wrapper);
	
	// Older DHS surveys may not have this data
	
	string target_filename_malesVisitCSW ="";
	string target_filename_ageMalesVisitCSWDistrib ="";
	string target_filename_nLifeSexPrtnrDistrib_f ="";
	string target_filename_nLifeSexPrtnrDistrib_m ="";
	
	if(DHSphase>4){
		
		target_filename_malesVisitCSW = _DIR_CALIB +  getParameterFromFile_string("target_filename_malesVisitCSW", target_filename_wrapper);
		target_filename_ageMalesVisitCSWDistrib = _DIR_CALIB +  getParameterFromFile_string("target_filename_ageMalesVisitCSWDistrib", target_filename_wrapper);
		target_filename_nLifeSexPrtnrDistrib_f = _DIR_CALIB +  getParameterFromFile_string("target_filename_nLifeSexPrtnrDistrib_f", target_filename_wrapper);
		target_filename_nLifeSexPrtnrDistrib_m = _DIR_CALIB +  getParameterFromFile_string("target_filename_nLifeSexPrtnrDistrib_m", target_filename_wrapper);
		
	}
	
	string target_filename_STIprevalences_names = _DIR_CALIB +  getParameterFromFile_string("target_filename_STIprevalences_names", target_filename_wrapper);
	string target_filename_STIprevalences_global = _DIR_CALIB +  getParameterFromFile_string("target_filename_STIprevalences_global", target_filename_wrapper);
	string target_filename_STIprevalences_matrix = _DIR_CALIB +  getParameterFromFile_string("target_filename_STIprevalences_matrix", target_filename_wrapper);
	
	
	string target_filename_HIV_prev_age_f = _DIR_CALIB + getParameterFromFile_string("target_filename_HIV_prev_age_f", target_filename_wrapper);
	string target_filename_HIV_prev_age_m = _DIR_CALIB + getParameterFromFile_string("target_filename_HIV_prev_age_m", target_filename_wrapper);
	
	// Set targets for Simulation object
	
	calibration_set_all_target(S,
							   target_filename_ageDistrib,
							   target_filename_ageGapDistrib,
							   target_filename_ageFirstSexDistrib_f,
							   target_filename_ageFirstSexDistrib_m,
							   target_filename_ageGapFirstSexSpouseDistrib_f,
							   target_filename_ageGapFirstSexSpouseDistrib_m,
							   target_filename_singleRatio,
							   target_filename_malesVisitCSW,
							   target_filename_ageMalesVisitCSWDistrib,
							   target_filename_nLifeSexPrtnrDistrib_f,
							   target_filename_nLifeSexPrtnrDistrib_m,
							   target_filename_STIprevalences_names,
							   target_filename_STIprevalences_global,
							   target_filename_STIprevalences_matrix,
							   target_filename_HIV_prev_age_f,
							   target_filename_HIV_prev_age_m);
}








dcDataFrame LHS_generated(vector<string> varname,
						  vector<double> Vmin,
						  vector<double> Vmax,
						  int samplingNb)
{
	/// GENERATES A DATA FRAME WITH
	/// SAMPLED VALUES FROM LATIN HYPERCUBE SAMPLING,
	/// WITH NAMED PARAMETERS
	
	dcMatrix M = LatinHypercubeSampling(Vmin,Vmax,samplingNb).transpose();
	
	dcDataFrame df(varname,M);
	
	return df;
}


dcDataFrame LHS_generate_all_samples2(string filename_limits,
									  int samplingNumber,
									  string filename_save)
{
	/// SAVE LHS SAMPLES INTO A CSV FILE
	
	/// COMPULSORY INPUT FILE FORMAT
	///
	/// FIRST COLUMN: PARAMETERS NAME
	/// SECOND COLUMN: MINIMUM VALUE
	/// THIRD COLUMN: MAXIMUM VALUE
	
	vector<string> varname;
	vector<double> vmin, vmax;
	vectorFromCSVfile_string(varname, filename_limits.c_str(), 1);
	vectorFromCSVfile(vmin, filename_limits.c_str(), 2);
	vectorFromCSVfile(vmax, filename_limits.c_str(), 3);
	
	dcDataFrame df = LHS_generated(varname,
								   vmin,
								   vmax,
								   samplingNumber);
	
	df.saveToCSV(filename_save, false);
	
	return(df);
}


void LHS_generate_all_samples(string filename_limits,
							  int samplingNumber,
							  string filename_save)
{
	/// SAVE LHS SAMPLES INTO A CSV FILE
	
	/// COMPULSORY INPUT FILE FORMAT
	///
	/// FIRST COLUMN: PARAMETERS NAME
	/// SECOND COLUMN: MINIMUM VALUE
	/// THIRD COLUMN: MAXIMUM VALUE
	
	vector<string> varname;
	vector<double> vmin, vmax;
	
	vectorFromCSVfile_string(varname, filename_limits.c_str(), 1);
	vectorFromCSVfile(vmin, filename_limits.c_str(), 2);
	vectorFromCSVfile(vmax, filename_limits.c_str(), 3);
	
	
	dcDataFrame df = LHS_generated(varname,
								   vmin,
								   vmax,
								   samplingNumber);
	
	df.saveToCSV(filename_save, false);
}


void split_LHS_limitsFile(string filename_limits, int n)
{
	/// SPLITS A FILE DEFINING THE LIMITS OF PARAMETER VALUES BEFORE LHS
	/// INTO SEVERAL FILES WITH LIMITS BEING A PARTITION OF THE ORIGINAL LIMIT
	/// i.e:
	/// [MIN;MAX]  -->  [MIN;A1] [A1;A2] .... [An-1;MAX]
	///
	
	
	// read original limits
	
	vector<string> varname;
	vector<double> vmin;
	vector<double> vmax;
	
	vectorFromCSVfile_string(varname, filename_limits.c_str(), 1);
	vectorFromCSVfile(vmin, filename_limits.c_str(), 2);
	vectorFromCSVfile(vmax, filename_limits.c_str(), 3);
	
	// calculate intermediate new limits
	// (store vectors into matrix)
	
	dcMatrix M(vmin);
	
	for (int i=1; i<n; i++)
	{
		vector<double> tmp;
		
		for (int k=0;k<vmin.size();k++)
			tmp.push_back(vmin[k]+(vmax[k]-vmin[k])*i/n);
		
		M.addColVector(tmp);
	}
	
	M.addColVector(vmax);
	
	// create new limit files
	
	string extname = "split";
	
	for (int c=0; c<M.nbCols-1; c++)
	{
		ofstream f(filename_limits+extname+int2string(c));
		for (int i=0; i<varname.size(); i++)
		{
			f<<varname[i]<<",";
			f<<M(i,c)<<",";
			f<<M(i,c+1)<<endl;
		}
		f.close();
	}
}




void calibration_LHS_fromFile(string target_filename_wrapper,
							  // Parameters (pre-sampled)
							  string param_LHS_DMG,
							  string param_LHS_FORM,
							  string param_LHS_SPOUSAL,
							  string param_LHS_DISSOL,
							  string param_LHS_SEXACT,
							  string param_LHS_STIFTR_V, // STI features - Viruses
							  string param_LHS_STIFTR_B, // STI features - Bacteria
							  string param_LHS_STIFTR_P, // STI features - Protozoa
							  string param_LHS_STIFTR_F, // STI features - Fungus
							  
							  int nRows,	//  nRows = nLHS/nJobs
							  int jobNumber,
							  Population initialPopulation,
							  double horizon, double timeStep,
							  int nMonteCarlo, int DHSphase)
{
	// ======================================================
	//
	/// SAVE THE DISTANCES (FROM ALL TARGETS)
	/// FOR EACH PARAMETERS COMBINATION.
	/// PARAMETERS SETS WERE PREVIOUSLY SAMPLED WITH LHS
	/// AND STORED IN FILES.
	
	/// THE FILE	WITH ALL DISTANCES (FOR EACH PARAM COMBINATION)
	/// WILL HAVE TO BE ANALYZED OUTSIDE C++ CODE (see calibration.R)
	
	/// Monte Carlo (average is calibrated)
	
	// ======================================================
	
	
	// Set up the simulation
	Simulation S(horizon, timeStep, initialPopulation, nMonteCarlo);
	
	cout << endl<<"Simulation object in 'calibration_LHS_fromFile':";
	S.displayInfo();// DEBUG
	
	// Retrieve targets from file and set them to Simulation object
	
	calibration_set_all_target_wrap(S, target_filename_wrapper,DHSphase);
	
	// ======================================
	// ====    SAMPLING PARAMETERS       ====
	// ======================================

	// File name where the calibration trials are saved
	
	string fname_DIST	= "SF_calib_distance_" ;
	string fname_AD		= "SF_ageDistribution_" ;
	string fname_AGD	= "SF_ageGapDistribution_";
	string fname_AFS_f	= "SF_ageFirstSexDistribution_f_";
	string fname_AFS_m	= "SF_ageFirstSexDistribution_m_";
	string fname_GFSS_f	= "SF_ageGapFirstSexSpouseDistribution_f_";
	string fname_GFSS_m	= "SF_ageGapFirstSexSpouseDistribution_m_";
	string fname_SR_f	= "SF_singleRatio_f_" ;
	string fname_SR_m	= "SF_singleRatio_m_" ;
	string fname_MVC	= "SF_malesVisitCSW_" ;
	string fname_AMVC	= "SF_ageMalesVisitCSW_" ;
	
	string fname_NLS_f	= "SF_nLifeSexPrtnr_f_";
	string fname_NLS_m	= "SF_nLifeSexPrtnr_m_";
	
	string fname_STI	= "SF_STI_prevalences_";
	
	
	
	// Sample the parameter space
	// within the lower and upper limits.
	// (Each row is a parameter set)
	
	dcMatrix LHS_DMG;		LHS_DMG.FromFile_Rows(param_LHS_DMG, nRows);
	dcMatrix LHS_FORM;	LHS_FORM.FromFile_Rows(param_LHS_FORM, nRows);
	dcMatrix LHS_SPOUSAL; LHS_SPOUSAL.FromFile_Rows(param_LHS_SPOUSAL, nRows);
	dcMatrix LHS_DISSOL;	LHS_DISSOL.FromFile_Rows(param_LHS_DISSOL, nRows);
	dcMatrix LHS_SEXACT;	LHS_SEXACT.FromFile_Rows(param_LHS_SEXACT, nRows);
	dcMatrix LHS_STIFTR_V;LHS_STIFTR_V.FromFile_Rows(param_LHS_STIFTR_V, nRows);
	dcMatrix LHS_STIFTR_B;LHS_STIFTR_B.FromFile_Rows(param_LHS_STIFTR_B, nRows);
	dcMatrix LHS_STIFTR_P;LHS_STIFTR_P.FromFile_Rows(param_LHS_STIFTR_P, nRows);
	dcMatrix LHS_STIFTR_F;LHS_STIFTR_F.FromFile_Rows(param_LHS_STIFTR_F, nRows);
	
	
	// Parameters used in each step of the loop
	int nParam_DMG		= LHS_DMG.getNbCols();
	int nParam_FORM		= LHS_FORM.getNbCols();
	int nParam_SPOUSAL	= LHS_SPOUSAL.getNbCols();
	int nParam_DISSOL	= LHS_DISSOL.getNbCols();
	int nParam_SEXACT	= LHS_SEXACT.getNbCols();
	int nParam_STIFTR_V	= LHS_STIFTR_V.getNbCols();
	int nParam_STIFTR_B	= LHS_STIFTR_B.getNbCols();
	int nParam_STIFTR_P	= LHS_STIFTR_P.getNbCols();
	int nParam_STIFTR_F	= LHS_STIFTR_F.getNbCols();
	
	vector<double> param_DMG(nParam_DMG);
	vector<double> param_FORM(nParam_FORM);
	vector<double> param_SPOUSAL(nParam_SPOUSAL);
	vector<double> param_DISSOL(nParam_DISSOL);
	vector<double> param_SEXACT(nParam_SEXACT);
	vector<double> param_STIFTR_V(nParam_STIFTR_V);
	vector<double> param_STIFTR_B(nParam_STIFTR_B);
	vector<double> param_STIFTR_P(nParam_STIFTR_P);
	vector<double> param_STIFTR_F(nParam_STIFTR_F);
	
	
	// Results of the objective function to minimize
	// at each step in the LHS loop
	vector<double> dist_k;
	dist_k.clear();
	
	
	// All response variables calculated
	// for each parameter set sampled will be
	// stored in these matrices (1 row=1 result)
	
	vector<double> SR_f(0);
	vector<double> SR_m(0);
	vector<double> MVC(0);
	
	dcMatrix AD(0);
	dcMatrix AGD(0);
	dcMatrix AFS_f(0);
	dcMatrix AFS_m(0);
	dcMatrix GFSS_f(0);
	dcMatrix GFSS_m(0);
	dcMatrix AMVC(0);
	dcMatrix NLS_f(0);
	dcMatrix NLS_m(0);
	
	dcMatrix stiPrev(0,S.get_population().get_maxRiskGroup()+2);
	
	// Retrieves bins from target (for histograms)
	vector<double> ageBreaks			= S.get_target_ageDistribution().extractColumn(0);
	vector<double> ageGapBreaks			= S.get_target_ageGapDistribution().extractColumn(0);
	vector<double> ageFirstSexBreaks_f	= S.get_target_ageFirstSexDistribution_f().extractColumn(0);
	vector<double> ageFirstSexBreaks_m	= S.get_target_ageFirstSexDistribution_m().extractColumn(0);
	vector<double> ageGapFirstSexSpouseBreaks_f	= S.get_target_ageGapFirstSexSpouseDistribution_f().extractColumn(0);
	vector<double> ageGapFirstSexSpouseBreaks_m	= S.get_target_ageGapFirstSexSpouseDistribution_m().extractColumn(0);
	vector<double> ageMalesVisitCSWBreaks		= S.get_target_ageMalesVisitCSWDistribution().extractColumn(0);
	
	vector<double> lifeSexPrtnrBreaks_f	= S.get_target_nLftSexPrtnrDistribution_f().extractColumn(0);
	vector<double> lifeSexPrtnrBreaks_m	= S.get_target_nLftSexPrtnrDistribution_m().extractColumn(0);
	
	
	ageBreaks.push_back(999.9); // makes sure the range is broad enough
	ageGapBreaks.push_back(999.9);
	ageFirstSexBreaks_f.push_back(999.9);
	ageFirstSexBreaks_m.push_back(999.9);
	ageGapFirstSexSpouseBreaks_f.push_back(999.9);
	ageGapFirstSexSpouseBreaks_m.push_back(999.9);
	ageMalesVisitCSWBreaks.push_back(999.9);
	lifeSexPrtnrBreaks_f.push_back(99999.99);
	lifeSexPrtnrBreaks_m.push_back(99999.99);
	
	
	// Loop through all parameters combinations
	// from the LHS
	
	for (int k=0; k< nRows; k++)
	{
		cout << endl << "Job #"<< jobNumber << " ; Calibration LHS: ";
		cout << k+1 <<"/"<<nRows <<endl;
		
		param_DMG		= LHS_DMG.extractRow(k);
		param_FORM		= LHS_FORM.extractRow(k);
		param_SPOUSAL	= LHS_SPOUSAL.extractRow(k);
		param_DISSOL	= LHS_DISSOL.extractRow(k);
		param_SEXACT	= LHS_SEXACT.extractRow(k);
		param_STIFTR_V	= LHS_STIFTR_V.extractRow(k);
		param_STIFTR_B	= LHS_STIFTR_B.extractRow(k);
		param_STIFTR_P	= LHS_STIFTR_P.extractRow(k);
		//param_STIFTR_F	= LHS_STIFTR_F.extractRow(k); // Not implemented yet
		
		// Reset to initial population size, etc
		S.set_population(initialPopulation);
		
		// Update new parameters
		S.calib_setParameters(param_DMG,
							  param_FORM,
							  param_SPOUSAL,
							  param_DISSOL);
		
		S.calib_setParameters_STI(param_SEXACT,
								  param_STIFTR_V,
								  param_STIFTR_B,
								  param_STIFTR_P,
								  param_STIFTR_F);
		
		// Run the simulation until 'horizon'
		// (Monte Carlo)
		dist_k.push_back(S.calib_distanceFromAllTargets_MC(nMonteCarlo));
		
		
		// Calculate the current parameters values
		// and store it in the matrix
		
		vector<double> tmp = S.get_population().census_ageDistribution(ageBreaks);
		AD.addRowVector(tmp);
		
		vector<double> tmp2 = S.get_population().census_ageGapDistribution(ageGapBreaks);
		AGD.addRowVector(tmp2);
		
		vector<double> tmp_afs_f = S.get_population().census_ageFirstSexDistribution(ageFirstSexBreaks_f,female);
		AFS_f.addRowVector(tmp_afs_f);
		
		vector<double> tmp_afs_m = S.get_population().census_ageFirstSexDistribution(ageFirstSexBreaks_m,male);
		AFS_m.addRowVector(tmp_afs_m);
		
		vector<double> tmp_gfss_f = S.get_population().census_ageGapFirstSexSpouseDistribution(ageGapFirstSexSpouseBreaks_f,female);
		GFSS_f.addRowVector(tmp_gfss_f);
		
		vector<double> tmp_gfss_m = S.get_population().census_ageGapFirstSexSpouseDistribution(ageGapFirstSexSpouseBreaks_m,male);
		GFSS_m.addRowVector(tmp_gfss_m);
		
		double maxDurationSinceLastVisit = 1.0;
		
		vector<double> tmp_amvc = S.get_population().census_ageMalesVisitCSWDistribution(ageMalesVisitCSWBreaks,
																						 maxDurationSinceLastVisit);
		AMVC.addRowVector(tmp_amvc);
		
		
		double tmp_SR_f = S.get_population().census_ratioSingles(female);
		SR_f.push_back(tmp_SR_f);
		
		double tmp_SR_m = S.get_population().census_ratioSingles(male);
		SR_m.push_back(tmp_SR_m);
		
		
		double tmp_MVC = S.get_population().census_maleVisitCSW(maxDurationSinceLastVisit);
		MVC.push_back(tmp_MVC);
		
		vector<double> tmp_nls_f = S.get_population().census_nLifeSexPrtnrDistrib(lifeSexPrtnrBreaks_f, female);
		NLS_f.addRowVector(tmp_nls_f);
		
		vector<double> tmp_nls_m = S.get_population().census_nLifeSexPrtnrDistrib(lifeSexPrtnrBreaks_m, male);
		NLS_m.addRowVector(tmp_nls_m);
		
		dcMatrix sti_tmp = S.get_population().STI_prevalence_by_riskGroup(S.get_target_STIprevalence_names());
		stiPrev = rowBind(stiPrev, sti_tmp);
		
	}
	
	// Save the results into files
	// (files will be analyzed by the serial farming reader)
	
	// Response variables calibrated on target
	
	AD.WriteToFileCSV(_DIR_CALIB + fname_AD + int2string(jobNumber) + ".out");
	AGD.WriteToFileCSV(_DIR_CALIB + fname_AGD + int2string(jobNumber) + ".out");
	AFS_f.WriteToFileCSV(_DIR_CALIB + fname_AFS_f + int2string(jobNumber) + ".out");
	AFS_m.WriteToFileCSV(_DIR_CALIB + fname_AFS_m + int2string(jobNumber) + ".out");
	
	GFSS_f.WriteToFileCSV(_DIR_CALIB + fname_GFSS_f + int2string(jobNumber) + ".out");
	GFSS_m.WriteToFileCSV(_DIR_CALIB + fname_GFSS_m + int2string(jobNumber) + ".out");
	
	vectorToFile(SR_f, _DIR_CALIB + fname_SR_f + int2string(jobNumber) + ".out");
	vectorToFile(SR_m, _DIR_CALIB + fname_SR_m + int2string(jobNumber) + ".out");
	
	vectorToFile(MVC, _DIR_CALIB + fname_MVC + int2string(jobNumber) + ".out");
	AMVC.WriteToFileCSV(_DIR_CALIB + fname_AMVC + int2string(jobNumber) + ".out");
	
	NLS_f.WriteToFileCSV(_DIR_CALIB + fname_NLS_f + int2string(jobNumber) + ".out");
	NLS_m.WriteToFileCSV(_DIR_CALIB + fname_NLS_m + int2string(jobNumber) + ".out");
	
	stiPrev.WriteToFileCSV(_DIR_CALIB + fname_STI + int2string(jobNumber) + ".out");
	
	// distance from target for each trial
	vectorToFile(dist_k, _DIR_CALIB + fname_DIST + int2string(jobNumber) + ".out");
}





void calib_LHS_SerialFarming_GenInputFiles(int nLHS, int nJobs, string fname,
										   // Demographic param limits
										   vector<double> param_DMG_LowerLim,
										   vector<double> param_DMG_UpperLim,
										   // Partnership formation param limits
										   vector<double> param_FORM_LowerLim,
										   vector<double> param_FORM_UpperLim,
										   // Spousal progression param limits
										   vector<double> param_SPOUSAL_LowerLim,
										   vector<double> param_SPOUSAL_UpperLim,
										   // Partnership dissolution param limits
										   vector<double> param_DISSOL_LowerLim,
										   vector<double> param_DISSOL_UpperLim,
										   // Sexual activity param limits
										   vector<double> param_SEXACT_LowerLim,
										   vector<double> param_SEXACT_UpperLim,
										   // STI virus param limits
										   vector<double> param_STIFTR_V_LowerLim,
										   vector<double> param_STIFTR_V_UpperLim,
										   // STI bacteria param limits
										   vector<double> param_STIFTR_B_LowerLim,
										   vector<double> param_STIFTR_B_UpperLim,
										   // STI protozoa param limits
										   vector<double> param_STIFTR_P_LowerLim,
										   vector<double> param_STIFTR_P_UpperLim,
										   // STI fungus param limits
										   vector<double> param_STIFTR_F_LowerLim,
										   vector<double> param_STIFTR_F_UpperLim)
{
	// GENERATES INPUT FILES FOR UNIT JOBS
	
	
	dcMatrix LHS_DMG		= LatinHypercubeSampling(param_DMG_LowerLim, param_DMG_UpperLim, nLHS);
	dcMatrix LHS_FORM		= LatinHypercubeSampling(param_FORM_LowerLim, param_FORM_UpperLim, nLHS);
	dcMatrix LHS_SPOUSAL	= LatinHypercubeSampling(param_SPOUSAL_LowerLim, param_SPOUSAL_UpperLim, nLHS);
	dcMatrix LHS_DISSOL	= LatinHypercubeSampling(param_DISSOL_LowerLim, param_DISSOL_UpperLim, nLHS);
	
	dcMatrix LHS_SEXACT	= LatinHypercubeSampling(param_SEXACT_LowerLim, param_SEXACT_UpperLim, nLHS);
	
	dcMatrix LHS_STIFTR_V	= LatinHypercubeSampling(param_STIFTR_V_LowerLim, param_STIFTR_V_UpperLim, nLHS);
	dcMatrix LHS_STIFTR_B	= LatinHypercubeSampling(param_STIFTR_B_LowerLim, param_STIFTR_B_UpperLim, nLHS);
	dcMatrix LHS_STIFTR_P	= LatinHypercubeSampling(param_STIFTR_P_LowerLim, param_STIFTR_P_UpperLim, nLHS);
	dcMatrix LHS_STIFTR_F	= LatinHypercubeSampling(param_STIFTR_F_LowerLim, param_STIFTR_F_UpperLim, nLHS);
	
	// DEBUG
	LHS_SPOUSAL.display();
	LHS_DISSOL.display();
	
	// number of single serial executions per jobs
	if (nLHS%nJobs)
	{
		cout << endl << "ERROR [calib_LHS_SerialFarming_GenInputFiles]: nJobs does not divide nLHS" <<endl;
		exit(1);
	}
	
	int ns = nLHS/nJobs;
	
	string cmdline = "rm -f " + fname + "*";
	system(cmdline.c_str());
	
	
	dcMatrix fDMG(0), fFORM(0), fSPOUSAL(0), fDISSOL(0);
	dcMatrix fSEXACT(0);
	dcMatrix fSTIFTR_V(0),fSTIFTR_B(0),fSTIFTR_P(0),fSTIFTR_F(0);
	int cnt = 1;
	
	for (int s=1; s<=nLHS; s++)
	{
		if (s%ns != 0)
		{
			// BUILD UP THE MATRIX FOR ONE JOB
			// MADE OF ns SERIAL EXECUTIONS
			fDMG.addRowVector(LHS_DMG.extractRow(s-1));
			fFORM.addRowVector(LHS_FORM.extractRow(s-1));
			fSPOUSAL.addRowVector(LHS_SPOUSAL.extractRow(s-1));
			fDISSOL.addRowVector(LHS_DISSOL.extractRow(s-1));
			
			fSEXACT.addRowVector(LHS_SEXACT.extractRow(s-1));
			
			fSTIFTR_V.addRowVector(LHS_STIFTR_V.extractRow(s-1));
			fSTIFTR_B.addRowVector(LHS_STIFTR_B.extractRow(s-1));
			fSTIFTR_P.addRowVector(LHS_STIFTR_P.extractRow(s-1));
			fSTIFTR_F.addRowVector(LHS_STIFTR_F.extractRow(s-1));
			
		}
		if (s%ns == 0)
		{
			fDMG.addRowVector(LHS_DMG.extractRow(s-1));
			fFORM.addRowVector(LHS_FORM.extractRow(s-1));
			fSPOUSAL.addRowVector(LHS_SPOUSAL.extractRow(s-1));
			fDISSOL.addRowVector(LHS_DISSOL.extractRow(s-1));
			
			fSEXACT.addRowVector(LHS_SEXACT.extractRow(s-1));
			
			fSTIFTR_V.addRowVector(LHS_STIFTR_V.extractRow(s-1));
			fSTIFTR_B.addRowVector(LHS_STIFTR_B.extractRow(s-1));
			fSTIFTR_P.addRowVector(LHS_STIFTR_P.extractRow(s-1));
			fSTIFTR_F.addRowVector(LHS_STIFTR_F.extractRow(s-1));
			
			string f_DMG	= fname + "DMG_" +  int2string(cnt) + ".csv";
			string f_FORM	= fname + "FORM_" + int2string(cnt) + ".csv";
			string f_SPOUSAL= fname + "SPOUSAL_" + int2string(cnt) + ".csv";
			string f_DISSOL	= fname + "DISSOL_" + int2string(cnt) + ".csv";
			
			string f_SEXACT	= fname + "SEXACT_" + int2string(cnt) + ".csv";
			
			string f_STRIFTR_V	= fname + "STIFTR_V_" + int2string(cnt) + ".csv";
			string f_STRIFTR_B	= fname + "STIFTR_B_" + int2string(cnt) + ".csv";
			string f_STRIFTR_P	= fname + "STIFTR_P_" + int2string(cnt) + ".csv";
			string f_STRIFTR_F	= fname + "STIFTR_F_" + int2string(cnt) + ".csv";
			
			cnt++;
			
			// The parameters for 'ns'
			// serial executions have been
			// stored in the matrix -> write to a file
			fDMG.WriteToFile(f_DMG);
			fFORM.WriteToFile(f_FORM);
			fSPOUSAL.WriteToFile(f_SPOUSAL);
			fDISSOL.WriteToFile(f_DISSOL);
			fSEXACT.WriteToFile(f_SEXACT);
			fSTIFTR_V.WriteToFile(f_STRIFTR_V);
			fSTIFTR_B.WriteToFile(f_STRIFTR_B);
			fSTIFTR_P.WriteToFile(f_STRIFTR_P);
			fSTIFTR_F.WriteToFile(f_STRIFTR_F);
			
			
			// re-initialize
			fDMG.resize(0);
			fFORM.resize(0);
			fSPOUSAL.resize(0);
			fDISSOL.resize(0);
			fSEXACT.resize(0);
			fSTIFTR_V.resize(0);
			fSTIFTR_B.resize(0);
			fSTIFTR_P.resize(0);
			fSTIFTR_F.resize(0);
		}
	}
}


dcMatrix calib_getTarget_Distribution(string filename)
{
	/// RETRIEVE TARGET DISTRIBUTION FOR CALIBRATION
	/// FORMAT EXPECTED: CSV
	/// 1ST COLUMN: BREAKS
	/// 2ND COLUMN: PROPORTIONS
	/// PROPORTION[i] for BREAK[i] <= x < BREAK[i+1]
	
	/// BREAKS MUST COVER THE FULL RANGE OF VALUES
	
	vector<double> breaks;
	vector<double> proportion;
	
	vectorFromCSVfile(breaks, filename.c_str(), 1);
	vectorFromCSVfile(proportion, filename.c_str(), 2);
	
	//cout<<endl<<"DEBUG calib_getTarget_Distribution";
	//displayVector(breaks);
	
	dcMatrix res(breaks);
	res.addColVector(proportion);
	
	return res;
}


dcMatrix calib_getTarget_AgeDistribution(string filename)
{
	/// RETRIEVE TARGET AGE DISTRIBUTION FOR CALIBRATION
	/// FORMAT EXPECTED: CSV
	/// 1ST COLUMN: AGE BUCKETS
	/// 2ND COLUMN: PROPORTION OF INDIV.
	/// PROPORTION[i] for AGE_BUCKET[i] <= age < AGE_BUCKET[i+1]
	
	vector<double> age_breaks;
	vector<double> age_proportion;
	
	vectorFromCSVfile(age_breaks, filename.c_str(), 1);
	vectorFromCSVfile(age_proportion, filename.c_str(), 2);
	
	//displayVector(age_breaks);
	
	dcMatrix res(age_breaks);
	res.addColVector(age_proportion);
	
	return res;
}



vector<string> calib_getTarget_STIprevalences_names(string filename)
{
	/// FIRST COLUMN OF FILE MUST BE STI NAMES
	/// THAT ARE CALIBRATED
	
	vector<string> stinames;
	vectorFromCSVfile_string(stinames, filename.c_str(), 1);
	return stinames;
}



vector<double> calib_getTarget_STIprevalences_values(string filename)
{
	/// FIRST COLUMN OF FILE MUST BE STI PREVALENCES
	/// THAT ARE CALIBRATED
	
	vector<double> prev;
	vectorFromCSVfile(prev, filename.c_str(), 1);
	return prev;
}



vector<vector<double> >	calib_get_parameters(string file_DMG,
											 string file_FORM,
											 string file_SPOUSAL,
											 string file_DISSOL)
{
	/// READ INITIAL VALUES OF PARAMETERS USED IN CALIBRATION
	/// * * * WARNING * * *
	/// THE ORDER IS IMPORTANT, AND *MUST* BE THE SAME AS IN
	/// FUNCTION 'Simulation::calib_setParameters'
	
	vector<vector<double> > calibparam(4); // size of vector same as number of files as input of this function
	
	// DEMOGRAPHIC PARAMETERS
	
	calibparam[0].resize(7); // size = as many parameters read just below
	// TO DO: do a 'push_back' instead of setting dimension beforehand????
	
	calibparam[0][0]	= getParameterFromFile("birthRate", file_DMG);
	calibparam[0][1]	= getParameterFromFile("infantMortality", file_DMG);
	calibparam[0][2]	= getParameterFromFile("childMortality", file_DMG);
	calibparam[0][3]	= getParameterFromFile("deathParam_Weibull_shape", file_DMG);
	calibparam[0][4]	= getParameterFromFile("deathParam_Weibull_scale", file_DMG);
	calibparam[0][5]	= getParameterFromFile("deathParam_Weibull_shape_hiv", file_DMG);
	calibparam[0][6]	= getParameterFromFile("deathParam_Weibull_scale_hiv", file_DMG);
	
	// FORMATION PARTNERSHIP PARAMETERS
	
	calibparam[1].clear();
	
	
	calibparam[1].push_back(getParameterFromFile("formation_MaxRate", file_FORM));
	
	calibparam[1].push_back(getParameterFromFile("formation_age_fullstart", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_age_pivot", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_age_shape", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_age_fmin", file_FORM));
	
	//	calibparam[1][1]	= getParameterFromFile("formation_meanAge_female", file_FORM);
	//	calibparam[1][2]	= getParameterFromFile("formation_varAge_female", file_FORM);
	calibparam[1].push_back(getParameterFromFile("formation_meanAgeGap", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_varAgeGap", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_correl_Age_AgeGap", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_shapeAge", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_RiskGroup_1", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_RiskGroup_2", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_PartnDeficit_1", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_STIsymptom_f", file_FORM));
	calibparam[1].push_back(getParameterFromFile("formation_STIsymptom_m", file_FORM));
	
	
	
	// SPOUSAL PROGRESSION
	
	calibparam[2].clear(); // size = as many parameters read just below
	
	
	calibparam[2].push_back(getParameterFromFile("spousalProgress_maxRate", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_meanAge_f", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_varAge_f", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_meanGap", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_varGap", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_durationK1", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_durationK2", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_meanDiffAgeGap", file_SPOUSAL));
	calibparam[2].push_back(getParameterFromFile("spousalProgress_varDiffAgeGap", file_SPOUSAL));
	
	// DISSOLUTION PARTNERSHIPS
	
	calibparam[3].resize(11); // size = as many parameters read just below
	
	int i=0;
	calibparam[3][i]	= getParameterFromFile("dissolution_MaxRate", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_spouse", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_RiskGroup_1", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_duration_1", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_duration_2", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_duration_3", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_age_mean", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_age_var", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_PartnerDeficit", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_ageConcPartn", file_DISSOL); i++;
	calibparam[3][i]	= getParameterFromFile("dissolution_STI_symptom", file_DISSOL); i++;
	
	return calibparam;
}


vector<vector<double> >	calib_get_parameters_STI(string file_SEXACT,
												 string file_STIFTR_V,
												 string file_STIFTR_B,
												 string file_STIFTR_P,
												 string file_STIFTR_F)
{
	// READ INITIAL VALUES OF PARAMETERS USED IN CALIBRATION
	// * * * WARNING * * *
	// THE ORDER IS IMPORTANT, AND *MUST* BE THE SAME AS IN
	// FUNCTION 'Simulation::calib_setParameters_STI'
	
	vector<vector<double> > calibparam(5); // size of vector same as number of files as input of this function

	int k = 0;
	
	// SEXUAL ACTIVITY
	
	calibparam[k].clear();
	
	calibparam[k].push_back(getParameterFromFile("_sexAct_maxRate_male", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_maxRate_female", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("prefSexWithSpouse", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_age_param_1", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_age_param_2", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_age_param_3", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_risk_param_1", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_STIsymptom_param_male", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_STIsymptom_param_female", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_reduce_nPartner_param", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_proba_distribute_partnerTypes_sexWorkParam_1", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_proba_distribute_partnerTypes_sexWorkParam_2", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_proba_sexWorker_param_1", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_proba_sexWorker_param_2", file_SEXACT));
	calibparam[k].push_back(getParameterFromFile("_sexAct_CostSexWork_reduction", file_SEXACT));
	
	
	// **** WARNING ****
	// ORDER IS IMPORTANT!
	// *******************
	
	// STI FEATURES - Virus
	
	k++;
	calibparam[k].clear();
	
	calibparam[k].push_back(getParameterFromFile("HIV_probaMaxSexTransm", file_STIFTR_V));
	calibparam[k].push_back(getParameterFromFile("HSV2_probaMaxSexTransm", file_STIFTR_V));
	calibparam[k].push_back(getParameterFromFile("HPV_probaMaxSexTransm", file_STIFTR_V));
	
	// STI FEATURES - Bacteria
	
	k++;
	calibparam[k].clear();
	
	calibparam[k].push_back(getParameterFromFile("Ct_probaMaxSexTransm", file_STIFTR_B));
	calibparam[k].push_back(getParameterFromFile("Ng_probaMaxSexTransm", file_STIFTR_B));
	calibparam[k].push_back(getParameterFromFile("Tp_probaMaxSexTransm", file_STIFTR_B));
	calibparam[k].push_back(getParameterFromFile("Hd_probaMaxSexTransm", file_STIFTR_B));
	calibparam[k].push_back(getParameterFromFile("Bv_probaMaxSexTransm", file_STIFTR_B));
	
	// STI FEATURES - PROTOZOA
	
	k++;
	calibparam[k].clear();
	
	calibparam[k].push_back(getParameterFromFile("Tv_probaMaxSexTransm", file_STIFTR_P));
	
	// STI FEATURES - Fungus : Not implemented yet
	
	return calibparam;
}



dcDataFrame LHS_explore(Population P,
						string limit_LHS_file_wrapper, //file name (wrapper) containing file names of LHS limit
						string file_init_STI,
						string file_intervention_wrapper,
						unsigned int nLHS,
						unsigned int nMC,
						unsigned int jobnum
						)
{
	// DEPRECATED
	dcDataFrame df0;
	return df0;
}

void calibration_LHS_new(Population P,
						 string limit_LHS_file_wrapper, //file name (wrapper) containing file names of LHS limit
						 string target_file_wrapper,	//file name (wrapper) containing file names of calibration targets
						 unsigned int nLHS,
						 unsigned int nMC,
						 unsigned int jobnum
						 )
{
	
	//	DEPRECATED
}
