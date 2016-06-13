/*
 *  calibration.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-10-06.
 *  Copyright 2013. All rights reserved.
 *
 */

#include "dcTools.h"
#include "simulation.h"
#include "globalVar.h"
#include "MCsimulation.h"




// === NEW STUFF ===


dcDataFrame LHS_explore(Population P,
						string limit_LHS_file_wrapper, //file name (wrapper) containing file names of LHS limit
						string file_init_STI,
						string file_intervention_wrapper,
						unsigned int nLHS,
						unsigned int nMC,
						unsigned int jobnum
						);


// === TARGET VALUES ===


dcMatrix calib_getTarget_AgeDistribution(string filename);
dcMatrix calib_getTarget_Distribution(string filename);


vector<string> calib_getTarget_STIprevalences_names(string filename);
vector<double> calib_getTarget_STIprevalences_values(string filename);


Population calib_AgeDistribution(string target_filename, 
								 Population initialPopulation,
								 double horizon, double timeStep,
								 int nMonteCarlo, int nLHS,
								 vector<double> param_initial,
								 vector<double> param_LowerLim,
								 vector<double> param_UpperLim,
								 bool doSex);

vector<vector<double> >	calib_get_parameters(string file_DMG,
											 string file_FORM,
											 string file_SPOUSAL,
											 string file_DISSOL);

vector<vector<double> >	calib_get_parameters_STI(string file_SEXACT,
												 string file_STIFTR_V,
												 string file_STIFTR_B,
												 string file_STIFTR_P,
												 string file_STIFTR_F);



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
											vector<double> param_STIFTR_F_UpperLim);


// Wrapper function to avoid call with numerous parameters
void calibration_set_all_target_wrap(Simulation& S,
									 string target_filename_wrapper,
									 int DHSphase);

// This function is usually called by the wrapper function above
// for convenience
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
								string target_filename_STIprevalences_matrix,
								string target_filename_HIV_prev_age_f,
								string target_filename_HIV_prev_age_m);



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
							  int nMonteCarlo, int DHSphase);


void calibration_LHS_fromFile2(string target_filename_wrapper,
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
							  int nMonteCarlo);


void calibration_LHS_fromFile___OLD(string target_filename_ageDistrib,
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
							  string target_filename_STIprevalences_matrix,
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
							  int nMonteCarlo);


void calibration_LHS_fromFile_wrap(string target_filename_wrapper,
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
								   int nMonteCarlo);


// ===========================
// NEW STUFF


dcDataFrame LHS_generated(vector<string> varname,
						  vector<double> Vmin,
						  vector<double> Vmax,
						  int samplingNb);

void LHS_generate_all_samples(string filename,int samplingNumber,string filename_save);
dcDataFrame LHS_generate_all_samples2(string filename,int samplingNumber,string filename_save);


void calibration_LHS_new(Population P,
						 string limit_LHS_file_wrapper, //file name (wrapper) containing file names of LHS limit
						 string target_file_wrapper,		//
						 unsigned int nLHS,
						 unsigned int nMC,
						 unsigned int nJobs);


void split_LHS_limitsFile(string filename_limits, int n);

