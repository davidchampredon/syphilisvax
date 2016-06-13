//
//  sensitivity.cpp
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 11/25/2013.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#include "sensitivity.h"









void sensitivity_distance_calib_from_files(string target_filename_wrapper,
										   Population initialPopulation,
										   double horizon, double timeStep,
										   int nMonteCarlo, double relativeBump,
										   int DHSphase)
{
	// ======================================
	// CALCULATE SENSITIVITY
	// OF THE DISTANCE TO CALIBRATION TARGETS
	// WITH RESPECT TO ALL CURRENT PARAMETERS
	// ======================================
	
	
	// Set up the simulation
	Simulation S(horizon, timeStep,
				 initialPopulation,
				 nMonteCarlo);
	
	
	// Retrieve from files and set targets
	calibration_set_all_target_wrap(S, target_filename_wrapper, DHSphase);
	
	
	// Set parameter values from files
	S.get_population().setAllParameters("../inputs/");
	
	
	// Store the baseline values of parameters in vectors
	vector<double> prm_DMG;
	vector<double> prm_FORM;
	vector<double> prm_SPOUSAL;
	vector<double> prm_DISSOL;
	// File names are hard coded because of consistency with 'setAllPameters()'
	vectorFromCSVfile(prm_DMG, "in_paramDMG.csv",2);
	vectorFromCSVfile(prm_FORM, "in_paramFORM.csv",2);
	vectorFromCSVfile(prm_SPOUSAL, "in_paramSPOUSAL.csv",2);
	vectorFromCSVfile(prm_DISSOL, "in_paramDISSOL.csv",2);
	
	vector<double> res = S.sensitivity_calib_distance(prm_DMG,
													  prm_FORM,
													  prm_SPOUSAL,
													  prm_DISSOL,
													  relativeBump,
													  nMonteCarlo);
	
	vectorToFile(res, _DIR_OUT+"sensi_calib.out");
	
	
}

