//
//  sensitivity.h
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 11/25/2013.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#ifndef __STIagent_AIR__sensitivity__
#define __STIagent_AIR__sensitivity__

#include "simulation.h"
#include "calibration.h"


void sensitivity_distance_calib_from_files(string target_filename_wrapper,
										   Population initialPopulation,
										   double horizon, double timeStep,
										   int nMonteCarlo, double relativeBump,
										   int DHSphase);



#endif /* defined(__STIagent_AIR__sensitivity__) */
