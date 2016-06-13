//
//  code_check.h
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-06-17.
//  Copyright (c) 2014 David CHAMPREDON. All rights reserved.
//

#ifndef __LocalSTI__code_check__
#define __LocalSTI__code_check__

#include "simulation.h"
#include "calibration.h"
#include "population.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <random>

#endif /* defined(__LocalSTI__code_check__) */




void CODECHECK_population_growth();

void CODECHECK_partnerships(Population &P);

void CODECHECK_sexActivity(Population &P);

void CODECHECK_STI_infectivity_curves();


void CODECHECK_mandatory();