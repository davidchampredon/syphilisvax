/*
 *  check_tools.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-10-04.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "simulation.h"
#include "calibration.h"
#include "population.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <random>

void	check_partnershipMatrix();
void	check_matrixOperations();

void	check_calib_get_parameters();


void	check_GSL_OOL();
void	check_GSL_multinomial();


void	check_speed_uniform();

void	check_speed_multinomial();
void    check_weibull();
