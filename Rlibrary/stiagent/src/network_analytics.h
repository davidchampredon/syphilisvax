//
//  network_analytics.h
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 12/3/2013.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#ifndef __STIagent__network_analytics__
#define __STIagent__network_analytics__

#include "dcTools.h"
#include "population.h"


vector<double>		degreeDistribution(Population &P);	// Calculate the degree distribution of partnerships

unsigned long		clusterSizes(Population P, unsigned long uid, unsigned long uid_exclude);	// Cluster (of partnerships) size of Individual #uid

#endif /* defined(__STIagent_AIR__network_analytics__) */



