//
//  network_analytics.cpp
//  STIagent_AIR
//
//  Created by David CHAMPREDON on 12/3/2013.
//  Copyright (c) 2013 David CHAMPREDON. All rights reserved.
//

#include "network_analytics.h"



vector<double> degreeDistribution(Population &P)
{

	/// Calculate the degree distribution of partnerships
	
	int N = P.get_size();
	unsigned int pmax = P.getMaxDegree();
	vector<unsigned long> countDegree(pmax+1,0);
	
	for (int i=0; i<N; i++)
	{
		if (P.getIndividual(i).isAlive())
		{
			int k = P.getIndividual(i).get_nCurrSexPartner();
			countDegree[k]++;
		}
	}
	
	
	// Total number of partnerships
	
	unsigned long totPartn = sumElements(countDegree);
	
	// Distribution in terms of proportions
	
	vector<double> distrib(pmax, 0.0);
	
	if (totPartn>0){
		for (int k=0; k<pmax; k++)
		{
			distrib[k] = (double)(countDegree[k])/totPartn;
		}
	}
	
	
	//	vector<unsigned long> data_breaks;
	//	for (unsigned long i=0;i<10;i++) data_breaks.push_back(i);
	//
	//	vector<double> distrib = distributionNormalized(countDegree, data_breaks);
	//
	return distrib;
}



unsigned long clusterSizes(Population P, unsigned long uid, unsigned long uid_exclude)
{
	// RETURNS THE SIZE OF PARTNERSHIP CLUSTER WHERE uid IS.
	// TO AVOID INFINITE LOOP, NEED TO SPECIFY uid_exclude
	// (FIRST CALL OF FUNCTION, INPUT A DUMMY uid_exclude=0)
	
	Individual I = P.getIndividual(uid);
	
	int np = I.get_nCurrSexPartner();
	
	unsigned long res=0;
	
	if (np==0) return 1;
	
	if (np==1 && uid_exclude>0) return 2;
	
	if (np==2 && uid_exclude>0) return 2;
	
	if (np>0)
	{
		res = np;
		
		vector<unsigned long> puid = I.getPartnerUID();
		
		for (int i=0; i<np; i++)
		{
			if (puid[i]!= uid_exclude)
				res += clusterSizes(P, puid[i], uid);
		}
	}
	
	cout << uid << " res = " << res << endl;
	
	return res;
}







