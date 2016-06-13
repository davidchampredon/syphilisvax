/*
 *  check_tools.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-10-04.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "check_tools.h"

void check_partnershipMatrix()
{
	Population P(0);
	
	P.initFromFile("startPopulation_small.csv",
				   "in_STI.csv",
				   "in_STI_SFincrease.csv",
				   "in_HIVrebound.csv");
	
	P.setAllParameters("../inputs");
	
	P.displayInfo(false);
	
	double timeStep = 1.00;
	double horizon = 5.0;

	bool save_trace_file = false;
	
	for (double t=0; t<=horizon; t+=timeStep)
	{
		if (P.totalNumberPartnerships()>0)
			P.dissolvePartnerships(timeStep,save_trace_file);
		
		P.formPartnerships(timeStep,save_trace_file);
		P.spousalScanAllCasual();
		P.updateAllDurations(timeStep);
		P.updateAllAges(timeStep);
	}
	
	P.displayInfo(true);
	
	int xx = 1;
	
	// Kill one partner
	dcMatrix puid = P.getPartnershipsUID();
	
	dcMatrix PUID = P.get_partnershipsMatrix();
	PUID.display();
	
	cout << "PARTNERSHIP position: " << P.findPositionIn_partnershipsMatrix(PUID(xx,0),PUID(xx,1))<<endl;

	P.kill(PUID(xx,0),save_trace_file);
	cout << "UID " << PUID(xx,0) << " is killed"<<endl;
	
	cout << " the line with UID "<< PUID(xx,0) << " should be removed:"<<endl;
	PUID = P.get_partnershipsMatrix();
	PUID.display();
	
	cout << "Spouses: "<<endl;
	dcMatrix SP = P.getSpousesUID();
	SP.display();
	
	P.displayInfo(true);
	P.saveToCSVFile(_DIR_OUT+ "pop_check.out");
	
	/*P.writeToFile_PartnerMatrix("P1.csv");
	
	P.writeToFile_all_UID_femaleUID("femUID.csv");
	P.writeToFile_all_UID_maleUID("malUID.csv");
	
	P.writeToFile_PartnerMatrix("P2.csv");
	
	P.writeToFile_PartnerMatrix("P3.csv");
	*/


	
}


void check_matrixOperations()
{
	dcMatrix A(5);
	A.RandomInit();
	A(0,0)=99;
	A.display();
	
	dcMatrix B(5); B=A;
	
	vector<double> xx(5,8);
	A.addRowVector(xx);
	A.display();
	
	
	
	cout <<"--- B --- "<<endl;
	
	B.display();
	
	B.resize(B.getNbRows()+1, B.getNbCols());
	B.display();
	B.setRowValues(5, xx);
	//B(5,4)=99;
	B.display();
	
	cout <<"--- C --- "<<endl;
	dcMatrix C(0,0);
	C.display();
	C.addRowVector(xx);
	C.display();		
	
	
	cout <<"--- Remove row --- "<<endl;
	
	B.display();
	B.removeRow(0);
	B.display();
	// ----- //
	
}




void	check_GSL_multinomial()
{
	const gsl_rng_type * T;
	gsl_rng * r;
	
	/* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */
	
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	
	/* print n random variates chosen from 
     the MULTINOMIAL distribution */

	// Set the parameters
	unsigned int N = 5;
	unsigned int K = 3;
	double *p = new double[K];
	unsigned int *res = new unsigned int[K];
	
	
	p[0] = 0.05;
	for (int k=1; k<K; k++) {
		p[k] = (1.0-p[0])/(K-1);
	}
	
	/*p[0] = 0.1;
	p[1] = 0.6;
	p[2] = 0.3;
	*/
	// Display one draw
	gsl_ran_multinomial(r,K,N,p,res);
	
	for (int i=0;i<K;i++)
	{
		cout << res[i]<<" ; " ;
	}
	cout<<endl;
	
	
	// Run multiple draws with GSL
	int H = 100000;
	
	timeval tim;
    gettimeofday(&tim, NULL);
    double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	for (int c=0; c<H; c++) {
		gsl_ran_multinomial(r,K,N,p,res);
		//cout << res[0]<<"|";
	}
	gettimeofday(&tim, NULL);
    double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	double time_GSL = (t2-t1);
	
	// Run multiple draws with RV.h
	
	vector<double> proba(K);
	for (int k=0; k<K; k++) {
		proba[k] = p[k];
	}
	
	gettimeofday(&tim, NULL);
    t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	for (int c=0; c<H; c++) {
		multinomial(N,proba);
		//cout << res[0]<<"|";
	}
	gettimeofday(&tim, NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	double time_RVh = (t2-t1);
	
	cout << endl<< "RATIO RV/GSL = " << time_RVh/time_GSL << endl;
	
	gsl_rng_free (r);
}



void check_speed_multinomial()
{
	int dim = 5;
	unsigned int N = 3;
	vector<double> proba(dim);
	
	proba[0] = 0.05;
	for (int k=1; k<dim; k++) {
		proba[k] = (1.0-proba[0])/(dim-1);
	}
	
	int H = 1e5; // number of iterations for speed test
    double t1,t2;
	timeval tim;
	
	// Run multiple draws with Direct method
	
    gettimeofday(&tim, NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	vector<unsigned int> resDirect(dim);
	for (int c=0; c<H; c++) 
	{
		resDirect = multinomial(N, proba);
	}
	gettimeofday(&tim, NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	double time_direct = (t2-t1);
	
	
	// Run with GSL (conditional binomial method)
	
	gettimeofday(&tim, NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	gsl_rng * r = GSL_generator(1234);
	
	vector<unsigned int> resGSL(dim);
	for (int c=0; c<H; c++) 
	{
		resGSL = multinomial_gsl(r, N, proba);
		//displayVector(resGSL);
	}
	gettimeofday(&tim, NULL);
	t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	double time_gsl = (t2-t1);	
	
	cout << endl<< "RATIO DIRECT/GSL = " << time_direct/time_gsl << endl;
	
}



void check_calib_get_parameters()
{
	vector<vector<double> > prm = calib_get_parameters("calib_param_DMG.csv",
													   "calib_param_FORM.csv",
													   "calib_param_SPOUSAL.csv",
													   "calib_param_DISSOL.csv");
	
	for (int i=0; i<prm.size(); i++) 
	{
		displayVector(prm[i]);
	}
}


void check_weibull()
{
    double shape = getParameterFromFile("death_Weibull_shape", "in_paramDMG.csv");
    double scale = getParameterFromFile("death_Weibull_scale", "in_paramDMG.csv");
    
    cout << "Weibull shape = " << shape <<endl;
    cout << "Weibull scale = " << scale <<endl;
    
    int n=100;
    double amin = 1.0;
    double da = 1.0;
    
    vector<double> age(n);
    vector<double> deathHazard(n);
    
    for (int i=0; i<n; i++)
    {
        age[i] = amin + i*da;
        deathHazard[i] = weibull_hazard(age[i],shape,scale);
        cout<<"deathHazard("<<age[i]<<") = "<<deathHazard[i]<<endl;
    }
}



void	check_speed_uniform()
{
	// COMPARE SPEED OF UNIFORM RANDOM DIST
	// BETWEEN rand() AND FUNCTIONS
	// FROM C++11 LIBRARY <random>
	
	
	int seedrand = 2;
	srand(seedrand);
	mt19937 eng(seedrand);
	
	unsigned int N=1e7;
	double x;
	
	timeval tim;
	gettimeofday(&tim, NULL);
	
	double z1=tim.tv_sec+(tim.tv_usec/1000000.0);
	
	for (int i=0; i<N; i++)
	{
		x=uniform(rand());
	}
	
	gettimeofday(&tim, NULL);
	double z2=tim.tv_sec+(tim.tv_usec/1000000.0) -z1;
	
	z1=tim.tv_sec+(tim.tv_usec/1000000.0);
	uniform_real_distribution<> dist(0,1);
	for (int i=0; i<N; i++)
	{
		x=dist(eng);
	}
	gettimeofday(&tim, NULL);
	double z3=tim.tv_sec+(tim.tv_usec/1000000.0)-z1;
	
	cout << "Uniform distribution generation:"<<endl;
	cout<<"Time ratio rand()/<random> = " << z2/z3 <<endl;
	
	x++;
	

}


/* ~~~ OOL DESACTIVATION ~~~ (uncomment to use)
 void check_GSL_OOL()
 {
	
	// Initialize empty Population object
	Population P(0);
	
	// Read starting population from a file
	// output of "generateIndividuals.R"
	P.initFromFile_withSTIprev("startPopulation_small.csv", "in_STI.csv","in_STI_SFincrease.csv");
	
	P.setAllParameters();
	
	P.displayInfo(false);
	
	// === Define simulation ===
	
	double horizon = 1.0;
	double timeStep = 10.0/365.0;
	int nMC_calibration = 1;
	Simulation S(horizon, timeStep, P, nMC_calibration);
 
	int dim = 2;
	vector<double> v0(dim);
	v0[0] = 0.60;
	v0[1] = 23.0;
	
	vector<double> LowerBound(dim);
	LowerBound[0] = 0.00;
	LowerBound[1] = 0.00;
	
	vector<double> UpperBound(dim);
	UpperBound[0] = 0.99;
	UpperBound[1] = 99.99;
	
	int maxIterations = 300;
	
	S.calib_minimization_example(v0,LowerBound,UpperBound,maxIterations);
	
 }
 */ // END OOL