/*
 *  simulation.cpp
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-08-30.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "simulation.h"
#include "network_analytics.h"
#include "globalVar.h"
#include "calibration.h"


Simulation::Simulation(double horizon,
					   double timeStep,
					   Population P,
					   int nMC_calibration)
{
	_horizon = horizon;
	_timeStep = timeStep;
	
	_schedule.clear();
	for (double t=0; t<_horizon; t+=_timeStep)
		_schedule.push_back(t);
	
	_population = P;
	_population_initial = P;
	
	_population_initial.STI_update_templateToAllIndividuals();
	_population.STI_update_templateToAllIndividuals();
	
	_nMC_calibration = nMC_calibration;
	
	_newborn_timestep = 0;
	
	_STI_incidence.resize(1,_population.get_STI().size());
	_STI_prevalence.resize(1,_population.get_STI().size());
	
	_save_trace_files = false;
}


// ===========================================================
// ===========================================================



vector<double> Simulation::get_STI_prevalence_final(){
	/// Return final prevalence for all STIs
	
	dcMatrix M = _STI_prevalence;
	unsigned int n = M.getNbRows();
	vector<double> x = M.extractRow(n-1);
	
	return x;
}

vector<double> Simulation::get_STI_cumIncidence_final(){
	// DEPRECATED
	return x;
}

void Simulation::STI_set_initial_prevalence(string filename){
	_population.STI_set_initial_prevalence(filename);
}

void Simulation::set_pregnant(unsigned long uid){
	_population.set_pregnant(uid);
}

void Simulation::increment_gestationDuration(unsigned long uid, double x){
	_population.increment_gestationDuration(uid, x);
}

void Simulation::set_gestationDuration(unsigned long uid, double x){
	_population.set_gestationDuration(uid, x);
}

void Simulation::increment_nChildBorn(unsigned long uid){
	_population.increment_nChildBorn(uid);
}


vector<bool> Simulation::MTCT(unsigned long uid)
{
	/// Returns a vector of STI transmission to child
	/// (vector size = number of STIs modelled)
	/// ('true' = transmission to child)
	
	Individual I = _population.getIndividual(uid);
	
	// MTCT at previous time step
	vector<bool> prev_mtct = I.get_STI_MTCT();
	
	// Loop through all STIs
	unsigned long nsti = _population.get_STI().size();
	vector<bool> mtct = prev_mtct;
	
	// DEBUG
	//displayVector(mtct);
	
	for (int i=0; i<nsti; i++)
	{
		STIname stiname = _population.get_STI()[i].get_name();
		
		if(!prev_mtct[i]){
			double stiduration = I.get_STIduration()[i];
			if(stiduration>0){
				// retrieve proba MTCT
				double proba = _population.STI_proba_MTCT(stiname, stiduration);
				// Draw if transmission occurs
				mtct[i] = (uniform01()<proba? true:false);
			}
		}
	}
	return mtct;
}



void Simulation::update_pregnancies(double timestep)
{
	/// Calculate the number of new pregnancies
	/// during a simulation time step.
	/// Loop through all females who had at least
	/// sex act and draw the number of newly
	/// pregnant women from a binomial distribution
	/// with probability _probaPregnantPerSexAct
	///
	/// Also increase the gestation duration by 'timestep'
	
	// --- NEW PREGNANCIES
	
	// Retrieve all female with potential to get pregnant
	vector<unsigned long> uid_pot_preg = _population.get_UID_pot_preg(); //_population.pregnantPotentialFemales();
	
	// DEBUG
//	unsigned int aa = uid_pot_preg.size();
//	unsigned int bb = _population.get_UID_pot_preg().size();
//	cout<< _simulationTime << " -- DEBUG  count= "<<aa<<" ; bookkeeping= "<<bb<<" ; diff = "<< bb-aa<< endl;
	
	// =====
	
	if (uid_pot_preg.size()>0)
	{
		// loop through all mothers-to-be
		// calculates the proba to become pregnant
		// based on number of sex acts (types 1 or 2).
		// Draw random variable for actual pregnancy
		
		double ppsex = _population.get_probaPregnantPerSexAct();
		
		for (int i=0; i<uid_pot_preg.size(); i++)
		{
			int nsex = _population.getIndividual(uid_pot_preg[i]).get_nSexActs_period();    //.nSexActType1or2(); <--- TO DO: WARNING condom not taken into account here!!!
			
			double proba_pregnant = 1 - pow(1 - ppsex, nsex);
			
			double u = uniform01();
			if (u<proba_pregnant)
				set_pregnant(uid_pot_preg[i]);
		}
	}
	
	// --- EXISTING PREGNANCIES
	
	vector<unsigned long> uid_preg = _population.census_pregnant_UID();
	
	_newborn_timestep = 0;  // reset at each time step
	
	for (int i=0; i<uid_preg.size(); i++)
	{
		double gest = _population.getIndividual(uid_preg[i]).get_gestationDuration();
		double age_mother = _population.getIndividual(uid_preg[i]).get_age();
		
		// Mother-to-child STI transmission.
		// Updated at each time step but
		// effective only at delivery
		vector<bool> sti_mtct = MTCT(uid_preg[i]);
		_population.set_STI_MTCT(uid_preg[i],sti_mtct);

		// increase pregnancy if birth not soon enough
		if (gest < 0.75-timestep)
			increment_gestationDuration(uid_preg[i],timestep);
		
		// trigger child birth if close to gestation period
		if (gest >= 0.75-timestep){
			Individual child;
			_nursery.add_child(child,
							   uid_preg[i],
							   age_mother,
							   _population.get_simulationTime(),
							   sti_mtct);
			
			// Birth book-keeping
			set_gestationDuration(uid_preg[i],0);
			increment_nChildBorn(uid_preg[i]);
			_population.update_STI_mtct_cumcount(sti_mtct);

			vector<bool> reset_mtct(_population.get_nSTImodelled(),false);
			_population.set_STI_MTCT(uid_preg[i],reset_mtct);
			
			_newborn_timestep++;
		}
		// DEBUG
		//cout<< "UID "<< uid_preg[i] << " : " << sti_mtct[0]<<" ; " << sti_mtct[1]<<endl;
		// ------
	}
}


double Simulation::STI_cumul_incidence(STIname s,int time_i)
{
	/// Calculates cumulative incidence up to i^th time point
	
	int sti_pos = positionSTIinVector(s, _population.get_STI());
	stopif(time_i>=_schedule.size(), "Cumulative incidence cannot be calculated beyond simulation horizon!");
	
	double sum = 0.0;
	for (int i=0; i<=time_i; i++)
		sum += _STI_incidence(i,sti_pos);
	
	return sum;
}


void Simulation::runAllEvents_timeStep(int numTimeStep,
									   bool doSex,
									   ofstream &f, ofstream &ff,
									   bool logIndivInfo)
{
	// DEPRECATED
}


void Simulation::runAllEvents_timeStep_obj(int numTimeStep,
										   bool doSex,
										   bool logIndivInfo)
{
	/// Execute all simulation events during a single time step
	
	bool debugflag = false;
	double t = _population.get_simulationTime();
	
	if (debugflag) cout << "SIZE = "<<_population.get_size()<<endl;
	
	// Youth arrivals
	if (debugflag) cout << "youth arrivals"<<endl;
	_population.youthArrivals(_timeStep,_save_trace_files);
	
	// Commercial sex workers (CSW)
	if (debugflag) cout << "CSW recruitment"<<endl;
	_population.CSWrecruitment(_timeStep);
	if (debugflag) cout << "CSW cessation"<<endl;
	_population.CSWcessation(_timeStep);
	
	// Resets variables at the start of this period
	if (doSex) _population.reset_sexActs();
	
	// Partnerships formation and dissolution
	if (debugflag) cout << "dissolve partnerships"<<endl;
	if (_population.totalNumberPartnerships()>0)
		_population.dissolvePartnerships(_timeStep,_save_trace_files);
	
	// Partnerships formations and dissolutions
	if (debugflag) cout << "spousalScanAllCasual"<<endl;
	_population.formPartnerships(_timeStep,_save_trace_files);
	_population.spousalScanAllCasual();
	
	// DEBUG
	_population.debug_check_partnerUID();
	
	// Sex acts
	if (debugflag) cout << "update_sexActs"<<endl;
	if (doSex) _population.update_sexActs(_timeStep, _save_trace_files);
	
	// STI transmissions
	if (debugflag) cout << "STI_transmissions"<<endl;
	vector<unsigned long> incidence;
	if (t>0 && doSex){
		incidence = _population.STI_transmissions(_timeStep,_save_trace_files);
	}
	
	// Pregnancies and MTCT
	if (debugflag) cout << "Pregnancies"<<endl;
	if (t>0 && doSex) update_pregnancies(_timeStep);
	
	// Duration, age updates
	if (debugflag) cout << "updateAllDurations_updateAllAges"<<endl;
	_population.updateAllDurations(_timeStep);
	_population.updateAllAges(_timeStep);
	
	// Natural clearance of all STIs
	if (debugflag) cout << "STI_update_naturalClearance"<<endl;
	_population.STI_update_naturalClearance();
	
	// STI Treatment and vaccine waning:
	for (int sti=0; sti<_population.get_STI().size(); sti++){
		STIname stiname = _population.get_STI()[sti].get_name();
		update_cure(stiname);
		update_vacc(stiname);
	}

	// Deaths
	if (debugflag) cout << "deathEvents"<<endl;
	_population.deathEvents(_timeStep,false /*_save_trace_files*/); //force not write trace file because huge!
	
	
	// update incidence
	if (doSex){
		if (numTimeStep==0){
			unsigned long n_sti = _population.get_STI().size();
			for (int j=0;j<n_sti; j++)
				_STI_incidence(0,j) = 0.0;
		}
		if (numTimeStep>0){
			_STI_incidence.addRowVector(incidence);
		}
	}
	
	// update reproductive numbers for all STIs
	if (doSex) {
		for (int sti=0; sti<_population.get_nSTImodelled(); sti++) {
			STIname stiname = _population.get_STI()[sti].get_name();
			_population.update_secondary_cases(stiname);
		}
	}
	
	vector<double> v;
	
	if(doSex)	// TO DO: do not impose no log file when no sex
	{
		//
		// * * * * WARNING * * * *
		//
		// ORDER MUST BE SAME AS DEFINED
		// BY HEADERS IN "runAllEvents_horizon_obj"
		
		v.push_back(t);
		v.push_back(_population.census_alive());
		v.push_back(_population.census_dead());
		v.push_back(_population.totalNumberPartnerships());
		v.push_back(_population.totalNumberSpousalPartneships());
		v.push_back(_population.census_Females());
		
		int nSTI = _population.get_nSTImodelled();
		for (int i=0; i<nSTI; i++)
			v.push_back(_population.census_STIinfected()[i]);
		
		v.push_back(_population.census_STIcoinfected(HIV, Tp));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,0));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,1));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,2));
		v.push_back(_population.census_STIcoinfected(HIV, Tp,9));

		v.push_back(_population.census_CSW().size());
		
		for (int r=0; r<= _population.get_maxRiskGroup(); r++)
			v.push_back(_population.census_riskGroup()[r]);
		
		v.push_back(_population.census_circum());
		v.push_back(_newborn_timestep);
		
		v.push_back(_population.census_pregnant(0));
		v.push_back(_population.census_pregnant(1));
		v.push_back(_population.census_pregnant(2));
		v.push_back(_population.census_pregnant(9));
		
		v.push_back(_population.get_STI_mtct_cumcount()[0]);
		v.push_back(_population.get_STI_mtct_cumcount()[1]);
		
		// HIV & Tp prevalence:
		v.push_back(_population.STI_prevalence(HIV));
		v.push_back(_population.STI_prevalence(Tp));
		
		// HIV prevalence by risk group
		v.push_back(_population.STI_prevalence(HIV, 0));
		v.push_back(_population.STI_prevalence(HIV, 1));
		v.push_back(_population.STI_prevalence(HIV, 2));
		v.push_back(_population.STI_prevalence(HIV, 9));
		
		// Tp prevalence by risk group
		v.push_back(_population.STI_prevalence(Tp, 0));
		v.push_back(_population.STI_prevalence(Tp, 1));
		v.push_back(_population.STI_prevalence(Tp, 2));
		v.push_back(_population.STI_prevalence(Tp, 9));
		
		// Reproductive numbers
		v.push_back(_population.Reff_cum_mean(HIV));
		v.push_back(_population.Reff_cum_mean(Tp));
		
		// Count of sex acts
		v.push_back(_population.count_sexActs(0));
		v.push_back(_population.count_sexActs(1));
		v.push_back(_population.count_sexActs(2));
		v.push_back(_population.count_sexActs(9));
		
		// * * * * WARNING * * * *
		//
		// ORDER MUST BE SAME AS DEFINED
		// BY HEADERS IN "runAllEvents_horizon_obj"
		// * * * * * * * * * * * *
		
		// Update data frame:
		_df_sim.addrow(to_string(numTimeStep), v);
	}
}


void Simulation::runAllEvents_horizon(bool doSex,
									  bool logIndivInfo,
									  bool traceNetwork,
									  int displayProgress,
									  int iter_mc)
{
	/// Run a simulation with all events until specified horizon
	
	
	// Trace files -----
	string filename1		= _DIR_OUT + "simul_mc" + int2string(iter_mc) + ".out";
	string filename2		= _DIR_OUT + "simul_indiv_mc" + int2string(iter_mc) +".out";
	string fileDegreeDist	= _DIR_OUT + "degreeDist_mc" + int2string(iter_mc) +".out";
	string filepop			= _DIR_OUT + "pop_mc"+ int2string(iter_mc);
	
	ofstream f(filename1.c_str(), ios::app);
	ofstream ff(filename2.c_str(), ios::app);
	ofstream fd(fileDegreeDist.c_str(),ios::app);
	// -----------------
	
	
	int nSTI = _population.get_nSTImodelled();
	_population.set_timeStep(_timeStep);
	_nursery.set_STI(_population.get_STI());
	
	// Number of interventions
	unsigned long n_intervention = _intervention.size();
	
	// DEBUG
	if ( _MC_trial_iter==1){
		cout << " -- Number of interventions:"<<n_intervention<<endl;
		for (int i=0; i<n_intervention; i++) {
			_intervention[i].displayInfo();
		}
	}
	
	
	bool tmp_save_trace_files = true; // <-- for debug purposes
	if (doSex && tmp_save_trace_files/*_save_trace_files*/)
	{
		// Headers of the log files
		// WARNING: headers must be consistent with values writtten in "runAllEvents_timeStep"
		f<<"time, nAlive, nDead, nPartn, nSp,nFemale,";
		for (int i=0; i<nSTI; i++)
		{
			f << STInameString(_population.get_STI()[i].get_name());
			//if (i<nSTI-1)
			f<<",";
		}
		f << "nHIVTp,";
		f << "nHIVTp0,";
		f << "nHIVTp1,";
		f << "nHIVTp2,";
		f << "nHIVTp9,";
		f << "nCSW,";
		
		for (int r=0; r<= _population.get_maxRiskGroup(); r++)
		{
			f << "nRskGrp"+int2string(r);
			f << ",";
		}
		
		f << "nCircum,nNewBorn";
		f << ",";
		
		f << "HIVprevRisk0,HIVprevRisk1,HIVprevRisk2,HIVprevRisk9,";
		f << "TpprevRisk0,TpprevRisk1,TpprevRisk2,TpprevRisk9";
		f << ",";
		
		f << "Reff_HIV, Reff_Tp";
		
		f <<endl;
	}
	
	if (logIndivInfo)
		ff << "time,UID,alive,UIDpartner0, nSexP0, HIVdur,HIVinf"<<endl;
	// ----
	
	
	
	// Loop through time until horizon
	
	int cnt = 0;
	int t_i = 0;
	
	int n_timesteps = (int)(_horizon/_timeStep);
	
	
	for (double t=0; t<_horizon; t+=_timeStep){
		// write network info
		if (traceNetwork){
			string file_cnt = filepop+ int2string(cnt) + ".out";
			_population.FileOutput(file_cnt);
			cnt++;
		}
		
		double prevSize = _population.get_size();
		
		
		// ---------------------------
		// Display simulation progress
		
		if (displayProgress==1){
			cout << " time: "<< t << " size: "<< prevSize;
			cout << " alive: " << _population.census_alive() <<endl;
		}
		
		if (displayProgress==11){
			if (t_i%(n_timesteps/10)==0)
			{
				cout << " time: "<< t << " size: "<< prevSize;
				cout << " alive: " << _population.census_alive();// <<endl;

				cout << " ; prev_HIV: " << _population.STI_prevalence(HIV);
				cout << " ; prev_Tp: " << _population.STI_prevalence(Tp) <<endl;
//				cout << " ; Reff_HIV: " << _population.Reff_cum_mean(HIV);
//				cout << " Reff_Tp: " << _population.Reff_cum_mean(Tp) << endl;
			}
		}
		
		if(displayProgress==2){
			cout << "alive: " <<	_population.census_alive() << endl;
			cout << "partner0: " << _population.census_Partnered(0) << endl;
			cout << "partner1: " << _population.census_Partnered(1) << endl;
			cout << "partner2: " << _population.census_Partnered(2) << endl;
			cout << "singleRatio="<<_population.census_ratioSingles()<<endl;
		}
		// ---------------------------
		
		
		// Records the time of simulation
		_simulationTime = t;
		_population.set_simulationTime(t);
		
		
		// ================================================
		// === Execute all events during the time step  ===
		// ================================================
		
		runAllEvents_timeStep(t_i,
							  doSex,
							  f,
							  ff,
							  logIndivInfo);
		
		// Record STI prevalences
		
		if (doSex)
		{
			vector<double> tmp;
			for (int s=0; s<_population.get_STI().size(); s++)
			{
				STIname stiname = _population.get_STI()[s].get_name();
				
				if (t==0) _STI_prevalence(t_i,s) = _population.STI_prevalence(stiname);
				if (t>0)
				{
					tmp.push_back(_population.STI_prevalence(stiname));
				}
			}
			
			if (t>0) _STI_prevalence.addRowVector(tmp);
			
			// DEBUG
			//_STI_prevalence.display();
		}
		
		
		// Records degree distribution of the network at each time step
		//
		// Headers for degree file
		// (Save in file just up to 3 concurrent partners)
		if (t<1E-5 && _save_trace_files) fd<<"time,P0,P1,P2,P3"<<endl;
		if(doSex && _save_trace_files)
		{
			vector<double> degDist = degreeDistribution(_population);
			vector<double> tmp;
			
			tmp.push_back(t);
			
			for (int i=0;i<=3;i++) tmp.push_back(degDist[i]);
			
			vectorToCSVFile_Row(tmp, fileDegreeDist);
		}
		
		
		// ========================
		// ====  CALIBRATION  =====
		// ========================
		
		if( doSex && isCalibrationTime(t) )
		{
			unsigned int calib_idx = whichCalibrationTime(t);
			
			calculate_calibDistance_all(calib_idx);
			
			// DEBUG
			//			cout<<" INTERIM CALIBRATION DISTANCE @ "<<t<<endl;
			//			for (int i=0;i<_calibrationDistances.size();i++) displayVector(_calibrationDistances[i]);
			
			// Dump outputs to be compared to calibration targets (outside C++)
			calibration_target_type_output_save(iter_mc,calib_idx);
		}		
		
		
		// ========================
		// ==== INTERVENTIONS =====
		// ========================
		
		
		// Intervention relevant when sex act occurs
		if (doSex && n_intervention>0)
		{
			// loop through all interventions
			for (int i=0; i<n_intervention; i++)
			{
				vector<double> sched = _intervention[i].get_schedule();
				
				// loop through the intervention schedule start/end dates
				for (int k=0; k<sched.size()/2; k++)
				{
					if (sched[2*k]<t && t<sched[2*k+1])
					{
						activate_intervention(i);
						//						//DEBUG
						//						cout << "time: "<<t<<" intervention #"<<i<<" activated (MC iter:"<<_MC_trial_iter<<")"<<endl;
					}
				}
			}
			
		}
		
		// Increase time step counter
		t_i++;
		
	} // end for loop on time
	
	if(_save_trace_files) _nursery.saveToCSVfile(_DIR_OUT + "nursery_"+int2string(_MC_trial_iter)+ ".out");
}


void Simulation::runAllEvents_horizon_obj(bool doSex,
										  bool logIndivInfo,
										  bool traceNetwork,
										  int displayProgress)
{
	/// Run a simulation with all events until specified horizon
	
	unsigned long nSTI = _population.get_nSTImodelled();
	_population.set_timeStep(_timeStep);
	_nursery.set_STI(_population.get_STI());
	
	// Number of interventions
	unsigned long n_intervention = _intervention.size();
	
//	// DEBUG
//	if ( _MC_trial_iter==1){
//		cout << " -- Number of interventions:"<<n_intervention<<endl;
//		for (int i=0; i<n_intervention; i++) {
//			_intervention[i].displayInfo();
//		}
//	}
	
	
	// Set-up data frame that will hold simulation outputs
	vector<string> colnames;
	
	if (doSex)
	{
		// Headers of the data frame
		//
		//  * * * WARNING * * * *
		// headers must be consistent with
		// values writtten in "runAllEvents_timeStep"
		//
		
		colnames.push_back("time");
		colnames.push_back("nAlive");
		colnames.push_back("nDead");
		colnames.push_back("nPartn");
		colnames.push_back("nSp");
		colnames.push_back("nFemale");

		for (int i=0; i<nSTI; i++)
			colnames.push_back(STInameString(_population.get_STI()[i].get_name()));
		
		colnames.push_back("nHIVTp");
		colnames.push_back("nHIVTp0");
		colnames.push_back("nHIVTp1");
		colnames.push_back("nHIVTp2");
		colnames.push_back("nHIVTp9");
		
		colnames.push_back("nCSW");
		
		for (int r=0; r<= _population.get_maxRiskGroup(); r++)
			colnames.push_back("nRskGrp"+int2string(r));
		
		colnames.push_back("nCircum");
		colnames.push_back("nNewBorn");
		
		colnames.push_back("nPregnantRisk0");
		colnames.push_back("nPregnantRisk1");
		colnames.push_back("nPregnantRisk2");
		colnames.push_back("nPregnantRisk9");
		
		colnames.push_back("mtctHIV");
		colnames.push_back("mtctTp");

		colnames.push_back("HIVprev");
		colnames.push_back("Tpprev");
		colnames.push_back("HIVprevRisk0");
		colnames.push_back("HIVprevRisk1");
		colnames.push_back("HIVprevRisk2");
		colnames.push_back("HIVprevRisk9");
		colnames.push_back("TpprevRisk0");
		colnames.push_back("TpprevRisk1");
		colnames.push_back("TpprevRisk2");
		colnames.push_back("TpprevRisk9");
		
		colnames.push_back("Reff_HIV");
		colnames.push_back("Reff_Tp");
		
		colnames.push_back("nSexActRisk0");
		colnames.push_back("nSexActRisk1");
		colnames.push_back("nSexActRisk2");
		colnames.push_back("nSexActRisk9");
		
		_df_sim.set_colname(colnames);
	}

	
	
//	if (logIndivInfo)
//		ff << "time,UID,alive,UIDpartner0, nSexP0, HIVdur,HIVinf"<<endl;
//	// ----
//	
	
	// ===================================================
	// ===    Loop through time until horizon
	// ===================================================
	
//	int cnt = 0;
	int t_i = 0;
	int n_timesteps = (int)(_horizon/_timeStep);
	

	for (double t=0; t<_horizon; t+=_timeStep)
	{
		// write network info
		if (traceNetwork){
//			string file_cnt = filepop+ int2string(cnt) + ".out";
//			_population.FileOutput(file_cnt);
//			cnt++;
		}

		// ---------------------------
		// Display simulation progress
		
		double prevSize = _population.get_size();

		if (displayProgress==1){
			cout << " time: "<< t << " size: "<< prevSize;
			cout << " alive: " << _population.census_alive() <<endl;
		}
		if (displayProgress==11){
			if (t_i%(n_timesteps/10)==0)
			{
				cout << " time: "<< t << " size: "<< prevSize;
				cout << " alive: " << _population.census_alive();// <<endl;
				
				cout << " ; prev_HIV: " << _population.STI_prevalence(HIV);
				cout << " ; prev_Tp: " << _population.STI_prevalence(Tp) <<endl;
				//				cout << " ; Reff_HIV: " << _population.Reff_cum_mean(HIV);
				//				cout << " Reff_Tp: " << _population.Reff_cum_mean(Tp) << endl;
			}
		}
		if(displayProgress==2){
			cout << "alive: " <<	_population.census_alive() << endl;
			cout << "partner0: " << _population.census_Partnered(0) << endl;
			cout << "partner1: " << _population.census_Partnered(1) << endl;
			cout << "partner2: " << _population.census_Partnered(2) << endl;
			cout << "singleRatio="<<_population.census_ratioSingles()<<endl;
		}
		// ---------------------------
		
		// Records the time of simulation
		_simulationTime = t;
		_population.set_simulationTime(t);
		
		
		// ================================================
		// === Execute all events during the time step  ===
		// ================================================
		
		runAllEvents_timeStep_obj(t_i,
								  doSex,
								  logIndivInfo);
		
		// Record STI prevalences
		if (doSex)
		{
			vector<double> tmp;
			for (int s=0; s<_population.get_STI().size(); s++)
			{
				STIname stiname = _population.get_STI()[s].get_name();
				
				if (t==0) _STI_prevalence(t_i,s) = _population.STI_prevalence(stiname);
				if (t>0) tmp.push_back(_population.STI_prevalence(stiname));
			}
			
			if (t>0) _STI_prevalence.addRowVector(tmp);
		}
		
		
		// Records degree distribution of the network at each time step
		//
		// Headers for degree file
		// (Save in file just up to 3 concurrent partners)
//		if (t<1E-5 && _save_trace_files) fd<<"time,P0,P1,P2,P3"<<endl;
		if(doSex && _save_trace_files)
		{
			vector<double> degDist = degreeDistribution(_population);
			vector<double> tmp;
			
			tmp.push_back(t);
			
			for (int i=0;i<=3;i++) tmp.push_back(degDist[i]);
		}
		
		
		// ========================
		// ====  CALIBRATION  =====
		// ========================
		
		if( doSex && isCalibrationTime(t) )
		{
			unsigned int calib_idx = whichCalibrationTime(t);
			
			calculate_calibDistance_all(calib_idx);
			
			// DEBUG
			//			cout<<" INTERIM CALIBRATION DISTANCE @ "<<t<<endl;
			//			for (int i=0;i<_calibrationDistances.size();i++) displayVector(_calibrationDistances[i]);
			
			// Dump outputs to be compared to calibration targets (outside C++)
//			calibration_target_type_output_save(iter_mc,calib_idx);
		}
		
		
		// ========================
		// ==== INTERVENTIONS =====
		// ========================
		
		
		// Intervention relevant when sex act occurs
		if (doSex && n_intervention>0)
		{
			// loop through all interventions
			for (int i=0; i<n_intervention; i++)
			{
				vector<double> sched = _intervention[i].get_schedule();
				
				// loop through the intervention schedule start/end dates
				for (int k=0; k<sched.size()/2; k++)
				{
					if (sched[2*k]<t && t<sched[2*k+1])
					{
						activate_intervention(i);
						//						//DEBUG
						//						cout << "time: "<<t<<" intervention #"<<i<<" activated (MC iter:"<<_MC_trial_iter<<")"<<endl;
					}
				}
			}
			
		}
		
		// Increase time step counter
		t_i++;
		
	} // end for loop on time
	
	if(_save_trace_files) {}
//		_nursery.saveToCSVfile(_DIR_OUT + "nursery_"+int2string(_MC_trial_iter)+ ".out");
}



/* ***************************************************/
/* ************** C A L I B R A T I O N **************/
/* ***************************************************/


void Simulation::set_calibration_schedule(string filename_calibrationTime,
										  string filename_all_calib_target,
										  string filename_weight_target)
{
	/// Set the calibration schedule
	/// which consists of:
	/// - a vector of times
	/// - for each of those times, a vector of
	/// strings identifying the filename of output
	/// to be calibrated (e.g. age distribution)
	
	/// Example to illustrate:
	///
	/// calibration times:
	/// ____
	/// 1.2
	/// 3.4
	/// 7.8
	///
	/// outputs to be calibrated at those times:
	/// _________________________________________________
	/// ageDist1.csv  |  ageDist2.csv   |  ageDist3.csv
	/// ageGap1.csv   |  ageGap2.csv    |  ageGap3.csv
	///               |  prevTp2.csv    |  prevTp3.csv
	///               |                 |  PrevHIV3.csv
	///
	/// Each file name must have a pre-specified name such that
	/// the _type_ of target is recognized and we can construct
	/// the associated table of calibration types:
	///
	/// _________________________________________________
	/// AGE           |  AGE            |  AGE
	/// AGE           |  GAP            |  GAP
	///               |  TP             |  TP
	///               |                 |  HIV
	///
	/// Finally, each calibration target has its weight
	/// in the global calibration score. For example:
	///
	/// _________________________________________________
	/// 1             |  1              |  1
	/// 1             |  1              |  1
	///               |  1              |  2
	///               |                 |  2
	
	
	// clean up
	_calibrationTime.clear();
	_calibrationTargetFile.clear();
	_calibrationTargetType.clear();
	_calibrationWeights.clear();
	_calibrationDistances.clear();
	
	// read calbration times from the file
	vectorFromCSVfile(_calibrationTime, filename_calibrationTime.c_str(), 1);
	
	// Integrity check
	unsigned int ncaltimes=getNumberColumns(filename_all_calib_target);
	stopif(_calibrationTime.size()!=ncaltimes,
		   "Calibration schedule files not consistent!");
	
	
	// Read all calibration target files
	
	vector< vector<string> > s_filename;
	s_filename.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		vectorFromCSVfile_string(s_filename[j],
								 filename_all_calib_target.c_str(),
								 j+1);
		s_filename[j] = trim(s_filename[j]);
		//displayVector(s_filename[j]);
	}
	_calibrationTargetFile = s_filename;
	
	
	// Convert filenames into calibration types
	
	vector< vector<string> > s_type;
	s_type.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		for(int i=0;i<s_filename[j].size(); i++)
			s_type[j].push_back(calibration_file_to_type(s_filename[j][i]));
		
		//displayVector(s_type[j]);
	}
	_calibrationTargetType = s_type;
	
	
	// Retrieve calibration weights
	
	_calibrationWeights.resize(ncaltimes);
	
	for(int j=0; j<ncaltimes; j++)
	{
		vectorFromCSVfile(_calibrationWeights[j],
						  filename_weight_target.c_str(),
						  j+1);
		
		//displayVector(_calibrationWeights[j]);
	}
	
	// Set correct sizes for calibration distances vectors.
	// Initialized values are set to obviously dummy ones (i.e. negative)
	// to facilitate bug detection
	
	_calibrationDistances.resize(ncaltimes);
	
	for (int i=0; i<ncaltimes; i++) {
		_calibrationDistances[i].resize(_calibrationTargetType[i].size(),-9E9);
	}
	
}

string calibration_file_to_type(string filename)
{
	/// converts a file name into a calibration type
	/// (the file name must obey the rules in this function)
	
	string res="";
	
	if(filename.substr(0,7)	=="ageDist")		res = "AGEDIST";
	if(filename.substr(0,10)=="ageGapDist")		res = "GAPDIST";
	if(filename.substr(0,6)	=="prevTp")			res = "TP";
	
	// USED?:
	if(filename.substr(0,7)	=="prevHIV")		res = "HIV";
	// -----
	
	if(filename.substr(0,12)=="HIVpos_age_f")	res = "HIVDAGEF";
	if(filename.substr(0,12)=="HIVpos_age_m")	res = "HIVDAGEM";
	if(filename.substr(0,14)=="HIV_prev_age_f")	res = "HIVPREVAGEF";
	if(filename.substr(0,14)=="HIV_prev_age_m")	res = "HIVPREVAGEM";
	
	if(filename.substr(0,7)	=="stiPrev")		res = "STIPREV";
	if(filename.substr(0,11)=="stiPrevRisk")	res = "STIPREVRISK";
	
	if(filename.substr(0,11)=="singleRatio")	res = "SINGLERATIO";
	
	if(filename.substr(0,9)=="age1sex_f")		res = "AGE1SEXF";
	if(filename.substr(0,9)=="age1sex_m")		res = "AGE1SEXM";
	
	if(filename.substr(0,7)=="lftNP_f")			res = "LFTNPF";
	if(filename.substr(0,7)=="lftNP_m")			res = "LFTNPM";
	
	if(filename.substr(0,8)=="visitCSW")		res = "VISITCSW";
	if(filename.substr(0,11)=="visitCSWage")	res = "VISITCSWAGE";
	if(filename.substr(0,9)=="condomCSW")		res = "CONDOMCSW";
	
	stopif(res=="","Cannot convert calibration file name ("+ filename+ ") to a calibration type!");
	return res;
}


bool Simulation::isCalibrationTime(double t)
{
	/// Check if this is a calibration time
	
	bool res = false;
	
	for(int i=0;i<_calibrationTime.size();i++)
	{
		if (fabs(t-_calibrationTime[i])<_timeStep/2) res=true;
	}
	return res;
}

unsigned int Simulation::whichCalibrationTime(double t)
{
	/// Given this is a calibration time, which set is it?
	
	unsigned int res = 0;
	
	for(int i=0;i<_calibrationTime.size();i++)
	{
		if (fabs(t-_calibrationTime[i])<_timeStep/2) res=i;
	}
	return res;
}

void Simulation::calculate_calibDistance_all(int calibtime_idx)
{
	/// Calculate distance of model from all targets
	
	for(int i=0; i<_calibrationTargetType[calibtime_idx].size(); i++)
	{
		calculate_calibDistance_one(calibtime_idx,i);
	}
}

void Simulation::calculate_calibDistance_one(int calibtime, int i)
{
	/// Calculate distance of model from ONE specific target
	
	bool type_found=false;
	double d = -999.999;
	
	string calibtype = _calibrationTargetType[calibtime][i];
	string targetfile = _DIR_CALIB + _calibrationTargetFile[calibtime][i];
	
	// Loop through all known calibration types
	
	
	// -- demographics
	
	if(calibtype=="AGEDIST")
	{
		// set age distribution found in file as a new target
		dcMatrix TAD = distributionFromFile(targetfile);
		set_target_ageDistribution(TAD);
		
		d = calib_distance_ageDistrib();
		type_found=true;
	}
	
	
	// -- partnerships
	
	if(calibtype=="GAPDIST")
	{
		// set age gap distribution found in file as a new target
		dcMatrix TAGD = distributionFromFile(targetfile);
		set_target_ageGapDistribution(TAGD);
		
		d = calib_distance_ageGapDistrib();
		type_found=true;
	}
	
	
	if(calibtype=="SINGLERATIO")
	{
		double TSR_f = getParameterFromFile("female", targetfile);
		double TSR_m = getParameterFromFile("male", targetfile);
		set_target_singleRatio_f(TSR_f);
		set_target_singleRatio_m(TSR_m);
		
		d = calib_distance_singleRatio(female)+calib_distance_singleRatio(male);
		
		type_found=true;
	}
	
	
	// -- sexual behaviour
	
	if(calibtype=="AGE1SEXF")
	{
		dcMatrix TAFS_f = calib_getTarget_Distribution(targetfile);
		set_target_ageFirstSexDistribution_f(TAFS_f);
		
		d = calib_distance_ageFirstSexDistrib(female);
		
		type_found=true;
	}
	
	if(calibtype=="AGE1SEXM")
	{
		dcMatrix TAFS_m = calib_getTarget_Distribution(targetfile);
		set_target_ageFirstSexDistribution_m(TAFS_m);
		
		d = calib_distance_ageFirstSexDistrib(male);
		
		type_found=true;
	}
	
	
	if(calibtype=="LFTNPF")
	{
		dcMatrix TNLS_f = calib_getTarget_Distribution(targetfile);
		set_target_nLftSexPrtnrDistribution_f(TNLS_f);
		
		d = calib_distance_nLifeSexPrtnrDistrib(female);
		
		type_found=true;
	}
	
	if(calibtype=="LFTNPM")
	{
		dcMatrix TNLS_m = calib_getTarget_Distribution(targetfile);
		set_target_nLftSexPrtnrDistribution_m(TNLS_m);
		
		d = calib_distance_nLifeSexPrtnrDistrib(male);
		
		type_found=true;
	}
	
	
	// -- epidemiological
	
	if(calibtype=="STIPREV")
	{
		// Retrieve STI names to be calibrated
		vector<string> stinames;
		vectorFromCSVfile_string(stinames, targetfile.c_str(), 1);
		set_target_STIprevalence_names(StringToSTIname_vector(stinames));
		
		// STI Prevalence global (across risk groups, genders,etc.)
		vector<double> stiprev;
		vectorFromCSVfile(stiprev, targetfile.c_str(), 2);
		set_target_STIprevalence(stiprev);
		
		d = calib_distance_STIprev();
		type_found=true;
	}
	
	
	if(calibtype=="STIPREVRISK")
	{
		// Retrieve STI names to be calibrated
		vector<string> stinames;
		vectorFromCSVfile_string(stinames, targetfile.c_str(), 1);
		set_target_STIprevalence_names(StringToSTIname_vector(stinames));
		
		// STI Prevalence by risk group
		int nRiskGroup = _population.get_maxRiskGroup()+2;
		
		dcMatrix STIprevTargetRiskGroup(stinames.size(),nRiskGroup);
		
		MatrixFromCSVfile(STIprevTargetRiskGroup,
						  targetfile.c_str(),
						  nRiskGroup);
		
		set_target_STIprevalence_by_riskGroup(STIprevTargetRiskGroup);
		
		d = calib_distance_STIprev_riskGroup();
		type_found=true;
	}
	
	
	if(calibtype=="HIVPREVAGEF")
	{
		dcMatrix THIVAGE_f = calib_getTarget_AgeDistribution(targetfile);
		set_target_HIV_prev_age_f(THIVAGE_f);
		
		d = calib_distance_HIV_prev_age(female);
		
		type_found=true;
	}
	
	if(calibtype=="HIVPREVAGEM")
	{
		dcMatrix THIVAGE_m = calib_getTarget_AgeDistribution(targetfile);
		set_target_HIV_prev_age_m(THIVAGE_m);
		
		d = calib_distance_HIV_prev_age(male);
		
		type_found=true;
	}
	
	
	
	// if not found, throw error message:
	string errmsg ="Distance from target of calibration" + calibtype+ " type not implemented!";
	stopif(!type_found,errmsg);
	
	// if found, update calibration distances:
	_calibrationDistances[calibtime][i] = d*_calibrationWeights[calibtime][i];
	
}

void Simulation::calibration_target_type_output_save(unsigned long iter_mc,int calibtime_idx)
{
	/// Save all output based on target types and calibration times
	
	
	string extfile = "_"+to_string(iter_mc)+".out";
	string folder = _DIR_OUT;
	
	int ct = calibtime_idx;
	
	bool type_found = false;
	
	// Loop through all calibration types at that date
	
	for(int i=0; i<_calibrationTargetType[ct].size(); i++)
	{
		string calibtype = _calibrationTargetType[ct][i];
		
		// -- demographics
		
		if(calibtype=="AGEDIST")
		{
			save_outputs_demog(folder,iter_mc,ct);
			type_found=true;
		}
		
		// -- partnerships
		
		if(calibtype=="GAPDIST" || calibtype=="SINGLERATIO")
		{
			save_outputs_prtnr(folder,iter_mc,ct);
			type_found=true;
		}
		
		// -- sexual behaviour
		
		if(calibtype=="AGE1SEXF" || calibtype=="AGE1SEXM" ||
		   calibtype=="LFTNPF" || calibtype=="LFTNPM")
		{
			save_outputs_sex(folder,iter_mc,ct);
			type_found=true;
		}
		
		// -- epidemiological
		
		if(calibtype=="STIPREV" ||calibtype=="STIPREVRISK" ||
		   calibtype=="HIVPREVAGEF" || calibtype=="HIVPREVAGEM")
		{
			save_outputs_epi(folder,iter_mc,ct);
			type_found=true;
		}
		
		// if not found, throw error message:
		string errmsg ="Calibration type" + calibtype+ " not implemented!";
		stopif(!type_found,errmsg);
	}
}

double Simulation::calibration_distance_targets()
{
	/// Returns the overall calibration distance from all targets
	/// (must be calculated after the last calibration date)
	
	double s=sumElements(_calibrationDistances);
	
	stopif(s<0,"Calibration distance cannot be negative!");
	
	return s;
}

double Simulation::calib_distance_ageDistrib()
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR AGE DISTRIBUTION
	
	// Extract targets: age
	vector<double> ageBreaks		= _target_ageDistribution.extractColumn(0);
	vector<double> targetProp		= _target_ageDistribution.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	ageBreaks.push_back(999.99);
	
	// Retrieve current age distribution
	vector<double> currentAgeDistrib = _population.census_ageDistribution(ageBreaks);
	
	// Distance from target
	int thepower = 2;
	
	double diff_age	= distanceLvector(thepower, currentAgeDistrib, targetProp);
	
	// Normalize distance
	double w_age	= normLvector(thepower,targetProp);
	
	return diff_age/w_age;
}

double Simulation::calib_distance_ageGapDistrib()
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR AGE GAP DISTRIBUTION
	
	// Extract targets: age gap
	vector<double> ageGapBreaks		= _target_ageGapDistribution.extractColumn(0);
	vector<double> targetPropGap	= _target_ageGapDistribution.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	ageGapBreaks.push_back(999.99);
	
	// Retrieve current age gap distribution
	vector<double> currentAgeGapDistrib	= _population.census_ageGapDistribution(ageGapBreaks);
	
	// Distance from target
	int thepower = 2;
	
	double diff_ageGaps	= distanceLvector(thepower, currentAgeGapDistrib, targetPropGap);
	
	// Normalize distance
	double w_ageG	= normLvector(thepower,targetPropGap);
	
	return diff_ageGaps/w_ageG;
}

double Simulation::calib_distance_ageFirstSexDistrib(Gender g)
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR AGE AT FIRST SEX DISTRIBUTION
	
	dcMatrix M = (g==female)?_target_ageFirstSexDistribution_f:_target_ageFirstSexDistribution_m;
	
	
	// Extract targets: age
	vector<double> ageBreaks		= M.extractColumn(0);
	vector<double> targetProp		= M.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	ageBreaks.push_back(999.99);
	
	// Retrieve current age distribution
	vector<double> currentAgeDistrib = _population.census_ageFirstSexDistribution(ageBreaks, g);
	
	// Distance from target
	int thepower = 2;
	
	double diff_age	= distanceLvector(thepower, currentAgeDistrib, targetProp);
	
	// Normalize distance
	double w_age	= normLvector(thepower,targetProp);
	
	return diff_age/w_age;
}

double Simulation::calib_distance_ageGapFirstSexSpouseDistrib(Gender g)
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR AGE GAP B/W FIRST SEX AND FIRST SPOUSE DISTRIBUTION
	
	dcMatrix M = (g==female)?_target_ageGapFirstSexSpouseDistribution_f:_target_ageGapFirstSexSpouseDistribution_m;
	
	
	// Extract targets: age
	vector<double> ageBreaks		= M.extractColumn(0);
	vector<double> targetProp		= M.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	ageBreaks.push_back(999.99);
	
	// Retrieve current age distribution
	vector<double> currentAgeDistrib = _population.census_ageGapFirstSexSpouseDistribution(ageBreaks, g);
	// Distance from target
	int thepower = 2;
	
	double diff_age	= distanceLvector(thepower, currentAgeDistrib, targetProp);
	
	// Normalize distance
	double w_age	= normLvector(thepower,targetProp);
	
	double res = diff_age/w_age;
	
	if (std::isnan(res))
		cout<<"DEBUG NAN"<<endl;
	
	return res;
}

double Simulation::calib_distance_singleRatio(Gender g)
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR SINGLE RATIO
	
	
	double currentSingleRatio = _population.census_ratioSingles(g);
	
	int thepower = 2;
	
	double target = (g==female)?_target_singleRatio_f:_target_singleRatio_m;
	
	double diff_singleR = pow(currentSingleRatio - target, thepower);
	
	return pow(diff_singleR/pow(target,thepower),1.0/thepower);
}

double Simulation::calib_distance_malesVisitCSW()
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR SINGLE RATIO
	
	
	double currentProp = _population.census_maleVisitCSW();
	
	int thepower = 2;
	
	
	double diff = pow(currentProp - _target_malesVisitCSW, thepower);
	
	return pow(diff/pow(_target_malesVisitCSW,thepower),1.0/thepower);
}

double Simulation::calib_distance_malesVisitCSW(double maxDurationSinceLastVisit)
{
	// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	// FOR SINGLE RATIO
	
	
	double currentProp = _population.census_maleVisitCSW(maxDurationSinceLastVisit);
	
	int thepower = 2;
	
	// DEBUG
	cout<<endl<<"currentProp="<<currentProp<<endl<<"_target_malesVisitCSW="<<_target_malesVisitCSW<<endl;
	
	
	double diff = pow(currentProp - _target_malesVisitCSW, thepower);
	
	return pow(diff/pow(_target_malesVisitCSW,thepower),1.0/thepower);
}


double Simulation::calib_distance_ageMalesVisitCSWDistrib()
{
	// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	// FOR AGE DISTRIBUTION
	
	// Extract targets: age
	vector<double> ageBreaks		= _target_ageDistribution.extractColumn(0);
	vector<double> targetProp		= _target_ageDistribution.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	ageBreaks.push_back(999.99);
	
	// Retrieve current age distribution
	double maxDurationSinceLastVisit = 1.0;
	vector<double> currentAgeDistrib = _population.census_ageMalesVisitCSWDistribution(ageBreaks,
																					   maxDurationSinceLastVisit);
	
	// Distance from target
	int thepower = 2;
	
	double diff_age	= distanceLvector(thepower, currentAgeDistrib, targetProp);
	
	// Normalize distance
	double w_age	= normLvector(thepower,targetProp);
	
	return diff_age/w_age;
}

double Simulation::calib_distance_nLifeSexPrtnrDistrib(Gender g)
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR LIFETIME NUMBER SEX PARTNERS DISTRIBUTION
	
	dcMatrix M = (g==female)?_target_nLftSexPrtnrDistribution_f:_target_nLftSexPrtnrDistribution_m;
	
	
	// Extract targets: age
	vector<double> Breaks		= M.extractColumn(0);
	vector<double> targetProp	= M.extractColumn(1);
	
	// add 'infinite' break points
	// to make sure the range is broad enough
	Breaks.push_back(99999.99);
	
	// Retrieve current age distribution
	vector<double> currentDistrib = _population.census_nLifeSexPrtnrDistrib(Breaks, g);
	// Distance from target
	int thepower = 2;
	
	double diff	= distanceLvector(thepower, currentDistrib, targetProp);
	
	// Normalize distance
	double w = normLvector(thepower,targetProp);
	
	return diff/w;
}


double Simulation::calib_distance_STIprev()
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR STI (GLOBAL) PREVALENCES
	
	// integrity checks
	int nSTI = _target_STIprevalence_names.size();
	stopif(nSTI==0,"STI names for target calibration not defined!");
	stopif(sumElements(_target_STIprevalence)<1E-10,"STI prevalences targets not defined!");
	
	// Current STI prevalences
	vector<double> currSTIprev;
	for( int i=0;i<nSTI;i++)
	{
		currSTIprev.push_back(_population.STI_prevalence(_target_STIprevalence_names[i]));
	}
	
	// Calculate distance from target
	double diff = distanceLvector(2.0, currSTIprev, _target_STIprevalence);
	double w = normLvector(2.0,_target_STIprevalence);
	
	return diff/w;
}


double Simulation::calib_distance_STIprev_riskGroup()
{
	/// RETURN NORMALIZED DISTANCE TO CALIBRATION TARGET
	/// FOR STI PREVALENCES BY RISK GROUPS
	
	dcMatrix currentSTIprevRiskGroup = _population.STI_prevalence_by_riskGroup(_target_STIprevalence_names);
	
	int thepower = 2;
	double diff_STIprev = distance_Matrix(currentSTIprevRiskGroup,
										  _target_STIprevalence_by_riskGroup,
										  thepower);
	
	// create matrix of 0
	// same dimension as target
	dcMatrix zeromatrix = _target_STIprevalence_by_riskGroup;
	zeromatrix.setAllValues(0.0);
	
	// Normalize distance
	double w_STI = distance_Matrix(_target_STIprevalence_by_riskGroup, zeromatrix, thepower);
	
	return diff_STIprev/w_STI;
}


double Simulation::calib_distance_HIV_prev_age(Gender g)
{
	/// Calculate the distance from the target
	/// for HIV prevalence by age and gender
	
	dcMatrix M;
	if(g==female)	M = _target_HIV_prev_age_f;
	if(g==male)		M = _target_HIV_prev_age_m;
	
	vector<double> curr_HIV_prev_g = _population.STI_prevalence_by_age(HIV, g, M.extractColumn(0)) ;
	
	// because the breaks must cover the range of the all the data
	// and convention in this code is 'left closed', the last row provide no information
	M.removeRow(M.getNbRows()-1);
	vector<double> target_prev =M.extractColumn(1);
	
	// Calculate and normalize distance
	int thepower = 2;
	double dist_g = distanceLvector(thepower, curr_HIV_prev_g, target_prev);
	double w = normLvector(thepower,target_prev);
	
	return dist_g/w;
}

double Simulation::calib_distance_HIV_prev_age()
{
	// For both females and males
	
	double d_f = calib_distance_HIV_prev_age(female);
	double d_m = calib_distance_HIV_prev_age(male);
	
	return d_f+d_m;
	
	//	vector<double> curr_HIV_prev_f = _population.STI_prevalence_by_age(HIV, female, _target_HIV_prev_age_f.extractColumn(0)) ;
	//	vector<double> curr_HIV_prev_m = _population.STI_prevalence_by_age(HIV, male, _target_HIV_prev_age_m.extractColumn(0)) ;
	//
	//	int thepower = 2;
	//	double dist_f = distanceLvector(thepower, curr_HIV_prev_f, _target_HIV_prev_age_f.extractColumn(1));
	//	double dist_m = distanceLvector(thepower, curr_HIV_prev_m, _target_HIV_prev_age_m.extractColumn(1));
	//
	//	// Normalize distance
	//	double w_f = normLvector(thepower,_target_HIV_prev_age_f.extractColumn(1));
	//	double w_m = normLvector(thepower,_target_HIV_prev_age_m.extractColumn(1));
	//
	//	return dist_f/w_f + dist_m/w_m;
}


double Simulation::calib_distanceFromAllTargets()
{
	/// DISTANCE BETWEEN CURRENT VALUES OF
	/// RESPONSE VARIABLES AND TARGETS
	
	// Run the simulation until horizon time
	bool	logIndivInfo	= false;
	bool	doSex			= true;
	bool	traceNetwork	= false;
	int		displayProgress = 0;
	
	runAllEvents_horizon(doSex, logIndivInfo,
						 traceNetwork, displayProgress);
	
	// Retrieve all distances b/w
	// response variables and their respective target:
	
	return calc_distanceFromAllTargets();
}


double Simulation::calc_distanceFromAllTargets()
{
	/// DISTANCE BETWEEN CURRENT VALUES OF
	/// RESPONSE VARIABLES AND TARGETS
	
	// Retrieve all distances b/w
	// response variables and their respective target:
	
	vector<double> dist_var;
	
	dist_var.push_back(calib_distance_ageDistrib());
	dist_var.push_back(calib_distance_ageGapDistrib());
	
	dist_var.push_back(calib_distance_ageFirstSexDistrib(female));
	dist_var.push_back(calib_distance_ageFirstSexDistrib(male));
	
	// not consistent measures
	//dist_var.push_back(calib_distance_ageGapFirstSexSpouseDistrib(female));
	//dist_var.push_back(calib_distance_ageGapFirstSexSpouseDistrib(male));
	
	dist_var.push_back(calib_distance_singleRatio(female));
	dist_var.push_back(calib_distance_singleRatio(male));
	
	bool DHSphaseRecent = false;
	if (DHSphaseRecent)
 {
	 // WHEN DATA SET DOES NOT INCLUDE THESE VARIABLES (DHSphase<=4)
	 
	 double maxDurationSinceLastVisit=1.0;
	 dist_var.push_back(calib_distance_malesVisitCSW(maxDurationSinceLastVisit));
	 
	 dist_var.push_back(calib_distance_ageMalesVisitCSWDistrib());
	 
	 dist_var.push_back(calib_distance_nLifeSexPrtnrDistrib(female));
	 dist_var.push_back(calib_distance_nLifeSexPrtnrDistrib(male));
 }
	
	dist_var.push_back(calib_distance_STIprev() );
	dist_var.push_back(calib_distance_STIprev_riskGroup());
	
	
	// Total distance
	
	// DEBUG
	cout << "calib_distanceFromAllTargets: "<<sumElements(dist_var);
	displayVector(dist_var);
	
	double res = sumElements(dist_var);
	
	return res;
}


double Simulation::calib_distanceFromAllTargets_MC(int nMC)
{
	/// MONTE-CARLO FOR EVALUATION OF MEAN DISTANCE FROM TARGET
	
	double s = 0.0;
	
	//cout <<"BIRTH="<<_population.get_birthRate()<<endl;
	
	// Initialize the population
	Population P0 = _population;
	
	for (int k=0; k<nMC; k++)
	{
		cout << "MC: "<<k+1<<"/"<<nMC << endl; // DEBUG
		
		_population = P0;
		
		s += calib_distanceFromAllTargets();
		
		// Clean up output files not needed when calibrating
		string cmd = "rm -f " + _DIR_OUT + "sexacts.out";
		system(cmd.c_str());
		cmd = "rm -f " + _DIR_OUT + "transmission_events.out";
		system(cmd.c_str());
		cmd = "rm -f " + _DIR_OUT + "sexDistribPartn.out";
		system(cmd.c_str());
		cmd = "rm -f " + _DIR_OUT + "partnership_tentatives.out";
		system(cmd.c_str());
		
	}
	
	return s/nMC;
}


double Simulation::calib_ageDistribution_get_diff(bool doSex)
{
	// Difference between current age distribution and target
	
	// Run the simulation until horizon time
	bool logIndivInfo = false;
	bool traceNetwork = false;
	int displayProgress = 1;
	
	runAllEvents_horizon(doSex,logIndivInfo, traceNetwork, displayProgress);
	
	// Extract and compare the age distribution to the target
	vector<double> ageBreaks = _target_ageDistribution.extractColumn(0);
	vector<double> targetProp = _target_ageDistribution.extractColumn(1);
	
	// Calculate distribution based on age breaks of target
	vector<double> currentAgeDistrib = _population.census_ageDistribution(ageBreaks);
	
	int thepower = 2;
	return distanceLvector(thepower, currentAgeDistrib, targetProp);
}


double Simulation::calib_ageDistribution_get_diff_MC(int nMC, bool doSex)
{
	// MONTE-CARLO FOR EVALUATION OF DIFF FCT
	// FOR CALIBRATION OF AGE DISTRIBUTION
	
	dcMatrix A(0); // DEBUG
	double s = 0.0;
	
	//cout <<"BIRTH="<<_population.get_birthRate()<<endl;
	
	// Initialize the population
	Population P0 = _population;
	
	
	for (int k=0; k<nMC; k++)
	{
		cout << "MC: "<<k+1<<"/"<<nMC << endl; // DEBUG
		
		_population = P0;
		
		s += calib_ageDistribution_get_diff(doSex);
		
		// DEBUG
		vector<double> tmp = _population.census_ageDistribution(_target_ageDistribution.extractColumn(0));
		A.addRowVector(tmp);
	}
	// DEBUG
	A.WriteToFileCSV(_DIR_OUT + "mc.out");
	
	return s/nMC;
}


dcMatrix Simulation::TMP_LHS(vector<double>loVal, vector<double> hiVal, int nLHS, int nMC)
{
	// TMP FCT: EXPLORE THE VALUES OF THE FUNCTION TO MINIMIZE
	dcMatrix M = LatinHypercubeSampling(loVal, hiVal, nLHS);
	M.display();
	
	vector<double> param(loVal.size());
	vector<double> res_k;
	
	_target_ageDistribution.WriteToFileCSV(_DIR_OUT + "calib_test1.out");
	dcMatrix AD(0);
	
	// Population for reset
	Population P0 = _population;
	
	for (int k=0; k< nLHS; k++)
	{
		cout << "LHS: "<< k <<endl;
		param = M.extractRow(k);
		//displayVector(param);
		
		// Reset to initial population size, etc
		_population = P0;
		
		// Update new parameters
		calib_ageDistribution_setParameters(param);
		bool doSex = false;
		res_k.push_back(calib_ageDistribution_get_diff_MC(nMC,doSex));
		
		vector<double> tmp = _population.census_ageDistribution(_target_ageDistribution.extractColumn(0));
		AD.addRowVector(tmp);
	}
	
	AD.WriteToFileCSV(_DIR_OUT + "calib_test2.out");
	
	M.addColVector(res_k);
	//displayVector(res_k);
	M.WriteToFileCSV(_DIR_OUT + "calib_test3.out");
	
	return M;
}

void Simulation::calib_ageDistribution_setParameters(vector<double> param)
{
	
	// ** WARNING **
	// ORDER IS IMPORTANT!
	
	int i=0;
	// birth rates
	//cout << "  I WANT = " <<param[i]<<endl; // DEBUG
	_population.set_birthRate(param[i]); i++;
	
	// death rates
	_population.set_infantMortality(param[i]); i++;
	_population.set_childMortality(param[i]); i++;
	
	_population.set_deathParam_Weibull_shape(param[i]); i++;
	_population.set_deathParam_Weibull_scale(param[i]); i++;
	_population.set_deathParam_Weibull_shape_hiv(param[i]); i++;
	_population.set_deathParam_Weibull_scale_hiv(param[i]); i++;
	
	// Redefine the new population as initial population
	_population_initial = _population;
	
	//_population.displayInfo(false); // DEBUG
}




void Simulation::calib_setParameters(vector<double> param_demographics,
									 vector<double> param_partner_formation,
									 vector<double> param_spousal_progression,
									 vector<double> param_partner_dissolution)
{
	// Log file
	ofstream logfile(_DIR_CALIB + "calib_setParameters_trace.out",ios::app);
	bool doTrace = true;
	
	if (doTrace) logfile<<endl<<"Simulation time: "<<_population.get_simulationTime()<<endl;
	
	// *********************************************
	// ** WARNING **
	// ORDER IS IMPORTANT!
	// MUST BE SAME ORDER AS IN calib_get_parameters
	// *********************************************
	
	//////////// DEMOGRAPHICS ////////////
	
	int i=0;
	
	// birth rates
	_population.set_birthRate(param_demographics[i]); i++;
	
	
	// death rates
	_population.set_infantMortality(param_demographics[i]); i++;
	_population.set_childMortality(param_demographics[i]); i++;
	
	_population.set_deathParam_Weibull_shape(param_demographics[i]); i++;
	_population.set_deathParam_Weibull_scale(param_demographics[i]); i++;
	_population.set_deathParam_Weibull_shape_hiv(param_demographics[i]); i++;
	_population.set_deathParam_Weibull_scale_hiv(param_demographics[i]); i++;
	
	if (i!=param_demographics.size())
	{
		cout << endl << "ERROR [calib_setParameters]: 'param_demographics' size no coherent "<<endl;
		exit(1);
	}
	
	if (doTrace)
	{
		logfile<<"birth rate: "<< _population.get_birthRate() <<endl;
		logfile<<"infantMortality: "<< _population.get_infantMortality() <<endl;
		logfile<<"childMortality: "<< _population.get_childMortality() <<endl;
		logfile<<"Weibull_shape: "<< _population.get_deathParam_Weibull_shape() <<endl;
		logfile<<"Weibull_scale: "<< _population.get_deathParam_Weibull_scale() <<endl;
		logfile<<"Weibull_shape_hiv: "<< _population.get_deathParam_Weibull_shape_hiv() <<endl;
		logfile<<"Weibull_scale_hiv: "<< _population.get_deathParam_Weibull_scale_hiv() <<endl;
	}
	
	// *********************************************
	// ** WARNING **
	// ORDER IS IMPORTANT!
	// MUST BE SAME ORDER AS IN calib_get_parameters
	// *********************************************
	
	
	//////////// PARTNERSHIPS ////////////
	
	// =================
	// === FORMATION ===
	// =================
	
	i = 0; // reset counter for new vector param_partner_formation
	
	// Formation -- Max rate
	_population.set_formation_MaxRate(param_partner_formation[i]); i++;
	
	// Formation -- Age
	// OLD: _population.set_formation_Age_female(param_partner_formation[i],param_partner_formation[i+1]); i=i+2;
	
	_population.set_formation_age_fullstart(param_partner_formation[i]); i++;
	_population.set_formation_age_pivot(param_partner_formation[i]); i++;
	_population.set_formation_age_shape(param_partner_formation[i]); i++;
	_population.set_formation_age_fmin(param_partner_formation[i]) ; i++;
	
	// Formation -- Age Gap
	_population.set_formation_agegap_mean(param_partner_formation[i]); i++;
	_population.set_formation_agegap_var(param_partner_formation[i]); i++;
	_population.set_formation_agegap_fmin(param_partner_formation[i]); i++;
	
	
	//	_population.set_formation_AgeGap_female(param_partner_formation[i],param_partner_formation[i+1]); i=i+2;
	//	_population.set_formation_Age_AgeGap_correl_shape(param_partner_formation[i],param_partner_formation[i+1]); i=i+2;
	
	// Formation -- Risk group
	vector<double> tmp;
	for (int k=0; k<_population.get_formation_RiskGroup().size(); k++)
	{
		tmp.push_back(param_partner_formation[i]);i++;
	}
	_population.set_formation_RiskGroup(tmp);
	
	// Formation -- Partner deficit
	tmp.clear();
	for (int k=0; k<_population.get_formation_PartnDeficit().size(); k++)
	{
		tmp.push_back(param_partner_formation[i]);i++;
	}
	_population.set_formation_PartnDeficit(tmp);
	
	// Formation -- STI symptoms
	_population.set_formation_STIsymptom_f(param_partner_formation[i]); i++;
	_population.set_formation_STIsymptom_m(param_partner_formation[i]); i++;
	
	
	
	if (i!=param_partner_formation.size())
	{
		cout << endl << "ERROR [calib_setParameters]: 'param_partner_formation' size no coherent "<<endl;
		exit(1);
	}
	
	
	if (doTrace)
	{
		logfile<<"formation_MaxRate: "<< _population.get_formation_MaxRate() <<endl;
		//		logfile<<"formation_meanAge_female: "<< _population.get_formation_meanAge_female() <<endl;
		//		logfile<<"formation_varAge_female: "<< _population.get_formation_varAge_female() <<endl;
		logfile<<"_formation_agegap_mean: "<< _population.get_formation_agegap_mean() <<endl;
		//	logfile<<"_formation_agegap_var: "<< _population.get_formation_agegap_var() <<endl;
		logfile<<"formation_RiskGroup: ";
		for (int a=0; a<_population.get_formation_RiskGroup().size(); a++) {
			logfile<<_population.get_formation_RiskGroup()[a]<<"  ";}
		logfile<<endl;
		logfile<<"formation_PartnDeficit: ";
		for (int a=0; a<_population.get_formation_PartnDeficit().size(); a++) {
			logfile<<_population.get_formation_PartnDeficit()[a]<<"  ";}
		logfile<<endl;
	}
	
	
	// *********************************************
	// ** WARNING **
	// ORDER IS IMPORTANT!
	// MUST BE SAME ORDER AS IN calib_get_parameters
	// *********************************************
	
	
	// ===================
	// === DISSOLUTION ===
	// ===================
	
	i = 0; // reset counter for new vector param_partner_dissolution
	unsigned int vecsize=0;
	
	// Dissolution -- Max rate
	_population.set_dissolution_MaxRate(param_partner_dissolution[i]); i++;
	
	// Dissolution -- Spousal reduction
	_population.set_dissolution_spouse(param_partner_dissolution[i]); i++;
	
	// Dissolution -- Risk group
	tmp.clear();
	vecsize =_population.get_dissolution_RiskGroup().size();
	
	for (int k=0; k<vecsize; k++)
	{
		tmp.push_back(param_partner_dissolution[i]); i++;
	}
	_population.set_dissolution_risk_param(tmp);
	
	
	// Dissolution -- Duration
	tmp.clear();
	vecsize =_population.get_dissolution_duration().size();
	
	for (int k=0; k<vecsize; k++)
	{
		tmp.push_back(param_partner_dissolution[i]); i++;
	}
	_population.set_dissolution_duration_param(tmp);
	
	
	// Dissolution -- Age
	tmp.clear();
	vecsize =_population.get_dissolution_age().size();
	
	for (int k=0; k<vecsize  ; k++)
	{
		tmp.push_back(param_partner_dissolution[i]); i++;
	}
	_population.set_dissolution_age_param(tmp);
	
	
	// Dissolution -- Partner deficit
	tmp.clear();
	vecsize =_population.get_dissolution_PartnerDeficit().size();
	for (int k=0; k< vecsize; k++)
	{
		tmp.push_back(param_partner_dissolution[i]); i++;
	}
	_population.set_dissolution_PartnerDeficit(tmp);
	
	
	// Dissolution -- Concurrent partnerships and age
	_population.set_dissolution_ageConcPartn(param_partner_dissolution[i]); i++;
	
	// Dissolution --STI symptoms
	_population.set_dissolution_STI_symptom(param_partner_dissolution[i]); i++;
	
	
	//int dummy =param_partner_dissolution.size();
	
	
	if (i!=param_partner_dissolution.size())
	{
		cout << endl << "ERROR [calib_setParameters]: 'param_partner_dissolution' size no coherent "<<endl;
		exit(1);
	}
	
	
	if (doTrace)
	{
		logfile<<"dissolution_MaxRate: "<< _population.get_dissolution_MaxRate() <<endl;
		
		logfile<<"dissolution_RiskGroup: ";
		for (int a=0; a<_population.get_dissolution_RiskGroup().size(); a++) {
			logfile<<_population.get_dissolution_RiskGroup()[a]<<"  ";}
		logfile<<endl;
		logfile<<"dissolution_duration: ";
		for (int a=0; a<_population.get_dissolution_duration().size(); a++) {
			logfile<<_population.get_dissolution_duration()[a]<<"  ";}
		logfile<<endl;
	}
	
	
	
	// === SPOUSAL PROGRESSION ===
	
	// *********************************************
	// ** WARNING **
	// ORDER IS IMPORTANT!
	// MUST BE SAME ORDER AS IN calib_get_parameters
	// *********************************************
	
	
	i = 0; // reset counter for new vector param_spousal_progression
	
	// Spousal progress -- Max Rate
	_population.set_spousalProgress_maxRate(param_spousal_progression[i]); i++;
	
	// Spousal progress -- female's age and age gap
	double m_af	= param_spousal_progression[i]; i++;
	double v_af	= param_spousal_progression[i]; i++;
	double m_g	= param_spousal_progression[i]; i++;
	double v_g	= param_spousal_progression[i]; i++;
	
	_population.set_spousalProgress_Age_param(m_af,v_af,m_g,v_g);
	
	// Spousal progress -- partnership duration
	double k1	= param_spousal_progression[i]; i++;
	double k2	= param_spousal_progression[i]; i++;
	
	_population.set_spousalProgress_Duration_param(k1,k2);
	
	// Spousal progress -- age gaps against other spouses
	double m_delta	= param_spousal_progression[i]; i++;
	double v_delta	= param_spousal_progression[i]; i++;
	
	_population.set_spousalProgress_DiffAgeGap_param(m_delta, v_delta);
	
	
	if (i!=param_spousal_progression.size())
	{
		cout << endl << "ERROR [calib_setParameters]: 'param_spousal_progression' size no coherent "<<endl;
		exit(1);
	}
	
	if (doTrace)
	{
		logfile<<"spousalProgress_maxRate: "<< _population.get_spousalProgress_maxRate() <<endl;
		logfile<<"spousalProgress_meanAge_f: "<< _population.get_spousalProgress_meanAge_f() <<endl;
		logfile<<"spousalProgress_varAge_f: "<< _population.get_spousalProgress_varAge_f() <<endl;
		logfile<<"spousalProgress_meanGap: "<< _population.get_spousalProgress_meanGap() <<endl;
		logfile<<"spousalProgress_varGap: "<< _population.get_spousalProgress_varGap() <<endl;
		logfile<<"spousalProgress_durationK1: "<< _population.get_spousalProgress_durationK1() <<endl;
		logfile<<"spousalProgress_durationK2: "<< _population.get_spousalProgress_durationK2() <<endl;
		logfile<<"spousalProgress_meanDiffAgeGap: "<< _population.get_spousalProgress_meanDiffAgeGap() <<endl;
		logfile<<"spousalProgress_varDiffAgeGap: "<< _population.get_spousalProgress_varDiffAgeGap() <<endl;
	}
	
	
	
	// Redefine the new population as initial population
	_population_initial = _population;
	
}


void Simulation::calib_setParameters_STI(vector<double> param_sexActivity,
										 vector<double> param_stiFeatures_virus,
										 vector<double> param_stiFeatures_bacteria,
										 vector<double> param_stiFeatures_protozoa,
										 vector<double> param_stiFeatures_fungus)
{
	//////////// SEXUAL ACTIVITY ////////////
	
	// **** WARNING ****
	// ORDER IS IMPORTANT!
	// *******************
	
	int i=0;
	
	//Sex acts rates BOTH male and then female
	_population.set_sexActMaxRate(param_sexActivity[i],param_sexActivity[i+1]); i=i+2;
	
	// Preference for sex act with spouse
	_population.set_sexAct_proba_distribute_partnerTypes_prefSpouse(param_sexActivity[i]); i++;
	
	// Sex acts reduction - age
	vector<double> tmp;
	tmp.push_back(param_sexActivity[i]); i++;
	tmp.push_back(param_sexActivity[i]); i++;
	tmp.push_back(param_sexActivity[i]); i++;
	
	_population.set_sexAct_reduce_age_param(tmp);
	tmp.clear();
	
	// Sex acts reduction - risk group
	
	tmp.push_back(param_sexActivity[i]); i++;
	
	_population.set_sexAct_reduce_risk_param(tmp);
	tmp.clear();
	
	// Sex acts reduction - symptoms
	
	tmp.push_back(param_sexActivity[i]); i++; // male
	tmp.push_back(param_sexActivity[i]); i++; // female
	
	_population.set_sexAct_reduce_STIsymptom_param(tmp);
	tmp.clear();
	
	// **** WARNING ****
	// ORDER IS IMPORTANT!
	// *******************
	
	// Sex acts reduction - number of current partners
	
	tmp.push_back(param_sexActivity[i]); i++;
	
	_population.set_sexAct_reduce_nPartner_param(tmp);
	tmp.clear();
	
	
	// Sex acts Distribution - preference for CSW
	
	tmp.push_back(param_sexActivity[i]); i++;
	tmp.push_back(param_sexActivity[i]); i++;
	
	_population.set_sexAct_proba_distribute_partnerTypes_sexWorkParam(tmp);
	tmp.clear();
	
	// Sex acts Distribution - proba sex with CSW
	
	tmp.push_back(param_sexActivity[i]); i++;
	tmp.push_back(param_sexActivity[i]); i++;
	
	_population.set_sexAct_proba_sexWorker_param(tmp);
	tmp.clear();
	
	// Sex acts Distribution - proba sex with CSW
	
	_population.set_sexAct_CostSexWork_reduction(param_sexActivity[i]);i++;
	
	
	
	//////////// STI TRANSMISSION ////////////
	
	
	// **** WARNING ****
	// ORDER IS IMPORTANT!
	// *******************
	
	
	/// ---- VIRUSES ---- ///
	
	int j=0;
	
	// ** ORDER IS IMPORTANT! **
	
	// === HIV ===
	_population.set_STI_probaMaxSexTransm(HIV, param_stiFeatures_virus[j]); j++;
	
	// === HSV2 ===
	_population.set_STI_probaMaxSexTransm(HSV2, param_stiFeatures_virus[j]); j++;
	
	// === HPV ===
	_population.set_STI_probaMaxSexTransm(HPV, param_stiFeatures_virus[j]); j++;
	
	
	
	/// ---- BACTERIA ---- ///
	
	j=0;	// reset counter for new parameters set
	
	// **** WARNING ****
	// ORDER IS IMPORTANT!
	// *******************
	
	// === Ct ===
	_population.set_STI_probaMaxSexTransm(Ct, param_stiFeatures_bacteria[j]); j++;
	
	// === Ng ===
	_population.set_STI_probaMaxSexTransm(Ng, param_stiFeatures_bacteria[j]); j++;
	
	// === Tp ===
	_population.set_STI_probaMaxSexTransm(Tp, param_stiFeatures_bacteria[j]); j++;
	
	// === Hd ===
	_population.set_STI_probaMaxSexTransm(Ng, param_stiFeatures_bacteria[j]); j++;
	
	// === Bv ===
	_population.set_STI_probaMaxSexTransm(Bv, param_stiFeatures_bacteria[j]); j++;
	
	
	/// ---- PROTOZOA ---- ///
	
	j=0;	// reset counter for new parameters set
	
	
	// === Tv ===
	
	_population.set_STI_probaMaxSexTransm(Tv, param_stiFeatures_protozoa[j]); j++;
	
	
	
	/// ---- Fungus ---- ///
	
	j=0;	// reset counter for new parameters set
	
	
	// === Cv ===
	
	// Not implemented yet
	
	
	/// APPLY ALL CHANGES OF STI FEATURES TO EVERY INDIVIDUALS
	_population.STI_update_templateToAllIndividuals();
	
	
	// Redefine the new population as initial population
	_population_initial = _population;
	
}


vector<double> Simulation::sensitivity_calib_distance(vector<double> param_DMG,
													  vector<double> param_FORM,
													  vector<double> param_SPOUSAL,
													  vector<double> param_DISSOL,
													  double relativeBump,int nMC)
{
	// CALCULATE SENSITIVITY
	// OF THE DISTANCE TO CALIBRATION TARGETS
	// WITH RESPECT TO ALL PARAMETERS
	
	
	// Baseline value
	force_seed_reset();	// Make sure the same seed is used
	cout<<endl<<"1st try:"<<endl;
	double baseline = calib_distanceFromAllTargets_MC(nMC);
	
	//	// DEBUG:
	cout<<endl<<"2nd try:"<<endl;
	force_seed_reset();
	_population = _population_initial;
	//double baseline2 = calib_distanceFromAllTargets_MC(nMC);
	
	
	// Number of input parameters
	vector<double> res;
	int n_DMG = param_DMG.size();
	int n_FORM = param_FORM.size();
	int n_SPOUSAL = param_SPOUSAL.size();
	int n_DISSOL = param_DISSOL.size();
	int n = n_DMG+n_FORM+n_SPOUSAL+n_DISSOL;
	
	for (int i=0; i<n; i++)
	{
		
		// Make sure the same seed is used
		force_seed_reset();
		
		// DEBUG
		
		cout << endl << "sensitivity_calib_distance param #"<<i<<"/"<<n<<endl;
		
		
		// Re-initialize the population at each step
		_population = _population_initial;
		
		// Bump the relevant parameter
		double tmp;
		bool check_case_found = false;
		
		if (i<n_DMG)
		{
			tmp = param_DMG[i];
			param_DMG[i] = (1.0+relativeBump)*param_DMG[i];
			check_case_found = true;
		}
		
		if (n_DMG<=i && i<(n_DMG+n_FORM))
		{
			tmp = param_FORM[i];
			param_FORM[i] = (1.0+relativeBump)*param_FORM[i];
			check_case_found = true;
		}
		
		if (n_DMG+n_FORM<=i && i<(n_DMG+n_FORM+n_SPOUSAL))
		{
			tmp = param_SPOUSAL[i];
			param_SPOUSAL[i] = (1.0+relativeBump)*param_SPOUSAL[i];
			check_case_found = true;
		}
		
		if (n_DMG+n_FORM+n_SPOUSAL<=i && i<(n_DMG+n_FORM+n_SPOUSAL+n_DISSOL))
		{
			tmp = param_DISSOL[i];
			param_DISSOL[i] = (1.0+relativeBump)*param_DISSOL[i];
			check_case_found = true;
		}
		
		if (!check_case_found)
		{
			cout << endl << "ERROR [sensitivity_calib_distance]: bump index not found!"<<endl;
			exit(1);
		}
		
		// Set the new parameter value
		calib_setParameters(param_DMG,param_FORM,param_SPOUSAL,param_DISSOL);
		
		
		//DEBUG:
		//_population.displayInfo(false);
		//_population_initial.displayInfo(false);
		_population =  _population_initial;
		
		
		// Calculate the function value with updated parameter
		double bumpedDistance = calib_distanceFromAllTargets_MC(nMC);
		res.push_back((bumpedDistance-baseline)/baseline);
		
		// back to initial value
		if (i<n_DMG) param_DMG[i] = tmp;
		if (n_DMG<=i && i<(n_DMG+n_FORM)) param_FORM[i] = tmp;
		if (n_DMG+n_FORM<=i && i<(n_DMG+n_FORM+n_SPOUSAL)) param_SPOUSAL[i] = tmp;
		if (n_DMG+n_FORM+n_SPOUSAL<=i && i<(n_DMG+n_FORM+n_SPOUSAL+n_DISSOL)) param_DISSOL[i] = tmp;
	}
	
	return res;
}


vector<double> Simulation::TMPsensi(vector<double> param, double relativeBump, int nMC)
{
	// CALCULATE SENSI
	
	vector<double> res;
	
	calib_ageDistribution_setParameters(param);
	bool doSex = false;
	
	double res0 = calib_ageDistribution_get_diff_MC(nMC,doSex);
	
	for (int k=0; k<param.size(); k++)
	{
		cout << "Sensi param #"<<k<<endl;
		
		_population = _population_initial;
		double tmp = param[k];
		
		param[k] = param[k]*(1.0+relativeBump); // bump
		
		calib_ageDistribution_setParameters(param);
		displayVector(param);
		
		bool doSex = false;
		res.push_back(calib_ageDistribution_get_diff_MC(nMC, doSex)-res0);
		
		param[k] = tmp; // back to initial value
	}
	vectorToFile(res, _DIR_OUT + "sensi.out");
	return res;
}


double Simulation::calib_ageDistribution_F(const gsl_vector *v, void *params)
{
	// Define the function to minimize here
	
	int dim = v->size;
	
	if (dim!=7)
	{
		cout << endl << "*** ERROR [Simulation::calib_ageDistribution_F]: ";
		cout << " wrong number of parameters. Expected 7, but have "<<dim<<endl;
		exit(1);
	}
	
	vector<double> param(dim);
	param = allocateGSLVector(v);
	
	// Reset population
	_population = _population_initial;
	
	// Update the values of the parameter to minimize
	calib_ageDistribution_setParameters(param);
	
	bool doSex = false;
	
	return calib_ageDistribution_get_diff_MC(_nMC_calibration, doSex);
}




// =================================================================
// =================================================================




void Simulation::runSimulation(bool logIndividualInfo)
{
	//
	// -- DEPRECATED --
	//
}



void Simulation::save_outputs_demog(string pathFolder,
									unsigned int iMC, unsigned int idate)
{
	get_population().save_outputs_demog(pathFolder,
										_target_ageDistribution.extractColumn(0),
										iMC,idate);
}


void Simulation::save_outputs_prtnr(string pathFolder,
									unsigned int iMC,unsigned int idate)
{
	get_population().save_outputs_prtnr(pathFolder,
										_target_ageGapDistribution.extractColumn(0),
										iMC,idate);
}


void Simulation::save_outputs_sex(string pathFolder,
								  unsigned int iMC,unsigned int idate)
{
	/// Save sexual behaviour output from model
	
	vector<double> dummy; for(int i=0;i<10;i++) dummy.push_back(i);
	
	vector<double> ab_gap_f(0);
	vector<double> ab_gap_m(0);
	
	if (_target_ageGapFirstSexSpouseDistribution_f.val.size()>0)
	{
		ab_gap_f = _target_ageGapFirstSexSpouseDistribution_f.extractColumn(0);
		ab_gap_m = _target_ageGapFirstSexSpouseDistribution_m.extractColumn(0);
	}
	
	get_population().save_outputs_sex(pathFolder,
									  _target_ageFirstSexDistribution_f.extractColumn(0),
									  _target_ageFirstSexDistribution_m.extractColumn(0),
									  ab_gap_f,
									  ab_gap_m,
									  dummy, //_target_nLftSexPrtnrDistribution_f.extractColumn(0),
									  dummy, //_target_nLftSexPrtnrDistribution_m.extractColumn(0),
									  dummy, //_target_ageMalesVisitCSWDistribution.extractColumn(0)
									  iMC,idate);
}


void Simulation::save_outputs_epi(string pathFolder,
								  unsigned int iMC,unsigned int idate)
{
	vector<double> ad_f(0);
	vector<double> ad_m(0);
	
	if (_target_HIV_prev_age_f.val.size()>0) ad_f = _target_HIV_prev_age_f.extractColumn(0);
	if (_target_HIV_prev_age_m.val.size()>0) ad_m = _target_HIV_prev_age_m.extractColumn(0);
	
	get_population().save_outputs_epi(pathFolder,
									  _target_STIprevalence_names,
									  ad_f,
									  ad_m,
									  iMC,idate);
}


void Simulation::save_incidence(unsigned int iter_mc){
	/// Save incidence time series to a file
	/// (must be used once the simulation has run until horizon)
	
	string filename_incidence = _DIR_OUT + "incid_mc"+to_string(iter_mc)+".out";
	vector<string> header;
	
	for(int s=0; s<get_population().get_nSTImodelled(); s++)
		header.push_back(STInameString(get_population().get_STI()[s].get_name()));
	header.push_back("time");
	
	dcMatrix tmp = get_STI_incidence();
	tmp.addColVector(vector_seq_by(0.0,_horizon,_timeStep));
	tmp.WriteToFileCSV(filename_incidence,header);
}


void Simulation::save_prevalence(unsigned int iter_mc){
	/// Save prevalence time series to a file
	/// (must be used once the simulation has run until horizon)
	
	string filename_prev = _DIR_OUT + "prev_mc"+to_string(iter_mc)+".out";
	vector<string> header;
	
	for(int s=0; s<get_population().get_nSTImodelled(); s++)
		header.push_back(STInameString(get_population().get_STI()[s].get_name()));
	header.push_back("time");
	
	dcMatrix tmp = get_STI_prevalence();
	tmp.addColVector(vector_seq_by(0.0,_horizon,_timeStep));
	tmp.WriteToFileCSV(filename_prev,header);
}




// ===================================================================
// ===================================================================
// =====  INTERVENTION  =====
// ===================================================================
// ===================================================================

void Simulation::activate_intervention(int i)
{
	/// ACTIVATE TREATMENT OR VACCINATION BASED ON INTERVENTION FEATURES
	/// (assumes simulation time is between start and end dates of intervention)
	/// God-like treatment: assumes we know everyone infected and
	/// everyone has access to intervention
	
	// Integrity checks
	stopif(_intervention[i].get_annCvgRate()<=0, "Proportion of intervention must be positive");
	
	// retrieve intervention features
	double target_dt	= _intervention[i].get_annCvgRate()*_timeStep;
	STIname sti			= _intervention[i].get_stiname();
	string interv_type	= _intervention[i].get_type();
	
	
	// * WARNING *
	// the proportion is understood as a rate and applied every time step
	// For example, propPerYear = 0.5 means 50% of infected will be treated over one year.
	// If the time step is 0.1 year, then the target proportion for that timestep
	// will be 0.1*0.5 = 0.05
	
	
	// check if sampling is actually needed
	// (when prop>1 it means everyone gets treatment during that time step)
	bool do_sample = (target_dt < 1.0);
	
	// count individuals (for log files)
	unsigned long cnt = 0;  // <-- number treated
	unsigned long cnt2 = 0; // <-- number targeted
	
	// scan the intervention type
	bool doTreat_mass		= (interv_type=="treatment_mass");
	bool doTreat_symptom	= (interv_type=="treatment_symptom");
	bool doVacc_mass		= (interv_type=="vaccination_mass");
	bool doVacc_symptom		= (interv_type=="vaccination_symptom");
	bool doVacc_female		= (interv_type=="vaccination_female");
	bool doVacc_femaleYoung	= (interv_type=="vaccination_femaleYoung");
	bool doVacc_hiRisk		= (interv_type=="vaccination_highRisk");
	
	bool doTreatment	= (interv_type.substr(0,9)=="treatment");
	bool doVaccination	= (interv_type.substr(0,11)=="vaccination");
	
	// Integrity checks
	bool type_known =	doTreat_mass || doTreat_symptom ||
						doVacc_mass || doVacc_symptom ||
						doVacc_female || doVacc_femaleYoung ||
						doVacc_hiRisk;
	stopif(!type_known, "Unknown intervention type!");
	
	unsigned int sti_i = positionSTIinVector(sti, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		
		// == Filter who is targeted by intervention ==
		
		Individual indiv = _population.getIndividual(uid);
		
		if(indiv.isAlive()){
			
			bool indivIsTargeted	= false;
			
			bool isSymptomatic		= indiv.get_STIsymptom()[sti_i];
			bool alreadyVacc		= indiv.get_STI_vacc()[sti_i];  
			
			if(doTreat_mass)	indivIsTargeted = indiv.get_STIduration()[sti_i]>0; // slow code: indiv.STI_infected(sti);
			if(doTreat_symptom) indivIsTargeted = indiv.get_STIsymptom()[sti_i]; // slow code: is_symptomatic(sti);
			if(doVacc_mass){
				// everyone is targeted
				// but exclude the ones already vaccinated
				indivIsTargeted =  !alreadyVacc;
			}
			if(doVacc_symptom){
				// only symptomatic (for THAT STI) individuals
				indivIsTargeted = isSymptomatic && !alreadyVacc;
			}
			if(doVacc_female){
				// only (all) females
				indivIsTargeted = (indiv.get_gender()==female) && !alreadyVacc;
			}
			if(doVacc_femaleYoung){
				// only _young_ females
				double young_age = 15.0;
				
				bool tmp1 = (indiv.get_gender()==female);
				bool tmp2 = (indiv.get_age()<young_age);
				bool tmp3 = !alreadyVacc;
				indivIsTargeted = tmp1 && tmp2 && tmp3;
			}
			if(doVacc_hiRisk){
				// only the highest risk group
				int mxrg = _population.get_maxRiskGroup();
				indivIsTargeted = (indiv.get_riskGroup()== mxrg) && !alreadyVacc;
			}
			
			// == Apply intervention on filtered individuals ==
			
			if(indivIsTargeted ){
				cnt2++;
				// no sampling because rate of intervention>1
				// (speed up code execution)
				if (!do_sample){
					if(doTreatment)		treat_indiv(uid, sti);
					if(doVaccination)	vaccinate_indiv(uid, sti);
					cnt++;
				}
				
				// random sample
				if (do_sample && (uniform01()<= target_dt) ){
					if(doTreatment)		treat_indiv(uid, sti);
					if(doVaccination)	vaccinate_indiv(uid, sti);
					cnt++;
				}
			}
		}
	} // end loop on individuals
	
	
	// Intervention information
	vector<double> v;
	vector<string> cnames;

	cnames.push_back("iMC");
	v.push_back(_MC_trial_iter);
	cnames.push_back("time");
	v.push_back(_population.get_simulationTime());
	cnames.push_back("n_reached");
	v.push_back(cnt);
	cnames.push_back("n_targeted");
	v.push_back(cnt2);
	
	string rowname = interv_type + "_" + STInameString(sti);
	_df_interv.addrow(rowname, v);
	if(_df_interv.get_nrows()==1) _df_interv.set_colname(cnames);
}




// ===================================================================
// ===================================================================
// =====  TREATMENT  =====
// ===================================================================
// ===================================================================


void Simulation::treat_indiv(unsigned long uid, STIname sti){
	_population.treat_indiv(uid, sti);
}


void Simulation::cure_indiv(unsigned long uid, STIname stiname)
{
	/// Cure and individual from a given STI
	/// (must be treated beforehand and Adherence>0.8)
	
	stopif(!_population.getIndividual(uid).STI_treated(stiname),
		   "Cannot cure without treatment!");
	stopif(!(_population.getIndividual(uid).get_STIduration(stiname)>0),
		   "Cannot cure inexistent STI!");
	
	// Cure only if:
	// - treated
	// - treatment is biologically successfull
	// - treatment duration larger than optimal duration
	// - STI is actually curable(!)
	
	
	int i_sti = positionSTIinVector(stiname, _population.get_STI());
	double optim_duration = _population.get_STI()[i_sti].get_optimalTreatmentDuration();
	
	if (_population.getIndividual(uid).get_STItreatDuration(stiname) > optim_duration &&
		_population.getIndividual(uid).get_STItreatTMS(stiname)==1)
	{
		double A = _population.getIndividual(uid).get_STItreatAdherence(stiname);
		if(A>0.8) {
			_population.cure_indiv(uid, stiname);
			//DEBUG
			//cout << "UID "<<uid<< " cured from "<<STInameString(stiname)<<endl;
		}
	}
}


void Simulation::update_cure(STIname sti)
{
	/// Update the cure for all individuals treated against that STI
	/// (mostly when treatment duration has reached optimal duration
	
	int sti_i = positionSTIinVector(sti, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		Individual tmp = _population.getIndividual(uid);
		if (tmp.isAlive() &&
			(tmp.get_STItreatDuration()[sti_i]>0) &&
			tmp.get_STIduration()[sti_i]>0) cure_indiv(uid,sti);
		
	}

// SLOW CODE: DELETE when sure (2015-10-06)
//	for (unsigned long uid=0; uid<_population.get_size(); uid++){
//		Individual tmp = _population.getIndividual(uid);
//		if (tmp.isAlive() &&
//			tmp.STI_treated(sti) &&
//			tmp.get_STIduration(sti)>0) cure_indiv(uid,sti);
}



// ===================================================================
// ===================================================================
// =====  VACCINATION  =====
// ===================================================================
// ===================================================================


void Simulation::vaccinate_indiv(unsigned long uid, STIname stiname)
{
	/// Vaccinate an individual against a given STI
	
	_population.vaccinate_indiv(uid, stiname);
}


void Simulation::update_vacc(STIname stiname){
	/// Update immunity provided by vaccination
	
	int sti_i = positionSTIinVector(stiname, _population.get_STI());
	
	for (unsigned long uid=0; uid<_population.get_size(); uid++){
		
		Individual tmp = _population.getIndividual(uid);
		
		if (tmp.isAlive() &&
			tmp.get_STI_vacc()[sti_i]){
			
			// immunity resulting from vaccination
			// (1 if success; 0 if failed)
			double imm_vax = tmp.get_STI_immunity(stiname);
			// time since vaccination:
			double tv = _simulationTime - tmp.get_STI_vacc_time()[sti_i];
			// waning rate:
			double w = _population.get_STI()[sti_i].get_vacc_waneRate();
			// immunity:
			double imm = imm_vax * exp(-tv * w);
			// update the value of waning immunity:
			_population.set_STI_immunity(uid, stiname, imm);
			
			// update susceptibility factor:
			double prev_sf = _population.get_STIsusceptFactor(stiname,uid);
			_population.set_STIsusceptFactor(uid, stiname, prev_sf*(1-imm));
		}
	}
}



// ===================================================================
// ===================================================================
// =====  MISCELLENAOUS  =====
// ===================================================================
// ===================================================================

void Simulation::displayInfo()
{
	coutline(80);
	cout<<"   * * * SIMULATION INFORMATION * * * "<<endl;
	coutline(80);
	
	cout << " Horizon:\t"<<_horizon<<endl;
	cout << " TimeStep:\t"<<_timeStep<<endl;
	cout << endl;
	
	cout << " Population information:";
	_population.displayInfo(false);
	
	coutline(80);
	
	cout << " Calibration Schedule"<<endl;
	
	cout << "Calibration Times:";
	displayVector(_calibrationTime);
	
	cout << "Output to be calibrated:"<<endl;
	
	//	for (int i=0; i<_calibrationTime.size(); i++)
	//	{
	//		cout <<endl<< "-- at time " << _calibrationTime[i]<<":";
	//		displayVector(_calibrationTargetFile[i]);
	//	}
	
	coutline(80);
}


