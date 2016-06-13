/*
 *  simulation.h
 *  LocalSTI
 *
 *  Created by David Champredon  on 13-08-30.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef simulation_h
#define simulation_h


#include "population.h"
#include "dcDataFrame.h"
#include "nursery.h"
#include "intervention.h"

class Simulation 
{
	
	// Monte-Carlo parameters
	double			_horizon;
	double			_timeStep;
	double			_simulationTime;
	vector<double>	_schedule;
	
	unsigned int	_MC_trial_iter;		// number of this Monte-Carlo trial (out of the total number of trials)
	
	// Population
	Population		_population_initial;	// store the initial population (used in calibration)
	Population		_population;
	
	// New borns
	unsigned long	_newborn_timestep;		// number of newborns (during a Simulation timestep)
	Nursery			_nursery;

	// STI related
	dcMatrix		_STI_incidence;		// rows=time ; columns = STI
	dcMatrix		_STI_prevalence;	// rows=time ; columns = STI
	
	// Interventions (Treatments, Vaccinations)
	vector<Intervention>	_intervention;

	// Data frames that record outputs
	dcDataFrame		_df_sim;		// time series
	dcDataFrame		_df_interv;		// intervention details
	
	
	// ===============================
	// === TARGETS FOR CALIBRATION ===
	// ===============================
	
	int		_nMC_calibration;	// Number of MC runs for the simulation; used for calibration purposes
	
	// Vector of times during the simulation
	// when a calibration is expected
	vector<double>				_calibrationTime;
	
	
	// For each calibration time, there can be
	// an arbitrary number of output to calibrate
	
	// filename of the associated target
	vector< vector<string> >	_calibrationTargetFile;
	// type of target (e.g., prevalence, age distribution, etc)
	vector< vector<string> >	_calibrationTargetType;

	// Relative weight in the glabal calibration score of each calibration target
	vector< vector<double> >	_calibrationWeights;
	
	// Weighted Distances model from target (needs to be calculated)
	vector< vector<double> >	_calibrationDistances;
	
	
	// Proportion of all sexually active indiv who have no partners
	// Segregated gender
	double			_target_singleRatio_f;
	double			_target_singleRatio_m;
	
	
	// FOR ALL AGE (AGE GAP) DISTRIBUTIONS:
	// 1st col = age breaks
	// 2nd col = proportion of population in associated age range age[i]<prop[i]<age[i+1]
	dcMatrix		_target_ageDistribution;
	dcMatrix		_target_ageGapDistribution;
	dcMatrix		_target_ageFirstSexDistribution_f;
	dcMatrix		_target_ageFirstSexDistribution_m;
	dcMatrix		_target_ageGapFirstSexSpouseDistribution_f;	// gap b/w 1st sex and 1st spouse
	dcMatrix		_target_ageGapFirstSexSpouseDistribution_m;	// gap b/w 1st sex and 1st spouse
	
	dcMatrix		_target_nLftSexPrtnrDistribution_f;	// Distribution of lifetime number of sex partners for females
	dcMatrix		_target_nLftSexPrtnrDistribution_m;	// Distribution of lifetime number of sex partners for males
	
	double			_target_malesVisitCSW;					// propotion of males ever visited CSW
	dcMatrix		_target_ageMalesVisitCSWDistribution;	// age distrib males ever visited CSW
	
	
	vector<STIname>	_target_STIprevalence_names;	// Names of the STI to be calibrated
	vector<double>	_target_STIprevalence;			// Prevalence targets of the STI to be calibrated
	dcMatrix		_target_STIprevalence_by_riskGroup; // Prevalence targets by risk group
	
	// TO DO: include in a more general variables (for all STIs)?
	// not sure, as we may rarely have STI prevalence by age...
	
	dcMatrix		_target_HIV_prev_age_f;		// HIV prevalence by age (for females)
	dcMatrix		_target_HIV_prev_age_m;		// HIV prevalence by age (for females)
	
	bool			_save_trace_files;
	
	
public:
	
	// === CONSTRUCTOR ===
	
	Simulation (){}
	
	Simulation (double horizon, double timeStep, Population P, 
				int nMC_calibration);
	
	void		set_population(Population P) {_population=P;}
	
		
	// ==== SET FUNCTIONS ====
	
	void		set_horizon(double x) {_horizon = x;}
	void		set_timeStep(double x) {_timeStep = x;}
	void		set_MC_trial_iter(unsigned int i){_MC_trial_iter=i;}
	
	void		set_newborn_timestep(unsigned long n) {_newborn_timestep=n;}
	
	void		set_intervention(vector<Intervention> v) {_intervention = v;}
	
	
	void		set_calibrationTime(vector<double> x) {_calibrationTime=x;}
	void		set_calibrationOutput(int i, vector<string> s) {_calibrationTargetFile[i] = s;}
	
	void		set_target_ageDistribution(dcMatrix M) { _target_ageDistribution = M;}
	void		set_target_ageGapDistribution(dcMatrix M) { _target_ageGapDistribution = M;}
	void		set_target_ageFirstSexDistribution_f(dcMatrix M) { _target_ageFirstSexDistribution_f = M;}
	void		set_target_ageFirstSexDistribution_m(dcMatrix M) { _target_ageFirstSexDistribution_m = M;}
	void		set_target_ageGapFirstSexSpouseDistribution_f(dcMatrix M) { _target_ageGapFirstSexSpouseDistribution_f = M;}
	void		set_target_ageGapFirstSexSpouseDistribution_m(dcMatrix M) { _target_ageGapFirstSexSpouseDistribution_m = M;}
	
	
	void		set_target_singleRatio_f(double r) {_target_singleRatio_f = r;}
	void		set_target_singleRatio_m(double r) {_target_singleRatio_m = r;}
	
	void		set_target_malesVisitCSW(double x) {_target_malesVisitCSW = x;}
	void		set_target_ageMalesVisitCSWDistribution(dcMatrix M) { _target_ageMalesVisitCSWDistribution = M;}
	
	void		set_target_nLftSexPrtnrDistribution_f(dcMatrix M) {_target_nLftSexPrtnrDistribution_f = M;}
	void		set_target_nLftSexPrtnrDistribution_m(dcMatrix M) {_target_nLftSexPrtnrDistribution_m = M;}
	
	
	void		set_target_STIprevalence_names(vector<STIname> stinames){_target_STIprevalence_names=stinames;}
	void		set_target_STIprevalence(vector<double> prev) {_target_STIprevalence=prev;}
	void		set_target_STIprevalence_by_riskGroup(dcMatrix M) {_target_STIprevalence_by_riskGroup=M;}
	
	void		set_target_HIV_prev_age_f(dcMatrix M) {_target_HIV_prev_age_f = M;}
	void		set_target_HIV_prev_age_m(dcMatrix M) {_target_HIV_prev_age_m = M;}
	
	void		set_save_trace_files(bool x) {_save_trace_files = x;}
	
	// =======================
	// =======================
	// ==== GET FUNCTIONS ====
	// =======================
	// =======================
	
	Population			get_population() {return _population;}
	Nursery				get_nursery() {return _nursery;}
	
	unsigned int		get_MC_trial_iter() {return _MC_trial_iter;}
	
	double				get_horizon(){return _horizon;}
	double				get_timeStep(){return _timeStep;}
	vector<double>		get_schedule(){return _schedule;}
	
	dcMatrix			get_STI_incidence(){return _STI_incidence;}
	dcMatrix			get_STI_prevalence(){return _STI_prevalence;}
	
	vector<double>		get_STI_prevalence_final();
	vector<double>		get_STI_cumIncidence_final();
	
	vector<STIname>		get_target_STIprevalence_names() {return _target_STIprevalence_names;}
	vector<double>		get_target_STIprevalence() {return _target_STIprevalence;}
	
	
	double				get_target_singleRatio_f() {return _target_singleRatio_f;}
	double				get_target_singleRatio_m() {return _target_singleRatio_m;}
	
	
	vector<double>			get_calibrationTime() {return _calibrationTime;}
	vector<vector<string> >	get_calibrationTargetFile() {return _calibrationTargetFile;}
	vector<vector<double> >	get_calibrationDistances() {return _calibrationDistances;}
	
	
	// FOR ALL AGE (AGE GAP) DISTRIBUTIONS:
	// 1st col = age breaks
	// 2nd col = proportion of population
	// in associated age range age[i]<prop[i]<age[i+1]
	
	dcMatrix			get_target_ageDistribution() {return _target_ageDistribution;}
	dcMatrix			get_target_ageGapDistribution() {return _target_ageGapDistribution;}
	dcMatrix			get_target_ageFirstSexDistribution_f() {return _target_ageFirstSexDistribution_f;}
	dcMatrix			get_target_ageFirstSexDistribution_m() {return _target_ageFirstSexDistribution_m;}
	dcMatrix			get_target_ageGapFirstSexSpouseDistribution_f() {return _target_ageGapFirstSexSpouseDistribution_f;}
	dcMatrix			get_target_ageGapFirstSexSpouseDistribution_m() {return _target_ageGapFirstSexSpouseDistribution_m;}
	dcMatrix			get_target_nLftSexPrtnrDistribution_f() {return _target_nLftSexPrtnrDistribution_f;}
	dcMatrix			get_target_nLftSexPrtnrDistribution_m() {return _target_nLftSexPrtnrDistribution_m;}
	double			get_target_malesVisitCSW() {return _target_malesVisitCSW;}
	dcMatrix			get_target_ageMalesVisitCSWDistribution() {return _target_ageMalesVisitCSWDistribution;}
	dcMatrix			get_target_STIprevalence_by_riskGroup() {return _target_STIprevalence_by_riskGroup;}

	dcMatrix			get_target_HIV_prev_age_f() {return _target_HIV_prev_age_f;}
	dcMatrix			get_target_HIV_prev_age_m() {return _target_HIV_prev_age_m;}
	
	dcDataFrame			get_df_sim(){return _df_sim;}
	dcDataFrame			get_df_interv(){return _df_interv;}
	
	// =======================
	// =======================
	// ==   UPDATE EVENTS   ==
	// =======================
	// =======================
	
	void			update_pregnancies(double timestep);
	void			set_pregnant(unsigned long uid);
	void			set_gestationDuration(unsigned long uid,double timestep);
	void			increment_gestationDuration(unsigned long uid,double timestep);
	void			increment_nChildBorn(unsigned long uid);

	
	// =======================
	// =======================
	// ====      STI      ====
	// =======================
	// =======================
	
	void			STI_set_initial_prevalence(string filename);
	vector<bool>	MTCT(unsigned long uid);
	double			STI_cumul_incidence(STIname s,int time_i);
	
	
	
	// =======================
	// =======================
	// ====  INTERVENTION ====
	// =======================
	// =======================
	
	void	activate_intervention(int i);
	
	
	// =======================
	// =======================
	// ====  TREATMENT    ====
	// =======================
	// =======================
	
	void	treat_indiv(unsigned long uid, STIname sti);
//	void	treat_all(STIname sti);
//	void	treat_any_proportion(STIname sti, double proportion);
//	void	treat_symptomatic(STIname sti);
	
	void	cure_indiv(unsigned long uid, STIname sti);
	void	update_cure(STIname sti);
	
	
	// =======================
	// =======================
	// ====  VACCINATION  ====
	// =======================
	// =======================
	
	void	vaccinate_indiv(unsigned long uid, STIname stiname);
	void	update_vacc(STIname stiname);

//	void	vaccinate_all(STIname stiname);
//
//	void	vaccinate_general(STIname stiname, double propPerYear);
	
	
	// =======================
	// =======================
	// ==== RUN SIMULATION ===
	// =======================
	// =======================
	
	// Run all events (births, partnerships, sex acts, deaths,...)
	// for a given time step
	void		runAllEvents_timeStep(int numTimeStep,
									  bool doSex,
									  ofstream &f, ofstream &ff,
									  bool logIndivInfo);
	
	void		runAllEvents_timeStep_obj(int numTimeStep,
										  bool doSex,
										  bool logIndivInfo);
	
	// Run all events until '_horizon'
	// returns output to be calibrated
	// One stochastic realization only
	void		runAllEvents_horizon(bool doSex,
									 bool logIndivInfo,
									 bool traceNetwork,
									 int displayProgress,
									 int iter_mc=1);
	
	void		runAllEvents_horizon_obj(bool doSex,
										 bool logIndivInfo,
										 bool traceNetwork,
										 int displayProgress);

	
	 // ??? DEPRECATED ???
	void		runSimulation(bool logIndividualsInfo); // ??? DEPRECATED ???
	

	
	// =====================
	// === MISCELLENAOUS ===
	// =====================
	
	void		save_outputs_demog(string pathFolder, unsigned int iMC=0, unsigned int idate=0);	// Save to file demographic outputs (e.g. age distribution) ; formatted for comparison w/ calibration targets
	void		save_outputs_prtnr(string pathFolder, unsigned int iMC=0, unsigned int idate=0);	// Save to file partnership outputs (e.g. age gap distribution, single ratio,...)
	void		save_outputs_sex(string pathFolder, unsigned int iMC=0, unsigned int idate=0);	// Save to file sexual behavior outputs (e.g. number of lifetime sex partners distribution)
	void		save_outputs_epi(string pathFolder, unsigned int iMC=0, unsigned int idate=0);	// Save to file epidemiological outputs (e.g. STI prevalences)

	void	save_incidence(unsigned int iMC=0);
	void	save_prevalence(unsigned int iter_mc=0);
	
	// =====================
	// ==== CALIBRATION ====
	// =====================
	
	
	void		set_calibration_schedule(string filename_calibrationTime,
										 string filename_all_calib_output,
										 string filename_weight_target);
	
	bool			isCalibrationTime(double t);
	unsigned int	whichCalibrationTime(double t);
	
	
	void		calculate_calibDistance_one(int calibtime, int i);
	void		calculate_calibDistance_all(int calibtime_idx);

	double		calibration_distance_targets();
	
	
	void		calibration_target_type_output_save(unsigned long iter_mc, int calibtime_idx);
	
	// -- Demographics
	double		calib_distance_ageDistrib();	// distance b/w current age distribution and calibration target
	double		calib_distance_ageGapDistrib();	// distance b/w current age gap distribution and calibration target

	// -- Sexual behaviour
	double		calib_distance_ageFirstSexDistrib(Gender g);	// distance b/w current age at 1st sex distribution and calibration target
	double		calib_distance_ageGapFirstSexSpouseDistrib(Gender g);// distance b/w current age gap 1st sex/1st spouse distrib and calib target
	
	double		calib_distance_singleRatio(Gender g);	// distance b/w current single ratio and calibration target

	double		calib_distance_malesVisitCSW();	// distance b/w current proportion males ever visited CSW and calibration target
	double		calib_distance_malesVisitCSW(double maxDurationSinceLastVisit);
	
	double		calib_distance_ageMalesVisitCSWDistrib();
	
	double		calib_distance_nLifeSexPrtnrDistrib(Gender g);	// number of lifetime sex partner distribution
	
	// -- STI prevalences
	double		calib_distance_STIprev();				// distance b/w current STI prevalences and calibration target
	double		calib_distance_STIprev_riskGroup();		// distance b/w current STI prevalences by risk group and calibration target

	double		calib_distance_HIV_prev_age(Gender g);			// distance b/w current HIV prevalences by age and calibration target (for both females and males)
	double		calib_distance_HIV_prev_age();			// distance b/w current HIV prevalences by age and calibration target (for both females and males)
	
	double		calib_distanceFromAllTargets();				// Distance b/w current simulation and target values
	double		calib_distanceFromAllTargets_MC(int nMC);	// Monte Carlo Mean distance from target
	
	
	// Calculate distance from _all_ targets
	// This is the summary measure assessing
	// the calibration quality
	double		calc_distanceFromAllTargets();
	

	// Set parameters
	void		calib_setParameters(vector<double> param_demographics,
										vector<double> param_partner_formation,
										vector<double> param_spousal_progression,
										vector<double> param_partner_dissolution);
	
	void		calib_setParameters_STI(vector<double> param_sexActivity,
											vector<double> param_stiFeatures_virus,
											vector<double> param_stiFeatures_bacteria,
											vector<double> param_stiFeatures_protozoa,
											vector<double> param_stiFeatures_fungus);
	

	void		calib_ageDistribution_setParameters(vector<double> param);
	double		calib_ageDistribution_get_diff_MC(int nMC, bool doSex);
	
	
	
	
	dcMatrix		TMP_LHS(vector<double>loVal, vector<double> hiVal, int nLHS, int nMC);
	
	
	// Age distributions
	double		calib_ageDistribution_get_diff(bool doSex);	// Difference between current age distribution and target
	
	
	// =======================
	// =======================
	// ==== SENSITIVITIES ====
	// =======================
	// =======================
	
	
	vector<double> sensitivity_calib_distance(vector<double> param_DMG,
											  vector<double> param_FORM,
											  vector<double> param_SPOUSAL,
											  vector<double> param_DISSOL,
											  double relativeBump,int nMC);
	
	/*DELETE?*/vector<double>	TMPsensi(vector<double> param, double relativeBump, int nMC);
	
	
		
	
	
	
	// OOL-GSL wrappers
	double			calib_ageDistribution_F(const gsl_vector *v, void *params);
	void			calib_ageDistribution_Gradient(const gsl_vector *v, void *params, gsl_vector *G);
	void			calib_ageDistribution_F_Gradient(const gsl_vector *v, void *params, double *f,gsl_vector *G);
	vector<double>	calib_ageDistribution_OOL(vector<double> v0,
											  vector<double> LowerBound,
											  vector<double> UpperBound,
											  int maxIterations,
											  bool displayInfo);
	
	static double   calib_ageDistribution_F_wrapper(const gsl_vector *v, void *params)
	{return static_cast<Simulation*>(params)->calib_ageDistribution_F(v,params);}
	static void     calib_ageDistribution_Gradient_wrapper(const gsl_vector *v, void *params, gsl_vector *G)
	{return static_cast<Simulation*>(params)->calib_ageDistribution_Gradient(v,params,G);}
	static void     calib_ageDistribution_F_Gradient_wrapper(const gsl_vector *v, void *params, double* f, gsl_vector *G)
	{return static_cast<Simulation*>(params)->calib_ageDistribution_F_Gradient(v,params,f,G);}
	
	

	
	
	// ==== Other =====
	
	
	void		displayInfo();
	
	void		longitudinal_Monitor(ofstream file, double t);	// Follow all individuals through a simulation
	
	
	
	
	/* ***************************************************************************/
	
	// *** EXAMPLE TEMPLATE: KEEP FOR FUTURE REFERENCE !!! *****
    // GSL WRAPPERS: Solve a purely technical issue of using C libraries in C++
    static double   example_fct_to_minimize_wrapper(const gsl_vector *v, void *params){ 
        return static_cast<Simulation*>(params)->example_fct_to_minimize(v,params);}
	
	static void     example_fct_to_minimize_Gradient_wrapper(const gsl_vector *v, void *params, gsl_vector *G){
        return static_cast<Simulation*>(params)->example_fct_to_minimize_Gradient(v,params,G);}	
	
	static void     example_fct_to_minimize_FctGradient_wrapper(const gsl_vector *v, void *params, double* f, gsl_vector *G){
        return static_cast<Simulation*>(params)->example_fct_to_minimize_FctGradient(v,params,f,G);}
	
	// GSL OOL Optimization
	double		example_fct_to_minimize(const gsl_vector *v, void *params); 
	void		example_fct_to_minimize_Gradient(const gsl_vector *v, void *params, gsl_vector *G);	
	void		example_fct_to_minimize_FctGradient(const gsl_vector *v, void *params, double *f,gsl_vector *G);
	
	double		dummy_c() 
	{
		double res = pow(_population.get_formation_MaxRate()-0.01234,2.0);
		res += pow(_population.get_dissolution_MaxRate()-9.876,2.0);
		return res;
	}
	
	vector<double>	calib_minimization_example(vector<double> v0,
											   vector<double> LowerBound,
											   vector<double> UpperBound,     
											   int maxIterations);
	// *** END EXAMPLE TEMPLATE ***
	
	/* ***************************************************************************/
	
};

string calibration_file_to_type(string filename);



#endif