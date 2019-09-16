import numpy as np
np.warnings.filterwarnings('ignore')
import Simulation
import Parameters
from ages import ages
import doses_distributed as dd
import vaccination_coverage as vc
from random import uniform
from scipy.optimize import basinhopping

## PVPWVals: is a list with zero and ones for each age class. 0 = vaccinate all in age class a with seasonal
## 1: vaccinate all in age class a with universal.

class optimization:
    DetailedobjectiveMap = {'Infections': 'infections',
                    'Deaths': 'deaths',
		     'Burden': 'DALY',
		     'Cost': 'Costs',
                    'Hospitalizations': 'hospitalizations'}
    
    objectiveMap = {'totalInfections': 'totalInfections',
                    'totalDeaths': 'totalDeaths',
		     'totalBurden': 'totalDALY',
		     'totalCost': 'totalCosts',
                    'totalHospitalizations': 'totalHospitalizations'}

    def __init__(self, objective = None, optimRuns = 1, season = None, proportion_universalVaccine_doses = 0,  paramValues = {}, index = None):

        self.optimRuns = optimRuns
	

        self.objective = objective
        self.proportional_universal = proportion_universalVaccine_doses
	self.season = season
	self.index = index
	#optimization is FAlse just to extract defaul Parameter values
	self.parameters = Parameters.Parameters(season, index, calibration = False, optimization=False)
	self.lowrisk_population = self.parameters.population_lowrisk[1:]
	self.highrisk_population = self.parameters.population_highrisk[1:]
	empirical_vax_coverage_lowrisk, empirical_vax_coverage_highrisk  = vc.age_specific_vaccination_coverage(self.season, self.parameters.population, self.parameters.population_lowrisk, self.parameters.population_highrisk)
	self.initial_vacDoses = dd.doses_applied_before_start_season(self.season)
	self.initial_universal_vacDoses = self.initial_vacDoses * self.proportional_universal
	self.initial_seasonal_vacDoses = self.initial_vacDoses - self.initial_universal_vacDoses
	total_doses_raw = (empirical_vax_coverage_lowrisk * self.parameters.population_lowrisk +  empirical_vax_coverage_highrisk * self.parameters.population_highrisk).sum()
	self.seasonal_lowrisk_vaxcoverage = ((self.initial_seasonal_vacDoses * empirical_vax_coverage_lowrisk)/total_doses_raw)[1:]
	self.seasonal_highrisk_vaxcoverage = ((self.initial_seasonal_vacDoses * empirical_vax_coverage_highrisk)/total_doses_raw)[1:]
	
	###########
	self.total_vacDoses = dd.total_seasonal_doses_for_season(self.season)
	self.universal_total_vacDoses = self.total_vacDoses * self.proportional_universal
	self.seasonal_total_vacDoses = self.total_vacDoses - self.universal_total_vacDoses
	self.seasonal_overall_lowrisk_vaxcoverage = ((self.seasonal_total_vacDoses * empirical_vax_coverage_lowrisk)/total_doses_raw)[1:]
	self.seasonal_overall_highrisk_vaxcoverage = ((self.seasonal_total_vacDoses * empirical_vax_coverage_highrisk)/total_doses_raw)[1:]
	self.seasonal_overall_coverage = list(self.seasonal_overall_lowrisk_vaxcoverage) + list(self.seasonal_overall_highrisk_vaxcoverage)
	


    def solve(self, universal_PVPWVals):
        # Only update for new PVPWVals
	self.s = Simulation.run_Simulation(season= self.season, proportion_universalVaccine_doses = self.proportional_universal, paramValues = {"PVuniversal": universal_PVPWVals}, index=self.index, optimization = True)
	   
	 	
    def evaluateObjective(self, universal_PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	if (universal_PVPWVals <0).any(): return np.inf
	self.solve(universal_PVPWVals)
	if (self.s.SUL <0).any() or (self.s.SUH <0).any() or (self.s.STL <0).any() or (self.s.STH <0).any() or (self.s.SNL <0).any() or (self.s.SNH <0).any() : return np.inf
	if np.isnan(universal_PVPWVals).any(): return np.inf
	#print ("--->"),  [round(num,2) for num in list(self.s.vaccine_doses_NL[1:]/ self.parameters.population_lowrisk[1:])], self.s.totalHospitalizations/1e3
	return getattr(self.s, self.objectiveMap[self.objective])
    
    
    def evaluateDetailedObjective(self, universal_PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	self.solve(universal_PVPWVals)
    
	return getattr(self.s, self.DetailedobjectiveMap[self.objective[5:]])

    def bounds(self):
	return [(0, 1.0 - self.seasonal_overall_coverage[i]) for i in range(self.VaccinatedLength)]
    
    def acceptance_criteria(self, **kwargs):
	print("in accept test")
	x_new = kwargs['x_new']
	test1 = [x <= 1 and x >=0 for x in x_new]
	test2 = [x_new[i] + self.seasonal_overall_coverage[i] <=1 for i in xrange(len(x_new))]
	return all(test1) and all(test2) 
    ###################################################################
    def print_fun(self, x, f, accepted):
	    print("x value accepted at"), list(x),f/1000., accepted
     ###################################################################
    def generate_PV0(self):
	
	trial = 0
	pseudo_coverage_limit = 1
	valid_condition = False
	while not valid_condition:
	    valid_doses_condition = False

	    while not valid_doses_condition:
		##generate 16+15 random numbers (one less than required)
		trial+=1
		PV0 = np.array([uniform(0, max(pseudo_coverage_limit-num,0)) for num in self.seasonal_overall_coverage])[:-1]
		doses_used = (PV0[:16] * self.parameters.population_lowrisk[1:]).sum() + (PV0[16:]*self.parameters.population_highrisk[1:-1]).sum()
		valid_doses_condition = doses_used < self.universal_total_vacDoses
		#print ("--->"), [round(num,2) for num in list(PV0)], pseudo_coverage_limit, trial
		if trial > 100:
		    pseudo_coverage_limit=pseudo_coverage_limit-0.01
		    trial=0

		
	    doses_left = self.universal_total_vacDoses - doses_used
	    if (self.seasonal_overall_coverage[-1 ] + (doses_left/ self.parameters.population_highrisk[-1])) < 1:
		condition1 =True
		PV0 = np.array(list(PV0) + [doses_left/ self.parameters.population_highrisk[-1]])

	    
		self.solve(PV0)
		
		valid_condition = condition1 and (self.s.SUL >=0).all() and (self.s.SUH >=0).all() and (self.s.STL >=0).all() and (self.s.STH >=0).all() and (self.s.SNL >=0).all() and (self.s.SNH >=0).all()
		print ("PV0 --->"),self.season, self.index, self.proportional_universal, valid_condition
		if not valid_condition: print ("check... "), self.season, self.index, self.proportional_universal, condition1 , (self.s.SUL >=0).all() , (self.s.SUH >=0).all() , (self.s.STL >=0).all() , (self.s.STH >=0).all() , (self.s.SNL >=0).all() , (self.s.SNH >=0).all()
		
	return PV0
	        
	    

    def optimize(self):
        from scipy.optimize import fmin_cobyla
	from scipy.optimize import minimize
	
	self.VaccinatedLength = (len(ages)-1)*2
	
        PV0 = self.generate_PV0()

	result = basinhopping(self.evaluateObjective, PV0, T=5, niter=100, stepsize=0.5, accept_test=self.acceptance_criteria, callback=self.print_fun, niter_success=10, disp=False, minimizer_kwargs= {"method":"L-BFGS-B", "bounds": self.bounds(), "tol": 0.5,  'options': {'eps':0.05,'maxfun': 100}})
	PVPWValsOpt=  result['x']

	self.PVBest = PVPWValsOpt
	


    def optimization_output(self):

	self.solve(self.PVBest)
	seasonal_vacDoses, universal_vacDoses, total_doses = self.s.doses_used()
	seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise = self.s.doses_used_agewise()
	return list(self.PVBest), seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, list(self.evaluateDetailedObjective(self.PVBest))



	
