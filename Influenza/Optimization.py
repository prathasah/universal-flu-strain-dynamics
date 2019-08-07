import numpy as np
np.warnings.filterwarnings('ignore')
import Simulation
import Parameters
from ages import ages
import doses_distributed as dd
import vaccination_coverage as vc
from random import shuffle
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


        self.PVUniversal= None


    def solve(self, PVPWVals):
        # Only update for new PVPWVals
	self.s = Simulation.run_Simulation(season= self.season, proportion_universalVaccine_doses = self.proportional_universal, paramValues = {"PVuniversal": PVPWVals}, index=self.index, optimization = True)
	   
	 	
    def evaluateObjective(self, PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	if (PVPWVals <0).any(): return np.inf
	self.solve(PVPWVals)
	if (self.s.SUL <0).any() or (self.s.SUH <0).any() or (self.s.STL <0).any() or (self.s.STH <0).any() or (self.s.SNL <0).any() or (self.s.SNH <0).any() : return np.inf
	if np.isnan(PVPWVals).any(): return np.inf
	print ("--->"), PVPWVals[:5], PVPWVals.sum(), self.s.totalHospitalizations/1e3
	return getattr(self.s, self.objectiveMap[self.objective])
    
    
    def evaluateDetailedObjective(self, PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	self.solve(PVPWVals)
    
	return getattr(self.s, self.DetailedobjectiveMap[self.objective[5:]])

    
    def totalVacsUsed(self):
	
	return lambda x:  1 - sum(x)

    
    def dosesL_less_than_popsize(self, i):
	
	##proportion vaccinated (both seasonal and universal vaccine) should be less than one
	return lambda PVPWVals: 1.0 - ((PVPWVals[i]*self.initial_universal_vacDoses)/self.lowrisk_population[i]) - self.seasonal_lowrisk_vaxcoverage[i]
    
    def dosesH_less_than_popsize(self, i):
	
	##proportion vaccinated should be less than one
	
	return lambda PVPWVals: 1.0 - ((PVPWVals[i]*self.initial_universal_vacDoses)/self.highrisk_population[i-self.proportionVaccinatedLength]) - self.seasonal_lowrisk_vaxcoverage[i-self.proportionVaccinatedLength]

    def lowerCondition(self, i):
	#the min values should be greater than zero
	return lambda PVPWVals: PVPWVals[i]

    def upperCondition(self, i):
	#1 - PVPWal should be greater than zero
        return lambda PVPWVals: 1.0 - PVPWVals[i]
    

    def optimize(self):
        from scipy.optimize import fmin_cobyla
	from scipy.optimize import minimize
	
	self.proportionVaccinatedLength = len(ages)-1
	
	conds = [self.totalVacsUsed()]

	conds.extend([self.lowerCondition(i) for i in range(self.proportionVaccinatedLength)])
	
        conds.extend([self.upperCondition(i) for i in range(self.proportionVaccinatedLength)])
	
	conds.extend([self.dosesL_less_than_popsize(i) for i in range(self.proportionVaccinatedLength)])
	conds.extend([self.dosesH_less_than_popsize(i) for i in range(self.proportionVaccinatedLength, self.proportionVaccinatedLength*2)])
	
        minObjective = None

        for i in range(self.optimRuns):

	    ## proportion of people in ageclass a that are vaccinated with universal vaccines (2).
	    ## initialize  with same proportion of universal doses to all age classes
            
	    valid_condition = False
	    while not valid_condition:
		condition1 = False
		while not condition1:
		    PV0_raw = np.random.rand(self.proportionVaccinatedLength*2)
		    # normalize so that they sum to 1
		    PV0_raw = PV0_raw/PV0_raw.sum()
		    #PV0_raw = sorted(PV0_raw, reverse = True)
		    PV0_raw_lowrisk = PV0_raw[:self.proportionVaccinatedLength]
		    PV0_raw_highrisk = PV0_raw[self.proportionVaccinatedLength:]
		    shuffle(PV0_raw_lowrisk)
		    PV0_raw_highrisk = sorted(PV0_raw_highrisk)
		    PV0 = np.array(list(PV0_raw_lowrisk) + list(PV0_raw_highrisk))
		    
		    
		    PV0_lowrisk = 1.0 - self.seasonal_lowrisk_vaxcoverage - np.array([((PV0[i]*self.initial_universal_vacDoses)/self.lowrisk_population[i])  for i in range(self.proportionVaccinatedLength)])
		    PV0_highrisk = 1.0 - self.seasonal_highrisk_vaxcoverage - np.array([((PV0[i]*self.initial_universal_vacDoses)/self.highrisk_population[i-self.proportionVaccinatedLength])  for i in range(self.proportionVaccinatedLength, self.proportionVaccinatedLength*2)])
		    condition1 = (PV0_lowrisk >=0).all() and  (PV0_highrisk >=0).all() 
		self.solve(PV0)
		condition2 = (self.s.SUL >=0).all() and (self.s.SUH >=0).all() and (self.s.STL >=0).all() and (self.s.STH >=0).all() and (self.s.SNL >=0).all() and (self.s.SNH >=0).all()
		
		
		valid_condition = condition1 and condition2
		print PV0[:5], condition1, condition2
		print ("all conds"), (self.s.SUL >=0).all(), (self.s.SUH >=0).all(),(self.s.STL >=0).all(), (self.s.STH >=0).all(),  (self.s.SNL >=0).all(), (self.s.SNH >=0).all()
		

	  

	    
            PVPWValsOpt = fmin_cobyla(self.evaluateObjective,
                                      PV0,
                                     conds,
                                      maxfun = 10000,
                                      rhobeg = 0.005,
	   			      rhoend=0.000001,
                                     disp = 0)
	    
	    
	   # minimizer_kwargs = dict(method="fmin_cobyla",constraints= conds)
	    #PVPWValsOpt =  basinhopping(self.evaluateObjective, PV0, minimizer_kwargs= minimizer_kwargs, niter=10)

	    print ("PVPWals"), PVPWValsOpt, self.evaluateObjective(PVPWValsOpt)
	    print ("PV0"), PV0, PV0
	    if (minObjective == None) \
                    or (self.evaluateObjective(PVPWValsOpt) < minObjective):
                
                minObjective = self.evaluateObjective(PVPWValsOpt)
                self.PVBest = PVPWValsOpt
	


    def optimization_output(self):

	self.solve(self.PVBest)
	seasonal_vacDoses, universal_vacDoses, total_doses = self.s.doses_used()
	seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise = self.s.doses_used_agewise()
	return list(self.PVBest), seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, list(self.evaluateDetailedObjective(self.PVBest))



	
