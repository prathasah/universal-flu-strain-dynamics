import numpy as np
np.warnings.filterwarnings('ignore')
import Simulation
from ages import ages

## PVPWVals: is a list with zero and ones for each age class. 0 = vaccinate all in age class a with seasonal
## 1: vaccinate all in age class a with universal.

class optimization:
    objectiveMap = {'totalInfections': 'totalInfections',
                    'totalDeaths': 'totalDeaths',
		     'totalBurden': 'totalDALY',
		     'totalCost': 'totalCosts',
                    'totalHospitalizations': 'totalHospitalizations',
		    'Infections': 'infections',
                    'Deaths': 'deaths',
		     'Burden': 'DALY',
		     'Cost': 'Costs',
                    'Hospitalizations': 'hospitalizations'}

    def __init__(self, objective = None, optimRuns = 1, season = None, proportion_universalVaccine_doses = 0,  paramValues = {}, index = None):

        self.optimRuns = optimRuns

        self.objective = objective
        self.proportional_universal = proportion_universalVaccine_doses
	self.season = season
	self.index = index

        #self.s = Simulation.run_Simulation(season= self.season, proportion_universalVaccine_doses = self.proportional_universal, paramValues = {}, index=self.index)
	

        self.PVUniversal= None

	self.totalvacsUsed = 0
	self.UniversalvacsUsed = 0
	self.SeasonalsvacsUsed = 0 

    def solve(self, PVPWVals):
        # Only update for new PVPWVals
        if np.any(PVPWVals != self.PVUniversal):
	    self.s = Simulation.run_Simulation(season= self.season, proportion_universalVaccine_doses = self.proportional_universal, paramValues = {"PVuniversal": PVPWVals}, index=self.index, optimization = True)
	    seasonal_vacDoses, universal_vacDoses, total_doses = self.s.doses_used()
	    self.totalvacsUsed = total_doses
	    self.UniversalvacsUsed = universal_vacDoses
	    self.SeasonalsvacsUsed = seasonal_vacDoses

    def evaluateObjective(self, PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	if (PVPWVals <0).any(): return np.inf
	self.solve(PVPWVals)
    
	return getattr(self.s, self.objectiveMap[self.objective])
    

    #def totalVacsConditions(self, PVPWVals):
    #	self.solve(PVPWVals)
    #	return (self.vacNumbers - self.vacsUsed)

    #def totalVacsCondition(self, i):
    #	return lambda PVPWVals: self.totalVacsConditions(PVPWVals)[i]
    
    def totalVacsUsed(self):
	return lambda x:  1 - sum(x)
    

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
	

	
        minObjective = None

        for i in range(self.optimRuns):

	    ## proportion of people in ageclass a that are vaccinated with universal vaccines (2).
	    ## initialize  with same proportion of universal doses to all age classes
            PV0 = np.full(self.proportionVaccinatedLength*2, 1./(self.proportionVaccinatedLength*2))
	   


            PVPWValsOpt = fmin_cobyla(self.evaluateObjective,
                                      PV0,
                                     conds,
                                      maxfun = 10000,
                                      rhobeg = 0.005,
	   			      rhoend=0.000001,
                                     disp = 0)

	    
	    if (minObjective == None) \
                    or (self.evaluateObjective(PVPWValsOpt) < minObjective):
                
                minObjective = self.evaluateObjective(PVPWValsOpt)
                self.PVBest = PVPWValsOpt
		
        
	self.solve(self.PVBest)
	bestPV = [round(num,2) for num in list(self.PVBest)]
	#print ("minimum objective"), minObjective, bestPV, sum(PVPWValsOpt), round(min(PVPWValsOpt),3), self.totalvacsUsed, self.UniversalvacsUsed, self.SeasonalsvacsUsed


    def short_optimization_output(self):

	return self.simulatedR0, list(self.vacsUsed)[0], self.evaluateObjective(self.PVBest), list(self.PVBest), list(self.Vaccinatedtotal), self.evaluateDetailedObjective(self.PVBest)


    def low_vaccine_optimization(self):
	return sum(list(self.infections)), sum(list(self.hospitalizations)), sum(list(self.deaths))
	
