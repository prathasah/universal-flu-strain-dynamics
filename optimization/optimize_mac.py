#!/usr/bin/python
import sys
sys.path.insert(0, r'../Influenza')
import numpy as np
import Simulation
import Optimization
import random
import subprocess


def optimize_universal_vaccine_distribution(season, proportion_universal_vacDoses, objective, index):
	#s = Simulation.run_Simulation(season = season, proportion_universalVaccine_doses = proportion_universal_vacDoses)
	o = Optimization.optimization(season = season, proportion_universalVaccine_doses = proportion_universal_vacDoses, objective = objective, index = index)
	o.optimize()
	
	PVbest, seasonal_vacDoses, universal_vacDoses, total_doses,seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise, outcome = o.optimization_output()
	
	return PVbest, seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise,list(outcome)
	
	
	
		
		
######################################################################################		
if __name__ == "__main__":
	
	for num in xrange(492,493):
		PVbest, seasonal_vacDoses, universal_vacDoses, total_doses, seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise,outcome = optimize_universal_vaccine_distribution('2016-17', 0.75, 'totalHospitalizations', num)
		#print PVbest, seasonal_vacDoses, universal_vacDoses, total_doses
    
    
	
		
	
		
