#!/usr/bin/python
import sys
sys.path.insert(0, r'../Influenza')
import numpy as np
import Simulation
import Optimization
import random
import subprocess


def optimize_universal_vaccine_distribution(season, proportion_universal_vacDoses, objective, index):
	o = Optimization.optimization(season = season, proportion_universalVaccine_doses = proportion_universal_vacDoses, objective = objective, index = index)
	o.optimize()
	
	seasonal_vacDoses, universal_vacDoses, total_doses,best_SV_agewise_doses, best_UV_agewise_doses, SV_coverage, UV_coverage,total_outcome, outcome_agewise = o.optimization_output()
	
	return seasonal_vacDoses, universal_vacDoses, total_doses,best_SV_agewise_doses, best_UV_agewise_doses, SV_coverage, UV_coverage,total_outcome, list(outcome_agewise)
	
		
######################################################################################		
if __name__ == "__main__":
	
	for num in xrange(830,831):
		
		seasonal_vacDoses, universal_vacDoses, total_doses,best_SV_agewise_doses, best_UV_agewise_doses, SV_coverage, UV_coverage,total_outcome, outcome_agewise = optimize_universal_vaccine_distribution('2012-13', 0.25, 'totalDeaths', num)
		print seasonal_vacDoses, universal_vacDoses, total_doses,total_outcome
    
    
	
		
	
		
