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
	
		
		
######################################################################################		
if __name__ == "__main__":
	
	for num in xrange(1):
		optimize_universal_vaccine_distribution('2011-12', 0.5, 'totalInfections', num)
    
    
	
		
	
		
