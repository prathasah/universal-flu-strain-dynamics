#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 

def run_efficacy_simulation(vacEfficacy, vacDoses, vaccination_type, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	if vaccination_type == "random": vaccineCoverage = s.comupte_random_vaccination(vacDoses)
	elif vaccination_type =="typical": vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy, vacDoses)
	infections_H1, infections_H3, infections_B, hospitalizations_H1, hospitalizations_H3, hospitalizations_B, deaths_H1, deaths_H3, deaths_B =  s.strain_output()
	return infections_H1, infections_H3, infections_B, hospitalizations_H1, hospitalizations_H3, hospitalizations_B, deaths_H1, deaths_H3, deaths_B 

######################################################################33
	