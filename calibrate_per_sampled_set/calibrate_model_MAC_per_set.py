#!/usr/bin/python
import sys
sys.path.insert(0, r'../Influenza')
import numpy as np
np.warnings.filterwarnings('ignore')
import pandas as pd
import Simulation
import csv
import itertools
import scipy.optimize
from scipy.optimize import fmin_cobyla
from scipy.optimize import basinhopping
import time
from copy import deepcopy
##########################################
def run_efficacy_simulation(season, sub_iter):
	s = Simulation.run_Simulation(season = season, index = sub_iter)
	seasonal_vacDoses, universal_vacDoses, total_doses = s.doses_used()
	#print ("check seasonal, universal doses "), seasonal_vacDoses, universal_vacDoses
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH
##########################################
result_list = []
def log_results(result):
    result_list.append(result)
################################################################
def universal_worker(input_pair):
    function, args0 = input_pair
    args = list(args0[0]) + [args0[1]]
    return function(*args)

def pool_args(function, *args):
    return zip(itertools.repeat(function), zip(itertools.repeat(args[:-1]) ,args[-1]))
   
###########################################################################

def run_simulation_incidence(paramList, year, sub_iter):
    
    vacEfficacy_universal = 0
    proportional_vacDoses_universal = 0
    start = time.time()
    ##set 1 for all hospitalization and death efficacy and scaling parameters. Calibrate these in step#2
    s = Simulation.run_Simulation(season = year, proportion_universalVaccine_doses = proportional_vacDoses_universal,paramValues = {"vacEfficacy_universal":vacEfficacy_universal, "betaList": paramList[0:3], "susceptibility_H1": paramList[3:7], "susceptibility_H3": paramList[7:11], "susceptibility_B": paramList[11:15],
 "vac_eff_hospitalization": 1, "vac_eff_mortality": 1, "prob_hosp_scaling": 1, "prob_death_scaling": 1}, index = sub_iter, calibration = True)
    seasonal_vacDoses, universal_vacDoses, total_doses = s.doses_used()
    print ("doses used  ===="), sub_iter, seasonal_vacDoses/1e6, universal_vacDoses/1e6, total_doses/1e6
    print ("time ==="), time.time() - start
    
    incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	
    incidence_raw = sum(incidenceL)+ sum(incidenceH)
    incidence = sum(incidenceL)/1e6 + sum(incidenceH)/1e6

    incidence_H1_0 = (100*sum(list(infections_H1)[0:2]))/(1.*incidence_raw)
    incidence_H1_5 = (100*sum(list(infections_H1)[2:6]))/(1.*incidence_raw)
    incidence_H1_25 = (100*sum(list(infections_H1)[6:14]))/(1.*incidence_raw)
    incidence_H1_65 = (100*sum(list(infections_H1)[14:]))/(1.*incidence_raw)
    
    incidence_H3_0 = (100*sum(list(infections_H3)[0:2]))/(1.*incidence_raw)
    incidence_H3_5 = (100*sum(list(infections_H3)[2:6]))/(1.*incidence_raw)
    incidence_H3_25 = (100*sum(list(infections_H3)[6:14]))/(1.*incidence_raw)
    incidence_H3_65 = (100*sum(list(infections_H3)[14:]))/(1.*incidence_raw)
    
    incidence_B_0 = (100*sum(list(infections_B)[0:2]))/(1.*incidence_raw)
    incidence_B_5 = (100*sum(list(infections_B)[2:6]))/(1.*incidence_raw)
    incidence_B_25 = (100*sum(list(infections_B)[6:14]))/(1.*incidence_raw)
    incidence_B_65 = (100*sum(list(infections_B)[14:]))/(1.*incidence_raw)
		
    return incidence, incidence_H1_0, incidence_H1_5, incidence_H1_25, incidence_H1_65, incidence_H3_0, incidence_H3_5, incidence_H3_25, incidence_H3_65, incidence_B_0, incidence_B_5, incidence_B_25, incidence_B_65

###########################################################################


def run_simulation_burden(paramList, betalist, susc_list, year,sub_iter):
	
	
    vacEfficacy_universal = 0
    proportional_vacDoses_universal = 0

    s = Simulation.run_Simulation(season = year, proportion_universalVaccine_doses = proportional_vacDoses_universal, paramValues = {"vacEfficacy_universal":vacEfficacy_universal,  "betaList": betalist, "susceptibility_H1": susc_list[0:4], "susceptibility_H3": susc_list[4:8], "susceptibility_B": susc_list[8:12], "vac_eff_hospitalization": paramList[0], "vac_eff_mortality": paramList[1], "prob_hosp_scaling": paramList[2], "prob_death_scaling": paramList[3]}, index = sub_iter, calibration = True)
    #total_seasonal_vacDoses, total_universal_vacDoses = s.simulateWithVaccine(year, vacEfficacy_universal, proportional_vacDoses_universal)
  
    incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
    
    hosp  = sum(hospitalizationsL)/1e3 + sum(hospitalizationsH)/1e3
    deaths = sum(deathsL)/1e3 + sum(deathsH)/1e3
	
	
    return hosp, deaths

#####################################################################
def optimize(cons1, cons2, data_list_incidence,  year, sub_iter):
    
    [iteration, data_infections,  data_hospitalization, data_mortality, data_incidence_H1_0, data_incidence_H1_5, data_incidence_H1_25, data_incidence_H1_65, data_incidence_H3_0, data_incidence_H3_5, data_incidence_H3_25, data_incidence_H3_65, data_incidence_B_0, data_incidence_B_5, data_incidence_B_25, data_incidence_B_65] = data_list_incidence
    
    beta_H1 = []
    beta_H3 = []
    beta_B=[]
    H1_0 = []
    H1_5 = []
    H1_25 = []
    H1_65 = []
    H3_0 = []
    H3_5 = []
    H3_25 = []
    H3_65 = []
    B_0= []
    B_5 = []
    B_25 = []
    B_65 = []
    
    #for sub_iter in range(sub_range[0], sub_range[1]):
    time1 = time.time()
    min_deviation = 1e7
    beta_opt = None
    
    for trial in xrange(2):
		bounds = [(0.00001,1)]*15
		beta_init = list(np.random.uniform(low=0.0015, high=0.0025, size=3))
		susc_init = list(np.random.uniform(low=0.4, high=1, size=12))
		BL0_incidence = beta_init +susc_init
		beta_opt_raw =  scipy.optimize.minimize(evaluateObjective_incidence, BL0_incidence, args=(data_list_incidence, year, sub_iter), method='TNC', jac=None, bounds=bounds, constraints = {cons1, cons2}, tol=None, callback=None, options={'disp': None, 'maxls': 20, 'iprint': -1, 'gtol': 1e-05, 'eps': 1e-08, 'maxiter': 1500000, 'ftol': 2.220446049250313e-09, 'maxcor': 10, 'maxfun': 15000})
		#beta_opt_raw = fmin_cobyla(evaluateObjective_incidence, BL0_incidence,  [cons1, cons2], args=  (data_incidence, sub_iter),  consargs=(),maxfun = 10e10, rhobeg = 0.001,rhoend=1e-7,catol= 0, disp = 0)
		deviation = beta_opt_raw['fun']
		print ("optimized = "), sub_iter, deviation, trial
		if deviation < min_deviation:
			min_deviation = deepcopy(deviation)
			beta_opt = [num for num in beta_opt_raw['x']]
			
		
	
    print ("final =="), sub_iter, min_deviation
    beta_H1.append(beta_opt[0])
    beta_H3.append(beta_opt[1])
    beta_B.append(beta_opt[2])
    H1_0.append(beta_opt[3])
    H1_5.append(beta_opt[4])
    H1_25.append(beta_opt[5])
    H1_65.append(beta_opt[6])
    H3_0.append(beta_opt[7])
    H3_5.append(beta_opt[8])
    H3_25.append(beta_opt[9])
    H3_65.append(beta_opt[10])
    B_0.append(beta_opt[11])
    B_5.append(beta_opt[12])
    B_25.append(beta_opt[13])
    B_65.append(beta_opt[14])
    
    
     
    eff_hosp =[]
    eff_death = []
    hosp_scaling  =[]
    death_scaling = []
   
    bounds = [(0.00000001,1)]*2 + [(0.00000001,1)]*2
    min_eff_deviation = 1e10
    for trial in xrange(1):
	BL0_burden = list(np.random.uniform(low=0.00000001, high=1, size=2)) + list(np.random.uniform(low=0.00000001, high=0.1, size=2))
	eff_opt_raw =  scipy.optimize.minimize(evaluateObjective_burden, BL0_burden, args=(beta_opt[0:3], beta_opt[3:], data_hospitalization[sub_iter], data_mortality[sub_iter], year, sub_iter), method='TNC', jac=None, bounds=bounds, constraints = {cons1, cons2}, tol=None, callback=None, options={'disp': None, 'maxls': 20, 'iprint': -1, 'gtol': 1e-08, 'eps': 1e-08, 'maxiter': 150000000, 'ftol': 2.220446049250313e-09, 'maxcor': 10, 'maxfun': 1500000})
	eff_deviation = eff_opt_raw['fun']
	#print ("optimized  efficacy = "), sub_iter, eff_deviation, trial
	if eff_deviation < min_eff_deviation:
	    min_eff_deviation = deepcopy(eff_deviation)
	    eff_opt = [num for num in eff_opt_raw['x']]
	    
	    
    
    print ("final eff =="), sub_iter, min_eff_deviation, eff_opt
    
    eff_hosp.append(eff_opt[0])
    eff_death.append(eff_opt[1])
    hosp_scaling.append(eff_opt[2])
    death_scaling.append(eff_opt[3])
    
    vax_hosp, vax_death = run_simulation_burden(eff_opt, beta_opt[0:3], beta_opt[3:], year,  sub_iter)
    print ("data vs calibration = "), sub_iter, data_hospitalization[sub_iter], data_mortality[sub_iter], vax_hosp, vax_death
    
     
    return beta_H1, beta_H3, beta_B, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65, eff_hosp, eff_death, hosp_scaling, death_scaling

###########################################################################
def evaluateObjective_incidence(paramList, data_list_incidence, year, sub_iter):
	
	
	time2 = time.time()
	[iteration, data_infections, data_hospitalization, data_mortality, data_incidence_H1_0, data_incidence_H1_5, data_incidence_H1_25, data_incidence_H1_65, data_incidence_H3_0, data_incidence_H3_5, data_incidence_H3_25, data_incidence_H3_65, data_incidence_B_0, data_incidence_B_5, data_incidence_B_25, data_incidence_B_65] = data_list_incidence
	
	vax_incidence, vax_incidence_H1_0, vax_incidence_H1_5, vax_incidence_H1_25, vax_incidence_H1_65, vax_incidence_H3_0, vax_incidence_H3_5, vax_incidence_H3_25, vax_incidence_H3_65, vax_incidence_B_0, vax_incidence_B_5, vax_incidence_B_25, vax_incidence_B_65 =run_simulation_incidence(paramList,year, sub_iter)
		
	
	## only incidence, H1 and H3 perc and not B because perc B = 100 - (perc_H1+H3)
	deviation = ((vax_incidence - data_infections[sub_iter])**2)**2 + \
	(vax_incidence_H1_0 - data_incidence_H1_0[sub_iter])**2  + \
	(vax_incidence_H1_5 - data_incidence_H1_5[sub_iter])**2 + \
	(vax_incidence_H1_25 - data_incidence_H1_25[sub_iter])**2 +\
	(vax_incidence_H1_65 - data_incidence_H1_65[sub_iter])**2 +\
	(vax_incidence_H3_0 - data_incidence_H3_0[sub_iter])**2  +\
	(vax_incidence_H3_5 - data_incidence_H3_5[sub_iter])**2 +\
	(vax_incidence_H3_25 - data_incidence_H3_25[sub_iter])**2 +\
	(vax_incidence_H3_65 - data_incidence_H3_65[sub_iter])**2 +\
	(vax_incidence_B_0 - data_incidence_B_0[sub_iter])**2  +\
	(vax_incidence_B_5 - data_incidence_B_5[sub_iter])**2 +\
	(vax_incidence_B_25 - data_incidence_B_25[sub_iter])**2 +\
	(vax_incidence_B_65 - data_incidence_B_65[sub_iter])**2 
	
	#print ("betalist"),  deviation, vax_incidence_H1_0 , data_incidence_H1_0[sub_iter]
	#print ("param fit --"), np.mean(vax_incidence), np.mean(vax_perc_H1), np.mean(vax_perc_H3), 100 - np.mean(vax_perc_H1) - np.mean(vax_perc_H3), np.mean(vax_incidence_H1_5), np.mean(vax_incidence_H3_5), np.mean(vax_incidence_B_5)
	#print ("data=="), data
	return deviation
    
###########################################################################
def evaluateObjective_burden(paramList, beta_list, susc_list,  data_hosp, data_mort, year, sub_iter):
    
	vax_hosp, vax_mortality = run_simulation_burden(paramList, beta_list, susc_list, year, sub_iter)
	
	deviation = (vax_hosp - data_hosp)**2 + (vax_mortality - data_mort)**2
	
	return deviation

###################################################

def cons1(x):
    return 1 - max(x)-0.001

def cons2(x):
    return min(x)-0.001

#######################################################################
def create_files():


	header = ["iter", "incidence(millions)",  "beta_H1", "beta_H3", "beta_B", "susceptibility_H1_0", "susceptibility_H1_5", "susceptibility_H1_25", "susceptibility_H1_65", "susceptibility_H3_0", "susceptibility_H3_5", "susceptibility_H3_25", "susceptibility_H3_65","susceptibility_B_0", "susceptibility_B_5", "susceptibility_B_25", "susceptibility_B_65", "vac_eff_hospitalization", "vac_eff_mortality", "prob_hosp_scaling", "prob_death_scaling"]
	
	myfile = open('MAC_results_calibrated_parameters.csv','wb')
	writer = csv.writer(myfile)
	writer.writerow(header)
	
	return myfile, writer
	
	
######################################################################33
if __name__ == "__main__":
	
	elements = {}
	
	yearlist = ['2011-12', '2012-13', '2013-14', '2014-15', '2015-16', '2016-17', '2017-18']
	yearlist = ['2012-13']
	for year in yearlist:
	    df  = pd.read_csv('sampled_parameter_1000_set_year_'+str(year)+'_10May2019.csv')
	    iteration = df['iter'].tolist()

	    data_infections = df['data_incidence'].tolist()
	    data_infections = [num/1e6 for num in data_infections]
	    data_hospitalization = df['data_hospitalizations'].tolist()
	    data_hospitalization = [num/1e3 for num in data_hospitalization]  
	    data_mortality = df['data_mortality'].tolist()
	    data_mortality = [num/1e3 for num in data_mortality]
	    data_incidence_H1_0 = df['data_H1_0'].tolist()
	    data_incidence_H1_5 = df['data_H1_5'].tolist()
	    data_incidence_H1_25 = df['data_H1_25'].tolist()
	    data_incidence_H1_65 = df['data_H1_65'].tolist()
	    data_incidence_H3_0 = df['data_H3_0'].tolist()
	    data_incidence_H3_5 = df['data_H3_5'].tolist()
	    data_incidence_H3_25 = df['data_H3_25'].tolist()
	    data_incidence_H3_65 = df['data_H3_65'].tolist()
	    data_incidence_B_0 = df['data_B_0'].tolist()
	    data_incidence_B_5 = df['data_B_5'].tolist()
	    data_incidence_B_25 = df['data_B_25'].tolist()
	    data_incidence_B_65 = df['data_B_65'].tolist()
	    
    
	    data_list_incidence = [iteration, data_infections, data_hospitalization, data_mortality, data_incidence_H1_0, data_incidence_H1_5, data_incidence_H1_25, data_incidence_H1_65, data_incidence_H3_0, data_incidence_H3_5, data_incidence_H3_25, data_incidence_H3_65, data_incidence_B_0, data_incidence_B_5, data_incidence_B_25, data_incidence_B_65]
	    
	    #myfile, writer=create_files()
	    for num in xrange(1):
		    print ("starting pool"), num
		    beta_H1, beta_H3, beta_B, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65, eff_hosp, eff_death, hosp_scaling, death_scaling = optimize(cons1, cons2, data_list_incidence, year, num)
		    incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH = run_efficacy_simulation(year, num)
		    
		    #elem1 = [num, data_infections[num],beta_H1[0], beta_H3[0], beta_B[0], H1_0[0], H1_5[0], H1_25[0], H1_65[0], H3_0[0], H3_5[0], H3_25[0], H3_65[0], B_0[0], B_5[0], B_25[0], B_65[0], eff_hosp[0], eff_death[0], hosp_scaling[0], death_scaling[0]]
		    #writer.writerow(elem1)
		    #myfile.flush()
	    

			
		
		
	