#!/usr/bin/python
import sys
sys.path.insert(0, r'./../Influenza')
import numpy as np
import Simulation
import pandas as pd

def run_efficacy_simulation(vacEfficacy_universal, vacEfficacy_seasonal, vacDoses, sub_iter):
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy_universal":vacEfficacy_universal, "vacEfficacy_seasonal": vacEfficacy_seasonal}, index = sub_iter)
	vaccineCoverage = s.compute_typical_vaccination(vacDoses)
	vacsUsedTypical, vacsUsedUniversal = s.simulateWithVaccine(vaccineCoverage, vacEfficacy_seasonal, vacEfficacy_universal, vacDoses)
	#print ("check seasonal, universal doses "), vacsUsedTypical, vacsUsedUniversal
	incidenceL, incidenceH, infections_H1, infections_H3, infections_B, perc_H1, perc_H3, perc_B,hospitalizationsL, hospitalizationsH, deathsL, deathsH = s.calibration_output()
	return incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH
		

######################################################################33
if __name__ == "__main__":
	
	
	df = pd.read_csv("sampled_parameters_with_calibration.csv")
	data_infections = df['data_incidence_x'].tolist()
	data_hosp = df['data_hospitalizations'].tolist()
	data_death = df['data_mortality'].tolist()
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
	
	
	incidence = []
	hosp = []
	total_doses = 141950000.0
	mort = []
	calib_incid = []
	calib_hosp = []
	calib_death = []
	incidence_H1_0 = []
	incidence_H1_5 = []
	incidence_H1_25 = []
	incidence_H1_65 = []
	incidence_H3_0 = []
	incidence_H3_5 = []
	incidence_H3_25 = []
	incidence_H3_65 = []
	incidence_B_0 = []
	incidence_B_5 = []
	incidence_B_25 = []
	incidence_B_65 = []
	
	for sub_index in xrange(5000):
		print sub_index
		scenario = 0
		dose2 = scenario * total_doses
		dose1 = total_doses - dose2
		doses = [dose1, dose2]
		#print ("doses"), doses, total_doses
		universal_efficacy = [0,0]
		seasonal_efficacy  = 0.46
		incidenceL, incidenceH, infections_H1, infections_H3, infections_B, hospitalizationsL, hospitalizationsH, deathsL, deathsH = run_efficacy_simulation(universal_efficacy, seasonal_efficacy,  doses, sub_index)
		
		
		hosp_raw = sum(hospitalizationsL) + sum(hospitalizationsH)
		calib_hosp.append(hosp_raw)
		death_raw = sum(deathsL) + sum(deathsH)
		calib_death.append(death_raw)
		hosp = (abs(hosp_raw - data_hosp[sub_index]))
		mort = (abs(death_raw - data_death[sub_index]))
		
		incidence_raw = sum(incidenceL) + sum(incidenceH)
		calib_incid.append(incidence_raw/1e6)
		
		incidence.append(abs(incidence_raw - data_infections[sub_index]))
		
		H1_0 = (100*sum(list(infections_H1)[0:2]))/(1.*incidence_raw)
                H1_5 = (100*sum(list(infections_H1)[2:6]))/(1.*incidence_raw)
                H1_25 = (100*sum(list(infections_H1)[6:14]))/(1.*incidence_raw)
                H1_65 = (100*sum(list(infections_H1)[14:]))/(1.*incidence_raw)
                        
                H3_0 = (100*sum(list(infections_H3)[0:2]))/(1.*incidence_raw)
                H3_5 = (100*sum(list(infections_H3)[2:6]))/(1.*incidence_raw)
                H3_25 = (100*sum(list(infections_H3)[6:14]))/(1.*incidence_raw)
                H3_65 = (100*sum(list(infections_H3)[14:]))/(1.*incidence_raw)
                        
                B_0 = (100*sum(list(infections_B)[0:2]))/(1.*incidence_raw)
                B_5 = (100*sum(list(infections_B)[2:6]))/(1.*incidence_raw)
                B_25 = (100*sum(list(infections_B)[6:14]))/(1.*incidence_raw)
                B_65 = (100*sum(list(infections_B)[14:]))/(1.*incidence_raw)
		
		incidence_H1_0.append(abs(H1_0 - data_incidence_H1_0[sub_index]))
		incidence_H1_5.append(abs(H1_0 - data_incidence_H1_5[sub_index]))
		incidence_H1_25.append(abs(H1_0 - data_incidence_H1_25[sub_index]))
		incidence_H1_65.append(abs(H1_0 - data_incidence_H1_65[sub_index]))
		
		incidence_H3_0.append(abs(H3_0 - data_incidence_H3_0[sub_index]))
		incidence_H3_5.append(abs(H3_0 - data_incidence_H3_5[sub_index]))
		incidence_H3_25.append(abs(H3_0 - data_incidence_H3_25[sub_index]))
		incidence_H3_65.append(abs(H3_0 - data_incidence_H3_65[sub_index]))
		
		incidence_B_0.append(abs(B_0 - data_incidence_B_0[sub_index]))
		incidence_B_5.append(abs(B_0 - data_incidence_B_5[sub_index]))
		incidence_B_25.append(abs(B_0 - data_incidence_B_25[sub_index]))
		incidence_B_65.append(abs(B_0 - data_incidence_B_65[sub_index]))
		
		
        #print ("----->"), incidence_H3_5
		
	print ("incidence"), np.mean(calib_incid), np.std(calib_incid)
	print ("hospitalizations"), np.mean(calib_hosp), np.std(calib_hosp)
	print ("deaths"), np.mean(calib_death), np.std(calib_death)
	print ("incidence error "), np.mean(incidence)/1e6, np.std(incidence)/1e6
	print ("hospitalizations"), np.mean(hosp)/1e3, np.std(hosp)/1e6
	print ("mortality"), np.mean(mort)/1e3, np.std(mort)/1e6
	print ("age specific H1 infections"), np.mean(incidence_H1_0), np.mean(incidence_H1_5), np.mean(incidence_H1_25),np.mean(incidence_H1_65)
	print ("age specific H3 infections"), np.mean(incidence_H3_0), np.mean(incidence_H3_5), np.mean(incidence_H3_25),np.mean(incidence_H3_65)
	print ("age specific B infections"), np.mean(incidence_B_0), np.mean(incidence_B_5), np.mean(incidence_B_25),np.mean(incidence_B_65)
	print ("std H1 infections"), np.std(incidence_H1_0), np.std(incidence_H1_5), np.std(incidence_H1_25),np.std(incidence_H1_65)
	print ("std H3 infections"), np.std(incidence_H3_0), np.std(incidence_H3_5), np.std(incidence_H3_25),np.std(incidence_H3_65)
	print ("std B infections"), np.std(incidence_B_0), np.std(incidence_B_5), np.std(incidence_B_25),np.std(incidence_B_65)
	
	
