import sys
sys.path.insert(0, r'../Influenza')
import Simulation
###########################################3
def run_efficacy_simulation(year, universal_coverage, sub_iter):
	s = Simulation.run_Simulation(season = year, proportion_universalVaccine_doses = universal_coverage, index = sub_iter)
	seasonal_vacDoses, universal_vacDoses, total_doses = s.doses_used()
	infections_H1, infections_H3, infections_B, hospitalizations_H1, hospitalizations_H3, hospitalizations_B, deaths_H1, deaths_H3, deaths_B, costs_H1, costs_H3, costs_B = s.detailed_strain_output()
	return infections_H1, infections_H3, infections_B, hospitalizations_H1, hospitalizations_H3, hospitalizations_B, deaths_H1, deaths_H3, deaths_B, costs_H1, costs_H3, costs_B
		
		
##########################################################

def compute_simulation(index, year, universal_coverage):


	infections_H1, infections_H3, infections_B, hospitalizations_H1, hospitalizations_H3, hospitalizations_B, deaths_H1, deaths_H3, deaths_B,costs_H1, costs_H3, costs_B =  run_efficacy_simulation(year, universal_coverage, index)
	print year, universal_coverage, index
	print ("total infections"), (sum(infections_H1)/1e6+ sum(infections_H3)/1e6 + sum(infections_B))/1e6
	print ("total deaths"), (sum(deaths_H1)+ sum(deaths_H3) + sum(deaths_B))
	
	

##############################################################################

## input param options
#yearlist = ['2011-12', '2012-13', '2013-14', '2014-15', '2015-16', '2016-17', '2017-18']
#universal_coverage_list = [0, 0.25, 0.5, 0.75, 1]
#index: 0 to 999

##input params = index, season, UV proportion
compute_simulation(565, "2012-13", 0.75)
	
