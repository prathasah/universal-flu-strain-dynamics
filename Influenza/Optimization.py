import numpy as np
np.warnings.filterwarnings('ignore')
from itertools import permutations
from itertools import combinations
import Simulation
import demography
from ages import ages
import doses_distributed as dd
import vaccination_coverage as vc
from operator import itemgetter

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

    def __init__(self, objective = None, season = None, proportion_universalVaccine_doses = 0,  paramValues = {}, index = None):


	self.agelist = [0,0.5,5,20,30,40,50,65]
	self.agelist_full = ages
	self.vaccinated_agelist = [0.5,5,20,30,40,50,65]
	self.full_agelist_dict = {0:[0], 0.5:[0.5], 5: [5,10,15], 20: [20,25],30: [30,35],40:[40,45],50:[50,55,60],65:[65,70,75]}
    
        self.objective = objective
        self.proportional_universal = proportion_universalVaccine_doses
	self.season = season
	self.index = index
	#optimization is FAlse just to extract defaul Parameter values
	self.population = demography.return_demography(season).full(self.agelist)
	self.population_full = demography.return_demography(season).full(self.agelist_full)
	
	
	self.total_vacDoses = dd.total_seasonal_doses_for_season(self.season)
	self.universal_total_vacDoses = self.total_vacDoses * self.proportional_universal
	self.seasonal_total_vacDoses = self.total_vacDoses - self.universal_total_vacDoses
	
	empirical_vaccination_coverage = vc.reduced_age_specific_vaccination_coverage(self.season, self.agelist)
	empirical_vaccination_coverage_full = vc.reduced_age_specific_vaccination_coverage(self.season, self.agelist_full)
	
	#raw dose uptake among age groups
	dosesVaccinated_raw =  {c: a*b for (a,b,c) in zip(empirical_vaccination_coverage, self.population, self.agelist)}
	self.dosesVaccinated = {age: (dosesVaccinated_raw[age]*self.total_vacDoses)/(1.*sum(dosesVaccinated_raw.values())) for age in self.agelist}
	
	dosesVaccinated_full_raw =  {c: a*b for (a,b,c) in zip(empirical_vaccination_coverage_full, self.population_full, self.agelist_full)}
	
	self.dosesVaccinated_full =  {age: (dosesVaccinated_full_raw[age]*self.total_vacDoses)/(1.*sum(dosesVaccinated_full_raw.values())) for age in self.agelist_full}
	
	
    def all_agegroup_receive_doses(self, vax_order):
	""" Do all the age groups chosen receive universal vaccines?"""
	doses_left = self.universal_total_vacDoses
	for group_id in vax_order:
	    popsize = self.dosesVaccinated[group_id]
	    doses_left -= popsize
	    ##check if doses have run out
	    if doses_left < 0 and group_id != vax_order[-1]: return False
	    
	return True
    
    
    def all_UVdoses_distributed(self, vax_order, debug = False):
	"""Are all the universal vaccine distributed? Ensure no vaccine wastage"""
	doses_left = self.universal_total_vacDoses
	if debug: print ("init universal doses"), doses_left
	for age in vax_order[:-1]:
	    ages_to_vax = self.full_agelist_dict[age]
	    for agegroup in ages_to_vax:
		popsize_vaxed = self.dosesVaccinated_full[agegroup]
		doses_left -= popsize_vaxed	    
	    if debug:
		print ("cond=="), agegroup, popsize_vaxed, doses_left
			
	#final age groups to vaccinate    
	age = vax_order[-1]
	ages_to_vax_raw = self.full_agelist_dict[age]
	ages_to_vax_raw = [(age, self.dosesVaccinated_full[age]) for age in ages_to_vax_raw]
	##sort by popsize to avoid vaccine wastage
	ages_to_vax = sorted(ages_to_vax_raw,key=itemgetter(1))
	num_agegroups = len(ages_to_vax)
	for agegroup, popsize_vaxed in ages_to_vax:
	    ##either vaccinate in equal portions or popsize, which ever is smaller
	    doses_left -= min(popsize_vaxed, (doses_left/(1. * num_agegroups)))
	    num_agegroups-=1
	    if debug: print ("cond=="), agegroup, popsize_vaxed, doses_left
	    
	##check if vaccines have been wasted 
	if doses_left >0: return False
	    
	return True
	  
	  
    def compute_full_agewise_doses(self, vax_order):
	
	UVdoses_full = {agegroup:0 for agegroup in self.agelist_full}
	SVdoses_full = {agegroup:self.dosesVaccinated_full[agegroup] for agegroup in self.agelist_full}
	
	#self.all_UVdoses_distributed(vax_order, debug = True)
	doses_left = self.universal_total_vacDoses
	##everone in age groups except last get vaccinated
	for age in vax_order[:-1]:
	    ages_to_vax = self.full_agelist_dict[age]
	    for agegroup in ages_to_vax:
		popsize_vaxed = self.dosesVaccinated_full[agegroup]
		UVdoses_full[agegroup] = popsize_vaxed
		SVdoses_full[agegroup] = 0
 		doses_left -= UVdoses_full[agegroup]
		
	#final age groups to vaccinate    
	age = vax_order[-1]
	ages_to_vax_raw = self.full_agelist_dict[age]
	ages_to_vax_raw = [(age, self.dosesVaccinated_full[age]) for age in ages_to_vax_raw]
	##sort by popsize to avoid vaccine wastage
	ages_to_vax = sorted(ages_to_vax_raw,key=itemgetter(1))
	num_agegroups = len(ages_to_vax)
	for agegroup, popsize_vaxed in ages_to_vax:
	    ##either vaccinate in equal portions or popsize, which ever is smaller
	    
	    UVdoses_full[agegroup] = min(popsize_vaxed, (doses_left/(1. * num_agegroups)))
	    SVdoses_full[agegroup] = popsize_vaxed - UVdoses_full[agegroup]
	    doses_left -= UVdoses_full[agegroup]
	    num_agegroups-=1
	    
	SV_agewise_doses_full = [SVdoses_full[age] for age in self.agelist_full]
	UV_agewise_doses_full = [UVdoses_full[age] for age in self.agelist_full]
	if doses_left>0: print ("WARNING-----> UV doses left"), self.proportional_universal , self.season ,self.index, vax_order, doses_left, sum(UV_agewise_doses_full), self.universal_total_vacDoses

	return np.array(SV_agewise_doses_full), np.array(UV_agewise_doses_full ) 
	    
    def compute_agewise_doses(self, vax_order):
	UVdoses = {agegroup:0 for agegroup in self.agelist}
	doses_left = self.universal_total_vacDoses
	for agegroup in vax_order:
	    popsize_vaxed = self.dosesVaccinated[agegroup]
	    UVdoses[agegroup] = min(popsize_vaxed, doses_left)
	    doses_left -= popsize_vaxed
	    
	UV_agewise_doses = [UVdoses[age] for age in self.agelist]
	return UV_agewise_doses
    
    
    def compute_seasonal_agewise_doses(self, UV_agewise_doses):
	
	total_doses = [self.dosesVaccinated[age] for age in self.agelist]
	return [a-b for (a,b) in zip(total_doses , UV_agewise_doses)]
    
    
    def compute_optimal_coverage(self,SV_agewise_doses, UV_agewise_doses):
	
	SV_coverage = [a/(1.*b) for (a,b) in zip(SV_agewise_doses, self.population)]
	UV_coverage = [a/(1.*b) for (a,b) in zip(UV_agewise_doses, self.population)]
	return SV_coverage, UV_coverage
	
    def solve(self, vax_order):
	
	SV_agewise_doses_full, UV_agewise_doses_full = self.compute_full_agewise_doses(vax_order)
	self.s = Simulation.run_Simulation(season= self.season, proportion_universalVaccine_doses = self.proportional_universal, paramValues = {"SV_doses": SV_agewise_doses_full, "UV_doses": UV_agewise_doses_full}, index=self.index, optimization = False)
	
	   
	 	
    def evaluateObjective(self, vax_order):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	self.solve(vax_order)
	return getattr(self.s, self.objectiveMap[self.objective])
    
    
    def evaluateDetailedObjective(self, vax_order):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""
	self.solve(vax_order)
    
	return getattr(self.s, self.DetailedobjectiveMap[self.objective[5:]])

 
    def optimize(self):

	# number of age groups to be vaccinated
	self.min_output = 10e100
	self.optimized_age_group_vaxed = []
	##evaluate increasing number of vaccinated age groups 
	for number_vax_groups in [1,2,3,4,5,6,7]:
	    if number_vax_groups > 2:
		##combinations where order of all but last is not important
		vax_options = combinations(self.vaccinated_agelist, number_vax_groups-1)
		vax_options = [tuple(list(num) + [a]) for num in vax_options for a in self.vaccinated_agelist if a not in num]
	    else:
		##order does matter when number of vaccinated age groups is<=2 
		vax_options = permutations(self.vaccinated_agelist, number_vax_groups)
	    for vax_order in vax_options:

		if self.all_agegroup_receive_doses(vax_order) and self.all_UVdoses_distributed(vax_order):
		    output =  self.evaluateObjective(list(vax_order))
		    ##debug
		    print vax_order, self.objective, output,output < self.min_output, self.seasonal_total_vacDoses/1e6, self.universal_total_vacDoses/1e6
		    UV_agewise_doses = self.compute_agewise_doses(vax_order)
		    SV_agewise_doses = self.compute_seasonal_agewise_doses(UV_agewise_doses)
		    seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise = self.s.doses_used_agewise()
		    if output < self.min_output:
			    self.min_output = output
			    self.optimized_age_group_vaxed = vax_order


    def optimization_output(self):

	self.solve(self.optimized_age_group_vaxed)
	best_UV_agewise_doses = self.compute_agewise_doses(self.optimized_age_group_vaxed)
	best_SV_agewise_doses = self.compute_seasonal_agewise_doses(best_UV_agewise_doses)
	SV_coverage, UV_coverage = self.compute_optimal_coverage(best_SV_agewise_doses, best_UV_agewise_doses)
	
	seasonal_vacDoses, universal_vacDoses, total_doses = self.s.doses_used()
	seasonal_vacDoses_agewise, universal_vacDoses_agewise, total_doses_agewise = self.s.doses_used_agewise()
	
	return seasonal_vacDoses, universal_vacDoses, total_doses, best_SV_agewise_doses, best_UV_agewise_doses, SV_coverage, UV_coverage, self.evaluateObjective(list(self.optimized_age_group_vaxed)),  list(self.evaluateDetailedObjective(self.optimized_age_group_vaxed))



	
