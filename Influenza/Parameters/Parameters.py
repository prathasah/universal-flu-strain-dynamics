from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages
import demography
import epidemiology
import costs
import os
import numpy as np
import types
import pandas as pd
#from .. import fileIO


import numpy

import UserDict

class ParamDict(UserDict.UserDict):
    def valueOrAttrFromOther(self, key, other):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return getattr(other, key)


    def epi_valueOrAttrFromOther(self, key, other, df, index):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return (getattr(other, key))(df, index)
    
class Parameters:
    def setPWAttr(self, namePW, value):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name.
        '''
        assert isinstance(value, PiecewiseAgeParameter)
        assert namePW.endswith('PW')
        name = namePW[ : -2]
        setattr(self, namePW, value)
        setattr(self, name, value.full(self.ages))
    
    def setPWAttrFromPassedOrOther(self, other, namePW):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.setPWAttr(namePW,self.passedParamValues.valueOrAttrFromOther(namePW, other))
	

    def setAttrFromPassedOrOther(self, other, name):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.valueOrAttrFromOther(name, other))
	
	
	
    def epi_setPWAttr(self, namePW, value):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name.
        '''

        assert isinstance(value, PiecewiseAgeParameter)
        assert namePW.endswith('PW')
        name = namePW[ : -2]
	
        setattr(self, namePW, value)
	
        setattr(self, name, value.full(self.ages))
    
    def epi_setPWAttrFromPassedOrOther(self, other, namePW, df, index):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.epi_setPWAttr(namePW,self.passedParamValues.epi_valueOrAttrFromOther(namePW, other, df, index))
	

    def epi_setAttrFromPassedOrOther(self, other, name, df, index):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.epi_valueOrAttrFromOther(name, other, df, index))

##################################################################

    def __init__(self, season, index, calibration = False, optimization = False, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)
	
	
	self.population = demography.return_demography(season).full(self.ages)

        # Load in parameters and expand as necessary
	# Go through each files
	for p in dir(costs):
	    #if module returns a numbers, then..
	    if isinstance(getattr(costs, p),(float, int)):
		self.setAttrFromPassedOrOther(costs, p)

	    ##if it is an agespecific parameter then..
	    elif isinstance(getattr(costs, p),
			    PiecewiseAgeParameter): 
		self.setPWAttrFromPassedOrOther(costs, p)
		
	if optimization:
	    if not "PVuniversal" in self.passedParamValues:
		print ("Warning! No PVuniversal value supplied.")
	    
	    self.PVuniversal = np.array(self.passedParamValues["PVuniversal"])
	    
		    
	## read file for epidemiology code
	if calibration:
	    df  = pd.read_csv("/Users/prathasah/Dropbox (Bansal Lab)/Git-files/universal-flu-strain-dynamics/calibrate_per_sampled_set/sampled_parameter_1000_set_year_"+season+"_10May2019.csv")
	else:
	    df = pd.read_csv("/Users/prathasah/Dropbox (Bansal Lab)/Git-files/universal-flu-strain-dynamics/calibrate_per_sampled_set/3.calibration_results_June25_2019/3.results_calibrated_parameters_year_"+season+"_COMBINED_June25_2019.csv")
        for p in dir(epidemiology):
	    #if module returns a numbers, then..
	    if calibration and p in ["prob_death_scaling", "prob_hosp_scaling", "beta_H1", "beta_H3", "beta_B", "vac_eff_hospitalization", "vac_eff_mortality", "susceptibility_H1PW", "susceptibility_H3PW", "susceptibility_BPW"]: continue
	    func = getattr(epidemiology, p)
	    if isinstance(func,types.FunctionType):
    
		if isinstance(func(df, index),(float, int)):
		    self.epi_setAttrFromPassedOrOther(epidemiology,p, df, index)
		##if it is an agespecific parameter then..
		elif isinstance(func(df, index), PiecewiseAgeParameter): 
			self.epi_setPWAttrFromPassedOrOther(epidemiology,p, df, index)

	self.population_highrisk = self.population * self.proportionHighRisk
        self.population_lowrisk = self.population - self.population_highrisk
		

        # Get contact matrix
        if 'contactMatrix' in paramValues:
            self.contactMatrix = paramValues.get('contactMatrix')
        else:
            from sys import modules
            import os.path
            import cPickle
            modulePath = os.path.dirname(modules[self.__module__].__file__)
            contactMatrixFile = os.path.join(modulePath, 'contactMatrix.p')
            self.contactMatrix = cPickle.load(open(contactMatrixFile))


   
	if "universalVac_efficacy" in self.passedParamValues:
	    self.Universalvaccine_efficacy = self.passedParamValues["universalVac_efficacy"]
	else: self.Universalvaccine_efficacy =  [0.75,0.75, 0.75]
	
	
	self.UniversalvaccineEfficacyVsInfection_H1 =np.array([min(1, num) for num in  (self.Universalvaccine_efficacy[0] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
	self.UniversalvaccineEfficacyVsInfection_H3 =np.array([min(1, num) for num in (self.Universalvaccine_efficacy[1] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
	self.UniversalvaccineEfficacyVsInfection_B = np.array([min(1, num) for num in (self.Universalvaccine_efficacy[2] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
	
	
	
	if calibration:
	    if "betaList" in self.passedParamValues:
		self.beta_H1 =  self.passedParamValues["betaList"][0]
		self.beta_H3 =  self.passedParamValues["betaList"][1]
		self.beta_B =   self.passedParamValues["betaList"][2]
	    
	    if "vac_eff_hospitalization" in self.passedParamValues: 
		self.vac_eff_hospitalization = self.passedParamValues["vac_eff_hospitalization"]
	    
	    if "vac_eff_mortality" in self.passedParamValues: 
		self.vac_eff_mortality = self.passedParamValues["vac_eff_mortality"]
		
	    if "susceptibility_H1" in self.passedParamValues:
		susceptibility_H1PW = PiecewiseAgeRate(self.passedParamValues["susceptibility_H1"], [0,5,25,65])
		susceptibility_H3PW = PiecewiseAgeRate(self.passedParamValues["susceptibility_H3"], [0,5,25,65])
		susceptibility_BPW = PiecewiseAgeRate(self.passedParamValues["susceptibility_B"], [0,5,25,65])
		setattr(self, "susceptibility_H1",susceptibility_H1PW.full(self.ages))
		setattr(self, "susceptibility_H3", susceptibility_H3PW.full(self.ages))
		setattr(self,"susceptibility_B", susceptibility_BPW.full(self.ages))
		
		
	    if "prob_hosp_scaling" in self.passedParamValues:
		self.prob_hosp_scaling = self.passedParamValues["prob_hosp_scaling"]
		
	    if "prob_death_scaling" in self.passedParamValues:
		self.prob_death_scaling = self.passedParamValues["prob_death_scaling"]

		
	
