from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages
import demography
import epidemiology
import costs
import os
import numpy as np
import types
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


    def epi_valueOrAttrFromOther(self, key, other, season, index):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return (getattr(other, key))(season, index)
    
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
    
    def epi_setPWAttrFromPassedOrOther(self, other, namePW, season, index):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.epi_setPWAttr(namePW,self.passedParamValues.epi_valueOrAttrFromOther(namePW, other, season, index))
	

    def epi_setAttrFromPassedOrOther(self, other, name, season, index):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.epi_valueOrAttrFromOther(name, other, season, index))

##################################################################

    def __init__(self, season, index, calibration, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)

        # Load in parameters and expand as necessary
	# Go through each files
        for m in (demography, costs):
	    # list all the modules
            for p in dir(m):
		#if module returns a numbers, then..
                if isinstance(getattr(m, p),(float, int)):
		    self.setAttrFromPassedOrOther(m, p)

		##if it is an agespecific parameter then..
                elif isinstance(getattr(m, p),
                                PiecewiseAgeParameter): 
		    self.setPWAttrFromPassedOrOther(m, p)

		    
	    
        for p in dir(epidemiology):
	    #if module returns a numbers, then..
	    func = getattr(epidemiology, p)
	    if isinstance(func,types.FunctionType):
    
		if isinstance(func(season, index),(float, int)):
		    self.epi_setAttrFromPassedOrOther(epidemiology,p, season, index)
		##if it is an agespecific parameter then..
		elif isinstance(func(season, index), PiecewiseAgeParameter): 
			self.epi_setPWAttrFromPassedOrOther(epidemiology,p, season, index)


	self.population_highrisk = self.population * self.proportionHighRisk
        self.population_lowrisk = self.population - self.population_highrisk
	
	self.UniversalvaccineEfficacyVsInfection_H1 =np.array([min(1, num) for num in  ([0.75] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
	self.UniversalvaccineEfficacyVsInfection_H3 =np.array([min(1, num) for num in ([0.75] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
	self.UniversalvaccineEfficacyVsInfection_B = np.array([min(1, num) for num in ([0.75] * self.age_specific_vaccineEfficacyVsInfection)/self.vaccineEfficacyVsInfection_all_ages])
		

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


   
    
	if calibration:
	    if "betaList" in self.passedParamValues:
		self.transmissionScaling_H1 =  self.passedParamValues["betaList"][0]
		self.transmissionScaling_H3 =  self.passedParamValues["betaList"][1]
		self.transmissionScaling_B =   self.passedParamValues["betaList"][2]
	    
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

		
	
