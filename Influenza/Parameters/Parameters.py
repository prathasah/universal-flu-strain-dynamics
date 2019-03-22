from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages, vaccinationAgesTypical, vaccinationAgesUniversal
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


    def epi_valueOrAttrFromOther(self, key, other, index):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	return (getattr(other, key))(index)
    
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
    
    def epi_setPWAttrFromPassedOrOther(self, other, namePW, index):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.epi_setPWAttr(namePW,self.passedParamValues.epi_valueOrAttrFromOther(namePW, other, index))
	

    def epi_setAttrFromPassedOrOther(self, other, name, index):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.epi_valueOrAttrFromOther(name, other, index))

##################################################################

    def __init__(self, index, calibration, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)
	self.vaccinationAgesTypical = numpy.array(vaccinationAgesTypical)
	self.vaccinationAgesUniversal = numpy.array(vaccinationAgesUniversal)


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
		if isinstance(func(index),(float, int)):
		    self.epi_setAttrFromPassedOrOther(epidemiology,p,index)
		##if it is an agespecific parameter then..
		elif isinstance(func(index),
				    PiecewiseAgeParameter): 
			self.epi_setPWAttrFromPassedOrOther(epidemiology,p, index)


	self.population_highrisk = [(a*b) for (a,b) in zip(self.population, self.proportionHighRisk)]
        self.population_lowrisk = self.population - self.population_highrisk
	
	self.TypicalvaccineEfficacyVsInfectionTypical_H1 = self.passedParamValues["vacEfficacy_seasonal"]  * self.relative_TypicalvaccineEfficacyVsInfection_H1
	self.TypicalvaccineEfficacyVsInfectionTypical_H3 = self.passedParamValues["vacEfficacy_seasonal"]  * self.relative_TypicalvaccineEfficacyVsInfection_H3
	self.TypicalvaccineEfficacyVsInfectionTypical_B = self.passedParamValues["vacEfficacy_seasonal"] * self.relative_TypicalvaccineEfficacyVsInfection_B
	
	self.UniversalvaccineEfficacyVsInfectionUniversal_H1 = self.passedParamValues["vacEfficacy_universal"] [0] * self.relative_UniversalvaccineEfficacyVsInfection_H1
	self.UniversalvaccineEfficacyVsInfectionUniversal_H3 = self.passedParamValues["vacEfficacy_universal"] [0] * self.relative_UniversalvaccineEfficacyVsInfection_H3
	self.UniversalvaccineEfficacyVsInfectionUniversal_B = self.passedParamValues["vacEfficacy_universal"] [1] * self.relative_UniversalvaccineEfficacyVsInfection_B

	# Set up proportion vaccinated vectors
        self.proportionVaccinatedTypicalLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedTypicalHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedUniversalLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
        self.proportionVaccinatedUniversalHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
	
	if len(vaccinationAgesTypical) != len(vaccinationAgesUniversal):
		raise ValueError, "The number of age group bins for low and high efficacy vaccine should be the same!"
	self.proportionVaccinatedTLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedNLPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
	self.proportionVaccinatedTHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesTypical),
            vaccinationAgesTypical)
	self.proportionVaccinatedNHPW = PiecewiseAgeRate([0.0] * len(vaccinationAgesUniversal),
            vaccinationAgesUniversal)
		
	
	
	
	
	self.proportionVaccinatedTL = self.proportionVaccinatedTypicalLPW.full(self.ages)
	self.proportionVaccinatedTH = self.proportionVaccinatedTypicalHPW.full(self.ages)
	self.proportionVaccinatedNL = self.proportionVaccinatedUniversalLPW.full(self.ages)
	self.proportionVaccinatedNH = self.proportionVaccinatedUniversalHPW.full(self.ages)
	self.proportionVaccinatedL = self.proportionVaccinatedTL + self.proportionVaccinatedNL
	self.proportionVaccinatedH = self.proportionVaccinatedTH + self.proportionVaccinatedNH
	self.proportionVaccinatedTypical = self.proportionVaccinatedTL + self.proportionVaccinatedTH
	self.proportionVaccinatedUniversal = self.proportionVaccinatedNL + self.proportionVaccinatedNH
	

	self.proportionVaccinatedLLength = len(vaccinationAgesTypical)
	self.proportionVaccinatedHLength = len(vaccinationAgesTypical)
        

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


        #self.transmissionScaling_H1 = 1.0
	#self.transmissionScaling_H1 *=  self.R0 / self.computeR0_H1()
	#self.transmissionScaling_H3 = 1.0
	#self.transmissionScaling_H3 *=  self.R0 / self.computeR0_H3()
	#self.transmissionScaling_B = 1.0
	#self.transmissionScaling_B *=  self.R0 / self.computeR0_B()
    
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

		
	
