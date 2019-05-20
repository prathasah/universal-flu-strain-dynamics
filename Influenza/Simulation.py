import random
import numpy
import sys
sys.path.insert(0, r'./Influenza/Parameters')
sys.path.insert(0, r'../Influenza/Parameters')
sys.path.insert(0, r'../../Influenza/Parameters')
import Parameters
import vaccination_coverage as vc
import doses_distributed as dd
        
class run_Simulation:
    def __init__(self, tMin=0 , tMax = 180, season = None, proportion_universalVaccine_doses = 0, paramValues = {}, index=None, calibration = False):
        self.tMax = tMax
	self.tMin = tMin
	self.season = season
	
	self.parameters = Parameters.Parameters(season, index, calibration, **paramValues)

        # Initial condition
        self.Y0 = numpy.zeros(46 * self.parameters.ages.size)
        self.hasSolution = False
	    
	#############################    
	self.resetSolution()
	
	self.proportional_universal = proportion_universalVaccine_doses
	
		   
        # Update initial condition for ODEs
        self.updateIC()
	
	tEnd = self.tMax
	tStart = self.tMin
	self.solve(tStart = tStart, tEnd = tEnd)
	
	import matplotlib.pyplot as plt
	times = [num for num in xrange(181)]
	plt.plot(times, self.IUL_H1.sum(axis=1))
	#plt.show()
	
        self.updateStats()
	
	

	#######################

    def getLastValues(self):
        return (self.SUL[-1, :], self.IUL_H1[-1, :], self.IUL_H3[-1, :], self.IUL_B[-1, :], self.RUL_H1[-1, :], self.RUL_H3[-1, :], self.RUL_B[-1, :], 
		self.SUH[-1, :], self.IUH_H1[-1, :], self.IUH_H3[-1, :], self.IUH_B[-1, :], self.RUH_H1[-1, :], self.RUH_H3[-1, :], self.RUH_B[-1, :], 
		self.STL[-1, :], self.ITL_H1[-1, :], self.ITL_H3[-1, :], self.ITL_B[-1, :], self.RTL_H1[-1, :], self.RTL_H3[-1, :], self.RTL_B[-1, :], 
		self.STH[-1, :], self.ITH_H1[-1, :], self.ITH_H3[-1, :], self.ITH_B[-1, :], self.RTH_H1[-1, :], self.RTH_H3[-1, :], self.RTH_B[-1, :],
		self.SNL[-1, :], self.INL_H1[-1, :], self.INL_H3[-1, :], self.INL_B[-1, :], self.RNL_H1[-1, :], self.RNL_H3[-1, :], self.RNL_B[-1, :],
		self.SNH[-1, :], self.INH_H1[-1, :], self.INH_H3[-1, :], self.INH_B[-1, :], self.RNH_H1[-1, :], self.RNH_H3[-1, :], self.RNH_B[-1, :],
		self.vacc_TL[-1, :], self.vacc_TH[-1,:],  self.vacc_NL[-1,:], self.vacc_NH[-1,:])
    
    ########################3
    def update_doses_distributed(self, seasonal_vacDoses, universal_vacDoses):

	
	empirical_vax_coverage_lowrisk, empirical_vax_coverage_highrisk  = vc.age_specific_vaccination_coverage(self.season)
	
	#raw dose uptake among age groups
	dosesVaccinatedL =  empirical_vax_coverage_lowrisk * self.parameters.population_lowrisk
	dosesVaccinatedH =  empirical_vax_coverage_highrisk * self.parameters.population_highrisk
	
	
	
	##coverage RELATIVE to each other (based is age-group 0.5 for Low risk age groups)
	relative_coverage_L = (1.*dosesVaccinatedL/dosesVaccinatedL[1]) 
	relative_coverage_H = (1.*dosesVaccinatedH/dosesVaccinatedL[1]) 
	
	#multiplication factor for final doses
	
	xFac_T = seasonal_vacDoses/(1.*(sum(relative_coverage_L)+ sum(relative_coverage_H)))
	xFac_N = universal_vacDoses/(1.*(sum(relative_coverage_L)+ sum(relative_coverage_H)))
	
	doses_TL = xFac_T *  relative_coverage_L
	doses_TH = xFac_T *  relative_coverage_H
	doses_NL = xFac_N *  relative_coverage_L
	doses_NH = xFac_N *  relative_coverage_H
	
        return doses_TL, doses_TH, doses_NL, doses_NH
    
    ###################################
    
    def updateIC(self):
        
	
	initial_vacDoses = dd.doses_applied_before_start_season(self.season)
	universal_vacDoses = initial_vacDoses * self.proportional_universal
	seasonal_vacDoses = initial_vacDoses - universal_vacDoses
	
	

	doses_TL, doses_TH, doses_NL, doses_NH = self.update_doses_distributed(seasonal_vacDoses, universal_vacDoses)
	
	proportionVaccinatedTL = doses_TL/(1.* self.parameters.population_lowrisk)
	proportionVaccinatedTH = doses_TH/(1.* self.parameters.population_highrisk)
	proportionVaccinatedNL = doses_NL/(1.* self.parameters.population_lowrisk)
	proportionVaccinatedNH = doses_NH/(1.* self.parameters.population_highrisk)
	

	if not self.hasSolution:
            # S
	     
	    ## SUL
            self.Y0[ 0: : 46] =  (1 - proportionVaccinatedTL -  proportionVaccinatedNL) * self.parameters.population_lowrisk 

	    ## SUH
            self.Y0[ 7: : 46] =  (1 - proportionVaccinatedTH -  proportionVaccinatedNH) * self.parameters.population_highrisk 
  
	    ## STL
            self.Y0[ 14: : 46] = proportionVaccinatedTL * self.parameters.population_lowrisk 
	    
	    ## STH
            self.Y0[ 21: : 46] =  proportionVaccinatedTH * self.parameters.population_highrisk

	    ## SNL
            self.Y0[ 28: : 46] = proportionVaccinatedNL * self.parameters.population_lowrisk 
	    ## SNH 
            self.Y0[ 35: : 46] = proportionVaccinatedNH  * self.parameters.population_highrisk
	
            # initially infection
	    initially_infected = 500
	    # IUL_H1, IUL_H3, IUL_B
	    self.Y0[ 1: : 46] = numpy.full(self.parameters.ages.size,initially_infected)
	    self.Y0[ 2: : 46] = numpy.full(self.parameters.ages.size,initially_infected)
	    self.Y0[ 3: : 46] = numpy.full(self.parameters.ages.size,initially_infected)
	    
	    # IUH_H1, IUH_H3, IUH_B
            self.Y0[ 8: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 9: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 10: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    
	    # ITL_H1, ITL_H3, ITL_B
	    self.Y0[ 15: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 16: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 17: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    
	    # ITH_H1, ITH_H3, ITH_B
	    self.Y0[ 22: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 23: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 24: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    
	    # INL_H1, INL_H3, INL_B
	    self.Y0[ 29: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 30: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 31: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    
	    # INH_H1, INH_H3, INH_B
	    self.Y0[ 36: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 37: : 46] = numpy.full(self.parameters.ages.size, initially_infected)
	    self.Y0[ 38: : 46] = numpy.full(self.parameters.ages.size, initially_infected)

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0: : 46] -= (self.Y0[ 1: : 46] + self.Y0[ 2: : 46] + self.Y0[ 3: : 46])
	    self.Y0[ 7: : 46] -= (self.Y0[ 8: : 46] + self.Y0[ 9: : 46] + self.Y0[ 10: : 46])
	    self.Y0[ 14: : 46] -= (self.Y0[ 15: : 46] + self.Y0[ 16: : 46] + self.Y0[ 17: : 46])
	    self.Y0[ 21: : 46] -= (self.Y0[ 22: : 46] + self.Y0[ 23: : 46] + self.Y0[ 24: : 46])
	    self.Y0[ 28: : 46] -= (self.Y0[ 29: : 46] + self.Y0[ 30: : 46] + self.Y0[ 31: : 46])
	    self.Y0[ 35: : 46] -= (self.Y0[ 32: : 46] + self.Y0[ 33: : 46] + self.Y0[ 34: : 46])
 
            # R
	    self.Y0[ 4:  : 46] = 0.
	    self.Y0[ 5:  : 46] = 0.
	    self.Y0[ 6:  : 46] = 0.
	    self.Y0[ 11:  : 46] = 0.
	    self.Y0[ 12:  : 46] = 0.
	    self.Y0[ 13:  : 46] = 0.
	    self.Y0[ 18:  : 46] = 0.
	    self.Y0[ 19:  : 46] = 0.
	    self.Y0[ 20:  : 46] = 0.
	    self.Y0[ 25:  : 46] = 0.
	    self.Y0[ 26:  : 46] = 0.
	    self.Y0[ 27:  : 46] = 0.
	    self.Y0[ 32:  : 46] = 0.
	    self.Y0[ 33:  : 46] = 0.
	    self.Y0[ 34:  : 46] = 0.
	    self.Y0[ 39:  : 46] = 0.
	    self.Y0[ 40:  : 46] = 0.
	    self.Y0[ 41:  : 46] = 0.
	    
	    ##doses
	    self.Y0[ 42:  : 46] = doses_TL
	    self.Y0[ 43:  : 46] = doses_TH
	    self.Y0[ 44:  : 46] = doses_NL
	    self.Y0[ 45:  : 46] = doses_NH
	    
        else:
	    SUL, IUL_H1, IUL_H3, IUL_B, RUL_H1, RUL_H3, RUL_B,
	    SUH, IUH_H1, IUH_H3, IUH_B, RUH_H1, RUH_H3, RUH_B,
	    STL, ITL_H1, ITL_H3, ITL_B, RTL_H1, RTL_H3, RTL_B,  
	    STH, ITH_H1, ITH_H3, ITH_B, RTH_H1, RTH_H3, RTH_B,
	    SNL, INL_H1, INL_H3, INL_B, RNL_H1, RNL_H3, RNL_B, 
	    SNH, INH_H1, INH_H3, INH_B, RNH_H1, RNH_H3, RNH_B,
	    vacc_TL ,vacc_TH, vacc_NL, vacc_NH = self.getLastValues()
	    
            self.Y0[ 0 : : 46] = (1 - proportionVaccinatedL) * SUL
	    self.Y0[ 7 : : 46] = (1 - proportionVaccinatedH) * SUH
	    self.Y0[ 14 : : 46] = STL + (1 - proportionVaccinatedTL) * SUL
	    self.Y0[ 21 : : 46] = STH + (1 - proportionVaccinatedTH) * SUH
	    self.Y0[ 28 : : 46] = SNL + (1 - proportionVaccinatedNL) * SUL
	    self.Y0[ 35 : : 46] = SNH + (1 - proportionVaccinatedNH) * SUH
	    
	    #I
	    self.Y0[ 1 : : 46] = IUL_H1
	    self.Y0[ 2 : : 46] = IUL_H3
	    self.Y0[ 3 : : 46] = IUL_B
	    
	    self.Y0[ 8 : : 46] = IUH_H1
	    self.Y0[ 9 : : 46] = IUH_H3
	    self.Y0[ 10 : : 46] = IUH_B
	    
	    self.Y0[ 15 : : 46] = ITL_H1
	    self.Y0[ 16 : : 46] = ITL_H3
	    self.Y0[ 17 : : 46] = ITL_B
	    
	    self.Y0[ 22 : : 46] = ITH_H1
	    self.Y0[ 23 : : 46] = ITH_H3
	    self.Y0[ 24 : : 46] = ITH_B
	    
	    self.Y0[ 29 : : 46] = INL_H1
	    self.Y0[ 30 : : 46] = INL_H3
	    self.Y0[ 31 : : 46] = INL_B
	    
	    self.Y0[ 36 : : 46] = INH_H1
	    self.Y0[ 37 : : 46] = INH_H3
	    self.Y0[ 38 : : 46] = INH_B
	    
	    #R class
	    self.Y0[ 4 : : 46] = RUL_H1
	    self.Y0[ 5 : : 46] = RUL_H3
	    self.Y0[ 6 : : 46] = RUL_B
	    
	    self.Y0[ 11 : : 46] = RUH_H1
	    self.Y0[ 12 : : 46] = RUH_H3
	    self.Y0[ 13 : : 46] = RUH_B
	    
	    self.Y0[ 18 : : 46] = RTL_H1
	    self.Y0[ 19 : : 46] = RTL_H3
	    self.Y0[ 20 : : 46] = RTL_B
	    
	    self.Y0[ 25 : : 46] = RTH_H1
	    self.Y0[ 26 : : 46] = RTH_H3
	    self.Y0[ 27 : : 46] = RTH_B
	    
	    self.Y0[ 32 : : 46] = RNL_H1
	    self.Y0[ 33 : : 46] = RNL_H3
	    self.Y0[ 34 : : 46] = RNL_B
	    
	    self.Y0[ 39 : : 46] = RNH_H1
	    self.Y0[ 40 : : 46] = RNH_H3
	    self.Y0[ 41 : : 46] = RNH_B
	    
	    ## vaccine doses
	    self.Y0[ 42 : : 46] = vacc_TL
	    self.Y0[ 43 : : 46] = vacc_TH
	    self.Y0[ 44 : : 46] = vacc_NL
	    self.Y0[ 45 : : 46] = vacc_NH
	    

    def RHS(self, Y, t):
        '''
        SIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL    = Y[ 0 : : 46]
        IUL_H1 = Y[ 1 : : 46]
        IUL_H3 = Y[ 2 : : 46]
        IUL_B  = Y[ 3 : : 46]
        RUL_H1 = Y[ 4 : : 46]
        RUL_H3 = Y[ 5 : : 46]
        RUL_B =  Y[ 6 : : 46]
            
        SUH    = Y[ 7 : : 46]
        IUH_H1 = Y[ 8 : : 46]
        IUH_H3 = Y[ 9 : : 46]
        IUH_B  = Y[ 10 : : 46]
        RUH_H1 = Y[ 11 : : 46]
        RUH_H3 = Y[ 12 : : 46]
        RUH_B =  Y[ 13 : : 46]
            
        STL    = Y[ 14 : : 46]
        ITL_H1 = Y[ 15 : : 46]
        ITL_H3 = Y[ 16 : : 46]
        ITL_B  = Y[ 17 : : 46]
        RTL_H1 = Y[ 18 : : 46]
        RTL_H3 = Y[ 19 : : 46]
        RTL_B =  Y[ 20 : : 46]
            
        STH    = Y[ 21 : : 46]
        ITH_H1 = Y[ 22 : : 46]
        ITH_H3 = Y[ 23 : : 46]
        ITH_B  = Y[ 24 : : 46]
        RTH_H1 = Y[ 25 : : 46]
        RTH_H3 = Y[ 26 : : 46]
        RTH_B =  Y[ 27 : : 46]
            
        SNL    = Y[ 28 : : 46]
        INL_H1 = Y[ 29 : : 46]
        INL_H3 = Y[ 30 : : 46]
        INL_B  = Y[ 31 : : 46]
        RNL_H1 = Y[ 32 : : 46]
        RNL_H3 = Y[ 33 : : 46]
        RNL_B =  Y[ 34 : : 46]
            
        SNH    = Y[ 35 : : 46]
        INH_H1 = Y[ 36 : : 46]
        INH_H3 = Y[ 37 : : 46]
        INH_B  = Y[ 38 : : 46]
        RNH_H1 = Y[ 39 : : 46]
        RNH_H3 = Y[ 40 : : 46]
        RNH_B =  Y[ 41 : : 46]
	
	vacc_TL = Y[42 : : 46]
	vacc_TH = Y[43 : : 46]
	vacc_NL = Y[44 : : 46]
	vacc_NH = Y[45 : : 46]
	
        N =  sum(SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1 + RUL_H3 + RUL_B + 
                SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH_H1 + RUH_H3 + RUH_B +
                STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL_H1 + RTL_H3 + RTL_B + 
                STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH_H1 + RTH_H3 + RTH_B + 
                SNL+ INL_H1+ INL_H3+ INL_B+ RNL_H1 + RNL_H3 + RNL_B + 
                SNH+ INH_H1+ INH_H3+ INH_B+ RNH_H1 + RNH_H3 + RNH_B )
            
        N_age_specific = SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1 + RUL_H3+ RUL_B + SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH_H1 + RUH_H3 + RUH_B +STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL_H1 + RTL_H3 + RTL_B +  STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH_H1 + RTH_H3 + RTH_B +  SNL+ INL_H1+ INL_H3+ INL_B+ RNL_H1 + RNL_H3 + RNL_B +  SNH+ INH_H1+ INH_H3+ INH_B+ RNH_H1 + RNH_H3 + RNH_B
            
	
        Lambda_H1 = self.parameters.transmissionScaling_H1 * self.parameters.susceptibility_H1\
		    * numpy.dot(self.parameters.contactMatrix, self.parameters.transmissibility * (IUL_H1 + IUH_H1 + ITL_H1 + ITH_H1+ INL_H1+ INH_H1)) / N_age_specific
	
	Lambda_H3 = self.parameters.transmissionScaling_H3 * self.parameters.susceptibility_H3 \
                 * numpy.dot(self.parameters.contactMatrix, 
                             self.parameters.transmissibility * (IUL_H3 + IUH_H3 + ITL_H3 + ITH_H3 + INL_H3+ INH_H3)) / N_age_specific
		
	Lambda_B = self.parameters.transmissionScaling_B * self.parameters.susceptibility_B \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL_B + IUH_B + ITL_B + ITH_B+ INL_B+ INH_B)) / N_age_specific
	
	
	vaccine_doses_t_raw = dd.doses_applied_per_time_function(self.season)(t)
	universal_vacDoses = vaccine_doses_t_raw* self.proportional_universal
	seasonal_vacDoses = vaccine_doses_t_raw - universal_vacDoses
	
	
	doses_TL, doses_TH, doses_NL, doses_NH = self.update_doses_distributed(seasonal_vacDoses, universal_vacDoses)
	
	proportionVaccinatedTL = doses_TL/(1.*SUL)
	proportionVaccinatedTH = doses_TH/(1.*SUH)
	proportionVaccinatedNL = doses_NL/(1.*SUL)
	proportionVaccinatedNH = doses_NH/(1.*SUH)
	
	
        # The right-hand sides
	
	#UL
        dSUL    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUL  -  (proportionVaccinatedTL +  proportionVaccinatedNL) * SUL 
        dIUL_H1 = (Lambda_H1 * SUL) - (self.parameters.recoveryRate ) * IUL_H1
	dIUL_H3 = (Lambda_H3 * SUL) - (self.parameters.recoveryRate ) * IUL_H3
	dIUL_B  = (Lambda_B * SUL) - (self.parameters.recoveryRate ) * IUL_B
	dRUL_H1    = self.parameters.recoveryRate * IUL_H1
	dRUL_H3    = self.parameters.recoveryRate * IUL_H3
	dRUL_B    = self.parameters.recoveryRate * IUL_B

	
	#UH
        dSUH    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUH -  (proportionVaccinatedTH +  proportionVaccinatedNH) * SUH 
        dIUH_H1 = (Lambda_H1 * SUH) - (self.parameters.recoveryRate ) * IUH_H1
	dIUH_H3 = (Lambda_H3 * SUH) - (self.parameters.recoveryRate ) * IUH_H3
	dIUH_B  = (Lambda_B * SUH) - (self.parameters.recoveryRate ) * IUH_B
	dRUH_H1    = self.parameters.recoveryRate * IUH_H1
	dRUH_H3    = self.parameters.recoveryRate * IUH_H3
	dRUH_B    = self.parameters.recoveryRate * IUH_B
	
	#TL
	dSTL = (proportionVaccinatedTL * SUL) - ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H1) * Lambda_H1 + (1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H3) * Lambda_H3+ (1 - self.parameters.SeasonalVaccineEfficacyVsInfection_B) * Lambda_B)  *STL

        dITL_H1 = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H1) * Lambda_H1 * STL) - (self.parameters.recoveryRate ) * ITL_H1
	dITL_H3 = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H3) * Lambda_H3 * STL) - (self.parameters.recoveryRate) * ITL_H3
	dITL_B  = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_B) * Lambda_B * STL) - (self.parameters.recoveryRate) * ITL_B
	dRTL_H1    = self.parameters.recoveryRate * ITL_H1
	dRTL_H3    = self.parameters.recoveryRate * ITL_H3
	dRTL_B    = self.parameters.recoveryRate * ITL_B
	
	#TH
	dSTH =(proportionVaccinatedTH * SUH) - ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H1) * Lambda_H1 + (1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H3) * Lambda_H3+ (1 - self.parameters.SeasonalVaccineEfficacyVsInfection_B) * Lambda_B)  *STH
	
        dITH_H1 = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H1) *Lambda_H1 * STH) - (self.parameters.recoveryRate) * ITH_H1
	dITH_H3 = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_H3) *Lambda_H3 * STH) - (self.parameters.recoveryRate) * ITH_H3
	dITH_B = ((1 - self.parameters.SeasonalVaccineEfficacyVsInfection_B) * Lambda_B * STH) - (self.parameters.recoveryRate) * ITH_B
	dRTH_H1    = self.parameters.recoveryRate * ITH_H1
	dRTH_H3    = self.parameters.recoveryRate * ITH_H3
	dRTH_B    = self.parameters.recoveryRate * ITH_B
	
	#NL
	dSNL = (proportionVaccinatedNL * SUL) - ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H1) * Lambda_H1 + (1 - self.parameters.UniversalvaccineEfficacyVsInfection_H3) * Lambda_H3+ (1 - self.parameters.UniversalvaccineEfficacyVsInfection_B) * Lambda_B)  *SNL
	
        dINL_H1 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H1) *Lambda_H1 * SNL) - (self.parameters.recoveryRate ) * INL_H1
	dINL_H3 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H3) *Lambda_H3 * SNL) - (self.parameters.recoveryRate ) * INL_H3
	dINL_B = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_B) * Lambda_B * SNL) - (self.parameters.recoveryRate ) * INL_B
	dRNL_H1    = self.parameters.recoveryRate * INL_H1
	dRNL_H3    = self.parameters.recoveryRate * INL_H3
	dRNL_B    = self.parameters.recoveryRate * INL_B
	
	#NH
	dSNH =  (proportionVaccinatedNH * SUH) - ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H1) * Lambda_H1 + (1 - self.parameters.UniversalvaccineEfficacyVsInfection_H3) * Lambda_H3+ (1 - self.parameters.UniversalvaccineEfficacyVsInfection_B) * Lambda_B)  *SNH
        dINH_H1 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H1) *Lambda_H1 * SNH) - (self.parameters.recoveryRate ) * INH_H1
	dINH_H3 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_H3) *Lambda_H3 * SNH) - (self.parameters.recoveryRate ) * INH_H3
	dINH_B = ((1 - self.parameters.UniversalvaccineEfficacyVsInfection_B) *Lambda_B * SNH) - (self.parameters.recoveryRate ) * INH_B
	dRNH_H1    = self.parameters.recoveryRate * INH_H1
	dRNH_H3    = self.parameters.recoveryRate * INH_H3
	dRNH_B    = self.parameters.recoveryRate * INH_B
	
	##vaccine doses distributed
	dvacc_TL  = doses_TL
	dvacc_TH = doses_TH
	dvacc_NL = doses_NL
	dvacc_NH = doses_NH
	
	
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 46] = dSUL
        dY[ 1 : : 46] = dIUL_H1
        dY[ 2 : : 46] = dIUL_H3
	dY[ 3 : : 46] = dIUL_B
	dY[ 4 : : 46] = dRUL_H1
	dY[ 5 : : 46] = dRUL_H3
	dY[ 6 : : 46] = dRUL_B
	
	dY[ 7 : : 46] = dSUH
        dY[ 8 : : 46] = dIUH_H1
        dY[ 9 : : 46] = dIUH_H3
	dY[ 10 : : 46] = dIUH_B
	dY[ 11 : : 46] = dRUH_H1
	dY[ 12 : : 46] = dRUH_H3
	dY[ 13 : : 46] = dRUH_B
	
	dY[ 14 : : 46] = dSTL
        dY[ 15 : : 46] = dITL_H1
        dY[ 16 : : 46] = dITL_H3
	dY[ 17 : : 46] = dITL_B
	dY[ 18 : : 46] = dRTL_H1
	dY[ 19 : : 46] = dRTL_H3
	dY[ 20 : : 46] = dRTL_B
	
	dY[ 21 : : 46] = dSTH
        dY[ 22 : : 46] = dITH_H1
        dY[ 23 : : 46] = dITH_H3
	dY[ 24 : : 46] = dITH_B
	dY[ 25 : : 46] = dRTH_H1
	dY[ 26 : : 46] = dRTH_H3
	dY[ 27 : : 46] = dRTH_B
	
	dY[ 28 : : 46] = dSNL
        dY[ 29 : : 46] = dINL_H1
        dY[ 30 : : 46] = dINL_H3
	dY[ 31 : : 46] = dINL_B
	dY[ 32 : : 46] = dRNL_H1
	dY[ 33 : : 46] = dRNL_H3
	dY[ 34 : : 46] = dRNL_B
	
	dY[ 35 : : 46] = dSNH
        dY[ 36 : : 46] = dINH_H1
        dY[ 37 : : 46] = dINH_H3
	dY[ 38 : : 46] = dINH_B
	dY[ 39 : : 46] = dRNH_H1
	dY[ 40 : : 46] = dRNH_H3
	dY[ 41 : : 46] = dRNH_B
	
	
	dY[42 : : 46] = dvacc_TL
	dY[43 : : 46] = dvacc_TH
	dY[44 : : 46] = dvacc_NL
	dY[45 : : 46] = dvacc_NH
	
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
	    
            SUL_Old = self.SUL.copy()
            IUL_H1_Old = self.IUL_H1.copy()
	    IUL_H3_Old = self.IUL_H3.copy()
	    IUL_B_Old = self.IUL_B.copy()
	    RUL_H1_Old = self.RUL_H1.copy()
	    RUL_H3_Old = self.RUL_H3.copy()
	    RUL_B_Old = self.RUL_B.copy()
	    
	    SUH_Old = self.SUH.copy()
            IUH_H1_Old = self.IUH_H1.copy()
	    IUH_H3_Old = self.IUH_H3.copy()
	    IUH_B_Old = self.IUH_B.copy()
	    RUH_H1_Old = self.RUH_H1.copy()
	    RUH_H3_Old = self.RUH_H3.copy()
	    RUH_B_Old = self.RUH_B.copy()
	    
	    STL_Old = self.STL.copy()
            ITL_H1_Old = self.ITL_H1.copy()
	    ITL_H3_Old = self.ITL_H3.copy()
	    ITL_B_Old = self.ITL_B.copy()
	    RTL_H1_Old = self.RTL_H1.copy()
	    RTL_H3_Old = self.RTL_H3.copy()
	    RTL_B_Old = self.RTL_B.copy()
	    
	    STH_Old = self.STH.copy()
            ITH_H1_Old = self.ITH_H1.copy()
	    ITH_H3_Old = self.ITH_H3.copy()
	    ITH_B_Old = self.ITH_B.copy()
	    RTH_H1_Old = self.RTH_H1.copy()
	    RTH_H3_Old = self.RTH_H3.copy()
	    RTH_B_Old = self.RTH_B.copy()
	    
	    SNL_Old = self.SNL.copy()
            INL_H1_Old = self.INL_H1.copy()
	    INL_H3_Old = self.INL_H3.copy()
	    INL_B_Old = self.INL_B.copy()
	    RNL_H1_Old = self.RNL_H1.copy()
	    RNL_H3_Old = self.RNL_H3.copy()
	    RNL_B_Old = self.RNL_B.copy()
	    
	    SNH_Old = self.SNH.copy()
            INH_H1_Old = self.INH_H1.copy()
	    INH_H3_Old = self.INH_H3.copy()
	    INH_B_Old = self.INH_B.copy()
	    RNH_H1_Old = self.RNH_H1.copy()
	    RNH_H3_Old = self.RNH_H3.copy()
	    RNH_B_Old = self.RNH_B.copy()
	    
	    vacc_TL_old = self.vacc_TL.copy()
	    vacc_TH_old = self.vacc_TH.copy()
	    vacc_NL_old = self.vacc_NL.copy()
	    vacc_NH_old = self.vacc_NH.copy()
	    
	
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))

        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y, infodict = odeint(self.RHS,self.Y0.copy(),self.T, mxstep = 1000, full_output =True)
	
	#print ("steps==="),  self.Y.shape, (numpy.unique(infodict['tcur']))
        Z = self.Y.copy()
	
	self.SUL    = Z[:, 0 : : 46]
        self.IUL_H1 = Z[:, 1 : : 46]
	self.IUL_H3 = Z[:, 2 : : 46]
	self.IUL_B  = Z[:, 3 : : 46]
	self.RUL_H1 = Z[:, 4 : : 46]
	self.RUL_H3 = Z[:, 5 : : 46]
	self.RUL_B  = Z[:, 6 : : 46]
	
	self.SUH    = Z[:, 7 : : 46]
        self.IUH_H1 = Z[:, 8 : : 46]
	self.IUH_H3 = Z[:, 9 : : 46]
	self.IUH_B  = Z[:, 10 : : 46]
	self.RUH_H1 = Z[:, 11 : : 46]
	self.RUH_H3 = Z[:, 12 : : 46]
	self.RUH_B  = Z[:, 13 : : 46]
	
	self.STL    = Z[:, 14 : : 46]
        self.ITL_H1 = Z[:, 15 : : 46]
	self.ITL_H3 = Z[:, 16 : : 46]
	self.ITL_B  = Z[:, 17 : : 46]
	self.RTL_H1 = Z[:, 18 : : 46]
	self.RTL_H3 = Z[:, 19 : : 46]
	self.RTL_B  = Z[:, 20 : : 46]
	
	self.STH    = Z[:, 21 : : 46]
        self.ITH_H1 = Z[:, 22 : : 46]
	self.ITH_H3 = Z[:, 23 : : 46]
	self.ITH_B  = Z[:, 24 : : 46]
	self.RTH_H1 = Z[:, 25 : : 46]
	self.RTH_H3 = Z[:, 26 : : 46]
	self.RTH_B  = Z[:, 27 : : 46]
	
	self.SNL    = Z[:, 28 : : 46]
        self.INL_H1 = Z[:, 29 : : 46]
	self.INL_H3 = Z[:, 30 : : 46]
	self.INL_B  = Z[:, 31 : : 46]
	self.RNL_H1 = Z[:, 32 : : 46]
	self.RNL_H3 = Z[:, 33 : : 46]
	self.RNL_B  = Z[:, 34 : : 46]
	
	self.SNH    = Z[:, 35 : : 46]
        self.INH_H1 = Z[:, 36 : : 46]
	self.INH_H3 = Z[:, 37 : : 46]
	self.INH_B  = Z[:, 38 : : 46]
	self.RNH_H1 = Z[:, 39 : : 46]
	self.RNH_H3 = Z[:, 40 : : 46]
	self.RNH_B  = Z[:, 41 : : 46]
	
	self.vacc_TL  = Z[:, 42 : : 46]
	self.vacc_TH  = Z[:, 43 : : 46]
	self.vacc_NL  = Z[:, 44 : : 46]
	self.vacc_NH  = Z[:, 45 : : 46]
	
	
	
        if self.hasSolution:
	   
            self.T = numpy.hstack((TOld, self.T))
	    #UL
            self.SUL = numpy.vstack((SUL_Old, self.SUL))
            self.IUL_H1 = numpy.vstack((IUL_H1_Old, self.IUL_H1))
	    self.IUL_H3 = numpy.vstack((IUL_H3_Old, self.IUL_H3))
	    self.IUL_B = numpy.vstack((IUL_B_Old, self.IUL_B))
	    self.RUL_H1 =  numpy.vstack((RUL_H1_Old, self.RUL_H1))
	    self.RUL_H3 =  numpy.vstack((RUL_H3_Old, self.RUL_H3))
	    self.RUL_B =  numpy.vstack((RUL_B_Old, self.RUL_B))
	    
	    #UH
	    self.SUH = numpy.vstack((SUH_Old, self.SUH))
            self.IUH_H1 = numpy.vstack((IUH_H1_Old, self.IUH_H1))
	    self.IUH_H3 = numpy.vstack((IUH_H3_Old, self.IUH_H3))
	    self.IUH_B = numpy.vstack((IUH_B_Old, self.IUH_B))
	    self.RUH_H1 =  numpy.vstack((RUH_H1_Old, self.RUH_H1))
	    self.RUH_H3 =  numpy.vstack((RUH_H3_Old, self.RUH_H3))
	    self.RUH_B =  numpy.vstack((RUH_B_Old, self.RUH_B))
	    
	    #TL
	    self.STL = numpy.vstack((STL_Old, self.STL))
            self.ITL_H1 = numpy.vstack((ITL_H1_Old, self.ITL_H1))
	    self.ITL_H3 = numpy.vstack((ITL_H3_Old, self.ITL_H3))
	    self.ITL_B = numpy.vstack((ITL_B_Old, self.ITL_B))
	    self.RTL_H1 =  numpy.vstack((RTL_H1_Old, self.RTL_H1))
	    self.RTL_H3 =  numpy.vstack((RTL_H3_Old, self.RTL_H3))
	    self.RTL_B =  numpy.vstack((RTL_B_Old, self.RTL_B))
	    
	    #TH
	    self.STH = numpy.vstack((STH_Old, self.STH))
            self.ITH_H1 = numpy.vstack((ITH_H1_Old, self.ITH_H1))
	    self.ITH_H3 = numpy.vstack((ITH_H3_Old, self.ITH_H3))
	    self.ITH_B = numpy.vstack((ITH_B_Old, self.ITH_B))
	    self.RTH_H1 =  numpy.vstack((RTH_H1_Old, self.RTH_H1))
	    self.RTH_H3 =  numpy.vstack((RTH_H3_Old, self.RTH_H3))
	    self.RTH_B =  numpy.vstack((RTH_B_Old, self.RTH_B))
	    
	    #NL
	    self.SNL = numpy.vstack((SNL_Old, self.SNL))
            self.INL_H1 = numpy.vstack((INL_H1_Old, self.INL_H1))
	    self.INL_H3 = numpy.vstack((INL_H3_Old, self.INL_H3))
	    self.INL_B = numpy.vstack((INL_B_Old, self.INL_B))
	    self.RNL_H1 =  numpy.vstack((RNL_H1_Old, self.RNL_H1))
	    self.RNL_H3 =  numpy.vstack((RNL_H3_Old, self.RNL_H3))
	    self.RNL_B =  numpy.vstack((RNL_B_Old, self.RNL_B))
	    
	    #NH
	    self.SNH = numpy.vstack((SNH_Old, self.SNH))
            self.INH_H1 = numpy.vstack((INH_H1_Old, self.INH_H1))
	    self.INH_H3 = numpy.vstack((INH_H3_Old, self.INH_H3))
	    self.INH_B = numpy.vstack((INH_B_Old, self.INH_B))
	    self.RNH_H1 =  numpy.vstack((RNH_H1_Old, self.RNH_H1))
	    self.RNH_H3 =  numpy.vstack((RNH_H3_Old, self.RNH_H3))
	    self.RNH_B =  numpy.vstack((RNH_B_Old, self.RNH_B))
	    
	    ## vaccine doses
	    self.vacc_NL = numpy.vstack((vacc_NL_old, self.vacc_NL))
	    self.vacc_NH = numpy.vstack((vacc_NL_old, self.vacc_NH))
	    self.vacc_TL = numpy.vstack((vacc_NL_old, self.vacc_TL))
	    self.vacc_TH = numpy.vstack((vacc_NL_old, self.vacc_TH))
	    
	    
        self.hasSolution = True

    def updateStats(self):
	    
	self.NUL = self.SUL +  self.IUL_H1 + self.IUL_H3 + self.IUL_B + self.RUL_H1 +  self.RUL_H3 +  self.RUL_B  
	self.NUH = self.SUH +  self.IUH_H1 + self.IUH_H3 + self.IUH_B + self.RUH_H1 +  self.RUH_H3 +  self.RUH_B
	self.NTL = self.STL +  self.ITL_H1 + self.ITL_H3 + self.ITL_B + self.RTL_H1 +  self.RTL_H3 +  self.RTL_B  
	self.NTH = self.STH +  self.ITH_H1 + self.ITH_H3 + self.ITH_B + self.RTH_H1 +  self.RTH_H3 +  self.RTH_B  
	self.NNL = self.SNL +  self.INL_H1 + self.INL_H3 + self.INL_B + self.RNL_H1 +  self.RNL_H3 +  self.RNL_B  
	self.NNH = self.SNH +  self.INH_H1 + self.INH_H3 + self.INH_B + self.RNH_H1 +  self.RNH_H3 +  self.RNH_B  
	
	self.NU = self.NUL + self.NUH
	self.NT = self.NTL + self.NTH
	self.NN = self.NNL + self.NNH
	
        self.N  = self.NU + self.NT + self.NN
	
	self.infectionsUL_H1 =  self.RUL_H1[-1,: ] +  self.IUL_H1[-1,: ]
	self.infectionsUL_H3 =  self.RUL_H3[-1,: ] +  self.IUL_H3[-1,: ]
	self.infectionsUL_B  =  self.RUL_B[-1,: ]  +  self.IUL_B[-1,: ]
	
	self.infectionsVL_H1 =  self.RTL_H1[-1,: ] + self.RNL_H1[-1,: ] + self.ITL_H1[-1,: ] + self.INL_H1[-1,: ]
	self.infectionsVL_H3 =  self.RTL_H3[-1,: ] + self.RNL_H3[-1,: ] + self.ITL_H3[-1,: ] + self.INL_H3[-1,: ]
	self.infectionsVL_B  =  self.RTL_B[-1,: ]  + self.RNL_B[-1,: ]  + self.ITL_B[-1,: ] + self.INL_B[-1,: ]
	
	self.infectionsUH_H1 =  self.RUH_H1[-1,: ] + self.IUH_H1[-1,: ]
	self.infectionsUH_H3 =  self.RUH_H3[-1,: ] + self.IUH_H3[-1,: ]
	self.infectionsUH_B  =  self.RUH_B[-1,: ]  + self.IUH_B[-1,: ]
	
	self.infectionsVH_H1 =  self.RTH_H1[-1,: ] + self.RNH_H1[-1,: ] + self.ITH_H1[-1,: ] + self.INH_H1[-1,: ]
	self.infectionsVH_H3 =  self.RTH_H3[-1,: ] + self.RNH_H3[-1,: ] + self.ITH_H3[-1,: ] + self.INH_H3[-1,: ]
	self.infectionsVH_B  =  self.RTH_B[-1,: ]  + self.RNH_B[-1,: ]  + self.ITH_B[-1,: ] + self.INH_B[-1,: ]


	self.infectionsL_H1 = self.infectionsUL_H1 + self.infectionsVL_H1
	self.infectionsH_H1 = self.infectionsUH_H1 + self.infectionsVH_H1
	self.infectionsL_H3 = self.infectionsUL_H3 + self.infectionsVL_H3
	self.infectionsH_H3 = self.infectionsUH_H3 + self.infectionsVH_H3
	self.infectionsL_B = self.infectionsUL_B + self.infectionsVL_B
	self.infectionsH_B = self.infectionsUH_B + self.infectionsVH_B
	
	self.infections_H1 = self.infectionsUL_H1 + self.infectionsUH_H1 + self.infectionsVL_H1 + self.infectionsVH_H1 
	self.infections_H3 = self.infectionsUL_H3 + self.infectionsUH_H3 + self.infectionsVL_H3 + self.infectionsVH_H3 
	self.infections_B = self.infectionsUL_B + self.infectionsUH_B + self.infectionsVL_B + self.infectionsVH_B 
	
	self.infectionsL  = self.infectionsUL_H1 + self.infectionsUL_H3 + self.infectionsUL_B  + self.infectionsVL_H1 + self.infectionsVL_H3 + self.infectionsVL_B
        self.infectionsH  = self.infectionsUH_H1 + self.infectionsUH_H3 + self.infectionsUH_B  + self.infectionsVH_H1 + self.infectionsVH_H3 + self.infectionsVH_B

	self.infectionsU =  self.infectionsUL_H1 + self.infectionsUL_H3 + self.infectionsUL_B  + self.infectionsUH_H1 + self.infectionsUH_H3 + self.infectionsUH_B
	self.infectionsV  = self.infectionsVL_H1 + self.infectionsVL_H3 + self.infectionsVL_B  + self.infectionsVH_H1 + self.infectionsVH_H3 + self.infectionsVH_B
	self.infectionsUL =  self.infectionsUL_H1 + self.infectionsUL_H3 + self.infectionsUL_B
	self.infectionsUH = self.infectionsUH_H1 + self.infectionsUH_H3 + self.infectionsUH_B
	self.infectionsVL =  self.infectionsVL_H1 + self.infectionsVL_H3 + self.infectionsVL_B
	self.infectionsVH = self.infectionsVH_H1 + self.infectionsVH_H3 + self.infectionsVH_B
	
	self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()
	
	self.vaccine_doses_T = (self.vacc_TL[-1,: ] + self.vacc_TH[-1,: ])
	self.vaccine_doses_N = (self.vacc_NL[-1,: ] + self.vacc_NH[-1,: ])
	self.vaccine_doses_total = self.vaccine_doses_T + self.vaccine_doses_N
	
	#################################
	self._RR_H3 = (self.parameters.lowRiskhospitalizationRate_H3 + self.parameters.highRiskhospitalizationRate_H3)/(1.*(self.parameters.lowRiskhospitalizationRate_H1 + self.parameters.highRiskhospitalizationRate_H1))
	    
	self._RR_B =  (self.parameters.lowRiskhospitalizationRate_B + self.parameters.highRiskhospitalizationRate_B)/(1.*(self.parameters.lowRiskhospitalizationRate_H1 + self.parameters.highRiskhospitalizationRate_H1))
	
	
	self._proportion_infections_H1 = self.infections_H1/(1.*(self.infections_H1+ self.infections_H3+ self.infections_B)) 
	self._proportion_infections_H3 = self.infections_H3/(1.*(self.infections_H1+ self.infections_H3+ self.infections_B)) 
	self._proportion_infections_B  = self.infections_B/(1.*(self.infections_H1+ self.infections_H3+ self.infections_B)) 
	
	
	self._prop_vaccinated = self.infectionsV/self.infections
	self._prop_unvaccinated = 1 - self._prop_vaccinated
	
	
	self.prob_hosp = self.parameters.prob_hosp_scaling * self.parameters.relative_prob_hosp
	
	self.prob_hospU_H1 = self.prob_hosp/(self._proportion_infections_H1+ self._RR_H3*self._proportion_infections_H3 + self._RR_B*self._proportion_infections_B)
	self.prob_hospU_H3 = self._RR_H3 * self.prob_hospU_H1
	self.prob_hospU_B = self._RR_B * self.prob_hospU_H1
	
	self.prob_hospV_H1 = self.prob_hospU_H1*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H1,1))
	self.prob_hospV_H3 = self.prob_hospU_H3*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H3,1)) 
	self.prob_hospV_B = self.prob_hospU_B*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_B,1))
	
	
	self.ratio_hosp_highriskU_H1 = self.parameters.ratio_hosp_highrisk_H1/(self._prop_unvaccinated + (1-numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H1,1))*self._prop_vaccinated) 
	self.ratio_hosp_highriskU_H3 = self.parameters.ratio_hosp_highrisk_H3/(self._prop_unvaccinated + (1 - numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H3,1))*self._prop_vaccinated) 
	self.ratio_hosp_highriskU_B = self.parameters.ratio_hosp_highrisk_B/(self._prop_unvaccinated + (1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_B,1))*self._prop_vaccinated) 
	
	self.ratio_hosp_highriskV_H1 = self.ratio_hosp_highriskU_H1*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H1,1)) 
	self.ratio_hosp_highriskV_H3 = self.ratio_hosp_highriskU_H3*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_H3,1)) 
	self.ratio_hosp_highriskV_B = self.ratio_hosp_highriskU_B*(1- numpy.minimum(self.parameters.vac_eff_hospitalization*self.parameters.relative_vaccineEfficacyVsHospitalization_B,1)) 
	

	self.case_hospitalizationUL_H1 = self.prob_hospU_H1
	self.case_hospitalizationUH_H1 = self.ratio_hosp_highriskU_H1 * self.case_hospitalizationUL_H1
	

	self.case_hospitalizationUL_H3 = self.prob_hospU_H3
	self.case_hospitalizationUH_H3 = self.ratio_hosp_highriskU_H3 * self.case_hospitalizationUL_H3
	
	self.case_hospitalizationUL_B = self.prob_hospU_B
	self.case_hospitalizationUH_B = self.ratio_hosp_highriskU_B * self.case_hospitalizationUL_B
	
	self.hospitalizationsUL_H1 = self.infectionsUL_H1 * self.case_hospitalizationUL_H1 
	self.hospitalizationsUL_H3 = self.infectionsUL_H3 * self.case_hospitalizationUL_H3 
	self.hospitalizationsUL_B = self.infectionsUL_B * self.case_hospitalizationUL_B 
	
	self.hospitalizationsUH_H1 = self.infectionsUH_H1 * self.case_hospitalizationUH_H1 
	self.hospitalizationsUH_H3 = self.infectionsUH_H3 * self.case_hospitalizationUH_H3 
	self.hospitalizationsUH_B = self.infectionsUH_B * self.case_hospitalizationUH_B
	
	################################################################################
	self.case_hospitalizationVL_H1 = self.prob_hospV_H1
	self.case_hospitalizationVH_H1 =  self.ratio_hosp_highriskV_H1 * self.case_hospitalizationVL_H1
	
	self.case_hospitalizationVL_H3 = self.prob_hospV_H3
	self.case_hospitalizationVH_H3 = self.ratio_hosp_highriskV_H3* self.case_hospitalizationVL_H3
	
	self.case_hospitalizationVL_B = self.prob_hospV_B
	self.case_hospitalizationVH_B =  self.ratio_hosp_highriskV_B * self.case_hospitalizationVL_B
	
	
	self.hospitalizationsVL_H1 =  self.infectionsVL_H1  * self.case_hospitalizationVL_H1
	self.hospitalizationsVL_H3 = self.infectionsVL_H3  * self.case_hospitalizationVL_H3
	self.hospitalizationsVL_B = self.infectionsVL_B  * self.case_hospitalizationVL_B
	
	self.hospitalizationsVH_H1 = self.infectionsVH_H1 * self.case_hospitalizationVH_H1
	self.hospitalizationsVH_H3 = self.infectionsVH_H3  * self.case_hospitalizationVH_H3
	self.hospitalizationsVH_B = self.infectionsVH_B * self.case_hospitalizationVH_B
	

	
	self.hospitalizationsL_H1 = self.hospitalizationsUL_H1 + self.hospitalizationsVL_H1
	self.hospitalizationsH_H1 = self.hospitalizationsUH_H1 + self.hospitalizationsVH_H1
	self.hospitalizationsL_H3 = self.hospitalizationsUL_H3 + self.hospitalizationsVL_H3
	self.hospitalizationsH_H3 = self.hospitalizationsUH_H3 + self.hospitalizationsVH_H3
	self.hospitalizationsL_B = self.hospitalizationsUL_B + self.hospitalizationsVL_B
	self.hospitalizationsH_B = self.hospitalizationsUH_B + self.hospitalizationsVH_B
	
	self.hospitalizations_H1 = self.hospitalizationsUL_H1 + self.hospitalizationsUH_H1 + self.hospitalizationsVL_H1 + self.hospitalizationsVH_H1
	self.hospitalizations_H3 = self.hospitalizationsUL_H3 + self.hospitalizationsUH_H3 + self.hospitalizationsVL_H3 + self.hospitalizationsVH_H3
	self.hospitalizations_B = self.hospitalizationsUL_B + self.hospitalizationsUH_B + self.hospitalizationsVL_B + self.hospitalizationsVH_B
	
	self.hospitalizationsL  = self.hospitalizationsUL_H1 + self.hospitalizationsUL_H3 + self.hospitalizationsUL_B + self.hospitalizationsVL_H1 + self.hospitalizationsVL_H3 + self.hospitalizationsVL_B
	self.hospitalizationsH  = self.hospitalizationsUH_H1 + self.hospitalizationsUH_H3 + self.hospitalizationsUH_B + self.hospitalizationsVH_H1 + self.hospitalizationsVH_H3 + self.hospitalizationsVH_B
	self.hospitalizationsU = self.hospitalizationsUL_H1 + self.hospitalizationsUH_H1 + self.hospitalizationsUL_H3 + self.hospitalizationsUH_H3 + self.hospitalizationsUL_B + self.hospitalizationsUH_B
	self.hospitalizationsV = self.hospitalizationsVL_H1 + self.hospitalizationsVH_H1 + self.hospitalizationsVL_H3 + self.hospitalizationsVH_H3 + self.hospitalizationsVL_B + self.hospitalizationsVH_B
	self.hospitalizationsUH =  self.hospitalizationsUH_H1 + self.hospitalizationsUH_H3 + self.hospitalizationsUH_B
	self.hospitalizationsUL =  self.hospitalizationsUL_H1 + self.hospitalizationsUL_H3 + self.hospitalizationsUL_B
	self.hospitalizationsVH =  self.hospitalizationsVH_H1 + self.hospitalizationsVH_H3 + self.hospitalizationsVH_B
	self.hospitalizationsVL =  self.hospitalizationsVL_H1 + self.hospitalizationsVL_H3 + self.hospitalizationsVL_B
	
	self.hospitalizations = self.hospitalizationsL + self.hospitalizationsH
	self.totalHospitalizations = self.hospitalizations.sum()
	
	#######################################################
	self.prob_death = self.parameters.prob_death_scaling * self.parameters.relative_prob_death
	
	self.prob_deathU_B = self.prob_death/(self._proportion_infections_B + self._proportion_infections_H1*self.parameters.ratio_death_strain_H1 + self._proportion_infections_H3*self.parameters.ratio_death_strain_H3) 
	self.prob_deathU_H1 = self.parameters.ratio_death_strain_H1 * self.prob_deathU_B
	self.prob_deathU_H3 = self.parameters.ratio_death_strain_H3 * self.prob_deathU_B
	
	self.prob_deathV_H1 = self.prob_deathU_H1 * (1- numpy.minimum(self.parameters.vac_eff_mortality * self.parameters.relative_vaccineEfficacyVsDeath_H1,1)) 
	self.prob_deathV_H3 = self.prob_deathU_H3 * (1- numpy.minimum(self.parameters.vac_eff_mortality * self.parameters.relative_vaccineEfficacyVsDeath_H3,1)) 
	self.prob_deathV_B = self.prob_deathU_B * (1- numpy.minimum(self.parameters.vac_eff_mortality * self.parameters.relative_vaccineEfficacyVsDeath_B,1))
	
	
	self.ratio_death_highriskU_H1 = self.parameters.ratio_death_highrisk_H1/(self._prop_unvaccinated + (1 - numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_H1,1))*self._prop_vaccinated) 
	self.ratio_death_highriskU_H3 = self.parameters.ratio_death_highrisk_H3/(self._prop_unvaccinated + (1 - numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_H3,1))*self._prop_vaccinated) 
	self.ratio_death_highriskU_B = self.parameters.ratio_death_highrisk_B/(self._prop_unvaccinated + (1 - numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_B,1))*self._prop_vaccinated)
	
	self.ratio_death_highriskV_H1 = self.ratio_death_highriskU_H1*(1- numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_H1,1)) 
	self.ratio_death_highriskV_H3 = self.ratio_death_highriskU_H3*(1- numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_H3,1)) 
	self.ratio_death_highriskV_B = self.ratio_death_highriskU_B*(1- numpy.minimum(self.parameters.vac_eff_mortality*self.parameters.relative_vaccineEfficacyVsDeath_B,1)) 
	
    	
        self.deathRateUL_H1 =  self.prob_deathU_H1/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskU_H1*self.parameters.proportionHighRisk) 
	self.deathRateUL_H3 =   self.prob_deathU_H3/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskU_H3*self.parameters.proportionHighRisk)
	self.deathRateUL_B =    self.prob_deathU_B/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskU_B*self.parameters.proportionHighRisk)
	
	##Death rate of high-risk unvaccinated individuals
        self.deathRateUH_H1 = self.ratio_death_highriskU_H1 *  self.deathRateUL_H1
	self.deathRateUH_H3 = self.ratio_death_highriskU_H3 *  self.deathRateUL_H3
	self.deathRateUH_B = self.ratio_death_highriskU_B *  self.deathRateUL_B
	
	self.deathsUL_H1 =   self.infectionsUL_H1 *  self.deathRateUL_H1  
	self.deathsUL_H3 =    self.infectionsUL_H3 *  self.deathRateUL_H3  
	self.deathsUL_B  =    self.infectionsUL_B *  self.deathRateUL_B
	
	self.deathsUH_H1 =    self.infectionsUH_H1 *  self.deathRateUH_H1 
	self.deathsUH_H3 =    self.infectionsUH_H3 *  self.deathRateUH_H3 
	self.deathsUH_B  =    self.infectionsUH_B *  self.deathRateUH_B
	
	###########################################################
	
	self.deathRateVL_H1 =  self.prob_deathV_H1/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskV_H1*self.parameters.proportionHighRisk) 
	self.deathRateVL_H3 =   self.prob_deathV_H3/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskV_H3*self.parameters.proportionHighRisk)
	self.deathRateVL_B =    self.prob_deathV_B/ ((1- self.parameters.proportionHighRisk) + self.ratio_death_highriskV_B*self.parameters.proportionHighRisk)
	
	##Death rate of high-risk unvaccinated individuals
        self.deathRateVH_H1 = self.ratio_death_highriskV_H1 *  self.deathRateVL_H1
	self.deathRateVH_H3 = self.ratio_death_highriskV_H3 *  self.deathRateVL_H3
	self.deathRateVH_B = self.ratio_death_highriskV_B *  self.deathRateVL_B
	
	self.deathsVL_H1 =   self.infectionsVL_H1 *  self.deathRateVL_H1  
	self.deathsVL_H3 =    self.infectionsVL_H3 *  self.deathRateVL_H3  
	self.deathsVL_B  =    self.infectionsVL_B *  self.deathRateVL_B
	
	self.deathsVH_H1 =    self.infectionsVH_H1 *  self.deathRateVH_H1 
	self.deathsVH_H3 =    self.infectionsVH_H3 *  self.deathRateVH_H3 
	self.deathsVH_B  =    self.infectionsVH_B *  self.deathRateVH_B
	
	
	self.deaths_H1 = self.deathsUL_H1 + self.deathsVL_H1 + self.deathsVL_H1 + self.deathsVH_H1
	self.deaths_H3 = self.deathsUL_H3 + self.deathsVL_H3 + self.deathsVL_H3 + self.deathsVH_H3
	self.deaths_B = self.deathsUL_B + self.deathsVL_B + self.deathsVL_B + self.deathsVH_B 
	self.deathsUL = self.deathsUL_H1 + self.deathsUL_H3 + self.deathsUL_B
	self.deathsUH = self.deathsUH_H1 + self.deathsUH_H3 + self.deathsUH_B
	self.deathsVL = self.deathsVL_H1 + self.deathsVL_H3 + self.deathsVL_B 
	self.deathsVH = self.deathsVH_H1 + self.deathsVH_H3 + self.deathsVH_B
	self.deathsV = self.deathsVL + self.deathsVH
	self.deathsU = self.deathsUL + self.deathsUH
	self.deathsL  = self.deathsUL + self.deathsVL
	self.deathsH  = self.deathsUH + self.deathsVH
        self.deaths   = self.deathsL + self.deathsH 
        self.totalDeaths = self.deaths.sum()
	
	self._lowrisk_outpatients_H1 = self.parameters.lowRiskOutpatientProb * (self.infectionsL_H1 - self.hospitalizationsL_H1)
	self._highrisk_outpatients_H1 = self.parameters.highRiskOutpatientProb * (self.infectionsH_H1 - self.hospitalizationsH_H1)
	self.outpatients_H1 = self._lowrisk_outpatients_H1+ self._highrisk_outpatients_H1
	# not medically attended = all infection cases - outpatients - hospitalized -deaths
	self._not_medically_attended_H1 = self.infections_H1 - self.outpatients_H1 - self.hospitalizations_H1 - self.deaths_H1
	self._cost_overcounterMeds_H1 = self._not_medically_attended_H1 * self.parameters.costOverCounterMeds
	self._cost_outpatient_H1 = self.parameters.costOutpatient * (self._lowrisk_outpatients_H1 + self._highrisk_outpatients_H1)
	self._cost_hospitalization_H1 = self.parameters.costHospitalization * self.hospitalizations_H1
	self.totalCosts_H1 = self._cost_overcounterMeds_H1 + self._cost_outpatient_H1 + self._cost_hospitalization_H1
	
	self._lowrisk_outpatients_H3 = self.parameters.lowRiskOutpatientProb * (self.infectionsL_H3 - self.hospitalizationsL_H3)
	self._highrisk_outpatients_H3 = self.parameters.highRiskOutpatientProb * (self.infectionsH_H3 - self.hospitalizationsH_H3)
	self.outpatients_H3 = self._lowrisk_outpatients_H3 + self._highrisk_outpatients_H3
	self._not_medically_attended_H3 = self.infections_H3 - self.outpatients_H3 - self.hospitalizations_H3 - self.deaths_H3
	self._cost_overcounterMeds_H3 = self._not_medically_attended_H3 * self.parameters.costOverCounterMeds
	self._cost_outpatient_H3 = self.parameters.costOutpatient * (self._lowrisk_outpatients_H3 + self._highrisk_outpatients_H3)
	self._cost_hospitalization_H3 = self.parameters.costHospitalization * self.hospitalizations_H3
	self.totalCosts_H3 = self._cost_overcounterMeds_H3 + self._cost_outpatient_H3 + self._cost_hospitalization_H3
	
	self._lowrisk_outpatients_B = self.parameters.lowRiskOutpatientProb * (self.infectionsL_B - self.hospitalizationsL_B)
	self._highrisk_outpatients_B = self.parameters.highRiskOutpatientProb * (self.infectionsH_B - self.hospitalizationsH_B)
	self.outpatients_B = self._lowrisk_outpatients_B + self._highrisk_outpatients_B
	self._not_medically_attended_B = self.infections_B - self.outpatients_B - self.hospitalizations_B - self.deaths_B
	self._cost_overcounterMeds_B = self._not_medically_attended_B * self.parameters.costOverCounterMeds
	self._cost_outpatient_B = self.parameters.costOutpatient * (self._lowrisk_outpatients_B + self._highrisk_outpatients_B)
	self._cost_hospitalization_B = self.parameters.costHospitalization * self.hospitalizations_B
	self.totalCosts_B = self._cost_overcounterMeds_B + self._cost_outpatient_B + self._cost_hospitalization_B
	
    
    
  
    def doses_used(self):
	return self.vaccine_doses_T.sum(), self.vaccine_doses_N.sum(), self.vaccine_doses_total.sum()
    
    def optimization_output(self):
	return self.parameters.proportionVaccinatedTypical ,  self.parameters.proportionVaccinatedUniversal

    def vaccinated_output(self):
        return list(self.parameters.proportionVaccinatedTL), list(self.parameters.proportionVaccinatedTH),list(self.parameters.proportionVaccinatedNL), list(self.parameters.proportionVaccinatedNH), [0]+ list(self.doses_TL), [0]+ list(self.doses_TH), [0]+ list(self.doses_NL), [0]+ list(self.doses_NH)

    def short_output(self):
	return list(self.infectionsL), list(self.infectionsH),  list(self.hospitalizationsL), list(self.hospitalizationsH), list(self.deathsL), list(self.deathsH)
    
    def strain_output(self):
	return sum(list(self.infections_H1)), sum(list(self.infections_H3)), sum(list(self.infections_B)), sum(list(self.hospitalizations_H1)), sum(list(self.hospitalizations_H3)), sum(list(self.hospitalizations_B)), sum(list(self.deaths_H1)), sum(list(self.deaths_H3)), sum(list(self.deaths_B))
    
    def detailed_strain_output(self):
	return list(self.infections_H1), list(self.infections_H3), list(self.infections_B), list(self.hospitalizations_H1), list(self.hospitalizations_H3), list(self.hospitalizations_B), list(self.deaths_H1), list(self.deaths_H3), list(self.deaths_B), list(self.totalCosts_H1), list(self.totalCosts_H3), list(self.totalCosts_B)
    
    def calibration_output(self):
	
	#self.updateStats()
	unvax  = ((self.NUL[0,:].sum() + self.NUH[0,:].sum() - self.SUL[-1,:].sum() - self.SUH[-1,:].sum()))/1e6
	typical = (self.NTL[0,:].sum() + self.NTH[0,:].sum() - self.STL[-1,:].sum() - self.STH[-1,:].sum())/1e6
	universal = ((self.NNL[0,:].sum() + self.NNH[0,:].sum() - self.SNL[-1,:].sum() - self.SNH[-1,:].sum()))/1e6
	incidence =  sum(list(self.infectionsL)) + sum(list(self.infectionsH)) 

	perc_H1 =  (sum(self.infections_H1) *100.)/(1.*incidence)
	perc_H3 =  (sum(self.infections_H3) *100.)/(1.*incidence)
	perc_B =  (sum(self.infections_B) *100.)/(1.*incidence)
	
	return list(self.infectionsL), list(self.infectionsH),  (self.infections_H1),  (self.infections_H3),  (self.infections_B), perc_H1, perc_H3, perc_B, list(self.hospitalizationsL), list(self.hospitalizationsH), list(self.deathsL), list(self.deathsH)
    

    def debug_info(self):
	return self.infectionsU.sum()

