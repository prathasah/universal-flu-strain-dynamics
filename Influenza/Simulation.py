import fileIO
import random
import numpy
import sys
sys.path.insert(0, r'./Influenza/Parameters')
sys.path.insert(0, r'../Influenza/Parameters')
sys.path.insert(0, r'../../Influenza/Parameters')
import Parameters
        
class run_Simulation:
    def __init__(self, options = None, tMin=0 , tMax = 365, paramValues = {}, index=None, calibration = False):
        self.tMax = tMax
	self.tMin = tMin

        if options != None:
            self.options = options
	   
        else:
            from getOptions import getOptions
	    

        # Must wait until after options, where RNG seed is set
       # Parameters = __import__('Parameters', globals())
	
        self.parameters = Parameters.Parameters(index, calibration, **paramValues)

        # Initial condition
        self.Y0 = numpy.zeros(42 * self.parameters.ages.size)

        self.hasSolution = False
	    
    def compute_typical_vaccination(self, doses):
	
    
	    
	    # source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
	    #------------------------------------------------------------
	    #Year	|  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
	    #-----------------------------------------------------------
	    #2012-13	|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
	    #2013-14	|	70.4|	55.3   |32.3	 |45.3     |	65
	    #2014-15	|	70.4|	55.8   |33.5	 |47	   |	66.7
	    #2015-16	|	70  |	55.9   |32.7	 |43.6	   | 	63.4
	    #2016-17	|	70  | 	55.6   |33.6	 |45.4	   |	65.3
	    #----------------------------------------------------------
	    
	    ## vaccination uptake in 2016 + adjustment to ensure total vaccines = 150million
	    
	    # length of vaccine coverage = 32. First 16 is for low risk people and last sixteen is for high risk people
	    ## vaccination rate of high risk groups is 1.34 times higher than low risk groups on an average (see data)
	    vaccination_coverage_low_risk = [ 0.70,  0.56,  0.56, 0.5152, 0.34, 0.34,0.34,0.34,0.34,0.34,
				   0.46,  0.46, 0.46, 0.653, 0.653, 0.653]
	    
	    vaccination_coverage_high_risk = [1.32* num for num in vaccination_coverage_low_risk]
	    
	    vaccination_coverage = vaccination_coverage_low_risk + vaccination_coverage_high_risk

	    # assume same vaccination rates for typical and universal vaccine
	    return vaccination_coverage + vaccination_coverage
	    

    #def computeR0(self):
    #    return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SUL[-1, :], self.IUL_H1[-1, :], self.IUL_H3[-1, :], self.IUL_B[-1, :], self.RUL_H1[-1, :], self.RUL_H3[-1, :], self.RUL_B[-1, :], 
		self.SUH[-1, :], self.IUH_H1[-1, :], self.IUH_H3[-1, :], self.IUH_B[-1, :], self.RUH_H1[-1, :], self.RUH_H3[-1, :], self.RUH_B[-1, :], 
		self.STL[-1, :], self.ITL_H1[-1, :], self.ITL_H3[-1, :], self.ITL_B[-1, :], self.RTL_H1[-1, :], self.RTL_H3[-1, :], self.RTL_B[-1, :], 
		self.STH[-1, :], self.ITH_H1[-1, :], self.ITH_H3[-1, :], self.ITH_B[-1, :], self.RTH_H1[-1, :], self.RTH_H3[-1, :], self.RTH_B[-1, :],
		self.SNL[-1, :], self.INL_H1[-1, :], self.INL_H3[-1, :], self.INL_B[-1, :], self.RNL_H1[-1, :], self.RNL_H3[-1, :], self.RNL_B[-1, :],
		self.SNH[-1, :], self.INH_H1[-1, :], self.INH_H3[-1, :], self.INH_B[-1, :], self.RNH_H1[-1, :], self.RNH_H3[-1, :], self.RNH_B[-1, :])
    
    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SUL
            self.Y0[ 0: : 42] =  (1 - self.parameters.proportionVaccinatedTL -  self.parameters.proportionVaccinatedNL) * self.parameters.population_lowrisk 

	    ## SUH
            self.Y0[ 7: : 42] =  (1 - self.parameters.proportionVaccinatedTH -  self.parameters.proportionVaccinatedNH) * self.parameters.population_highrisk 
  
	    ## STL
            self.Y0[ 14: : 42] = self.parameters.proportionVaccinatedTL * self.parameters.population_lowrisk 
	    
	    ## STH
            self.Y0[ 21: : 42] =  self.parameters.proportionVaccinatedTH * self.parameters.population_highrisk

	    ## SNL
            self.Y0[ 28: : 42] = self.parameters.proportionVaccinatedNL * self.parameters.population_lowrisk 
	    ## SNH 
            self.Y0[ 35: : 42] = self.parameters.proportionVaccinatedNH  * self.parameters.population_highrisk
	
            # I: Add a single infectious person in each age
	    # IUL_H1, IUL_H3, IUL_B
	    self.Y0[ 1: : 42] = numpy.full(self.parameters.ages.size,10)
	    self.Y0[ 2: : 42] = numpy.full(self.parameters.ages.size,10)
	    self.Y0[ 3: : 42] = numpy.full(self.parameters.ages.size,10)
	    
	    # IUH_H1, IUH_H3, IUH_B
            self.Y0[ 8: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 9: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 10: : 42] = numpy.full(self.parameters.ages.size, 10)
	    
	    # ITL_H1, ITL_H3, ITL_B
	    self.Y0[ 15: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 16: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 17: : 42] = numpy.full(self.parameters.ages.size, 10)
	    
	    # ITH_H1, ITH_H3, ITH_B
	    self.Y0[ 22: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 23: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 24: : 42] = numpy.full(self.parameters.ages.size, 10)
	    
	    # INL_H1, INL_H3, INL_B
	    self.Y0[ 29: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 30: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 31: : 42] = numpy.full(self.parameters.ages.size, 10)
	    
	    # INH_H1, INH_H3, INH_B
	    self.Y0[ 36: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 37: : 42] = numpy.full(self.parameters.ages.size, 10)
	    self.Y0[ 38: : 42] = numpy.full(self.parameters.ages.size, 10)

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0: : 42] -= (self.Y0[ 1: : 42] + self.Y0[ 2: : 42] + self.Y0[ 3: : 42])
	    self.Y0[ 7: : 42] -= (self.Y0[ 8: : 42] + self.Y0[ 9: : 42] + self.Y0[ 10: : 42])
	    self.Y0[ 14: : 42] -= (self.Y0[ 15: : 42] + self.Y0[ 16: : 42] + self.Y0[ 17: : 42])
	    self.Y0[ 21: : 42] -= (self.Y0[ 22: : 42] + self.Y0[ 23: : 42] + self.Y0[ 24: : 42])
	    self.Y0[ 28: : 42] -= (self.Y0[ 29: : 42] + self.Y0[ 30: : 42] + self.Y0[ 31: : 42])
	    self.Y0[ 35: : 42] -= (self.Y0[ 32: : 42] + self.Y0[ 33: : 42] + self.Y0[ 34: : 42])
 
            # R
	    self.Y0[ 4:  : 42] = 0.
	    self.Y0[ 5:  : 42] = 0.
	    self.Y0[ 6:  : 42] = 0.
	    self.Y0[ 11:  : 42] = 0.
	    self.Y0[ 12:  : 42] = 0.
	    self.Y0[ 13:  : 42] = 0.
	    self.Y0[ 18:  : 42] = 0.
	    self.Y0[ 19:  : 42] = 0.
	    self.Y0[ 20:  : 42] = 0.
	    self.Y0[ 25:  : 42] = 0.
	    self.Y0[ 26:  : 42] = 0.
	    self.Y0[ 27:  : 42] = 0.
	    self.Y0[ 32:  : 42] = 0.
	    self.Y0[ 33:  : 42] = 0.
	    self.Y0[ 34:  : 42] = 0.
	    self.Y0[ 39:  : 42] = 0.
	    self.Y0[ 40:  : 42] = 0.
	    self.Y0[ 41:  : 42] = 0.
	    
	    
        else:
	    SUL, IUL_H1, IUL_H3, IUL_B, RUL_H1, RUL_H3, RUL_B,
	    SUH, IUH_H1, IUH_H3, IUH_B, RUH_H1, RUH_H3, RUH_B,
	    STL, ITL_H1, ITL_H3, ITL_B, RTL_H1, RTL_H3, RTL_B,  
	    STH, ITH_H1, ITH_H3, ITH_B, RTH_H1, RTH_H3, RTH_B,
	    SNL, INL_H1, INL_H3, INL_B, RNL_H1, RNL_H3, RNL_B, 
	    SNH, INH_H1, INH_H3, INH_B, RNH_H1, RNH_H3, RNH_B = self.getLastValues()
	    
            self.Y0[ 0 : : 42] = (1 - self.parameters.proportionVaccinatedL) * SUL
	    self.Y0[ 7 : : 42] = (1 - self.parameters.proportionVaccinatedH) * SUH
	    self.Y0[ 14 : : 42] = STL + (1 - self.parameters.proportionVaccinatedTL) * SUL
	    self.Y0[ 21 : : 42] = STH + (1 - self.parameters.proportionVaccinatedTH) * SUH
	    self.Y0[ 28 : : 42] = SNL + (1 - self.parameters.proportionVaccinatedNL) * SUL
	    self.Y0[ 35 : : 42] = SNH + (1 - self.parameters.proportionVaccinatedNH) * SUH
	    
	    #I
	    self.Y0[ 1 : : 42] = IUL_H1
	    self.Y0[ 2 : : 42] = IUL_H3
	    self.Y0[ 3 : : 42] = IUL_B
	    
	    self.Y0[ 8 : : 42] = IUH_H1
	    self.Y0[ 9 : : 42] = IUH_H3
	    self.Y0[ 10 : : 42] = IUH_B
	    
	    self.Y0[ 15 : : 42] = ITL_H1
	    self.Y0[ 16 : : 42] = ITL_H3
	    self.Y0[ 17 : : 42] = ITL_B
	    
	    self.Y0[ 22 : : 42] = ITH_H1
	    self.Y0[ 23 : : 42] = ITH_H3
	    self.Y0[ 24 : : 42] = ITH_B
	    
	    self.Y0[ 29 : : 42] = INL_H1
	    self.Y0[ 30 : : 42] = INL_H3
	    self.Y0[ 31 : : 42] = INL_B
	    
	    self.Y0[ 36 : : 42] = INH_H1
	    self.Y0[ 37 : : 42] = INH_H3
	    self.Y0[ 38 : : 42] = INH_B
	    
	    #R class
	    self.Y0[ 4 : : 42] = RUL_H1
	    self.Y0[ 5 : : 42] = RUL_H3
	    self.Y0[ 6 : : 42] = RUL_B
	    
	    self.Y0[ 11 : : 42] = RUH_H1
	    self.Y0[ 12 : : 42] = RUH_H3
	    self.Y0[ 13 : : 42] = RUH_B
	    
	    self.Y0[ 18 : : 42] = RTL_H1
	    self.Y0[ 19 : : 42] = RTL_H3
	    self.Y0[ 20 : : 42] = RTL_B
	    
	    self.Y0[ 25 : : 42] = RTH_H1
	    self.Y0[ 26 : : 42] = RTH_H3
	    self.Y0[ 27 : : 42] = RTH_B
	    
	    self.Y0[ 32 : : 42] = RNL_H1
	    self.Y0[ 33 : : 42] = RNL_H3
	    self.Y0[ 34 : : 42] = RNL_B
	    
	    self.Y0[ 39 : : 42] = RNH_H1
	    self.Y0[ 40 : : 42] = RNH_H3
	    self.Y0[ 41 : : 42] = RNH_B
	    
	    

    def RHS(self, Y, t):
        '''
        SIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors

	SUL    = Y[ 0 : : 42]
        IUL_H1 = Y[ 1 : : 42]
        IUL_H3 = Y[ 2 : : 42]
        IUL_B  = Y[ 3 : : 42]
        RUL_H1 = Y[ 4 : : 42]
        RUL_H3 = Y[ 5 : : 42]
        RUL_B =  Y[ 6 : : 42]
            
        SUH    = Y[ 7 : : 42]
        IUH_H1 = Y[ 8 : : 42]
        IUH_H3 = Y[ 9 : : 42]
        IUH_B  = Y[ 10 : : 42]
        RUH_H1 = Y[ 11 : : 42]
        RUH_H3 = Y[ 12 : : 42]
        RUH_B =  Y[ 13 : : 42]
            
        STL    = Y[ 14 : : 42]
        ITL_H1 = Y[ 15 : : 42]
        ITL_H3 = Y[ 16 : : 42]
        ITL_B  = Y[ 17 : : 42]
        RTL_H1 = Y[ 18 : : 42]
        RTL_H3 = Y[ 19 : : 42]
        RTL_B =  Y[ 20 : : 42]
            
        STH    = Y[ 21 : : 42]
        ITH_H1 = Y[ 22 : : 42]
        ITH_H3 = Y[ 23 : : 42]
        ITH_B  = Y[ 24 : : 42]
        RTH_H1 = Y[ 25 : : 42]
        RTH_H3 = Y[ 26 : : 42]
        RTH_B =  Y[ 27 : : 42]
            
        SNL    = Y[ 28 : : 42]
        INL_H1 = Y[ 29 : : 42]
        INL_H3 = Y[ 30 : : 42]
        INL_B  = Y[ 31 : : 42]
        RNL_H1 = Y[ 32 : : 42]
        RNL_H3 = Y[ 33 : : 42]
        RNL_B =  Y[ 34 : : 42]
            
        SNH    = Y[ 35 : : 42]
        INH_H1 = Y[ 36 : : 42]
        INH_H3 = Y[ 37 : : 42]
        INH_B  = Y[ 38 : : 42]
        RNH_H1 = Y[ 39 : : 42]
        RNH_H3 = Y[ 40 : : 42]
        RNH_B =  Y[ 41 : : 42]
	
	
        N =  sum(SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1 + RUL_H3 + RUL_B + 
                SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH_H1 + RUH_H3 + RUH_B +
                STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL_H1 + RTL_H3 + RTL_B + 
                STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH_H1 + RTH_H3 + RTH_B + 
                SNL+ INL_H1+ INL_H3+ INL_B+ RNL_H1 + RNL_H3 + RNL_B + 
                SNH+ INH_H1+ INH_H3+ INH_B+ RNH_H1 + RNH_H3 + RNH_B )
            
        N_age_specific = SUL+ IUL_H1+ IUL_H3+ IUL_B+ RUL_H1 + RUL_H3+ RUL_B + SUH+ IUH_H1+ IUH_H3+ IUH_B+ RUH_H1 + RUH_H3 + RUH_B +STL+ ITL_H1+ ITL_H3+ ITL_B+ RTL_H1 + RTL_H3 + RTL_B +  STH+ ITH_H1+ ITH_H3+ ITH_B+ RTH_H1 + RTH_H3 + RTH_B +  SNL+ INL_H1+ INL_H3+ INL_B+ RNL_H1 + RNL_H3 + RNL_B +  SNH+ INH_H1+ INH_H3+ INH_B+ RNH_H1 + RNH_H3 + RNH_B
            
          
        # The force of infection
	#print ("check!!!"), self.parameters.transmissionScaling_H1, self.parameters.transmissionScaling_H3 ,self.parameters.transmissionScaling_B 
	#print self.parameters.susceptibility_H1
	#print self.parameters.susceptibility_H3
	#print self.parameters.susceptibility_B
	
        Lambda_H1 = self.parameters.transmissionScaling_H1 * self.parameters.susceptibility_H1\
		    * numpy.dot(self.parameters.contactMatrix, self.parameters.transmissibility * (IUL_H1 + IUH_H1 + ITL_H1 + ITH_H1+ INL_H1+ INH_H1)) / N_age_specific
	
	Lambda_H3 = self.parameters.transmissionScaling_H3 * self.parameters.susceptibility_H3 \
                 * numpy.dot(self.parameters.contactMatrix, 
                             self.parameters.transmissibility * (IUL_H3 + IUH_H3 + ITL_H3 + ITH_H3 + INL_H3+ INH_H3)) / N_age_specific
		
	Lambda_B = self.parameters.transmissionScaling_B * self.parameters.susceptibility_B \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IUL_B + IUH_B + ITL_B + ITH_B+ INL_B+ INH_B)) / N_age_specific
	

	
        
        # The right-hand sides
	
	#UL
        dSUL    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUL 
        dIUL_H1 = (Lambda_H1 * SUL) - (self.parameters.recoveryRate ) * IUL_H1
	dIUL_H3 = (Lambda_H3 * SUL) - (self.parameters.recoveryRate ) * IUL_H3
	dIUL_B  = (Lambda_B * SUL) - (self.parameters.recoveryRate ) * IUL_B
	dRUL_H1    = self.parameters.recoveryRate * IUL_H1
	dRUL_H3    = self.parameters.recoveryRate * IUL_H3
	dRUL_B    = self.parameters.recoveryRate * IUL_B

	
	#UH
        dSUH    = - (Lambda_H1 + Lambda_H3 + Lambda_B) * SUH 
        dIUH_H1 = (Lambda_H1 * SUH) - (self.parameters.recoveryRate ) * IUH_H1
	dIUH_H3 = (Lambda_H3 * SUH) - (self.parameters.recoveryRate ) * IUH_H3
	dIUH_B  = (Lambda_B * SUH) - (self.parameters.recoveryRate ) * IUH_B
	dRUH_H1    = self.parameters.recoveryRate * IUH_H1
	dRUH_H3    = self.parameters.recoveryRate * IUH_H3
	dRUH_B    = self.parameters.recoveryRate * IUH_B
	
	#TL
	dSTL = - ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 + (1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H3) * Lambda_H3+ (1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_B) * Lambda_B)  *STL

        dITL_H1 = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 * STL) - (self.parameters.recoveryRate ) * ITL_H1
	dITL_H3 = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H3) * Lambda_H3 * STL) - (self.parameters.recoveryRate) * ITL_H3
	dITL_B  = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_B) * Lambda_B * STL) - (self.parameters.recoveryRate) * ITL_B
	dRTL_H1    = self.parameters.recoveryRate * ITL_H1
	dRTL_H3    = self.parameters.recoveryRate * ITL_H3
	dRTL_B    = self.parameters.recoveryRate * ITL_B
	
	#TH
	dSTH =- ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H1) * Lambda_H1 + (1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H3) * Lambda_H3+ (1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_B) * Lambda_B)  *STH
	
        dITH_H1 = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H1) *Lambda_H1 * STH) - (self.parameters.recoveryRate) * ITH_H1
	dITH_H3 = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_H3) *Lambda_H3 * STH) - (self.parameters.recoveryRate) * ITH_H3
	dITH_B = ((1 - self.parameters.TypicalvaccineEfficacyVsInfectionTypical_B) * Lambda_B * STH) - (self.parameters.recoveryRate) * ITH_B
	dRTH_H1    = self.parameters.recoveryRate * ITH_H1
	dRTH_H3    = self.parameters.recoveryRate * ITH_H3
	dRTH_B    = self.parameters.recoveryRate * ITH_B
	
	#NL
	dSNL = - ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 + (1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3+ (1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_B) * Lambda_B)  *SNL
	
        dINL_H1 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNL) - (self.parameters.recoveryRate ) * INL_H1
	dINL_H3 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNL) - (self.parameters.recoveryRate ) * INL_H3
	dINL_B = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_B) * Lambda_B * SNL) - (self.parameters.recoveryRate ) * INL_B
	dRNL_H1    = self.parameters.recoveryRate * INL_H1
	dRNL_H3    = self.parameters.recoveryRate * INL_H3
	dRNL_B    = self.parameters.recoveryRate * INL_B
	
	#NH
	dSNH =  - ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H1) * Lambda_H1 + (1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H3) * Lambda_H3+ (1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_B) * Lambda_B)  *SNH
        dINH_H1 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H1) *Lambda_H1 * SNH) - (self.parameters.recoveryRate ) * INH_H1
	dINH_H3 = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_H3) *Lambda_H3 * SNH) - (self.parameters.recoveryRate ) * INH_H3
	dINH_B = ((1 - self.parameters.UniversalvaccineEfficacyVsInfectionUniversal_B) *Lambda_B * SNH) - (self.parameters.recoveryRate ) * INH_B
	dRNH_H1    = self.parameters.recoveryRate * INH_H1
	dRNH_H3    = self.parameters.recoveryRate * INH_H3
	dRNH_B    = self.parameters.recoveryRate * INH_B
	
	
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 42] = dSUL
        dY[ 1 : : 42] = dIUL_H1
        dY[ 2 : : 42] = dIUL_H3
	dY[ 3 : : 42] = dIUL_B
	dY[ 4 : : 42] = dRUL_H1
	dY[ 5 : : 42] = dRUL_H3
	dY[ 6 : : 42] = dRUL_B
	
	dY[ 7 : : 42] = dSUH
        dY[ 8 : : 42] = dIUH_H1
        dY[ 9 : : 42] = dIUH_H3
	dY[ 10 : : 42] = dIUH_B
	dY[ 11 : : 42] = dRUH_H1
	dY[ 12 : : 42] = dRUH_H3
	dY[ 13 : : 42] = dRUH_B
	
	dY[ 14 : : 42] = dSTL
        dY[ 15 : : 42] = dITL_H1
        dY[ 16 : : 42] = dITL_H3
	dY[ 17 : : 42] = dITL_B
	dY[ 18 : : 42] = dRTL_H1
	dY[ 19 : : 42] = dRTL_H3
	dY[ 20 : : 42] = dRTL_B
	
	dY[ 21 : : 42] = dSTH
        dY[ 22 : : 42] = dITH_H1
        dY[ 23 : : 42] = dITH_H3
	dY[ 24 : : 42] = dITH_B
	dY[ 25 : : 42] = dRTH_H1
	dY[ 26 : : 42] = dRTH_H3
	dY[ 27 : : 42] = dRTH_B
	
	dY[ 28 : : 42] = dSNL
        dY[ 29 : : 42] = dINL_H1
        dY[ 30 : : 42] = dINL_H3
	dY[ 31 : : 42] = dINL_B
	dY[ 32 : : 42] = dRNL_H1
	dY[ 33 : : 42] = dRNL_H3
	dY[ 34 : : 42] = dRNL_B
	
	dY[ 35 : : 42] = dSNH
        dY[ 36 : : 42] = dINH_H1
        dY[ 37 : : 42] = dINH_H3
	dY[ 38 : : 42] = dINH_B
	dY[ 39 : : 42] = dRNH_H1
	dY[ 40 : : 42] = dRNH_H3
	dY[ 41 : : 42] = dRNH_B
	
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
	
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
	
	self.SUL    = Z[:, 0 : : 42]
        self.IUL_H1 = Z[:, 1 : : 42]
	self.IUL_H3 = Z[:, 2 : : 42]
	self.IUL_B  = Z[:, 3 : : 42]
	self.RUL_H1 = Z[:, 4 : : 42]
	self.RUL_H3 = Z[:, 5 : : 42]
	self.RUL_B  = Z[:, 6 : : 42]
	
	self.SUH    = Z[:, 7 : : 42]
        self.IUH_H1 = Z[:, 8 : : 42]
	self.IUH_H3 = Z[:, 9 : : 42]
	self.IUH_B  = Z[:, 10 : : 42]
	self.RUH_H1 = Z[:, 11 : : 42]
	self.RUH_H3 = Z[:, 12 : : 42]
	self.RUH_B  = Z[:, 13 : : 42]
	
	self.STL    = Z[:, 14 : : 42]
        self.ITL_H1 = Z[:, 15 : : 42]
	self.ITL_H3 = Z[:, 16 : : 42]
	self.ITL_B  = Z[:, 17 : : 42]
	self.RTL_H1 = Z[:, 18 : : 42]
	self.RTL_H3 = Z[:, 19 : : 42]
	self.RTL_B  = Z[:, 20 : : 42]
	
	self.STH    = Z[:, 21 : : 42]
        self.ITH_H1 = Z[:, 22 : : 42]
	self.ITH_H3 = Z[:, 23 : : 42]
	self.ITH_B  = Z[:, 24 : : 42]
	self.RTH_H1 = Z[:, 25 : : 42]
	self.RTH_H3 = Z[:, 26 : : 42]
	self.RTH_B  = Z[:, 27 : : 42]
	
	self.SNL    = Z[:, 28 : : 42]
        self.INL_H1 = Z[:, 29 : : 42]
	self.INL_H3 = Z[:, 30 : : 42]
	self.INL_B  = Z[:, 31 : : 42]
	self.RNL_H1 = Z[:, 32 : : 42]
	self.RNL_H3 = Z[:, 33 : : 42]
	self.RNL_B  = Z[:, 34 : : 42]
	
	self.SNH    = Z[:, 35 : : 42]
        self.INH_H1 = Z[:, 36 : : 42]
	self.INH_H3 = Z[:, 37 : : 42]
	self.INH_B  = Z[:, 38 : : 42]
	self.RNH_H1 = Z[:, 39 : : 42]
	self.RNH_H3 = Z[:, 40 : : 42]
	self.RNH_B  = Z[:, 41 : : 42]
	
	
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

        #self.infectionsUL = self.NUL[0,:] - self.SUL[-1,:]
	#self.infectionsUH = self.NUH[0,:] - self.SUH[-1,:]
	#self.infectionsTL = self.NTL[0,:] - self.STL[-1,:]
	#self.infectionsTH = self.NTH[0,:] - self.STH[-1,:]
	#self.infectionsNL = self.NNL[0,:] - self.SNL[-1,:]
	#self.infectionsNH = self.NNH[0,:] - self.SNH[-1,:]
	
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
	
	#################################
	self._RR_H3 = [(a+b)/(1.*(c+d)) for (a,b,c,d) in zip(self.parameters.lowRiskhospitalizationRate_H3,self.parameters.highRiskhospitalizationRate_H3, self.parameters.lowRiskhospitalizationRate_H1,self.parameters.highRiskhospitalizationRate_H1)]
	self._RR_B = [(a+b)/(1.*(c+d)) for (a,b,c,d) in zip(self.parameters.lowRiskhospitalizationRate_B,self.parameters.highRiskhospitalizationRate_B, self.parameters.lowRiskhospitalizationRate_H1,self.parameters.highRiskhospitalizationRate_H1)]
	self._proportion_infections_H1 = [a/(1.*(a+b+c)) for (a,b,c) in zip(self.infections_H1, self.infections_H3, self.infections_B)]
	self._proportion_infections_H3 = [b/(1.*(a+b+c)) for (a,b,c) in zip(self.infections_H1, self.infections_H3, self.infections_B)]
	self._proportion_infections_B = [c/(1.*(a+b+c)) for (a,b,c) in zip(self.infections_H1, self.infections_H3, self.infections_B)]
	
	
	self._prop_vaccinated = self.infectionsV/self.infections
	self._prop_unvaccinated = 1 - self._prop_vaccinated
	
	
	self.prob_hosp = self.parameters.prob_hosp_scaling * self.parameters.relative_prob_hosp
	
	self.prob_hospU_H1 = [a/(b+ c*d + e*f) for (a,b,c,d,e,f) in zip(self.prob_hosp, self._proportion_infections_H1,self._RR_H3, self._proportion_infections_H3, self._RR_B,self._proportion_infections_B)]
	self.prob_hospU_H3 = [(a*b) for (a,b) in zip(self._RR_H3, self.prob_hospU_H1)]
	self.prob_hospU_B = [(a*b) for (a,b) in zip(self._RR_B, self.prob_hospU_H1)]
	
	self.prob_hospV_H1 = [a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.prob_hospU_H1, self.parameters.relative_vaccineEfficacyVsHospitalization_H1)]
	self.prob_hospV_H3 = [a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.prob_hospU_H3,  self.parameters.relative_vaccineEfficacyVsHospitalization_H3)]
	self.prob_hospV_B =[a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.prob_hospU_B, self.parameters.relative_vaccineEfficacyVsHospitalization_B)]
	
	
	self.ratio_hosp_highriskU_H1 = [a/(b + (1-min(self.parameters.vac_eff_hospitalization*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_hosp_highrisk_H1, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsHospitalization_H1,self._prop_vaccinated)]
	self.ratio_hosp_highriskU_H3 = [a/(b + (1-min(self.parameters.vac_eff_hospitalization*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_hosp_highrisk_H3, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsHospitalization_H3,self._prop_vaccinated)]
	self.ratio_hosp_highriskU_B = [a/(b + (1-min(self.parameters.vac_eff_hospitalization*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_hosp_highrisk_B, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsHospitalization_B,self._prop_vaccinated)]
	
	self.ratio_hosp_highriskV_H1 = [a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.ratio_hosp_highriskU_H1, self.parameters.relative_vaccineEfficacyVsHospitalization_H1)]
	self.ratio_hosp_highriskV_H3 = [a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.ratio_hosp_highriskU_H3, self.parameters.relative_vaccineEfficacyVsHospitalization_H3)]
	self.ratio_hosp_highriskV_B = [a*(1- min(self.parameters.vac_eff_hospitalization*c,1)) for (a,c) in zip(self.ratio_hosp_highriskU_B, self.parameters.relative_vaccineEfficacyVsHospitalization_B)]
	

	self.case_hospitalizationUL_H1 = self.prob_hospU_H1
	self.case_hospitalizationUH_H1 =  [(a*b) for (a,b) in zip(self.ratio_hosp_highriskU_H1 , self.case_hospitalizationUL_H1)]
	

	self.case_hospitalizationUL_H3 = self.prob_hospU_H3
	self.case_hospitalizationUH_H3 =   [(a*b) for (a,b) in zip(self.ratio_hosp_highriskU_H3 , self.case_hospitalizationUL_H3)]
	
	self.case_hospitalizationUL_B = self.prob_hospU_B
	self.case_hospitalizationUH_B =   [(a*b) for (a,b) in zip(self.ratio_hosp_highriskU_B , self.case_hospitalizationUL_B)]
	
	self.hospitalizationsUL_H1 = self.infectionsUL_H1 * self.case_hospitalizationUL_H1 
	self.hospitalizationsUL_H3 = self.infectionsUL_H3 * self.case_hospitalizationUL_H3 
	self.hospitalizationsUL_B = self.infectionsUL_B * self.case_hospitalizationUL_B 
	
	self.hospitalizationsUH_H1 = self.infectionsUH_H1 * self.case_hospitalizationUH_H1 
	self.hospitalizationsUH_H3 = self.infectionsUH_H3 * self.case_hospitalizationUH_H3 
	self.hospitalizationsUH_B = self.infectionsUH_B * self.case_hospitalizationUH_B
	
	################################################################################
	self.case_hospitalizationVL_H1 = self.prob_hospV_H1
	self.case_hospitalizationVH_H1 =    [(a*b) for (a,b) in zip(self.ratio_hosp_highriskV_H1, self.case_hospitalizationVL_H1)]
	
	self.case_hospitalizationVL_H3 = self.prob_hospV_H3
	self.case_hospitalizationVH_H3 =   [(a*b) for (a,b) in zip(self.ratio_hosp_highriskV_H3, self.case_hospitalizationVL_H3)]
	
	self.case_hospitalizationVL_B = self.prob_hospV_B
	self.case_hospitalizationVH_B =   [(a*b) for (a,b) in zip(self.ratio_hosp_highriskV_B, self.case_hospitalizationVL_B)]
	
	
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
	
	#print ("!!!!!!"), sum(self.hospitalizationsV)/1e3 <  sum(self.hospitalizationsU)/1e3
	#######################################################
	self.prob_death = self.parameters.prob_death_scaling * self.parameters.relative_prob_death
	
	self.prob_deathU_B = [a/(b+ c*d + e*f) for (a,b,c,d,e,f) in zip(self.prob_death, self._proportion_infections_B ,self._proportion_infections_H1,self.parameters.ratio_death_strain_H1, self._proportion_infections_H3,self.parameters.ratio_death_strain_H3)]
	self.prob_deathU_H1 = [(a*b) for (a,b) in zip(self.parameters.ratio_death_strain_H1, self.prob_deathU_B)]
	self.prob_deathU_H3 = [(a*b) for (a,b) in zip(self.parameters.ratio_death_strain_H3, self.prob_deathU_B)]
	
	self.prob_deathV_H1 = [a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.prob_deathU_H1, self.parameters.relative_vaccineEfficacyVsDeath_H1)]
	self.prob_deathV_H3 = [a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.prob_deathU_H3,  self.parameters.relative_vaccineEfficacyVsDeath_H3)]
	self.prob_deathV_B =[a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.prob_deathU_B, self.parameters.relative_vaccineEfficacyVsDeath_B)]
	
	
	self.ratio_death_highriskU_H1 = [a/(b + (1-min(self.parameters.vac_eff_mortality*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_death_highrisk_H1, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsDeath_H1,self._prop_vaccinated)]
	self.ratio_death_highriskU_H3 = [a/(b + (1-min(self.parameters.vac_eff_mortality*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_death_highrisk_H3, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsDeath_H3,self._prop_vaccinated)]
	self.ratio_death_highriskU_B = [a/(b + (1-min(self.parameters.vac_eff_mortality*d,1))*e) for (a,b,d,e) in zip(self.parameters.ratio_death_highrisk_B, self._prop_unvaccinated, self.parameters.relative_vaccineEfficacyVsDeath_B,self._prop_vaccinated)]
	
	self.ratio_death_highriskV_H1 = [a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.ratio_death_highriskU_H1, self.parameters.relative_vaccineEfficacyVsDeath_H1)]
	self.ratio_death_highriskV_H3 = [a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.ratio_death_highriskU_H3, self.parameters.relative_vaccineEfficacyVsDeath_H3)]
	self.ratio_death_highriskV_B = [a*(1- min(self.parameters.vac_eff_mortality*c,1)) for (a,c) in zip(self.ratio_death_highriskU_B, self.parameters.relative_vaccineEfficacyVsDeath_B)]
	
    	
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
	
	
	"""
	self.YLL_L = numpy.multiply(self.parameters.expectationOfLife, self.deathsL)
	self.YLL_H = numpy.multiply(self.parameters.expectationOfLife, self.deathsH)
	self.YLL = numpy.multiply(self.parameters.expectationOfLife, self.deaths)

	#Years lived with disability
	################################################################################
	### Complications cases that are hospitalized. NOTE DALY to be updated with low and high risk groups
	self._ARDS = numpy.multiply(self.infections, self.parameters.caseARDSfraction)
	self.YLD_ARDS = self.parameters.disabilityWeightARDS * (numpy.multiply(self._ARDS, self.parameters.expectationOfLife))
	
	self._pneumonia = numpy.multiply(self.infections, self.parameters.casePneumoniafraction)
	self.YLD_pneumonia = self.parameters.disabilityWeightPneumonia * (numpy.multiply(self._pneumonia, self.parameters.durationPneumonia))

	self._HospitalizedComplicated_cases = self._ARDS + self._pneumonia
	self.YLD_HospitalizedComplicated = self.YLD_ARDS + self.YLD_pneumonia
	################################################################################
	## Complications that are not hospitalized
	
	self._otitis = numpy.multiply(self.infections, self.parameters.caseOtitisfraction)
	self.YLD_otitis = self.parameters.disabilityWeightOtitis * (numpy.multiply(self._otitis, self.parameters.durationOtitis))
	
	
	self._deafness = numpy.multiply(self.infections, self.parameters.caseDeafnessfraction)
	self.YLD_deafness = numpy.multiply(self.parameters.disabilityWeightDeafness, (numpy.multiply(self._deafness, self.parameters.expectationOfLife)))
	
	self._NonhospitalizedComplicated_cases =  self._otitis + self._deafness
	self.YLD_NonhospitalizedComplicated = self.YLD_otitis + self.YLD_deafness
	################################################################################
	## Hospitalizations with no complications
	
	self._HospitalizedUncomplicated_cases = self.hospitalizations - self._HospitalizedComplicated_cases
	self.YLD_HospitalizedUncomplicated =  self.parameters.disabilityWeightHospitalizedUncomplicated * numpy.multiply(self._HospitalizedUncomplicated_cases, (1./(self.parameters.recoveryRate*365)))
	
	################################################################################
	## No hospitalizations with no complications 
	self._NonhospitalizedUncomplicated_cases = self.infections - (self._HospitalizedUncomplicated_cases + self._HospitalizedComplicated_cases + self._NonhospitalizedComplicated_cases)	
	self.YLD_NonhospitalizedUncomplicated = self.parameters.disabilityWeightUncomplicated * numpy.multiply(self._NonhospitalizedUncomplicated_cases, (1./(self.parameters.recoveryRate*365)))
	
	########################################################################
	##total YLD
	self.YLD = self.YLD_NonhospitalizedUncomplicated + self.YLD_HospitalizedUncomplicated + self.YLD_HospitalizedComplicated + self.YLD_NonhospitalizedComplicated
	

	self.DALY = (self.YLL + self.YLD)
	self.totalDALY = self.DALY.sum()
	"""

        
        
    def simulate(self):
        self.updateIC()
        self.solve()
        self.updateStats()

    def updateProportionVaccinated(self, PVPWVal, nVacTypes, vacDoses):
        # Update propotion vaccinated

	 # Convert flat vector to 2-D array
        if numpy.ndim(PVPWVal) != 2:
	  
            PVPWVal = numpy.asarray(PVPWVal).reshape(
                (nVacTypes,
                 self.parameters.proportionVaccinatedLLength + self.parameters.proportionVaccinatedHLength))

	
	## raw coverage among age groups
	proportionVaccinatedTL = PVPWVal[0][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedTH = list(PVPWVal[0][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedTL)
	proportionVaccinatedNL = PVPWVal[1][: self.parameters.proportionVaccinatedLLength]
	proportionVaccinatedNH = list(PVPWVal[1][self.parameters.proportionVaccinatedLLength:])*len(proportionVaccinatedNL)
	
	
	#raw dose uptake among age groups
	dosesVaccinatedTL =  [(a*b) for a, b in zip(proportionVaccinatedTL, self.parameters.population_lowrisk[1:])]
	dosesVaccinatedTH =  [(a*b) for a, b in zip(proportionVaccinatedTH, self.parameters.population_highrisk[1:])]
	dosesVaccinatedNL =  [(a*b) for a, b in zip(proportionVaccinatedNL, self.parameters.population_lowrisk[1:])]
	dosesVaccinatedNH =  [(a*b) for a, b in zip(proportionVaccinatedNH, self.parameters.population_highrisk[1:])]
	
	
	##coverage RELATIVE to each other (based is age-group 0 for Low risk age groups)
	relative_coverage_TL = [num/dosesVaccinatedTL[0] for num in dosesVaccinatedTL]
	relative_coverage_TH = [num/dosesVaccinatedTL[0] for num in dosesVaccinatedTH]
	relative_coverage_NL = [num/dosesVaccinatedNL[0] for num in dosesVaccinatedNL]
	relative_coverage_NH = [num/dosesVaccinatedNL[0] for num in dosesVaccinatedNH]
	
	#multiplication factor for final doses
	xFac_T = vacDoses[0]/(1.*(sum(relative_coverage_TL)+ sum(relative_coverage_TH)))
	xFac_N = vacDoses[1]/(1.*(sum(relative_coverage_NL)+ sum(relative_coverage_NH)))
	
	    
	self.doses_TL = [xFac_T * num for num in relative_coverage_TL]
	self.doses_TH = [xFac_T * num for num in relative_coverage_TH]
	self.doses_NL = [xFac_N * num for num in relative_coverage_NL]
	self.doses_NH = [xFac_N * num for num in relative_coverage_NH]
	
	    
        self.parameters.proportionVaccinatedTLPW.values = [(a/(1.*b)) for a,b in zip(self.doses_TL,  self.parameters.population_lowrisk[1:])]
	self.parameters.proportionVaccinatedTHPW.values = [(a/(1.*b)) for a,b in zip(self.doses_TH,  self.parameters.population_highrisk[1:])]
	self.parameters.proportionVaccinatedNLPW.values = [(a/(1.*b)) for a,b in zip(self.doses_NL,  self.parameters.population_lowrisk[1:])]
	self.parameters.proportionVaccinatedNHPW.values = [(a/(1.*b)) for a,b in zip(self.doses_NH,  self.parameters.population_highrisk[1:])]

	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.142
	
	
	self.parameters.proportionVaccinatedTL = self.parameters.proportionVaccinatedTLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedTH = self.parameters.proportionVaccinatedTHPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNL = self.parameters.proportionVaccinatedNLPW.full(self.parameters.ages)
	self.parameters.proportionVaccinatedNH = self.parameters.proportionVaccinatedNHPW.full(self.parameters.ages)
	
	

	vacsUsedTypical = sum(self.doses_TL + self.doses_TH)
	vacsUsedUniversal = sum(self.doses_NL + self.doses_NH)

	   
        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsedTypical, vacsUsedUniversal
    
        
    def simulateWithVaccine(self, PVPWVals, vacEfficacy_seasonal, vacEfficacy_universal, vacDoses):
	
        nVacTypes = 2

        self.resetSolution()

        # Vaccinate the population
	vacsUsedTypical, vacsUsedUniversal = self.updateProportionVaccinated(PVPWVals, nVacTypes, vacDoses)
	
	if vacsUsedUniversal <0 and min(list(PVPWVals)) >0 : print ("check!!!!"), PVPWVals, nVacTypes, vacsUsedUniversal
        tEnd = self.tMax
	tStart = self.tMin
	self.solve(tStart = tStart, tEnd = tEnd)
	
	self.updateStats()
	
	#import matplotlib.pyplot as plt
	#times = [num for num in xrange(tEnd+1)]
	#plt.plot(times, (self.IUL_H1).sum(axis=1), color = "red", linestyle= "-")
	#plt.plot(times, (self.IUL_H3).sum(axis=1), color = "blue", linestyle= "-")
	#plt.plot(times, (self.IUL_B).sum(axis=1), color = "green", linestyle = "-")	 
	#plt.show()
	
        

	return vacsUsedTypical, vacsUsedUniversal

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

