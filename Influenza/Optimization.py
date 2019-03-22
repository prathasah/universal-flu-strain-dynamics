import fileIO
import numpy
numpy.warnings.filterwarnings('ignore')

class Optimization:
    objectiveMap = {'Infections': 'totalInfections',
                    'Deaths': 'totalDeaths',
		     'Burden': 'totalDALY',
                    'Hospitalizations': 'totalHospitalizations'}
    
    detailed_objectiveMap = {'Infections': 'infections',
                    'Deaths': 'deaths',
		     'Burden': 'DALY',
                    'Hospitalizations': 'hospitalizations'}

    def __init__(self, options = None, optimRuns = 1,
                 *args, **kwargs):
        self.optimRuns = optimRuns

        if options != None:
            self.options = options
        else:
            from getOptions import getOptions
            self.options = getOptions('Optimization')
        
	
        self.vacNumbers = numpy.array([v[0] for v in self.options.vacSchedule])
	self.vacEfficacy = numpy.array([v[1] for v in self.options.vacSchedule])
        self.nVacTypes = len(self.vacNumbers)

	if list(self.vacEfficacy) != list(sorted(self.vacEfficacy)):
		raise ValueError, "First vaccine should be of lower efficacy>!"
	

        self.objective = self.options.objective 

        from Simulation import run_Simulation

        self.s = run_Simulation(options = self.options, paramValues = {"vacEfficacy":self.vacEfficacy}, *args, **kwargs)
	
        self.PVUsed = None

	self.vacsUsed = numpy.array([0] * self.nVacTypes)

    def solve(self, PVPWVals):
        # Only update for new PVPWVals
	
	if numpy.any(PVPWVals != self.PVUsed):
	    self.vacsUsedLow, self.vacsUsedHigh, self.proportionVaccinatedLow, self.proportionVaccinatedHigh = self.s.simulateWithVaccine(PVPWVals, self.vacEfficacy)
	    self.PVUsed = PVPWVals.copy()

    def evaluateObjective(self, PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (totalinfections or ....)"""

	
	self.solve(PVPWVals)
	return getattr(self.s, self.objectiveMap[self.objective])
    
    def evaluateDetailedObjective(self, PVPWVals):
	""" main objective function to minimize. Returns infection simulation instance and the objective (infections or ....)"""
	self.solve(PVPWVals)
	return getattr(self.s, self.detailed_objectiveMap[self.objective])

    def totalVacsConditions(self, PVPWVals):
	#if numpy.any(PVPWVals[:self.s.parameters.proportionVaccinatedLength]+PVPWVals[self.s.parameters.proportionVaccinatedLength:] >1.01):
	#    print PVPWVals, PVPWVals[:self.s.parameters.proportionVaccinatedLength]+PVPWVals[self.s.parameters.proportionVaccinatedLength:]
	self.solve(PVPWVals)
	## all these values should be more than zero. For convenience, return the min value
	return [self.vacNumbers[0] -  self.vacsUsedLow.sum(), self.vacNumbers[1] - self.vacsUsedHigh.sum()] 
    
    def totalVacsCondition(self):
	return [lambda PVPWVals: self.totalVacsConditions(PVPWVals)[0], lambda PVPWVals: self.totalVacsConditions(PVPWVals)[1]]
    

    
    

    def lowerCondition(self, i):
	#the min values should be greater than zero
	## return smallest numbers out of the two prop. vaccinated for the same age grup

	return lambda PVPWVals: min(PVPWVals[i],PVPWVals[i+self.s.parameters.proportionVaccinatedLength])

    def upperCondition(self, i):
	#1 - PVPWal across low and high efficacy vaccine should be greater than zero
        return lambda PVPWVals: 1.0 - (PVPWVals[i] + PVPWVals[i+self.s.parameters.proportionVaccinatedLength])

    def optimize(self):
        from scipy.optimize import fmin_cobyla
	from scipy.optimize import minimize
	

	##returns [# vax remaining, current prop. vax. for age group 1, 2, 3,4,5]
	##, (1- cureent prop. vax. for age groups 1,2,3,4,5]
	conds = []
        conds.extend(self.totalVacsCondition())
	
	conds.extend([self.lowerCondition(i) for i in
                      range(self.s.parameters.proportionVaccinatedLength)])
	
	
        conds.extend([self.upperCondition(i) for i in
                      range(self.s.parameters.proportionVaccinatedLength)])


	
        minObjective = None

        for i in range(self.optimRuns):
            if self.options.debug:
                print 'Run %d' % (i + 1)
	    ## pick random vaccination levels between 0 and 0.5
	    ## shape of array is 1D
            PV0 = (0.4 * numpy.random.rand(self.nVacTypes*
                self.s.parameters.proportionVaccinatedLength))
	    
	    #print PV0.shape, len(conds), conds, self.evaluateObjective
            PVPWValsOpt = fmin_cobyla(self.evaluateObjective,
                                      PV0,
                                     conds,
                                      maxfun = 100000,
                                      rhobeg = 0.01,
	   			      rhoend=0.0001,
				      catol= 0,
                                     disp = 0)

	    if (minObjective == None) \
                    or (self.evaluateObjective(PVPWValsOpt) < minObjective):
                if self.options.debug:
                    print 'Optimum improved.', PVPWValsOpt
                
                minObjective = self.evaluateObjective(PVPWValsOpt)
                self.PVBest = PVPWValsOpt
		self.solve(self.PVBest)
		self.totalVacsUsedLowBest = self.vacsUsedLow.sum()
		self.totalVacsUsedHighBest = self.vacsUsedHigh.sum()
		self.totalVaccinatedLowBest = self.vacsUsedLow
		self.totalVaccinatedHighBest = self.vacsUsedHigh
		self.propVaccinatedLowBest = self.proportionVaccinatedLow 
		self.propVaccinatedHighBest =  self.proportionVaccinatedHigh 
		self.objectiveBest = self.evaluateDetailedObjective(self.PVBest)
	
	self.solve(self.PVBest)

        #if self.options.write:
        #    self.dump(self.openDumpFile())

    def outputInfo(self):
	print '\nValues:\t\t',
        for PVI in self.PVBest:
            print ' %g' % round(PVI,3),
        print '\n'
        
        print 'Doses:\t\t\t %s' % ', '.join(['%g' % v for v in self.vacsUsed])
        print 'Objective:\t\t %g' % self.evaluateObjective(self.PVBest)
  
        self.s.outputInfo()

    def short_output(self):
	return self.totalVacsUsedLowBest, self.totalVacsUsedHighBest, list(self.propVaccinatedLowBest), list(self.propVaccinatedHighBest),list(self.totalVaccinatedLowBest),list(self.totalVaccinatedHighBest), list(self.objectiveBest)
	
	    
    
    def encodeVacSchedule(self):
        return '_'.join(['%g,%d,%g' % (t, int(v), float(e)) for (t, v,e)
                         in self.options.vacSchedule])

