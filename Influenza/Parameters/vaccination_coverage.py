from PiecewiseAgeParameter import PiecewiseAgeNumber, PiecewiseAgeRate
import ages
import Parameters

"""
##calculation for high risk vs low risk
D(all) = D(lr) + D(hr)
cov(all)*N(all) = cov(lr)N(lr) + cov(hr)N(hr)
cov(all)*N(all) = cov(lr)N(lr) + 1.34*cov(lr)N(hr)
cov(all)*N(all) = cov(lr)[N(lr) + 1.34*N(hr)]
cov(lr) = cov(all)*N(all)  / [N(lr) + 1.34*N(hr)]

##calculation for high risk vs general population

D(all) = D(lr) + D(hr)
cov(all)*N(all) = cov(lr)N(lr) + cov(hr)N(hr)
cov(all)*N(all) = cov(lr)N(lr) + 1.34*cov(all)N(hr)
cov(lr)N(lr)  = cov*all(N(all) - 1.34*N(hr))
cov(lr)  = [cov*all( N(all) - 1.34N(hr))]/N(lr)
"""

def age_specific_vaccination_coverage(season, population, population_lowrisk, population_highrisk):
            
        
                
                # source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
                #------------------------------------------------------------
                #Year	         |  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
                #-----------------------------------------------------------'
                #2010-11        |       63.6|   46.3   | 30.5    |44.5     |    66.6
                #2011-12        |       67.6|   45.4   | 28.6    |42.7     |    64.9
                #2012-13	|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
                #2013-14	|	70.4|	55.3   |32.3	 |45.3     |	65
                #2014-15	|	70.4|	55.8   |33.5	 |47	   |	66.7
                #2015-16	|	70  |	55.9   |32.7	 |43.6	   | 	63.4
                #2016-17	|	70  | 	55.6   |33.6	 |45.4	   |	65.3
                #2017-18        |       67.8|   54.8   | 26.9    | 39.7    |    59.6
                #average        |       68.8|   53.8   | 31.5    | 44.4    |    64.7
               #----------------------------------------------------------
   
                
                empirical_vax_coverage = {'2010-11': [0, 0.64, 0.46, 0.30, 0.44, 0.67],'2011-12': [0, 0.68, 0.45,0.29, 0.43,0.65], '2012-13': [0, 0.70, 0.53, 0.31, 0.45, 0.66], '2013-14':[0, 0.70, 0.55, 0.32, 0.45, 0.65], '2014-15': [0, 0.70, 0.56, 0.33, 0.47, 0.67], '2015-16': [0,0.70, 0.56, 0.33, 0.44, 0.63], '2016-17': [0, 0.70, 0.56, 0.34, 0.45, 0.65], '2017-18': [0, 0.68, 0.55, 0.27, 0.40, 0.60]}
                
                ##2010-11, 2011-12 and and 2012-13 ratios are for high risk vs general population
                empirical_high_vs_low_coverage = {'2010-11': 1.38,'2011-12': 1.36, '2012-13': 1.32, '2013-14':1.37, '2014-15': 1.35, '2015-16': 1.37, '2016-17': 1.33, '2017-18': 1.27}
                vax_coverage_all = PiecewiseAgeRate(empirical_vax_coverage[season],[0, 0.5, 5, 18,50,65]) 
                
                vaccination_coverage_all = vax_coverage_all.full(ages.ages)
                if season in ['2010-11', '2011-12', '2012-13']:
                        #cov(lr)  = [cov*all( N(all) - 1.34N(hr))]/N(lr)
                        vaccination_coverage_low_risk = (vaccination_coverage_all*(population - empirical_high_vs_low_coverage[season]*population_highrisk))/ population_lowrisk
                else:
                        #cov(lr) = cov(all)*N(all)  / [N(lr) + 1.34*N(hr)]
                        vaccination_coverage_low_risk =  (vaccination_coverage_all *population)/(population_lowrisk + empirical_high_vs_low_coverage[season]*population_highrisk)
                        
                vaccination_coverage_high_risk = empirical_high_vs_low_coverage[season]* vaccination_coverage_low_risk
                
                
                return vaccination_coverage_low_risk, vaccination_coverage_high_risk
        
        

def reduced_age_specific_vaccination_coverage(season, agelist):
            
        
                
                # source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
                #------------------------------------------------------------
                #Year	         |  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
                #-----------------------------------------------------------'
                #2010-11        |       63.6|   46.3   | 30.5    |44.5     |    66.6
                #2011-12        |       67.6|   45.4   | 28.6    |42.7     |    64.9
                #2012-13	|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
                #2013-14	|	70.4|	55.3   |32.3	 |45.3     |	65
                #2014-15	|	70.4|	55.8   |33.5	 |47	   |	66.7
                #2015-16	|	70  |	55.9   |32.7	 |43.6	   | 	63.4
                #2016-17	|	70  | 	55.6   |33.6	 |45.4	   |	65.3
                #2017-18        |       67.8|   54.8   | 26.9    | 39.7    |    59.6
                #average        |       68.8|   53.8   | 31.5    | 44.4    |    64.7
               #----------------------------------------------------------
   
                
                empirical_vax_coverage = {'2010-11': [0, 0.64, 0.46, 0.30, 0.44, 0.67],'2011-12': [0, 0.68, 0.45,0.29, 0.43,0.65], '2012-13': [0, 0.70, 0.53, 0.31, 0.45, 0.66], '2013-14':[0, 0.70, 0.55, 0.32, 0.45, 0.65], '2014-15': [0, 0.70, 0.56, 0.33, 0.47, 0.67], '2015-16': [0,0.70, 0.56, 0.33, 0.44, 0.63], '2016-17': [0, 0.70, 0.56, 0.34, 0.45, 0.65], '2017-18': [0, 0.68, 0.55, 0.27, 0.40, 0.60]}
                
 
                vax_coverage_all = PiecewiseAgeRate(empirical_vax_coverage[season],[0, 0.5, 5, 18,50,65]) 
                
                vaccination_coverage_all = vax_coverage_all.full(agelist)
                
                return vaccination_coverage_all


######################################################################33
if __name__ == "__main__":
    
    print age_specific_vaccination_coverage("2011-12")