def age_specific_vaccination_coverage(season):
    

        
        # source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
        #------------------------------------------------------------
        #Year	         |  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
        #-----------------------------------------------------------'
        #2010-11        |       63.6|   46.3   | 30.5    |44.5     |    66.6
        #2012-13	|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
        #2013-14	|	70.4|	55.3   |32.3	 |45.3     |	65
        #2014-15	|	70.4|	55.8   |33.5	 |47	   |	66.7
        #2015-16	|	70  |	55.9   |32.7	 |43.6	   | 	63.4
        #2016-17	|	70  | 	55.6   |33.6	 |45.4	   |	65.3
        #2017-18        |       67.8|   54.8   | 26.9    | 39.7    |    59.6
        #average        |       68.8|   53.8   | 31.5    | 44.4    |    64.7
       #----------------------------------------------------------
       
        empirical_vax_coverage = {'2010-11': [0, 0.64, 0.46, 0.30, 0.44, 0.67], '2012-13': [0, 0.70, 0.53, 0.31, 0.45, 0.66], '2013-14':[0, 0.70, 0.55, 0.32, 0.45, 0.65], '2014-15': [0, 0.70, 0.56, 0.33, 0.47, 0.67], '2015-16': [0,0.70, 0.56, 0.33, 0.44, 0.63], '2016-17': [0, 0.70, 0.56, 0.34, 0.45, 0.65], '2017-18': [0, 0.68, 0.55, 0.27, 0.40, 0.60]}
        
        vax_coverage = PiecewiseAgeRate(empirical_vax_coverage[season],[0, 0.5, 5, 18,50,65]) 
        
        vaccination_coverage_low_risk = list(vax_coverage.full(ages.vaccinationAges))
        vaccination_coverage_high_risk = [1.34* num for num in vaccination_coverage_low_risk]
        
        vaccination_coverage = vaccination_coverage_low_risk + vaccination_coverage_high_risk
        vaccination_coverage = [round(num,2) for num in vaccination_coverage]
        
        # assume same vaccination rates for typical and universal vaccine
        return vaccination_coverage + vaccination_coverage


######################################################################33
if __name__ == "__main__":
    
    print age_specific_vaccination_coverage()