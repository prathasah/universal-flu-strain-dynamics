import sys
sys.path.insert(0, r'./vaccine_efficacy_data/')
import csv
import numpy
import scipy.stats as stats
import random
import pandas as pd
import vaccine_effectiveness as ve
#####################################
def gamma_random_sample(mean, variance):
    """Yields a list of random numbers following a gamma distribution defined by mean and variance"""
    g_alpha = mean*mean/variance
    g_beta = mean/variance
    return random.gammavariate(g_alpha,1/g_beta)
#####################################
def get_incidence_data(year):
    
    incidence_data = {}

    incidence_data['2010-11'] = numpy.random.triangular(20e6, 21e6 , 25e6)
    incidence_data['2011-12'] = numpy.random.triangular(8.7e6, 9.3e6, 12e6)
    incidence_data['2012-13'] = numpy.random.triangular(32e6, 34e6, 38e6)
    incidence_data['2013-14'] = numpy.random.triangular(28e6, 30e6, 33e6)
    incidence_data['2014-15'] = numpy.random.triangular(29e6, 30e6, 33e6)
    incidence_data['2015-16'] = numpy.random.triangular(24e6, 25e6,28e6)
    incidence_data['2016-17'] = numpy.random.triangular(28e6, 30e6,32e6)
    incidence_data['2017-18']= numpy.random.triangular(46e6, 49e6,53e6)
    incidence_data['2018-19'] = numpy.random.uniform(25.5e6,29.3e6)
    
    return incidence_data[year]
#####################################
def get_hospitalization_data(year):
    
    
    hospitalization_data = {}
    hospitalization_data['2010-11'] = numpy.random.triangular(270e3, 290e3, 370e3)
    hospitalization_data['2011-12'] = numpy.random.triangular(130e3, 140e3, 190e3)
    hospitalization_data['2012-13'] = numpy.random.triangular(530e3, 570e3,680e3)
    hospitalization_data['2013-14'] = numpy.random.triangular(320e3, 350e3, 390e3)
    hospitalization_data['2014-15'] = numpy.random.triangular(540e3, 590e3, 680e3)
    hospitalization_data['2015-16'] = numpy.random.triangular(290e3, 310e3, 340e3)
    hospitalization_data['2016-17'] = numpy.random.triangular(520e3,580e3,660e3)
    hospitalization_data['2017-18'] = numpy.random.triangular(870e3,960e3, 1100e3)
    hospitalization_data['2018-19'] = numpy.random.uniform(327e3,394e3)
    
    return hospitalization_data[year]
    
#####################################
def get_mortality_data(year):
    
    
    mortality_data = {}
    mortality_data['2010-11'] = numpy.random.triangular(32e3,37e3,51e3)
    mortality_data['2011-12'] = numpy.random.triangular(11e3, 12e3,23e3)
    mortality_data['2012-13'] = numpy.random.triangular(37e3,43e3, 57e3)
    mortality_data['2013-14'] = numpy.random.triangular(33e3,38e3, 50e3)
    mortality_data['2014-15'] = numpy.random.triangular(44e3,51e3,64e3)
    mortality_data['2015-16'] = numpy.random.triangular(21e3,25e3,31e3)
    mortality_data['2016-17'] = numpy.random.triangular(44e3,51e3,64e3)
    mortality_data['2017-18'] = numpy.random.triangular(69e3,79e3,99e3)
    mortality_data['2018-19'] = numpy.random.uniform(21.5,35.5)
    
    return mortality_data[year]
#######################################
def get_age_virologic_profile(year):
    
    df = pd.read_csv("age_virologic_profile.csv")
    H1_0 = df.loc[(df['Season'] == year), "H1_0"].iloc[0]
    H1_5 =df.loc[(df['Season'] == year), "H1_5"].iloc[0]
    H1_25 =df.loc[(df['Season'] == year), "H1_25"].iloc[0]
    H1_65 =df.loc[(df['Season'] == year), "H1_65"].iloc[0]
    H3_0 =df.loc[(df['Season'] == year), "H3_0"].iloc[0]
    H3_5 =df.loc[(df['Season'] == year), "H3_5"].iloc[0]
    H3_25 =df.loc[(df['Season'] == year), "H3_25"].iloc[0]
    H3_65 =df.loc[(df['Season'] == year), "H3_65"].iloc[0]
    B_0 =df.loc[(df['Season'] == year), "B_0"].iloc[0]
    B_5 =df.loc[(df['Season'] == year), "B_5"].iloc[0]
    B_25 =df.loc[(df['Season'] == year), "B_25"].iloc[0]
    B_65 = df.loc[(df['Season'] == year), "B_65"].iloc[0]
    
    
    return H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65
######################################
def get_vacEfficacy(year,num):
    
    df = ve.ve_data[year]
    efficacy = {}
    for strain in ["H1N1", "H3N2", "B"]:
        efficacy[strain] = {}
    
        for age_group in ['0', '0.5-8', '9-17','18-49','50-64','65+']:
            #print ("check!!!!"), strain, age_group,len(df.loc[df['Subtype'] == strain]), df.loc[df['Subtype'] == strain, 'Age group'].iloc[0]
            #print ("check!!!!"), year, strain, age_group, strain not in df['Subtype'].unique()
                
            
            
            
            ## assume ve=0 for all age-groups if strain not present in the data-frame
            if strain not in df['Subtype'].unique():
                efficacy[strain][age_group] = 0
                
            elif (len(df.loc[df['Subtype'] == strain])==1) and (df.loc[df['Subtype'] == strain, 'Age group'].iloc[0] =="Any"):
                mean_ve = df.loc[(df['Subtype'] == strain) & (df['Age group']== "Any"), "Mean"].iloc[0]
                left_ci = max(df.loc[(df['Subtype'] == strain) & (df['Age group']== "Any"), "Low 95% CI"].iloc[0],0)
                right_ci = df.loc[(df['Subtype'] == strain) & (df['Age group']== "Any"), "High 95% CI"].iloc[0]
                efficacy[strain][age_group] = numpy.random.triangular(left_ci, mean_ve, right_ci)
                
            elif age_group=='0': efficacy[strain][age_group] = 0


            else:
                mean_ve = df.loc[(df['Subtype'] == strain) & (df['Age group']==age_group), "Mean"].iloc[0]
                left_ci = max(df.loc[(df['Subtype'] == strain) & (df['Age group']==age_group), "Low 95% CI"].iloc[0],0)
                right_ci = df.loc[(df['Subtype'] == strain) & (df['Age group']==age_group), "High 95% CI"].iloc[0]
                
                efficacy[strain][age_group] = numpy.random.triangular(left_ci, mean_ve, right_ci)
             
            
            
            #print year, strain, age_group, efficacy[strain][age_group]
                

    return efficacy

######################################
def get_vacDoses(year):
    
    doses = {}
    doses['2010-11'] = 155.1e6
    doses['2011-12'] = 132e6
    doses['2012-13'] = 134.9e6
    doses['2013-14'] = 134.5e6
    doses['2014-15'] = 147.8e6
    doses['2015-16'] = 146.4e6
    doses['2016-17'] = 145.9e6
    doses['2017-18'] = 155.3e6
    doses['2018-19'] = 169.1e6
    
    return doses[year]
######################################################################################      
if __name__ == "__main__":
    
    yearlist = ['2010-11', '2011-12', '2012-13', '2013-14', '2014-15', '2015-16', '2016-17', '2017-18']
    
    for year in yearlist:
        print ("computing....."), year
    
        header = ["iter", "year","data_incidence", "data_hospitalizations", "data_mortality", "data_vacDoses", "data_H1_0", "data_H1_5", "data_H1_25", "data_H1_65", "data_H3_0", "data_H3_5", "data_H3_25", "data_H3_65", "data_B_0", "data_B_5", "data_B_25", "data_B_65", "infectious_period_0", "infectious_period_15", "proportionHighRisk_0", "proportionHighRisk_2","proportionHighRisk_5","proportionHighRisk_19", "proportionHighRisk_25", "proportionHighRisk_50","proportionHighRisk_65", 
                  "seasonal_vaccineEfficacy_H1_0", "seasonal_vaccineEfficacy_H1_0.5",  "seasonal_vaccineEfficacy_H1_9","seasonal_vaccineEfficacy_H1_18","seasonal_vaccineEfficacy_H1_50", "seasonal_vaccineEfficacy_H1_65",
                  "seasonal_vaccineEfficacy_H3_0", "seasonal_vaccineEfficacy_H3_0.5",  "seasonal_vaccineEfficacy_H3_9","seasonal_vaccineEfficacy_H3_18","seasonal_vaccineEfficacy_H3_50", "seasonal_vaccineEfficacy_H1_65",
                  "seasonal_vaccineEfficacy_B_0", "seasonal_vaccineEfficacy_B_0.5",  "seasonal_vaccineEfficacy_B_9","seasonal_vaccineEfficacy_B_18","seasonal_vaccineEfficacy_B_50" , "seasonal_vaccineEfficacy_H1_65",
                  "age_specific_vaccineEfficacyVsInfection_0", "age_specific_vaccineEfficacyVsInfection_0.5",  "age_specific_vaccineEfficacyVsInfection_5","age_specific_vaccineEfficacyVsInfection_18","age_specific_vaccineEfficacyVsInfection_50",
              "vaccineEfficacyVsInfection_all_ages",
                   "relative_vaccineEfficacyVsHospitalization_H1_0", "relative_vaccineEfficacyVsHospitalization_H1_0.5", "relative_vaccineEfficacyVsHospitalization_H1_16", "relative_vaccineEfficacyVsHospitalization_H1_65",
                  "relative_vaccineEfficacyVsHospitalization_H3_0", "relative_vaccineEfficacyVsHospitalization_H3_0.5", "relative_vaccineEfficacyVsHospitalization_H3_16", "relative_vaccineEfficacyVsHospitalization_H3_65",
                  "relative_vaccineEfficacyVsHospitalization_B_0", "relative_vaccineEfficacyVsHospitalization_B_0.5", "relative_vaccineEfficacyVsHospitalization_B_16", "relative_vaccineEfficacyVsHospitalization_B_65",
                  "relative_vaccineEfficacyVsDeath_H1_0", "relative_vaccineEfficacyVsDeath_H1_0.5", "relative_vaccineEfficacyVsDeath_H1_18", "relative_vaccineEfficacyVsDeath_H1_65",
                  "relative_vaccineEfficacyVsDeath_H3_0", "relative_vaccineEfficacyVsDeath_H3_0.5", "relative_vaccineEfficacyVsDeath_H3_18", "relative_vaccineEfficacyVsDeath_H3_65",
                  "relative_vaccineEfficacyVsDeath_B_0", "relative_vaccineEfficacyVsDeath_B_0.5", "relative_vaccineEfficacyVsDeath_B_18", "relative_vaccineEfficacyVsDeath_B_65",
                  "relative_highRiskvaccineEfficacyVsDeath_H1_0", "relative_highRiskvaccineEfficacyVsDeath_H1_0.5", "relative_highRiskvaccineEfficacyVsDeath_H1_18", "relative_highRiskvaccineEfficacyVsDeath_H1_65",
                  "relative_highRiskvaccineEfficacyVsDeath_H3_0", "relative_highRiskvaccineEfficacyVsDeath_H3_0.5", "relative_highRiskvaccineEfficacyVsDeath_H3_18", "relative_highRiskvaccineEfficacyVsDeath_H3_65",
                  "relative_highRiskvaccineEfficacyVsDeath_B_0", "relative_highRiskvaccineEfficacyVsDeath_B_0.5", "relative_highRiskvaccineEfficacyVsDeath_B_18", "relative_highRiskvaccineEfficacyVsDeath_B_65",
                    "relative_prob_death_0", "relative_prob_death_5", "relative_prob_death_18", "relative_prob_death_50", "relative_prob_death_65",
                    "ratio_death_strain_H1_0", "ratio_death_strain_H1_5", "ratio_death_strain_H1_18", "ratio_death_strain_H1_50", "ratio_death_strain_H1_65", "ratio_death_strain_H1_75",
                    "ratio_death_strain_H3_0", "ratio_death_strain_H3_5", "ratio_death_strain_H3_18", "ratio_death_strain_H3_50", "ratio_death_strain_H3_65", "ratio_death_strain_H3_75",
                   "ratio_death_highrisk_H1_0","ratio_death_highrisk_H1_5","ratio_death_highrisk_H1_18","ratio_death_highrisk_H1_50","ratio_death_highrisk_H1_65","ratio_death_highrisk_H1_75",
                   "ratio_death_highrisk_H3_0","ratio_death_highrisk_H3_5","ratio_death_highrisk_H3_18","ratio_death_highrisk_H3_50","ratio_death_highrisk_H3_65","ratio_death_highrisk_H3_75",  
                        "ratio_death_highrisk_B_0","ratio_death_highrisk_B_5","ratio_death_highrisk_B_18","ratio_death_highrisk_B_50","ratio_death_highrisk_B_65", "ratio_death_highrisk_B_75",
                  "relative_prob_hosp_0", "relative_prob_hosp_5", "relative_prob_hosp_18", "relative_prob_hosp_50", "relative_prob_hosp_65", 
                  "ratio_hosp_highrisk_H1_0", "ratio_hosp_highrisk_H1_5", "ratio_hosp_highrisk_H1_18", "ratio_hosp_highrisk_H1_50", "ratio_hosp_highrisk_H1_65","ratio_hosp_highrisk_H1_75",
                  "ratio_hosp_highrisk_H3_0", "ratio_hosp_highrisk_H3_5", "ratio_hosp_highrisk_H3_18", "ratio_hosp_highrisk_H3_50", "ratio_hosp_highrisk_H3_65","ratio_hosp_highrisk_H3_75",
                  "ratio_hosp_highrisk_B_0", "ratio_hosp_highrisk_B_5", "ratio_hosp_highrisk_B_18", "ratio_hosp_highrisk_B_50", "ratio_hosp_highrisk_B_65","ratio_hosp_highrisk_B_75",
                  "lowRiskhospitalizationRate_H1_0","lowRiskhospitalizationRate_H1_5","lowRiskhospitalizationRate_H1_18","lowRiskhospitalizationRate_H1_50","lowRiskhospitalizationRate_H1_65", "lowRiskhospitalizationRate_H1_75",
                  "lowRiskhospitalizationRate_H3_0","lowRiskhospitalizationRate_H3_5","lowRiskhospitalizationRate_H3_18","lowRiskhospitalizationRate_H3_50","lowRiskhospitalizationRate_H3_65", "lowRiskhospitalizationRate_H3_75",
                  "lowRiskhospitalizationRate_B_0","lowRiskhospitalizationRate_B_5","lowRiskhospitalizationRate_B_18","lowRiskhospitalizationRate_B_50","lowRiskhospitalizationRate_B_65", "lowRiskhospitalizationRate_B_75",
                  "highRiskhospitalizationRate_H1_0","highRiskhospitalizationRate_H1_5","highRiskhospitalizationRate_H1_18","highRiskhospitalizationRate_H1_50","highRiskhospitalizationRate_H1_65", "highRiskhospitalizationRate_H1_75",
                  "highRiskhospitalizationRate_H3_0","highRiskhospitalizationRate_H3_5","highRiskhospitalizationRate_H3_18","highRiskhospitalizationRate_H3_50","highRiskhospitalizationRate_H3_65", "highRiskhospitalizationRate_H3_75",
                  "highRiskhospitalizationRate_B_0","highRiskhospitalizationRate_B_5","highRiskhospitalizationRate_B_18","highRiskhospitalizationRate_B_50","highRiskhospitalizationRate_B_65", "highRiskhospitalizationRate_B_75", "vac_eff_hospitalization", "vac_eff_mortality",
                  "prob_outpatient_lowrisk_0","prob_outpatient_lowrisk_5","prob_outpatient_lowrisk_18" ,"prob_outpatient_lowrisk_65",
                        "prob_outpatient_highrisk_0" ,"prob_outpatient_highrisk_5" ,"prob_outpatient_highrisk_18" ,"prob_outpatient_highrisk_65"]
        writer = csv.writer(open('sampled_parameter_1000_set_year_'+str(year)+'_1April2019.csv','wb'))
        writer.writerow(header)
         
        for num in xrange(1000):
            print year, num
            #############################
            ##data
             
            #see supplement or data folder for raw data from CDC
            incidence = get_incidence_data(year)
            hospitalizations = get_hospitalization_data(year)
            mortality = get_mortality_data(year)
            seasonal_vacEfficacy = get_vacEfficacy(year, num)
            seasonal_vacDoses = get_vacDoses(year)
            H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65 = get_age_virologic_profile(year)
     
            ######################
            # infectious period from https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.1912
            infectious_period_0 = numpy.random.triangular(2.3,3.6,5.2)
            infectious_period_15 = numpy.random.triangular(3.2,3.9,4.9)
            
            # high risk proportion from Zimmerman et al., 2010
            prop_high_risk_0 = numpy.random.normal(0.0415, 0.0044)
            prop_high_risk_2 = numpy.random.normal(0.0883, 0.0051)
            prop_high_risk_5 = numpy.random.normal(0.1168, 0.0030)
            prop_high_risk_19 = numpy.random.normal(0.1235, 0.0055)
            prop_high_risk_25 = numpy.random.normal(0.1570, 0.0027)
            prop_high_risk_50 = numpy.random.normal(0.3056, 0.0044)
            prop_high_risk_65 = numpy.random.normal(0.4701, 0.0050)
             


            ##################3
            ##vaccine efficacy against infection.
            ##ref Table 2 (complete dataset information) of Assessment of influenza vaccine effectiveness in a sentinel surveillance network 2010-13, United States
            ##Benjamin J. Cowlinga, Shuo Fenga, Lyn Finelli, Andrea Steffens, Ashley Fowlkes
             
            vac_eff_inf_0 = 0
            vac_eff_inf_6mo = numpy.random.triangular(0.5, 0.6, 0.68)
            vac_eff_inf_5 = numpy.random.triangular(0.35, 0.46, 0.56)
            vac_eff_inf_18 = numpy.random.triangular(0.26, 0.39, 0.50)
            vac_eff_inf_50 = numpy.random.triangular(0.08, 0.33, 0.51)
             
            vac_eff_inf_all_ages = numpy.random.triangular(0.43, 0.49, 0.54)
             
            #############################
            ##vaccine efficacy against hospitalization
            ## for 6mo - 15years: Table 3 of Vaccine effectiveness against laboratoryconfirmed influenza
            ## hospitalizations among young children during the 2010-11 to 2013-14 influenza seasons in Ontario, Canada
             
            ## for 16-64 and 65+ : Table 2 of
            ## Effectiveness of influenza vaccines in preventing severe influenza illness among adults: A systematic review and
            ##meta-analysis of test-negative design case-control studies
             
            relative_vac_eff_hosp_H1_0 = 0
            relative_vac_eff_hosp_H1_6mo = 1
            relative_vac_eff_hosp_H1_16 = numpy.random.triangular(0.34, 0.55, 0.76)/numpy.random.triangular(0.273, 0.821,0.956)
            relative_vac_eff_hosp_H1_65 = numpy.random.triangular(0.26, 0.54, 0.82)/numpy.random.triangular(0.273, 0.821,0.956)
             
            relative_vac_eff_hosp_H3_0 = 0
            relative_vac_eff_hosp_H3_6mo = 1
            relative_vac_eff_hosp_H3_16 = numpy.random.triangular(0.38, 0.50, 0.62)/numpy.random.triangular(0.273, 0.821,0.956)
            relative_vac_eff_hosp_H3_65 = numpy.random.triangular(0.21, 0.33, 0.45)/numpy.random.triangular(0.273, 0.821,0.956)
             
             
            relative_vac_eff_hosp_B_0 = 0
            relative_vac_eff_hosp_B_6mo = 1
            relative_vac_eff_hosp_B_16 = numpy.random.triangular(0.08, 0.45, 0.81)/numpy.random.triangular(0.273, 0.821,0.956)
            relative_vac_eff_hosp_B_65 = numpy.random.triangular(0.11, 0.31, 0.51)/numpy.random.triangular(0.273, 0.821,0.956)
            ###########################
             
            #################################
            ## vaccine efficacy against death
            ##for 6mo-17 years: Table 3 of nfluenza Vaccine Effectiveness Against Pediatric Deaths: 2010-2014; Flanerry et al.
            # for elderly: Table 4 of A Cohort Study of the Effectiveness of Influenza Vaccine in Older People, Performed Using the United Kingdom General Practice Research Database
            # for healthy adults: text in Prioritization of Influenza Pandemic Vaccination to Minimize Years of Life Lost
            # for high risk adults: Table 3 of Clinical Effectiveness of Influenza Vaccination in Persons Younger Than 65 Years With High-Risk Medical Conditions
             
            relative_vac_eff_death_H1_0 = 0
            relative_vac_eff_death_H1_6mo = 1
            relative_vac_eff_death_H1_18 = numpy.random.uniform(0.7, 0.9)/ numpy.random.triangular(0.31, 0.59,0.77)
            relative_vac_eff_death_H1_65 = numpy.random.triangular(0.11, 0.21, 0.29)/ numpy.random.triangular(0.31, 0.59,0.77)
             
            relative_vac_eff_death_H3_0 = 0
            relative_vac_eff_death_H3_6mo = 1
            relative_vac_eff_death_H3_18 = numpy.random.uniform(0.7, 0.9)/ numpy.random.triangular(0.31, 0.59,0.77)
            relative_vac_eff_death_H3_65 = numpy.random.triangular(0.11, 0.21, 0.29)/ numpy.random.triangular(0.31, 0.59,0.77)
             
            relative_vac_eff_death_B_0 = 0
            relative_vac_eff_death_B_6mo = 1
            relative_vac_eff_death_B_18 = numpy.random.uniform(0.7, 0.9)/numpy.random.triangular(0.43, 0.71,0.87)
            relative_vac_eff_death_B_65 = numpy.random.triangular(0.11, 0.21, 0.29)/numpy.random.triangular(0.43, 0.71,0.87)
             
            relative_high_risk_vac_eff_death_H1_0 = 0
            relative_high_risk_vac_eff_death_H1_6mo = 1
            relative_high_risk_vac_eff_death_H1_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.35, 0.59,0.74)
            relative_high_risk_vac_eff_death_H1_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.35, 0.59,0.74)
             
            relative_high_risk_vac_eff_death_H3_0 = 0
            relative_high_risk_vac_eff_death_H3_6mo = 1
            relative_high_risk_vac_eff_death_H3_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.35, 0.59,0.74)
            relative_high_risk_vac_eff_death_H3_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.35, 0.59,0.74)
             
            relative_high_risk_vac_eff_death_B_0 = 0
            relative_high_risk_vac_eff_death_B_6mo = 1
            relative_high_risk_vac_eff_death_B_18 = numpy.random.triangular(0.39, 0.78,0.92)/numpy.random.triangular(0.13, 0.35,0.63)
            relative_high_risk_vac_eff_death_B_65 = numpy.random.triangular(0.22, 0.29, 0.34)/numpy.random.triangular(0.13, 0.35,0.63)
            ######################
            # probability of death from https://www.sciencedirect.com/science/article/pii/S0264410X07003854?via%3Dihub and
            #https://www.sciencedirect.com/science/article/pii/S1098301516305034
            relative_prob_death_0 = 1   
            relative_prob_death_5 =  numpy.random.triangular(0.000008, 0.00001 ,0.000012)/ numpy.random.triangular(0.00002, 0.00004 ,0.00006)
            relative_prob_death_18 = numpy.random.triangular(0.00003, 0.00009 ,0.00015)/ numpy.random.triangular(0.00002, 0.00004 ,0.00006)
            relative_prob_death_50 = numpy.random.triangular(0.000458, 0.00134 ,0.00222)/ numpy.random.triangular(0.00002, 0.00004 ,0.00006)
            relative_prob_death_65 = numpy.random.triangular(0.00406, 0.0117 ,0.0193)/ numpy.random.triangular(0.00002, 0.00004 ,0.00006)
            ###########################################################
             #################
            ##case mortality
            ## rate calculated from Table 5 of
            ##Estimates of mortality attributable to influenza and RSV in the United States during 1997-2009 by influenza type or subtype, age, cause of death, and risk status
            ## Goncalo Matias,Robert Taylor,Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig, Vivek Shinde
            ##see mortality rate data sheet for calculations
             
            ratio_death_strain_H1_0 = 0.444
            ratio_death_strain_H1_5 = 0.5882
            ratio_death_strain_H1_18 = 0.1245
            ratio_death_strain_H1_50 = 0.0277
            ratio_death_strain_H1_65 = 0.0016
            ratio_death_strain_H1_75 = 0
             
            ratio_death_strain_H3_0 = 1.333
            ratio_death_strain_H3_5 = 0.8823
            ratio_death_strain_H3_18 = 1.5897
            ratio_death_strain_H3_50 = 2.7944
            ratio_death_strain_H3_65 = 3.1693
            ratio_death_strain_H3_75 = 1.4822
             
            
            ratio_death_highrisk_H1_0 = 5/11.
            ratio_death_highrisk_H1_5 = 8/12.
            ratio_death_highrisk_H1_18 = 19/15.
            ratio_death_highrisk_H1_50 = 12.
            ratio_death_highrisk_H1_65 = 1.
            ratio_death_highrisk_H1_75 = 1.
            
            
            ratio_death_highrisk_H3_0 = 12/36.
            ratio_death_highrisk_H3_5 = 11/19.
            ratio_death_highrisk_H3_18 = 265/169.
            ratio_death_highrisk_H3_50 = 967/243.
            ratio_death_highrisk_H3_65 = 1709/294.
            ratio_death_highrisk_H3_75 = 7166/2695.
            
            
            ratio_death_highrisk_B_0 = 11/25.
            ratio_death_highrisk_B_5 = 25./9.
            ratio_death_highrisk_B_18 = 154/119.
            ratio_death_highrisk_B_50 = 326/107.
            ratio_death_highrisk_B_65 = 541/91.
            ratio_death_highrisk_B_75 = 2968/990.
             
            ##################
            # probability of hospitalization from  Molinari 2007
            # https://www.sciencedirect.com/science/article/pii/S0264410X07003854 and
            # de 2016 https://www.sciencedirect.com/science/article/pii/S1098301516305034?via%3Dihub        
             
            relative_prob_hosp_0 = 1
            relative_prob_hosp_5 =  numpy.random.triangular(0.0002,0.0006,0.001)/numpy.random.triangular(0.0049, 0.0141,0.0233)
            relative_prob_hosp_18 =  numpy.random.triangular(0.0015,0.0042,0.0069)/numpy.random.triangular(0.0049, 0.0141,0.0233)
            relative_prob_hosp_50 =  numpy.random.triangular(0.00676,0.0193,0.0318)/numpy.random.triangular(0.0049, 0.0141,0.0233)
            relative_prob_hosp_65 = numpy.random.triangular(0.0147,0.0421,0.0695)/numpy.random.triangular(0.0049, 0.0141,0.0233)
             
     
            ############################################################
            #case hospitalization
            ##ref Table 4 of Estimates of hospitalization attributable to influenza and RSV in the US during 1997-2009, by age and risk status
            # Goncalo Matias, Robert Taylor, Francois Haguinet, Cynthia Schuck-Paim, Roger Lustig and Vivek Shinde
            #######################################
             
             
            low_risk_hosp_rate_H1_0 = numpy.random.triangular(0,8,23)
            low_risk_hosp_rate_H1_5 = numpy.random.triangular(0,3,9)
            low_risk_hosp_rate_H1_18 = numpy.random.triangular(0,1,3)
            low_risk_hosp_rate_H1_50 =numpy.random.triangular(0,1,1)
            low_risk_hosp_rate_H1_65 = 0
            low_risk_hosp_rate_H1_75 = 0
             
            low_risk_hosp_rate_H3_0 = numpy.random.triangular(1,61,124)
            low_risk_hosp_rate_H3_5 = numpy.random.triangular(0,8,17)
            low_risk_hosp_rate_H3_18 =numpy.random.triangular(0,10,20)
            low_risk_hosp_rate_H3_50 = numpy.random.triangular(0,17,35)
            low_risk_hosp_rate_H3_65 = numpy.random.triangular(1,39,80)
            low_risk_hosp_rate_H3_75 =numpy.random.triangular(3,117,235)
             
            low_risk_hosp_rate_B_0 = numpy.random.triangular(2,41,86)
            low_risk_hosp_rate_B_5 = numpy.random.triangular(0,8,15)
            low_risk_hosp_rate_B_18 = numpy.random.triangular(0,6,11)
            low_risk_hosp_rate_B_50 = numpy.random.triangular(0,7,15)
            low_risk_hosp_rate_B_65 = numpy.random.triangular(0,9,18)
            low_risk_hosp_rate_B_75 =numpy.random.triangular(1,40,80)
             
            high_risk_hosp_rate_H1_0 = 0
            high_risk_hosp_rate_H1_5 = numpy.random.triangular(0,3,9)
            high_risk_hosp_rate_H1_18 = numpy.random.triangular(0,9,30)
            high_risk_hosp_rate_H1_50 = numpy.random.triangular(0,10,33)
            high_risk_hosp_rate_H1_65 = numpy.random.triangular(0,20,81)
            high_risk_hosp_rate_H1_75 = numpy.random.triangular(0,21,86)
            ## the parameter realistically cannot be zero, so compute value based on highrisk/lowrisk ratio
            #low_risk_hosp_rate_H1_65 = high_risk_hosp_rate_H1_65/7.9
            #low_risk_hosp_rate_H1_75 = high_risk_hosp_rate_H1_65/4.9
             
            high_risk_hosp_rate_H3_0 = numpy.random.triangular(0,11,26)
            high_risk_hosp_rate_H3_5 = numpy.random.triangular(0,4,8)
            high_risk_hosp_rate_H3_18 = numpy.random.triangular(1,52,110)
            high_risk_hosp_rate_H3_50 =numpy.random.triangular(3,149,313)
            high_risk_hosp_rate_H3_65 = numpy.random.triangular(6,286,591)
            high_risk_hosp_rate_H3_75 = numpy.random.triangular(12,587,1198)
             
            high_risk_hosp_rate_B_0 =numpy.random.triangular(0,10,24)
            high_risk_hosp_rate_B_5 = numpy.random.triangular(0,3,7)
            high_risk_hosp_rate_B_18 =numpy.random.triangular(1,30,64)
            high_risk_hosp_rate_B_50 = numpy.random.triangular(2,60,134)
            high_risk_hosp_rate_B_65 =numpy.random.triangular(1,44,105)
            high_risk_hosp_rate_B_75 = numpy.random.triangular(4,161,380)
            
            ratio_hosp_highrisk_H1_0 =  high_risk_hosp_rate_H1_0 /min(low_risk_hosp_rate_H1_0,1)
            ratio_hosp_highrisk_H1_5 =  high_risk_hosp_rate_H1_5 /min(low_risk_hosp_rate_H1_5,1)
            ratio_hosp_highrisk_H1_18 =  high_risk_hosp_rate_H1_18 /min(low_risk_hosp_rate_H1_18,1)
            ratio_hosp_highrisk_H1_50 =  high_risk_hosp_rate_H1_50 /min(low_risk_hosp_rate_H1_50,1)
            ratio_hosp_highrisk_H1_65 =  high_risk_hosp_rate_H1_65 /max(low_risk_hosp_rate_H1_65,1)
            ratio_hosp_highrisk_H1_75 =  high_risk_hosp_rate_H1_75 /max(low_risk_hosp_rate_H1_75,1)
            
            
            ratio_hosp_highrisk_H3_0 =  high_risk_hosp_rate_H3_0 /min(low_risk_hosp_rate_H3_0,1)
            ratio_hosp_highrisk_H3_5 =  high_risk_hosp_rate_H3_5 /min(low_risk_hosp_rate_H3_5,1)
            ratio_hosp_highrisk_H3_18 =  high_risk_hosp_rate_H3_18 /min(low_risk_hosp_rate_H3_18,1)
            ratio_hosp_highrisk_H3_50 =  high_risk_hosp_rate_H3_50 /min(low_risk_hosp_rate_H3_50,1)
            ratio_hosp_highrisk_H3_65 =  high_risk_hosp_rate_H3_65 /min(low_risk_hosp_rate_H3_65,1)
            ratio_hosp_highrisk_H3_75 =  high_risk_hosp_rate_H3_75 /min(low_risk_hosp_rate_H3_75,1)
            
            
            ratio_hosp_highrisk_B_0 =  high_risk_hosp_rate_B_0 /min(low_risk_hosp_rate_B_0,1)
            ratio_hosp_highrisk_B_5 =  high_risk_hosp_rate_B_5 /min(low_risk_hosp_rate_B_5,1)
            ratio_hosp_highrisk_B_18 =  high_risk_hosp_rate_B_18 /min(low_risk_hosp_rate_B_18,1)
            ratio_hosp_highrisk_B_50 =  high_risk_hosp_rate_B_50 /min(low_risk_hosp_rate_B_50,1)
            ratio_hosp_highrisk_B_65 =  high_risk_hosp_rate_B_65 /min(low_risk_hosp_rate_B_65,1)
            ratio_hosp_highrisk_B_75 =  high_risk_hosp_rate_B_75 /min(low_risk_hosp_rate_B_75,1)
             
            vac_eff_hospitalization = 0
            vac_eff_mortality = 0
            
            #####################
            ## economic costs
            
            #source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3253949/pdf/IRV-6-167.pdf
            prob_outpatient_lowrisk_0 = numpy.random.normal(0.455, 0.098)
            prob_outpatient_lowrisk_5 = numpy.random.normal(0.318, 0.061)
            prob_outpatient_lowrisk_18 = numpy.random.normal(0.313, 0.014)
            prob_outpatient_lowrisk_65 = numpy.random.normal(0.620, 0.027)
            
            prob_outpatient_highrisk_0 = numpy.random.normal(0.910, 0.250)
            prob_outpatient_highrisk_5 = numpy.random.normal(0.635, 0.167)
            prob_outpatient_highrisk_18 = numpy.random.normal(0.625, 0.118)
            prob_outpatient_highrisk_65 = numpy.random.normal(0.850, 0.093)
 
            ##################
     
             
            elements = [num, year, incidence, hospitalizations, mortality, seasonal_vacDoses, H1_0, H1_5, H1_25, H1_65, H3_0, H3_5, H3_25, H3_65, B_0, B_5, B_25, B_65, infectious_period_0, infectious_period_15,  prop_high_risk_0, prop_high_risk_2,prop_high_risk_5,prop_high_risk_19, prop_high_risk_25, prop_high_risk_50, prop_high_risk_65,
                        seasonal_vacEfficacy['H1N1']['0'], seasonal_vacEfficacy['H1N1']['0.5-8'],seasonal_vacEfficacy['H1N1']['9-17'],seasonal_vacEfficacy['H1N1']['18-49'],seasonal_vacEfficacy['H1N1']['50-64'], seasonal_vacEfficacy['H1N1']['65+'],
                         seasonal_vacEfficacy['H3N2']['0'], seasonal_vacEfficacy['H3N2']['0.5-8'],seasonal_vacEfficacy['H3N2']['9-17'],seasonal_vacEfficacy['H3N2']['18-49'],seasonal_vacEfficacy['H3N2']['50-64'], seasonal_vacEfficacy['H3N2']['65+'],
                          seasonal_vacEfficacy['B']['0'], seasonal_vacEfficacy['B']['0.5-8'],seasonal_vacEfficacy['B']['9-17'],seasonal_vacEfficacy['B']['18-49'],seasonal_vacEfficacy['B']['50-64'], seasonal_vacEfficacy['B']['65+'],
                     vac_eff_inf_0, vac_eff_inf_6mo, vac_eff_inf_5, vac_eff_inf_18, vac_eff_inf_50,
                    vac_eff_inf_all_ages,
                        relative_vac_eff_hosp_H1_0, relative_vac_eff_hosp_H1_6mo, relative_vac_eff_hosp_H1_16, relative_vac_eff_hosp_H1_65,
                        relative_vac_eff_hosp_H3_0, relative_vac_eff_hosp_H3_6mo, relative_vac_eff_hosp_H3_16, relative_vac_eff_hosp_H3_65,
                        relative_vac_eff_hosp_B_0, relative_vac_eff_hosp_B_6mo, relative_vac_eff_hosp_B_16, relative_vac_eff_hosp_B_65,
     
                         
                        relative_vac_eff_death_H1_0, relative_vac_eff_death_H1_6mo, relative_vac_eff_death_H1_18, relative_vac_eff_death_H1_65,
                        relative_vac_eff_death_H3_0, relative_vac_eff_death_H3_6mo, relative_vac_eff_death_H3_18, relative_vac_eff_death_H3_65,
                        relative_vac_eff_death_B_0, relative_vac_eff_death_B_6mo, relative_vac_eff_death_B_18, relative_vac_eff_death_B_65,
     
                        relative_high_risk_vac_eff_death_H1_0, relative_high_risk_vac_eff_death_H1_6mo, relative_high_risk_vac_eff_death_H1_18, relative_high_risk_vac_eff_death_H1_65,
                        relative_high_risk_vac_eff_death_H3_0, relative_high_risk_vac_eff_death_H3_6mo, relative_high_risk_vac_eff_death_H3_18, relative_high_risk_vac_eff_death_H3_65,
                        relative_high_risk_vac_eff_death_B_0, relative_high_risk_vac_eff_death_B_6mo, relative_high_risk_vac_eff_death_B_18, relative_high_risk_vac_eff_death_B_65,
                         
                        relative_prob_death_0, relative_prob_death_5, relative_prob_death_18, relative_prob_death_50, relative_prob_death_65,
                    ratio_death_strain_H1_0, ratio_death_strain_H1_5, ratio_death_strain_H1_18, ratio_death_strain_H1_50, ratio_death_strain_H1_65, ratio_death_strain_H1_75,
                    ratio_death_strain_H3_0, ratio_death_strain_H3_5, ratio_death_strain_H3_18, ratio_death_strain_H3_50, ratio_death_strain_H3_65, ratio_death_strain_H3_75,
                    
                    ratio_death_highrisk_H1_0, ratio_death_highrisk_H1_5, ratio_death_highrisk_H1_18, ratio_death_highrisk_H1_50, ratio_death_highrisk_H1_65, ratio_death_highrisk_H1_75,
                    ratio_death_highrisk_H3_0, ratio_death_highrisk_H3_5, ratio_death_highrisk_H3_18, ratio_death_highrisk_H3_50, ratio_death_highrisk_H3_65, ratio_death_highrisk_H3_75,  
                         ratio_death_highrisk_B_0, ratio_death_highrisk_B_5, ratio_death_highrisk_B_18, ratio_death_highrisk_B_50, ratio_death_highrisk_B_65, ratio_death_highrisk_B_75,
                         
                         
                        relative_prob_hosp_0, relative_prob_hosp_5,  relative_prob_hosp_18,  relative_prob_hosp_50,  relative_prob_hosp_65, 
                        ratio_hosp_highrisk_H1_0, ratio_hosp_highrisk_H1_5, ratio_hosp_highrisk_H1_18, ratio_hosp_highrisk_H1_50, ratio_hosp_highrisk_H1_65,ratio_hosp_highrisk_H1_75,
                        ratio_hosp_highrisk_H3_0, ratio_hosp_highrisk_H3_5, ratio_hosp_highrisk_H3_18, ratio_hosp_highrisk_H3_50, ratio_hosp_highrisk_H3_65,ratio_hosp_highrisk_H3_75,
                        ratio_hosp_highrisk_B_0, ratio_hosp_highrisk_B_5, ratio_hosp_highrisk_B_18, ratio_hosp_highrisk_B_50, ratio_hosp_highrisk_B_65,ratio_hosp_highrisk_B_75,
                        
                        low_risk_hosp_rate_H1_0, low_risk_hosp_rate_H1_5, low_risk_hosp_rate_H1_18, low_risk_hosp_rate_H1_50, low_risk_hosp_rate_H1_65, low_risk_hosp_rate_H1_75,
                        low_risk_hosp_rate_H3_0, low_risk_hosp_rate_H3_5, low_risk_hosp_rate_H3_18, low_risk_hosp_rate_H3_50, low_risk_hosp_rate_H3_65, low_risk_hosp_rate_H3_75,
                        low_risk_hosp_rate_B_0, low_risk_hosp_rate_B_5, low_risk_hosp_rate_B_18, low_risk_hosp_rate_B_50, low_risk_hosp_rate_B_65, low_risk_hosp_rate_B_75,
                        high_risk_hosp_rate_H1_0, high_risk_hosp_rate_H1_5, high_risk_hosp_rate_H1_18, high_risk_hosp_rate_H1_50, high_risk_hosp_rate_H1_65, high_risk_hosp_rate_H1_75,
                        high_risk_hosp_rate_H3_0, high_risk_hosp_rate_H3_5, high_risk_hosp_rate_H3_18, high_risk_hosp_rate_H3_50, high_risk_hosp_rate_H3_65, high_risk_hosp_rate_H3_75,
                        high_risk_hosp_rate_B_0, high_risk_hosp_rate_B_5, high_risk_hosp_rate_B_18, high_risk_hosp_rate_B_50, high_risk_hosp_rate_B_65, high_risk_hosp_rate_B_75, vac_eff_hospitalization, vac_eff_mortality,
                        prob_outpatient_lowrisk_0 ,prob_outpatient_lowrisk_5 ,prob_outpatient_lowrisk_18 ,prob_outpatient_lowrisk_65,
                        prob_outpatient_highrisk_0 ,prob_outpatient_highrisk_5 ,prob_outpatient_highrisk_18 ,prob_outpatient_highrisk_65]
            writer.writerow(elements)
             
        