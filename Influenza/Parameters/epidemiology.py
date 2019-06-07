# -*- coding: iso-8859-1 -*-
#
# Epidemiological parameter values
#

from PiecewiseAgeParameter import PiecewiseAgeRate

def recoveryRatePW(df, index):
    
    return PiecewiseAgeRate(
    [1/(1.*df.at[index, "infectious_period_0"]),
     1/(1.*df.at[index, "infectious_period_15"])],
    [0, 15])


def proportionHighRiskPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "proportionHighRisk_0"],
     df.at[index, "proportionHighRisk_2"],
     df.at[index, "proportionHighRisk_5"],
     df.at[index, "proportionHighRisk_19"],
     df.at[index, "proportionHighRisk_25"],
     df.at[index, "proportionHighRisk_50"],
     df.at[index, "proportionHighRisk_65"]],
    [0, 2, 5, 19, 25, 50, 65])

def SeasonalVaccineEfficacyVsInfection_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "seasonal_vaccineEfficacy_H1_0"],
     df.at[index, "seasonal_vaccineEfficacy_H1_0.5"],
      df.at[index, "seasonal_vaccineEfficacy_H1_9"],
     df.at[index, "seasonal_vaccineEfficacy_H1_18"],
     df.at[index, "seasonal_vaccineEfficacy_H1_50"],
     df.at[index, "seasonal_vaccineEfficacy_H1_65"]],
    [0, 0.5,  5, 18, 50, 65])

def SeasonalVaccineEfficacyVsInfection_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "seasonal_vaccineEfficacy_H3_0"],
     df.at[index, "seasonal_vaccineEfficacy_H3_0.5"],
      df.at[index, "seasonal_vaccineEfficacy_H3_9"],
     df.at[index, "seasonal_vaccineEfficacy_H3_18"],
     df.at[index, "seasonal_vaccineEfficacy_H3_50"],
     df.at[index, "seasonal_vaccineEfficacy_H3_65"]],
    [0, 0.5,  5, 18, 50, 65])

def age_specific_vaccineEfficacyVsInfectionPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "age_specific_vaccineEfficacyVsInfection_0"],
     df.at[index, "age_specific_vaccineEfficacyVsInfection_0.5"],
      df.at[index, "age_specific_vaccineEfficacyVsInfection_5"],
     df.at[index, "age_specific_vaccineEfficacyVsInfection_18"],
     df.at[index, "age_specific_vaccineEfficacyVsInfection_50"]],
    [0, 0.5,  5, 18, 50])


def vaccineEfficacyVsInfection_all_ages(df, index):
    
    return df.at[index, "vaccineEfficacyVsInfection_all_ages"]

def SeasonalVaccineEfficacyVsInfection_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "seasonal_vaccineEfficacy_B_0"],
     df.at[index, "seasonal_vaccineEfficacy_B_0.5"],
      df.at[index, "seasonal_vaccineEfficacy_B_9"],
     df.at[index, "seasonal_vaccineEfficacy_B_18"],
     df.at[index, "seasonal_vaccineEfficacy_B_50"],
     df.at[index, "seasonal_vaccineEfficacy_B_65"]],
    [0, 0.5,  5, 18, 50, 65])

def relative_vaccineEfficacyVsHospitalization_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsHospitalization_H1_0"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H1_0.5"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H1_16"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H1_65"]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsHospitalization_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsHospitalization_H3_0"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H3_0.5"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H3_16"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_H3_65"]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsHospitalization_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsHospitalization_B_0"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_B_0.5"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_B_16"],
     df.at[index, "relative_vaccineEfficacyVsHospitalization_B_65"]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsDeath_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsDeath_H1_0"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H1_0.5"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H1_18"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H1_65"]],
    [0, 0.5, 18, 65])

def relative_vaccineEfficacyVsDeath_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsDeath_H3_0"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H3_0.5"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H3_18"],
     df.at[index, "relative_vaccineEfficacyVsDeath_H3_65"]],
    [0, 0.5,18, 65])

def relative_vaccineEfficacyVsDeath_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsDeath_B_0"],
     df.at[index, "relative_vaccineEfficacyVsDeath_B_0.5"],
     df.at[index, "relative_vaccineEfficacyVsDeath_B_18"],
     df.at[index, "relative_vaccineEfficacyVsDeath_B_65"]],
    [0, 0.5,18,65])


def relative_prob_deathPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_prob_death_0"],
     df.at[index, "relative_prob_death_5"],
     df.at[index, "relative_prob_death_18"],
     df.at[index, "relative_prob_death_50"],
    df.at[index, "relative_prob_death_65"]],
    [0, 5, 18, 50, 65])

def ratio_death_strain_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_death_strain_H1_0"],
     df.at[index, "ratio_death_strain_H1_5"],
     df.at[index, "ratio_death_strain_H1_18"],
     df.at[index, "ratio_death_strain_H1_50"],
    df.at[index, "ratio_death_strain_H1_65"]],
    [0, 5, 18, 50, 65])

def ratio_death_strain_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_death_strain_H3_0"],
     df.at[index, "ratio_death_strain_H3_5"],
     df.at[index, "ratio_death_strain_H3_18"],
     df.at[index, "ratio_death_strain_H3_50"],
    df.at[index, "ratio_death_strain_H3_65"]],
    [0, 5, 18, 50, 65])


def ratio_death_highrisk_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_death_highrisk_H1_0"],
     df.at[index, "ratio_death_highrisk_H1_5"],
     df.at[index, "ratio_death_highrisk_H1_18"],
     df.at[index, "ratio_death_highrisk_H1_50"],
    df.at[index, "ratio_death_highrisk_H1_65"]],
    [0, 5, 18, 50, 65])

def ratio_death_highrisk_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_death_highrisk_H3_0"],
     df.at[index, "ratio_death_highrisk_H3_5"],
     df.at[index, "ratio_death_highrisk_H3_18"],
     df.at[index, "ratio_death_highrisk_H3_50"],
    df.at[index, "ratio_death_highrisk_H3_65"]],
    [0, 5, 18, 50, 65])


def ratio_death_highrisk_BPW(df, index):
        
    return PiecewiseAgeRate(
    [df.at[index, "ratio_death_highrisk_B_0"],
     df.at[index, "ratio_death_highrisk_B_5"],
     df.at[index, "ratio_death_highrisk_B_18"],
     df.at[index, "ratio_death_highrisk_B_50"],
    df.at[index, "ratio_death_highrisk_B_65"]],
    [0, 5, 18, 50, 65])


def relative_prob_hospPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "relative_prob_hosp_0"],
     df.at[index, "relative_prob_hosp_5"],
     df.at[index, "relative_prob_hosp_18"],
     df.at[index, "relative_prob_hosp_50"],
    df.at[index, "relative_prob_hosp_65"]],
    [0, 5, 18, 50,65])

def ratio_hosp_highrisk_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_hosp_highrisk_H1_0"],
     df.at[index, "ratio_hosp_highrisk_H1_5"],
     df.at[index, "ratio_hosp_highrisk_H1_18"],
     df.at[index, "ratio_hosp_highrisk_H1_50"],
    df.at[index, "ratio_hosp_highrisk_H1_65"],
    df.at[index, "ratio_hosp_highrisk_H1_75"]],
    [0, 5, 18, 50, 65, 75])


def ratio_hosp_highrisk_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_hosp_highrisk_H3_0"],
     df.at[index, "ratio_hosp_highrisk_H3_5"],
     df.at[index, "ratio_hosp_highrisk_H3_18"],
     df.at[index, "ratio_hosp_highrisk_H3_50"],
    df.at[index, "ratio_hosp_highrisk_H3_65"],
    df.at[index, "ratio_hosp_highrisk_H3_75"]],
    [0, 5, 18, 50, 65, 75])



def ratio_hosp_highrisk_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "ratio_hosp_highrisk_B_0"],
     df.at[index, "ratio_hosp_highrisk_B_5"],
     df.at[index, "ratio_hosp_highrisk_B_18"],
     df.at[index, "ratio_hosp_highrisk_B_50"],
    df.at[index, "ratio_hosp_highrisk_B_65"],
    df.at[index, "ratio_hosp_highrisk_B_75"]],
    [0, 5, 18, 50, 65, 75])


def lowRiskhospitalizationRate_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskhospitalizationRate_H1_0"],
     df.at[index, "lowRiskhospitalizationRate_H1_5"],
     df.at[index, "lowRiskhospitalizationRate_H1_18"],
     df.at[index, "lowRiskhospitalizationRate_H1_50"],
    df.at[index, "lowRiskhospitalizationRate_H1_65"],
    df.at[index, "lowRiskhospitalizationRate_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskhospitalizationRate_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskhospitalizationRate_H3_0"],
     df.at[index, "lowRiskhospitalizationRate_H3_5"],
     df.at[index, "lowRiskhospitalizationRate_H3_18"],
     df.at[index, "lowRiskhospitalizationRate_H3_50"],
    df.at[index, "lowRiskhospitalizationRate_H3_65"],
    df.at[index, "lowRiskhospitalizationRate_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def lowRiskhospitalizationRate_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "lowRiskhospitalizationRate_B_0"],
     df.at[index, "lowRiskhospitalizationRate_B_5"],
     df.at[index, "lowRiskhospitalizationRate_B_18"],
     df.at[index, "lowRiskhospitalizationRate_B_50"],
    df.at[index, "lowRiskhospitalizationRate_B_65"],
    df.at[index, "lowRiskhospitalizationRate_B_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "highRiskhospitalizationRate_H1_0"],
     df.at[index, "highRiskhospitalizationRate_H1_5"],
     df.at[index, "highRiskhospitalizationRate_H1_18"],
     df.at[index, "highRiskhospitalizationRate_H1_50"],
    df.at[index, "highRiskhospitalizationRate_H1_65"],
    df.at[index, "highRiskhospitalizationRate_H1_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "highRiskhospitalizationRate_H3_0"],
     df.at[index, "highRiskhospitalizationRate_H3_5"],
     df.at[index, "highRiskhospitalizationRate_H3_18"],
     df.at[index, "highRiskhospitalizationRate_H3_50"],
    df.at[index, "highRiskhospitalizationRate_H3_65"],
    df.at[index, "highRiskhospitalizationRate_H3_75"]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "highRiskhospitalizationRate_B_0"],
     df.at[index, "highRiskhospitalizationRate_B_5"],
     df.at[index, "highRiskhospitalizationRate_B_18"],
     df.at[index, "highRiskhospitalizationRate_B_50"],
    df.at[index, "highRiskhospitalizationRate_B_65"],
    df.at[index, "highRiskhospitalizationRate_B_75"]],
    [0, 5, 18, 50, 65, 75])


def lowRiskOutpatientProbPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "prob_outpatient_lowrisk_0"],
     df.at[index, "prob_outpatient_lowrisk_5"],
     df.at[index, "prob_outpatient_lowrisk_18"],
    df.at[index, "prob_outpatient_lowrisk_65"]],
    [0, 5, 18, 65])

def highRiskOutpatientProbPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "prob_outpatient_highrisk_0"],
     df.at[index, "prob_outpatient_highrisk_5"],
     df.at[index, "prob_outpatient_highrisk_18"],
    df.at[index, "prob_outpatient_highrisk_65"]],
    [0, 5, 18, 65])

def susceptibility_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_H1_0"],
     df.at[index, "susceptibility_H1_5"],
     df.at[index, "susceptibility_H1_25"],
     df.at[index, "susceptibility_H1_65"]],
    [0, 5, 25, 65])

def susceptibility_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_H3_0"],
     df.at[index, "susceptibility_H3_5"],
     df.at[index, "susceptibility_H3_25"],
     df.at[index, "susceptibility_H3_65"]],
    [0, 5, 25, 65])

def susceptibility_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_B_0"],
     df.at[index, "susceptibility_B_5"],
     df.at[index, "susceptibility_B_25"],
     df.at[index, "susceptibility_B_65"]],
    [0, 5, 25, 65])

def transmissionScaling_H1(df, index):
    
    return df.at[index, "beta_H1"]

def transmissionScaling_H3(df, index):
    
    return df.at[index, "beta_H3"]

def transmissionScaling_B(df, index):
    
    return df.at[index, "beta_B"]

def vac_eff_hospitalization(df, index):
    
    return df.at[index, "vac_eff_hospitalization"]

def vac_eff_mortality(df, index):
    
    return df.at[index, "vac_eff_mortality"]

def prob_hosp_scaling(df, index):
    
    return df.at[index, "prob_hosp_scaling"]

def prob_death_scaling(df, index):

    return df.at[index, "prob_death_scaling"]

