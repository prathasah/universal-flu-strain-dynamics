# -*- coding: iso-8859-1 -*-
#
# Epidemiological parameter values
#

from PiecewiseAgeParameter import PiecewiseAgeRate

def recoveryRatePW(df, index):
    
    return PiecewiseAgeRate(
    [1/(1.*df.loc[df['iter'] == index, "infectious_period_0"].iloc[0]),
     1/(1.*df.loc[df['iter'] == index, "infectious_period_15"].iloc[0])],
    [0, 15])


def proportionHighRiskPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "proportionHighRisk_0"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_2"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_5"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_19"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_25"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_50"].iloc[0],
     df.loc[df['iter'] == index, "proportionHighRisk_65"].iloc[0]],
    [0, 2, 5, 19, 25, 50, 65])

def SeasonalVaccineEfficacyVsInfection_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_0.5"].iloc[0],
      df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_9"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_50"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H1_65"].iloc[0]],
    [0, 0.5,  5, 18, 50, 65])

def SeasonalVaccineEfficacyVsInfection_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_0.5"].iloc[0],
      df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_9"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_50"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_H3_65"].iloc[0]],
    [0, 0.5,  5, 18, 50, 65])

def age_specific_vaccineEfficacyVsInfectionPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index,  "age_specific_vaccineEfficacyVsInfection_0"].iloc[0],
      df.loc[df['iter'] == index, "age_specific_vaccineEfficacyVsInfection_0.5"].iloc[0],
      df.loc[df['iter'] == index, "age_specific_vaccineEfficacyVsInfection_5"].iloc[0],
     df.loc[df['iter'] == index, "age_specific_vaccineEfficacyVsInfection_18"].iloc[0],
     df.loc[df['iter'] == index, "age_specific_vaccineEfficacyVsInfection_50"].iloc[0]],
    [0, 0.5,  5, 18, 50])


def vaccineEfficacyVsInfection_all_ages(df, index):
    
    return df.loc[df['iter'] == index, "vaccineEfficacyVsInfection_all_ages"].iloc[0]

def SeasonalVaccineEfficacyVsInfection_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_0"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_0.5"].iloc[0],
      df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_9"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_18"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_50"].iloc[0],
     df.loc[df['iter'] == index, "seasonal_vaccineEfficacy_B_65"].iloc[0]],
    [0, 0.5,  5, 18, 50, 65])

def relative_vaccineEfficacyVsHospitalization_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H1_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H1_16"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H1_65"].iloc[0]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsHospitalization_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H3_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H3_16"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_H3_65"].iloc[0]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsHospitalization_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_B_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_B_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_B_16"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsHospitalization_B_65"].iloc[0]],
    [0, 0.5, 16, 65])

def relative_vaccineEfficacyVsDeath_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H1_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H1_65"].iloc[0]],
    [0, 0.5, 18, 65])

def relative_vaccineEfficacyVsDeath_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H3_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_H3_65"].iloc[0]],
    [0, 0.5,18, 65])

def relative_vaccineEfficacyVsDeath_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_B_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_B_0.5"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_B_18"].iloc[0],
     df.loc[df['iter'] == index, "relative_vaccineEfficacyVsDeath_B_65"].iloc[0]],
    [0, 0.5,18,65])


def relative_prob_deathPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_prob_death_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_death_5"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_death_18"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_death_50"].iloc[0],
    df.loc[df['iter'] == index, "relative_prob_death_65"].iloc[0]],
    [0, 5, 18, 50, 65])

def ratio_death_strain_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_death_strain_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H1_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_death_strain_H1_65"].iloc[0]],
    [0, 5, 18, 50, 65])

def ratio_death_strain_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_death_strain_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_strain_H3_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_death_strain_H3_65"].iloc[0]],
    [0, 5, 18, 50, 65])


def ratio_death_highrisk_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_death_highrisk_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H1_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_death_highrisk_H1_65"].iloc[0]],
    [0, 5, 18, 50, 65])

def ratio_death_highrisk_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_death_highrisk_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_H3_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_death_highrisk_H3_65"].iloc[0]],
    [0, 5, 18, 50, 65])


def ratio_death_highrisk_BPW(df, index):
        
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_death_highrisk_B_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_B_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_B_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_death_highrisk_B_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_death_highrisk_B_65"].iloc[0]],
    [0, 5, 18, 50, 65])


def relative_prob_hospPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "relative_prob_hosp_0"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_hosp_5"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_hosp_18"].iloc[0],
     df.loc[df['iter'] == index, "relative_prob_hosp_50"].iloc[0],
    df.loc[df['iter'] == index, "relative_prob_hosp_65"].iloc[0]],
    [0, 5, 18, 50,65])

def ratio_hosp_highrisk_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_65"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_H1_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])


def ratio_hosp_highrisk_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_65"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_H3_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])



def ratio_hosp_highrisk_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_0"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_5"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_18"].iloc[0],
     df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_50"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_65"].iloc[0],
    df.loc[df['iter'] == index, "ratio_hosp_highrisk_B_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])


def lowRiskhospitalizationRate_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_50"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_65"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H1_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])

def lowRiskhospitalizationRate_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_50"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_65"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_H3_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])

def lowRiskhospitalizationRate_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_0"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_5"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_18"].iloc[0],
     df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_50"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_65"].iloc[0],
    df.loc[df['iter'] == index, "lowRiskhospitalizationRate_B_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_18"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_50"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_65"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_H1_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_18"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_50"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_65"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_H3_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])

def highRiskhospitalizationRate_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_0"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_5"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_18"].iloc[0],
     df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_50"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_65"].iloc[0],
    df.loc[df['iter'] == index, "highRiskhospitalizationRate_B_75"].iloc[0]],
    [0, 5, 18, 50, 65, 75])


def lowRiskOutpatientProbPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "prob_outpatient_lowrisk_0"].iloc[0],
     df.loc[df['iter'] == index, "prob_outpatient_lowrisk_5"].iloc[0],
     df.loc[df['iter'] == index, "prob_outpatient_lowrisk_18"].iloc[0],
    df.loc[df['iter'] == index, "prob_outpatient_lowrisk_65"].iloc[0]],
    [0, 5, 18, 65])

def highRiskOutpatientProbPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "prob_outpatient_highrisk_0"].iloc[0],
     df.loc[df['iter'] == index, "prob_outpatient_highrisk_5"].iloc[0],
     df.loc[df['iter'] == index, "prob_outpatient_highrisk_18"].iloc[0],
    df.loc[df['iter'] == index, "prob_outpatient_highrisk_65"].iloc[0]],
    [0, 5, 18, 65])

def susceptibility_H1PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "susceptibility_H1_0"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H1_5"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H1_25"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H1_65"].iloc[0]],
    [0, 5, 25, 65])

def susceptibility_H3PW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "susceptibility_H3_0"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H3_5"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H3_25"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_H3_65"].iloc[0]],
    [0, 5, 25, 65])

def susceptibility_BPW(df, index):
    
    return PiecewiseAgeRate(
    [df.loc[df['iter'] == index, "susceptibility_B_0"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_B_5"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_B_25"].iloc[0],
     df.loc[df['iter'] == index, "susceptibility_B_65"].iloc[0]],
    [0, 5, 25, 65])

def beta_H1(df, index):
    
    return df.loc[df['iter'] == index, "beta_H1"].iloc[0]

def beta_H3(df, index):
    
    return df.loc[df['iter'] == index, "beta_H3"].iloc[0]

def beta_B(df, index):
    
    return df.loc[df['iter'] == index, "beta_B"].iloc[0]

def vac_eff_hospitalization(df, index):
    
    return df.loc[df['iter'] == index, "vac_eff_hospitalization"].iloc[0]

def vac_eff_mortality(df, index):
    
    return df.loc[df['iter'] == index, "vac_eff_mortality"].iloc[0]

def prob_hosp_scaling(df, index):
    
    return df.loc[df['iter'] == index, "prob_hosp_scaling"].iloc[0]

def prob_death_scaling(df, index):

    return df.loc[df['iter'] == index, "prob_death_scaling"].iloc[0]

