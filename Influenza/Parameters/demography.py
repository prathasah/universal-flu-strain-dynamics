from PiecewiseAgeParameter import PiecewiseAgeNumber
import pandas as pd
import os

def return_demography(season):
    ##consider the initial year for demography
    year = season[:4]
    __file__ = '.'

    df = pd.read_csv(os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "../"
                     "Influenza/Parameters/data/demography_yearwise.csv")))
    agelist = df['Age'].tolist()
    popsize = df[year].tolist()
    return  PiecewiseAgeNumber(popsize, agelist)

######################################################################33
if __name__ == "__main__":
    demography("2011-12")