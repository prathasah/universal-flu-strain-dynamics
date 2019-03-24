import os
import pandas as pd

if '__file__' not in globals():
    __file__ = '.'

location_data = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "./"))
"""Location of directory that contains data for the model."""


ve_data = {}
"""Dictionary with data frames containing the vaccine effectiveness of
each influenza season. Seasons are the keys of the dictionary.

"""

ve_data['2010-11'] = pd.read_csv(
    os.path.join(location_data, "treanor2012.csv"))

ve_data['2011-12'] = pd.read_csv(
    os.path.join(location_data, "ohnmit2013.csv"))

ve_data['2012-13'] = pd.read_csv(
    os.path.join(location_data, "mclean2014.csv"))

ve_data['2013-14'] = pd.read_csv(
    os.path.join(location_data, "gaglani2016.csv"))

ve_data['2014-15'] = pd.read_csv(
    os.path.join(location_data, "zimmerman2016.csv"))

ve_data['2015-16'] = pd.read_csv(
    os.path.join(location_data, "jackson2017.csv"))

ve_data['2016-17'] = pd.read_csv(
    os.path.join(location_data, "unpublished2018.csv"))

ve_data['2017-18'] = pd.read_csv(
    os.path.join(location_data, "rolfes2019.csv"))
