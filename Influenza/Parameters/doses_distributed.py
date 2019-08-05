"""Provide continuous functions for the vaccination doses distributed
and applied per unit of time.

This library provides two functions to use in an ODE system that
models influenza transmission dynamics:

    - :func:`doses_applied_per_time_function`
    - :func:`doses_applied_before_start_season`

Both of these functions should be evaluated at a season, which can be
chosen to be any of the following ones:

    - '2011-12'
    - '2012-13'
    - '2013-14'
    - '2014-15'
    - '2015-16'
    - '2016-17'
    - '2017-18'

Configure these three variables in the Settings section within the
code:

    - :const:`start_of_season`
    - :const:`offset_distributed_to_applied`
    - :const:`time_unit`

"""

import os
import pandas as pd

import numpy as np
from scipy.integrate import odeint

import matplotlib.pyplot as plt

## * Settings

start_of_season = {'month':10, 'day':1}
"""Start of the season. A dictionary with the month and the day."""


offset_distributed_to_applied = 15
"""Time delay in days from distribution of doses to actual
vaccination.

"""

time_unit = 'day'
"""Time unit for the model. Could be either 'year' or 'day'."""



time_units_in_a_day = 1


## * Import data

if '__file__' not in globals():
    __file__ = '.'

df = pd.read_csv(
    os.path.abspath(
        os.path.join(os.path.dirname(__file__),
                     "../"
                     "Parameters/data/"
                     "doses_distributed.csv")),
    header=0,
    names=["date distributed", "doses"])

df['date distributed'] = pd.to_datetime(df['date distributed'],
                                        format="%m/%d/%Y")
df['doses'] = df['doses']*10**6

df['date applied'] = (df['date distributed']
                      + pd.to_timedelta(offset_distributed_to_applied,
                                        unit='d'))

df['season year'] = df['date applied'].apply(lambda x: (x.year - 1
                                                        if x.month < 7
                                                        else x.year))

df['season day'] = (df['date applied'] - pd.to_datetime(dict({'year':df['season year']}.items() + start_of_season.items()))).dt.days

df['season'] = df['season year'].apply(lambda x: "{}-{}".format(x, x+1 - 2000))

df['extrapolated'] = False



## * Extrapolate day 0

def extrapolate_day_0():
    """Fill the data frame with extrapolation of the doses value at day 0
    for each season.

    """
    global df
    for season in set(df.season):
        df_season = df[df['season'] == season]
        ## if day 0 is already there we are done, otherwise:
        if df_season[df_season['season day'] == 0].empty:
            if df_season[df_season['season day'] < 0].empty:
                t_minus_1 = 0
                y_minus_1 = 0
            else:
                t_minus_1 = max(
                    df_season['season day'][df_season['season day'] < 0])
                y_minus_1 = max(
                    df_season['doses'][df_season['season day'] < 0])

            t1 = min(df_season['season day'][df_season['season day'] >= 0])
            y1 = min(df_season['doses'][df_season['season day'] >= 0])

            y0 = y_minus_1 - (y1 - y_minus_1)/(t1 - t_minus_1) * t_minus_1

            season_year = int(season[:-3])

            date_day_0 = pd.to_datetime(
                dict({'year':df['season year']}.items() + start_of_season.items()))[0]
               

            date_distributed_0 = (date_day_0 - pd.to_timedelta(
                offset_distributed_to_applied,
                unit='d'))

            to_add = {'date distributed':[date_distributed_0],
                      'doses':[y0],
                      'date applied':[date_day_0],
                      'season year':[season_year],
                      'season day':[0],
                      'season':[season],
                      'extrapolated':[True]}

            df = df.append(pd.DataFrame(to_add), ignore_index=True, sort=True)


extrapolate_day_0()


## * Compute continuous functions of doses applied per time

def doses_applied_per_day_function_1(season):
    """Return a piecewise function with the number of doses applied per
    day in the `season`.

    :type season: str
    :param season: Influenza season.

    :rtype: function
    :return: Continuous function that can be evaluated at any time t,
        to give the doses applied per :const:`time_unit` in the
        `season`.
    """
    df_season = df[df['season'] == season]

    def doses_applied_per_day(t):
        if (df_season[df_season['season day'] >= t].empty or
            df_season[df_season['season day'] < t].empty):
            return 0
        else:
            y2 = min(df_season['doses'][df_season['season day'] >= t])
            y1 = max(df_season['doses'][df_season['season day'] < t])
            t2 = min(df_season['season day'][df_season['season day'] >= t])
            t1 = max(df_season['season day'][df_season['season day'] < t])
            return (y2-y1)/(t2-t1)

    return doses_applied_per_day


def doses_applied_per_day_function_2(season):
    """Return a piecewise function with the number of doses applied per
    day in the `season`. It seems like this function is faster than
    :func:`doses_applied_per_day_function_1`, but less exact.

    :type season: str
    :param season: Influenza season.

    :rtype: function
    :return: Continuous function that can be evaluated at any time t,
        to give the doses applied per :const:`time_unit` in the
        `season`.

    """
    df_season = df[df['season'] == season]

    t1s_t2s_y1s_y2s = [(t1, t2, y1, y2) for t1, t2, y1, y2
                       in zip(df_season['season day'][1:],
                              df_season['season day'][:-1],
                              df_season['doses'][1:],
                              df_season['doses'][:-1])
                       if t1 >= 0]

    # return t1s_t2s_y1s_y2s

    def doses_applied_per_day(t):
        return np.piecewise(t,
                            [t < t2
                             for t1, t2, y1, y2 in t1s_t2s_y1s_y2s],
                            [(y2-y1)/(t2-t1)
                             for t1, t2, y1, y2 in t1s_t2s_y1s_y2s])

    return doses_applied_per_day


ts_each_day = np.arange(0, 370)

doses_each_day = {season:[doses_applied_per_day_function_1(season)(t)
                          for t in ts_each_day]
                  for season in set(df.season)}

def doses_applied_per_time_function(season):
    """Return a piecewise function with the number of doses applied per
    :const:`time_unit` in the `season`.

    :type season: str
    :param season: Influenza season.

    :rtype: function
    :return: Continuous function that can be evaluated at any time t,
        to give the doses applied per :const:`time_unit` in the
        `season`.

    """
    return lambda t: (doses_each_day[season][int(t)+1])


def doses_applied_before_start_season(season):
    """Return the number of doses applied before the day 0 in the
    `season`. Day 0 is determined by :const:`start_of_season`.

    :type season: str
    :param season: Influenza season.

    :rtype: int
    :return: Number of doses applied before the start of
        `season`. Interpolated from the closest days available before
        and after day 0.

    """
    doses = list(df['doses'][(df['season day'] == 0) &
                             (df['season'] == season)])
    
    if len(doses) != 1:
        raise Exception(
            "Data contains more than one day 0)")
    else:
        return int(round(doses[0]))
    


## * Test function

if __name__ == '__main__':
    
    print sum([doses_applied_per_time_function("2011-12")(t) for t in xrange(180)])
    """
    for season in set(df['season']):
        df_season = df[df['season'] == season]

        df_season_non_extrapolated = df_season[~df_season['extrapolated']]
        fig, ax = plt.subplots()
        ax.scatter(
            df_season_non_extrapolated['date distributed'],
            df_season_non_extrapolated['doses'],
            color='Pink', label='Distributed (data)')
        ax.scatter(
            df_season_non_extrapolated['date applied'],
            df_season_non_extrapolated['doses'],
            color='Red', label='Applied (data)')

        y0 = doses_applied_before_start_season(season)

        season_duration_days = 181
        year = int(season[:-3])
        start_date = pd.to_datetime(dict({'year':df['season year']}.items() + start_of_season.items()))[0]

        real_days = pd.date_range(start_date, freq='D',
                                  periods=season_duration_days)

        days = np.arange(0, season_duration_days)

        ts = days*time_units_in_a_day
        for num in xrange(10):
            ys, infodict = odeint(lambda y, t: doses_applied_per_time_function(season)(t), y0, ts, full_output =True)
            if season =="2011-12":
                #print ("checkkkk==="),  doses_applied_per_time_function("2011-12")(60)
                print season, num, max(ys), [(round(t/2),round(doses_applied_per_time_function("2011-12")(t),2)) for t in xrange(210)]
                print ("total==="), sum([doses_applied_per_time_function("2011-12")(t) for t in sorted(ts)])
   
        ax.plot(real_days, ys, label='Applied (computed)')

        ylim_upper = 160*10**6
        ytick_step = 50*10**6
        ytick_minor_step = 10*10**6
        yticks = np.arange(0, ylim_upper+ytick_minor_step, ytick_step)
        yticks_minor = np.arange(0,
                                 ylim_upper+ytick_minor_step,
                                 ytick_minor_step)
        ax.set_yticks(yticks)
        ax.set_yticks(yticks_minor, minor=True)
        ax.set_yticklabels(['{:n}'.format(y/10**6) for y in yticks])
        ax.set_ylim(0, ylim_upper)
        ax.set_ylabel("Vaccination doses (millions)")
        ax.vlines([start_date], [0], [ylim_upper],
                  linestyles='dashed', color='Gray', label="Season start")

        xticks = pd.date_range(
            start_date,
            end=start_date + pd.to_timedelta(season_duration_days, unit='d'),
            freq='MS')
        ax.set_xticks(xticks)
        ax.set_xticklabels([x.strftime("%m/%d/%y") for x in xticks],
                           fontsize=11)
        ax.set_xlabel("Date")

        ax.legend(loc='center right')
        

        # ax.plot(real_days,
        #         [doses_applied_per_time_function(season)(t) for t in ts])
        
        if save_plots:
            fig.savefig("doses_function_{}.pdf".format(season))
            plt.close()
        else:
            ax.set_title("Season: {}".format(season))
            plt.show()
        """