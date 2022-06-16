# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:53:03 2022

@author: alexia
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def lignes(name):
    """
    Calculate the amount of ligns in a file

    Parameters
    ----------
    name : str
        File name.

    Returns
    -------
    nligne : int
        Number of ligns.

    """
    fichier = open(name)
    nligne = 0
    while fichier.readline() != "":
        nligne = nligne +1
    fichier.close()
    return nligne

def colonnes(name):
    """
    Calculate the number of columns in a file

    Parameters
    ----------
    name : str
        File name.

    Returns
    -------
    int
        Amount of columns.

    """
    fichier = open(name)
    f = fichier.readline().split(" ")
    return len(f)

def colonnes2(name):
    """
    Calculate the number of columns in a file

    Parameters
    ----------
    name : str
        File name.

    Returns
    -------
    int
        Amount of columns.

    """
    fichier = open(name)
    f = fichier.readline().split("\t")
    return len(f)

def data(name, ncol, nligne):
    """
    Extract the data from the file into a matrix (used for SIE and SI volume)

    Parameters
    ----------
    name : str
        File name.
    ncol : int
        amount of columns.
    nligne : int
        amount of rows.

    Returns
    -------
    A : list
        Data written in the file.

    """
    fichier = open(name)
    A = [[0.0]*ncol for i in range (nligne)]
    for i in range (nligne):
        f = fichier.readline().split(" ")
        for k in range(ncol):
            A[i][k] = float(f[k])
    fichier.close()
    return A

def data2(name, ncol, nligne):
    """
    Extract the data from the file into a matrix (used for SI thickness)

    Parameters
    ----------
    name : str
        File name.
    ncol : int
        amount of columns.
    nligne : int
        amount of rows.

    Returns
    -------
    A : list
        Data extracted from the file.

    """
    fichier = open(name)
    A = [[0.0]*ncol for i in range (nligne)]
    for i in range (nligne):
        f = fichier.readline().split("\t")
        for k in range (ncol):
            A[i][k] = float(f[k])
    fichier.close()
    return A

def extract_year(fichier, nlig, sie_all_months):
    """
    Extract the first column corresponding to the years from a file

    Parameters
    ----------
    fichier : str
        File name.
    nlig : int
        Number of rows in the file.
    sie_all_months : list
        SIE extracted from the file for every month and every year,
        for example from the function data: data(fichier, ncol, nlig)

    Returns
    -------
    year : list
        List of years.

    """
    year = []
    for i in range(0,nlig):
        if sie_all_months[i][2] == 12.0:
            year.append(sie_all_months[i][1])
    return year

def extract_sie(fichier, nlig, ncol, sie_all_months, month):
    """
    

    Parameters
    ----------
    fichier : str
        file name.
    nlig : int
        amount of rows.
    ncol : int
        amount of columns.
    sie_all_months : float
        sea ice extent.
    month : TYPE
        DESCRIPTION.

    Returns
    -------
    sie : float
        sea ice extent from the chosen month for all years.

    """
    sie = []
    for i in range(0,nlig):
        if sie_all_months[i][2] == month:
            sie.append(sie_all_months[i][4])
    return sie

def extract_thickness(fichier, nlig, ncol, month):
    thickness_all = data2(fichier, ncol, nlig)
    thickness = []
    for i in range(nlig):
        if thickness_all[i][1] == month:
            thickness.append(thickness_all[i][2])
    return thickness

def extract_volume(fichier, nlig, ncol, month):
    volume_all = data(fichier, ncol, nlig)
    volume = []
    for i in range(nlig):
        volume.append(volume_all[i][month])
    return volume

def regression_line(sie, year):
    sie_mean = sum(sie)/len(sie)
    year_mean = sum(year)/len(year)
    xy = 0.0
    xx = 0.0
    for i in range(len(year)):
        xy += (year[i] - year_mean) * (sie[i] - sie_mean)
        xx += (year[i] - year_mean)**2
    a = xy/xx
    b = sie_mean - a * year_mean
    print(a)
    sie_regression = [0.0 for i in range(len(year))]
    for i in range(len(year)):
        sie_regression[i] = b + a * year[i]
    sie_reg_2022 = b + a * 2022
    return sie_regression, sie_reg_2022

def event_forecast(year, event, mu, sigma, sie_obs):
    """
    To calculate event occurrence probability and corresponding observed frequencies
    """
    p_event = np.zeros(len(year)-2)
    obs_freq =  np.zeros(len(year)-2)
    event_occ = np.zeros(len(year)-2)
    proba_color = []
    for i in range(len(year)-2):
        p_event[i] = norm.cdf((event[i+1] - mu[i]) / sigma[i])
        if sie_obs[i+2] < event[i+1]:
            event_occ[i] = 1
            proba_color.append("#50D51F")
        else:
            proba_color.append("#FF5733")
        for j in range(i):
            if sie_obs[j+2] < event[i+1]:
                obs_freq[i] += 1
        obs_freq[i] /= (i+1)
    return p_event, obs_freq, event_occ, proba_color

def stat_forecast(sie_sept, sie_may, year):
    mu = [0 for i in range(len(year)-2)]
    var = [0 for i in range(len(year)-2)]
    sd = [0 for i in range(len(year)-2)]
    sd_min = [0 for i in range(len(year)-2)]
    sd_max = [0 for i in range(len(year)-2)]
    for i in range(2,len(year)):
        mu[i-2] = np.mean(sie_sept[0:i]) + (sie_may[i] - np.mean(sie_may[0:i]))
        var[i-2] = (np.var(sie_sept[0:i]) + np.var(sie_may[0:i])) / i
        sd[i-2] = var[i-2]**0.5
        sd_min[i-2] = mu[i-2] - 2 * sd[i-2]
        sd_max[i-2] = mu[i-2] + 2 * sd[i-2]
    return mu, sd, sd_min, sd_max

def stat_forecast_all(sie_sept, sie_oct, sie_nov, sie_dec, sie_jan, sie_feb, sie_mar, sie_apr, sie_may, year):
    mu = [0 for i in range(len(year)-2)]
    var = [0 for i in range(len(year)-2)]
    sd = [0 for i in range(len(year)-2)]
    sd_min = [0 for i in range(len(year)-2)]
    sd_max = [0 for i in range(len(year)-2)]
    for i in range(2,len(year)):
        mu[i-2] = np.mean(sie_sept[0:i]) + ((sie_may[i] - np.mean(sie_may[0:i])) \
            + (sie_oct[i-1] - np.mean(sie_oct[0:i-1])) \
               + (sie_nov[i-1] - np.mean(sie_nov[0:i-1])) \
                   + (sie_dec[i-1] - np.mean(sie_dec[0:i-1])) \
                       + (sie_jan[i] - np.mean(sie_jan[0:i])) \
                           + (sie_feb[i] - np.mean(sie_feb[0:i])) \
                               + (sie_mar[i] - np.mean(sie_mar[0:i])) \
                                   + (sie_apr[i] - np.mean(sie_apr[0:i]))) / 8
        var[i-2] = (np.var(sie_sept[0:i]) + np.var(sie_may[0:i]) \
                    + np.var(sie_oct[0:i-1]) \
                        + np.var(sie_nov[0:i-1]) \
                            + np.var(sie_dec[0:i-1]) \
                                + np.var(sie_jan[0:i]) \
                                    + np.var(sie_feb[0:i]) \
                                        + np.var(sie_mar[0:i]) \
                                            + np.var(sie_apr[0:i]) ) / i
        sd[i-2] = var[i-2]**0.5
        sd_min[i-2] = mu[i-2] - 2 * sd[i-2]
        sd_max[i-2] = mu[i-2] + 2 * sd[i-2]
    return mu, sd, sd_min, sd_max


def verif_forecast(p, o, event_occ, n):
    bs = 0 # Brier score of the forecasts over 1981-2021
    bs_ref = 0 # Reference Brier score
    o_average = 0
    for i in range(n):
        bs += (p[i] - event_occ[i])**2
        #bs_ref += (o[i] - event_occ[i])**2
        o_average += event_occ[i]
    o_average /= n
    for i in range(n):
        bs_ref += (o_average - event_occ[i])**2
    bs /= n
    bs_ref /= n
    #p_climato = sum(event_occ)/np.size(event_occ)
    #bs_ref = p_climato * (1 - p_climato)
    bss = (bs - bs_ref) / (0 - bs_ref) # Brier skill score
    return bs, bs_ref, bss


#%%

"""
Forecast based on SIE observations.

NB: As no value is available for May 1986, this year isn't taken into account for any part of this forecast.
"""

f = "osisaf_nh_sie_monthly.txt" # File containing the SIE data monthly, coming from https://osisaf-hl.met.no/v2p1-sea-ice-index

# Extract data
year = extract_year(f, lignes(f), data(f, colonnes(f), lignes(f))) # All years except 1986
sie_sept = extract_sie(f, lignes(f), colonnes(f), data(f, colonnes(f), lignes(f)), 9.0) # September SIE for all years except 1986
sie_may = extract_sie(f, lignes(f), colonnes(f), data(f, colonnes(f), lignes(f)), 5.0) # May SIE for all years except 1986

# Regression line
sie_regression, sie_reg_2022 = regression_line(sie_sept, year)

#---------------------------------------------------------------------------------------------------------------

#Event definition: September Arctic sea ice extent will be less than previous year
event = [] # Event
event_year = [] # Years for which the event occurred
for i in range(1,len(year)):
    event.append(sie_sept[i-1])
    if sie_sept[i] < sie_sept[i-1]:
        event_year.append(year[i])

# Statistical forecasting system
mu, sd, sd_min, sd_max = stat_forecast(sie_sept, sie_may, year)

# Regression line associated to the forecast
forecast_trend, forecast_2022_extrapolate = regression_line(mu, year[2:])    

# Post-processing: removing the trend bias to the forecast
mu_pp = np.empty(len(year)-2)
for i in range(len(year)-2):
    mu_pp[i] = mu[i] - (forecast_trend[i] - sie_regression[i+2])

# Plot observed september SIE, regression line and forecast before and after post-processing
"""plt.figure()
plt.plot(year,sie_sept, color='red', label="September SIE")
#plt.plot(year, sie_regression, color='black', label = "Regression line")
plt.plot(2022, sie_reg_2022,'*', color='black', label="September 2022 SIE extrapolation forecast")
plt.plot(year[2:], mu, color='blue', label="September SIE forecast")
plt.fill_between(year[2:], sd_min, sd_max, alpha = 0.3, color='blue', label="Standard deviation")
plt.plot(year[2:], forecast_trend, '--', color='blue',label="Regression line")
plt.plot(year[2:], mu_pp, color='purple', label="September SIE forecast without trend bias")
plt.legend(fontsize="8")
plt.xlabel("Year")
plt.ylabel("SIE [km$^2$]")
plt.title("Arctic Sea Ice Extent: September time series")
plt.show()"""

# Retrospective probabilistic forecast of the event
p_event, obs_freq, event_occ, proba_color = event_forecast(year, event, mu, sd, sie_sept)

# Verification of retrospective forecast
bs, bs_ref, bss = verif_forecast(p_event, obs_freq, event_occ, len(year)-2)
#print("Forecast based on SIE \n Before post-processing : \n BS =", bs,"\n BS ref =", bs_ref,"\n BSS =", bss)

# Verification of retrospective forecast after bias correction
p_event_pp, obs_freq_pp, event_occ_pp, proba_color_pp = event_forecast(year, event, mu_pp, sd, sie_sept)
bs_pp, bs_ref_pp, bss_pp = verif_forecast(p_event_pp, obs_freq_pp, event_occ_pp, len(year)-2)
print("With sea ice extent \n After post-processing : \n BS =", bs_pp,"\n BS ref =", bs_ref_pp,"\n BSS =", bss_pp)

"""plt.figure()
plt.bar(year[2:], p_event_3_pp, color = proba_color_3_pp, alpha=0.5, label="Probability of event occurrence")
plt.plot(year[2:], obs_freq_3_pp,'-o', label="Observed frequencies")
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
plt.ylabel("Frequency")
plt.xlabel("Year")
plt.title("Probability of event occurrence and corresponding observed frequencies")
plt.show()"""

# Scatter plot
"""plt.figure()
plt.plot(sie_sept[2:], mu, 'o')
plt.plot([4*1e6,9*1e6], [4*1e6,9*1e6], color='black')
plt.title("Forecast (y-axis) vs observations (x-axis)")
plt.show()"""

#%%

"""
Forecast based on sea ice volume and thickness

NB: Here, every year is associated to data so they are all used for forecasting.
"""

# Extract SIE observations
f_2 = "osisaf_nh_sie_monthly_1.txt"

year_2 = extract_year(f_2, lignes(f_2), data(f_2, colonnes(f_2), lignes(f_2))) # All years of observation
sie_sept_all = extract_sie(f_2, lignes(f_2), colonnes(f_2), data(f_2, colonnes(f_2), lignes(f_2)), 9.0) # Observed September SIE for all years  

# Regression line based on SIE observations
print("reg")
sie_regression_2, sie_reg_2022_2 = regression_line(sie_sept_all, year_2)
print("stop reg")

# Event definition for all years: SIE extent lower than previous year
event_year_2 = [] # Years for which the event occurred
event_2 = [] # September SIE defined for all years based on event definition
for i in range(1,len(year_2)):
    event_2.append(sie_sept_all[i-1])
    if sie_sept_all[i] < sie_sept_all[i-1]:
        event_year_2.append(year_2[i])


# Extract sea ice thickness observations from September to May
f_thickness = "thickness_monthly.txt"
lig_thickness = lignes(f_thickness)
col_thickness = colonnes2(f_thickness)
thick_may = extract_thickness(f_thickness, lig_thickness, col_thickness, 5) # Sea ice thickness in May [m]
thick_sept = extract_thickness(f_thickness, lig_thickness, col_thickness, 9) # Sea ice thickness in September [m]
thick_oct = extract_thickness(f_thickness, lig_thickness, col_thickness, 10)
thick_nov = extract_thickness(f_thickness, lig_thickness, col_thickness, 11)
thick_dec = extract_thickness(f_thickness, lig_thickness, col_thickness, 12)
thick_jan = extract_thickness(f_thickness, lig_thickness, col_thickness, 1)
thick_feb = extract_thickness(f_thickness, lig_thickness, col_thickness, 2)
thick_mar = extract_thickness(f_thickness, lig_thickness, col_thickness, 3)
thick_apr = extract_thickness(f_thickness, lig_thickness, col_thickness, 4)

# Extract sea ice volume observations for May and September
f_volume = "PIOMAS_2sst_monthly_Current_v2_1.txt"
lig_volume = lignes(f_volume)
col_volume = colonnes(f_volume)
vol_may = extract_volume(f_volume, lig_volume, col_volume, 5) # Sea ice volume in May [10^3 km^3]
vol_sept = extract_volume(f_volume, lig_volume, col_volume, 9) # Sea ice volume in September [10^3 km^3]
vol_oct = extract_volume(f_volume, lig_volume, col_volume, 10)
vol_nov = extract_volume(f_volume, lig_volume, col_volume, 11)
vol_dec = extract_volume(f_volume, lig_volume, col_volume, 12)
vol_jan = extract_volume(f_volume, lig_volume, col_volume, 1)
vol_feb = extract_volume(f_volume, lig_volume, col_volume, 2)
vol_mar = extract_volume(f_volume, lig_volume, col_volume, 3)
vol_apr = extract_volume(f_volume, lig_volume, col_volume, 4)

# Calculate SIE based on sea ice volume and thickness
sie_may_2 = np.empty(len(year_2)) 
sie_sept_2 = np.empty(len(year_2)) 
sie_oct_2 = np.empty(len(year_2))
sie_nov_2 = np.empty(len(year_2))
sie_dec_2 = np.empty(len(year_2))
sie_jan_2 = np.empty(len(year_2))
sie_feb_2 = np.empty(len(year_2))
sie_mar_2 = np.empty(len(year_2))
sie_apr_2 = np.empty(len(year_2))
for i in range(len(year_2)):
    sie_may_2[i] = vol_may[i] * 1e3 / (thick_may[i] / 1e3)
    sie_sept_2[i] = vol_sept[i] * 1e3 / (thick_sept[i] / 1e3)
    sie_oct_2[i] = vol_oct[i] * 1e3 / (thick_oct[i] / 1e3)
    sie_nov_2[i] = vol_nov[i] * 1e3 / (thick_nov[i] / 1e3)
    sie_dec_2[i] = vol_dec[i] * 1e3 / (thick_dec[i] / 1e3)
    sie_jan_2[i] = vol_jan[i] * 1e3 / (thick_jan[i] / 1e3)
    sie_feb_2[i] = vol_feb[i] * 1e3 / (thick_feb[i] / 1e3)
    sie_mar_2[i] = vol_mar[i] * 1e3 / (thick_mar[i] / 1e3)
    sie_apr_2[i] = vol_apr[i] * 1e3 / (thick_apr[i] / 1e3)

# Forecast September SIE based on May and September calculated values of SIE
mu_2, sd_2, sd_min_2, sd_max_2  = stat_forecast_all(sie_sept_2, sie_oct_2, sie_nov_2, sie_dec_2, \
                                                    sie_jan_2, sie_feb_2, sie_mar_2, sie_apr_2, \
                                                        sie_may_2, year_2)

# Trend associated to the forecast
forecast_trend_2, forecast_2022_extrapolate_2 = regression_line(mu_2, year_2[2:])    

# Removing the trend bias
mu_pp_2 = np.empty(len(year_2)-2)
sd_min_pp_2 = np.empty(len(year_2)-2)
sd_max_pp_2 = np.empty(len(year_2)-2)
for i in range(len(year_2)-2):
    mu_pp_2[i] = mu_2[i] - (forecast_trend_2[i] - sie_regression_2[i+2])
sd_min_pp_2 = mu_pp_2 - 2 * np.array(sd_2)
sd_max_pp_2 = mu_pp_2 + 2 * np.array(sd_2)

# Probabilistic forecast based on SIE observations
p_event_pp_2, obs_freq_pp_2, event_occ_pp_2, proba_color_pp_2 = event_forecast(year_2, event_2, mu_pp_2, sd_2, sie_sept_all)
bs_pp_2, bs_ref_pp_2, bss_pp_2 = verif_forecast(p_event_pp_2, obs_freq_pp_2, event_occ_pp_2, len(year_2)-2)
print("Forecast based on SIV and thickness \n After post-processing : \n BS =", bs_pp_2,"\n BS ref =", bs_ref_pp_2,"\n BSS =", bss_pp_2)

# Plot September SIE
plt.figure()
plt.plot(year_2, sie_sept_all, '--', label="Observations", color="black", linewidth=0.75)
#plt.plot(year[2:], mu_pp, label="SIE forecast after pp")
#plt.plot(year_2[2:], mu_2, label="SIE volume-based forecast")
#plt.plot(year_2[2:], forecast_trend_2, label= "SIE volume-based forecast trend")
plt.plot(year_2[2:], mu_pp_2, label="SIE forecast", color="steelblue")
plt.fill_between(year_2[2:], sd_min_pp_2, sd_max_pp_2, alpha = 0.3, color='steelblue', label="95% CI")
plt.legend(fontsize="8")
plt.xlabel("Time [Years]")
plt.ylabel("SIE [km$^2$]")
#plt.title("Arctic Sea Ice Extent: September time series")
#plt.savefig("forecast.png", transparent=True)
plt.show()

# Plot probabilistic forecast
plt.figure()
plt.bar(year_2[2:], p_event_pp_2, color = proba_color_pp_2, alpha=0.7, label="Probability of event occurrence")
plt.plot(year_2[2:], obs_freq_pp_2,'-o', label="Observed frequencies")
#plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
plt.ylabel("Frequency")
plt.xlabel("Year")
#plt.title("Probability of event occurrence and corresponding observed frequencies for ice-based volume forecast")
plt.show()

# Scatter plot
"""plt.figure()
plt.plot(sie_sept_all[2:], mu_pp_2, 'o')
plt.plot([4*1e6,9*1e6], [4*1e6,9*1e6], color='black')
plt.title("Forecast (y-axis) vs observations (x-axis) for thickness")
plt.show()"""

#%%
"""
Forecast for September 2022
"""

# Extract sea ice thickness for May 2022, based on May mean thickness and 2022 April anomaly
mean_thick_may = np.mean(thick_may)
thick_april = extract_thickness(f_thickness, lig_thickness, col_thickness, 4)
thick_april_2022 = 1.778
thick_anomaly_april = (thick_april_2022 - np.mean(thick_april)) #/ np.std(thick_april)
thick_may_2022 = mean_thick_may + thick_anomaly_april

# Extract sea ice volume for May 2022, based on May mean volume and 2022 April anomaly
mean_vol_may = np.mean(vol_may)
vol_april = extract_volume(f_volume, lig_volume, col_volume, 4) # Sea ice volume in September [10^3 km^3]
vol_april_2022 = 22.997
vol_anomaly_april = (vol_april_2022 - np.mean(vol_april)) #/ np.std(vol_april)
vol_may_2022 = mean_vol_may + vol_anomaly_april

# Sea ice volume calculated for January, February, March, April and May 2022
for i in range(len(year_2)):
    sie_jan_2022 = 16.999 * 1e3 / (1.339 / 1e3)
    sie_feb_2022 = 19.745 * 1e3 / (1.479 / 1e3)
    sie_mar_2022 = 21.734 * 1e3 / (1.632 / 1e3)
    sie_apr_2022 = 22.997 * 1e3 / (1.778 / 1e3)
    sie_may_2022 = vol_may_2022 * 1e3 / (thick_may_2022 / 1e3)


# Forecast
mu_2022 = np.mean(sie_sept_2) + ((sie_may_2022 - np.mean(sie_may_2)) \
    + (sie_jan_2022 - np.mean(sie_jan_2)) \
        + (sie_feb_2022 - np.mean(sie_feb_2)) \
            + (sie_mar_2022 - np.mean(sie_mar_2)) \
                + (sie_apr_2022 - np.mean(sie_apr_2)) \
                    + (sie_oct_2[-1] - np.mean(sie_oct_2[0:-1])) \
                        + (sie_nov_2[-1] - np.mean(sie_nov_2[0:-1])) \
                            + (sie_dec_2[-1] - np.mean(sie_dec_2[0:-1]))) / 8

var_2022 = (np.var(sie_sept_2) + np.var(sie_may_2) \
            + np.var(sie_jan_2) \
                + np.var(sie_feb_2) \
                    + np.var(sie_mar_2) \
                        + np.var(sie_apr_2) \
                            + np.var(sie_oct_2[0:-1]) \
                                + np.var(sie_nov_2[0:-1]) \
                                    + np.var(sie_dec_2[0:-1])) / 44 # 44 = len(year_2) + 1
sd_2022 = var_2022**0.5
sd_min_2022 = mu_2022 - 2 * sd_2022
sd_max_2022 = mu_2022 + 2 * sd_2022

# Removing the trend bias
mu_pp_2022 = mu_2022 - (forecast_2022_extrapolate_2 - sie_reg_2022_2)
sd_min_pp_2022 = mu_pp_2022 - 2 * sd_2022
sd_max_pp_2022 = mu_pp_2022 + 2 * sd_2022

# Adding 2022 forecast to data from previous years
sd_min_3 = np.append(sd_min_pp_2, sd_min_pp_2022)
sd_max_3 = np.append(sd_max_pp_2, sd_max_pp_2022)
year_3 = np.append(year_2, 2022)
mu_3 = np.append(mu_pp_2, mu_pp_2022)

# Plot September SIE
plt.figure()
plt.plot(year_2, sie_sept_all, '--', label="Observations", color="black", linewidth=0.75)
#plt.plot(year_2, sie_regression_2, color='black', label = "Trend")
#plt.plot(2022, sie_reg_2022_2, 'x', color='black', markersize=3)
#plt.plot(year[2:], mu_pp, label="SIE forecast after pp")
#plt.plot(year_2[2:], mu_2, label="SIE volume-based forecast")
#plt.plot(year_2[2:], forecast_trend_2, label= "SIE volume-based forecast trend")
plt.plot(year_3[2:], mu_3, label="SIE forecast", color="steelblue")
#plt.plot(2022, mu_pp_2022, '*', color='steelblue')
plt.fill_between(year_3[2:], sd_min_3, sd_max_3, alpha = 0.3, color='steelblue', label="95% CI")
#plt.fill_between(2022, [sd_min_pp_2022, sd_max_pp_2022], alpha = 0.3, color='steelblue')
plt.legend()
plt.savefig("forecast_final.png", transparent=True)
plt.show()

# Probability for event occurrence in September 2022
p_2022 = norm.cdf((sie_sept_all[-1] - mu_pp_2022) / sd_2022)

# Plot probabilistic forecast
plt.figure()
plt.bar(year_2[2:], p_event_pp_2, color = proba_color_pp_2, alpha=0.7, label="Probability of event occurrence")
plt.bar(2022, p_2022, color = "steelblue", alpha=0.5)
#plt.plot(year_2[2:], obs_freq_pp_2,'-o', label="Observed frequencies")
#plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
plt.ylabel("Frequency")
plt.xlabel("Year")
#plt.title("Probability of event occurrence and corresponding observed frequencies for ice-based volume forecast")
plt.savefig("probability.png", transparent=True)
plt.show()
