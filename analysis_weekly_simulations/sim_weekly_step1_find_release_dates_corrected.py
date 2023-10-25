#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
import sss_toolbox           as to
from scipy                   import interpolate
from scipy                   import stats
from matplotlib              import dates as mdates
from datetime                import datetime 
import shutil
from os                      import path
import sim_toolbox           as sto
import pickle
import glob
"""
Open SMAP data and select weekly dates to release 
massive 7-day backward simulations.

written by Bàrbara Barceló-Llull on 15-January-2020 at Mallorca
  
  
"""

def convert_smapdateformat2py(time_orig):  
      # first date
      d_ini_or    = mdates.date2num(datetime(2000, 1, 1, 0, 0, 0))  
  
      # Return seconds as days
      time_days = np.zeros(time_orig.shape)
  
      for ind, tt in enumerate(time_orig):
         time_days[ind] = mdates.seconds(tt) # Return seconds as days.
    
      # Sum these days to d_ini_or
      time = d_ini_or + time_days
      
      time = time[0]
      
      print('')
      print('SSS - center of the product time interval...') 
      print(mdates.num2date(time).strftime("%Y-%m-%d %H:%M"))  
      
      return time
  
if __name__ == '__main__':
  
  plt.close('all')
  
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  dir_savedic = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_weekly_all/'

  #/year/4 days before/file    
  files    = sorted(glob.glob(dir_SSS + '/*/*/*.nc'))
  
  dates_all = []
  
  for file in files: 
      
      ''' Open date of the file '''
      
      nc   = netcdf.Dataset(file, 'r')
      date_orig = nc.variables['time'][:].data # seconds since 2000-01-01T00:00:00Z
      nc.close() 

      ''' convert time to python format '''
     
      date_py = convert_smapdateformat2py(date_orig)
     
      dates_all = np.append(dates_all, date_py)

  ''' 
  There is a gap in SMAP data between
  19 June - 22 July, 2019  days: 171-203 '''
  
  # search the beggining of the gap
  
  arg_gap = np.where(np.diff(dates_all)>1)
  
  # There are two gaps!!
  
  arg_gap1 = arg_gap[0][0]
  arg_gap2 = arg_gap[0][1]
  
  ''' The simulation release dates are weekly dates, start with the index 7
  because the simulations are runned backwards 7 days'''
  
  
  
  plt.figure()
  plt.plot_date(dates_all, dates_all, 'ob')
  plt.plot_date(dates_all[arg_gap1], dates_all[arg_gap1], 'xr')
  plt.plot_date(dates_all[arg_gap2], dates_all[arg_gap2], 'xg')
  
  
  # dates release will be:
  
  # from the beggining to the first gap
  dates_sim_1 = dates_all[7:arg_gap1+1:7]
  dates_sim_2 = dates_all[arg_gap1+1+7:arg_gap2+1:7]
  dates_sim_3 = dates_all[arg_gap2+1+7::7]
  
  plt.figure()
  plt.plot_date(dates_all, dates_all, 'ob')  
  plt.plot_date(dates_sim_1, dates_sim_1, 'xy')  
  plt.plot_date(dates_sim_2, dates_sim_2, 'xr')
  plt.plot_date(dates_sim_3, dates_sim_3, 'xc')  
  
  
  dates_sim = np.append(np.append(dates_sim_1, dates_sim_2), dates_sim_3)

  plt.figure()
  plt.plot_date(dates_all, dates_all, 'ob')  
  plt.plot_date(dates_sim, dates_sim, 'xy')  
      
  # save dates
  dic_tsg_transects = {'date_release' : dates_sim}

  f = open(dir_savedic + 'sim_weekly_relase_dates_corrected_gaps.pkl','wb')
  pickle.dump(dic_tsg_transects,f)
  f.close()   