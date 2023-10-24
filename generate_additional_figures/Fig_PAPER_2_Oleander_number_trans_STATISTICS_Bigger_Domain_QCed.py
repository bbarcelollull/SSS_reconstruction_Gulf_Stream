#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
import pickle
from matplotlib              import dates as mdates
from datetime                import datetime, timedelta
from scipy                   import interpolate
import sim_toolbox           as sto
import geopy.distance
import rt_anatools_3         as rt_anatools
import matplotlib.gridspec   as gridspec
from matplotlib.ticker import FormatStrFormatter

"""
Figure temporal coverage and spatial distribution, 
and distribution along years and seasons
of the original transects.

In this code we do statistics from the file 
'Oleander_TSG_per_transect_QCed_BiggerDomain.pkl'
that is the file where Oleander transects are isolated
in a Big Domain including Bermuda and quality controlled.
This is the final file to use.

>> FIGURE FOR PAPER <<

written by Bàrbara Barceló-Llull on 30-March-2020 at Mallorca.
"""


def distance2apoint(lon, lat, lonp, latp):
  dist = np.zeros(lon.shape)
  for k in np.arange(lon.shape[0]):
    coords_1 = (lon[k], lat[k])
    coords_2 = (lonp, latp)
    #dist[k]  = geopy.distance.vincenty(coords_1, coords_2).km  
    dist[k]  = geopy.distance.geodesic(coords_1, coords_2).km  
    
  return dist


if __name__ == '__main__':

  dirfig     = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_paper/'
   
  
  ''' File with original TSG data '''  
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'
  file_dic   = 'Oleander_TSG_per_transect_QCed_BiggerDomain.pkl'

  f = open(dir_dic + file_dic, 'rb')
  dict_tsg_ori = pickle.load(f)
  f.close() 

  dates_strings_pre   = list(dict_tsg_ori.keys())
  dates_strings_ori   = dates_strings_pre.copy()
  dates_strings_ori.sort(key=int)
  print(dates_strings_ori)
  
  print('')
  #print('STATISTICS')
  print('Number of transects...', len(dates_strings_ori))
  

  ''' Loop for each transect '''    

  num_2016  = 0
  num_2017  = 0
  num_2018  = 0
  num_2019  = 0
  save_month = []

  for date in dates_strings_ori:
      
      lat_tran  = dict_tsg_ori[date]['lat']
      lon_tran  = dict_tsg_ori[date]['lon']
      time_tran = dict_tsg_ori[date]['time']
      sal_tran  = dict_tsg_ori[date]['sal']
      tem_tran  = dict_tsg_ori[date]['tem']
      
      date_tran = datetime.strptime(date, '%Y%m%d')
      
          
      save_month = np.append(save_month, date_tran.month)
          
      if date_tran.year == 2016:
              num_2016 = num_2016 + 1

      elif date_tran.year == 2017:
              num_2017 = num_2017 + 1
              
      elif date_tran.year == 2018:
              num_2018 = num_2018 + 1
              
      elif date_tran.year == 2019:
              num_2019 = num_2019 + 1   
    
          
  print('Number of valid in 2016...', num_2016)     
  print('Number of valid in 2017...', num_2017)       
  print('Number of valid in 2018...', num_2018) 
  print('Number of valid in 2019...', num_2019)       
  
  years_analyzed = [2016, 2017, 2018, 2019]
  Ntrans_year = [num_2016, num_2017, num_2018, num_2019] 

  # Transecs per season
  i_winter = np.where(np.logical_or(np.logical_or(save_month==1, 
                        save_month==2), save_month==3))
  i_spring = np.where(np.logical_or(np.logical_or(save_month==4, 
                        save_month==5), save_month==6))
  i_summer = np.where(np.logical_or(np.logical_or(save_month==7, 
                        save_month==8), save_month==9))
  i_autumn = np.where(np.logical_or(np.logical_or(save_month==10, 
                        save_month==11), save_month==12))  
  
  num_winter = len(i_winter[0])
  num_spring = len(i_spring[0])
  num_summer = len(i_summer[0])
  num_autumn = len(i_autumn[0])
  
  Ntrans_season = [num_winter, num_spring, num_summer, num_autumn]
  
  print('Number of valid in winter...', num_winter)  
  print('Number of valid in spring...', num_spring)  
  print('Number of valid in summer...', num_summer)  
  print('Number of valid in autumn...', num_autumn)  
  

  fig = plt.figure(figsize=(10,10))
  #gs   = gridspec.GridSpec(2, 2, width_ratios=[18,1, 12]) 
  gs   = gridspec.GridSpec(2, 2, width_ratios=[6,4]) 
  tz = 14
  
  ax1  = plt.subplot(gs[:,0])
  #axcb  = plt.subplot(gs[:,1])
  #ax2  = plt.subplot(gs[0,2])
  #ax3  = plt.subplot(gs[1,2])
  ax2  = plt.subplot(gs[0,1])
  ax3  = plt.subplot(gs[1,1]) 
  
  for date in dates_strings_ori:
      
      lat_tran  = dict_tsg_ori[date]['lat']
      lon_tran  = dict_tsg_ori[date]['lon']
      time_tran = dict_tsg_ori[date]['time']
      sal_tran  = dict_tsg_ori[date]['sal']
      
      if len(sal_tran) >= 1000: # minimum number of data points

        sc= ax1.scatter(lat_tran, time_tran, c=sal_tran, vmin=32.5, vmax=36.5, s=1)
    
    
  date_form = mdates.DateFormatter("%b-%Y")
  ax1.yaxis.set_major_formatter(date_form)
  ax1.yaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,6)))
  
  ax1.set_xlabel('Latitude [$^{\circ}$N]', fontsize=tz)
  ax1.invert_xaxis()
  ax1.tick_params(labelsize=tz-2)  
  
  #cb = plt.colorbar(sc, cax=axcb)  
  cb = plt.colorbar(sc, ax=ax1)  
  cb.ax.tick_params(labelsize=tz-2) 
  cb.ax.set_title('[PSU]', fontsize=tz-2)
  
  ax2.bar(years_analyzed, Ntrans_year, color='teal',align='center')
  ax2.xaxis.set_major_formatter(FormatStrFormatter('%i'))
  ax2.set_xlabel('Year', fontsize=tz)
  ax2.set_ylabel('Number of transects', fontsize=tz)
  ax2.set_xticks(years_analyzed, years_analyzed)
  ax2.tick_params(labelsize=tz-2)  
  #ax2.tick_params(labelsize=tz-2, bottom='off', top='off')
  
  ax3.bar(np.arange(4), Ntrans_season, color='teal',align='center')
  ax3.xaxis.set_major_formatter(FormatStrFormatter('%i'))
  ax3.set_xlabel('Season', fontsize=tz)
  ax3.set_ylabel('Number of transects', fontsize=tz)
  ax3.set_xticks(np.arange(4))
  ax3.set_xticklabels(['Winter', 'Spring', 'Summer', 'Autumn'], fontsize=11)
  ax3.tick_params(labelsize=tz-2, bottom='off', top='off')

  fig.text(0.015, 0.97, '(a)' , fontsize = tz+1, zorder=200)
  fig.text(0.61, 0.97, '(b)' , fontsize = tz+1, zorder=200)
  fig.text(0.61 , 0.5, '(c)' , fontsize = tz+1, zorder=200)

  plt.tight_layout()
  fig.savefig(dirfig + 'TSG_Ol_statistics_original_lat_BiggerDomain.png', dpi=200)

 
