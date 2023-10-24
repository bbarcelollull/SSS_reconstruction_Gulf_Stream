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
from scipy                   import interpolate

"""
Open interpolated data and plot S_adv, S_SMAP and S_TSG 
on '20170212' - Figure for paper. Not used!!

written by Bàrbara Barceló-Llull on 13 April 2020 at Mallorca
"""

if __name__ == '__main__':
    
  plt.close('all')
  
  dir_dic = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all_Bigger_Domain/post-processing/'
                
  file_dic = 'SimOl_B_alt_BD_sal_int_to_TSG_final.pkl'

  dirfig_ppr      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/' + \
                'fig_Oleander_BiggerRegion/forpaper/'

  # FINAL BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46
  
  ''' Open interpolated data '''
  
  f = open(dir_dic + file_dic, 'rb')
  data_all = pickle.load(f)
  f.close()  

  date_to_plot = '20170212'
  
  data_date = data_all[date_to_plot]['07']
  
  
  lon_tsg  = data_date['lon_tsg']
  lat_tsg  = data_date['lat_tsg']
  sal_tsg  = data_date['sal_tsg']
  sal_adv  = data_date['s_advi']
  sal_smap = data_date['s_smapi']
  time_tsg = data_date['time_tsg']
  
  ''' Plot data '''
#  fsize = 14
#  lw = 2
#  fig = plt.figure(figsize=(10,5))
#
#  plt.plot(lon_tsg, sal_tsg, linestyle='-', color='k' ,
#                   linewidth=lw,label=r'S$_{TSG}$', zorder=1)
#
#  plt.plot(lon_tsg, sal_adv, linestyle='--', color='steelblue', 
#             linewidth=lw, 
#             label=r'S$_{adv}$', zorder=3)
#  
#  plt.plot(lon_tsg, sal_smap, linestyle='dashdot', color='pink' ,
#             linewidth=lw, 
#             label=r'S$_{SMAP}$', zorder=2)  
#      
#
#
#      
#  plt.legend(loc=4,prop={'size': fsize})
#     
#  #plt.xlabel('Longitude', fontsize=fsize)
#  plt.ylabel('Salinity [PSU]', fontsize=fsize)
#        
#  plt.tick_params(axis='both', labelsize=fsize)
#  plt.tight_layout()
#  plt.xlim(lon_tsg.min()-0.5, lon_tsg.max()+0.5)      
#  plt.xticks([-72, -70, -68, -66], 
#             ['72$^{\circ}$W', '70$^{\circ}$W', '68$^{\circ}$W', '66$^{\circ}$W'])
#             
#  plt.savefig(dirfig_ppr + 'sim_Oleander_back_alt_BD_int_' + date_to_plot + '.png', dpi=200)
#          
#
  fsize = 13
  lw = 2
  fig = plt.figure(figsize=(10,5))

  plt.plot(lon_tsg, sal_tsg, linestyle='-', color='k' ,
                   linewidth=lw,label=r'S$_{TSG}$', zorder=1)

  plt.plot(lon_tsg, sal_adv, linestyle='--', color='darkviolet', 
             linewidth=lw, 
             label=r'S$_{adv}$', zorder=3)
  
  plt.plot(lon_tsg, sal_smap, linestyle='dashdot', color='khaki' ,
             linewidth=lw, 
             label=r'S$_{SMAP}$', zorder=2)  
      


      
  plt.legend(loc=4,prop={'size': fsize})
     
  #plt.xlabel('Longitude', fontsize=fsize)
  plt.ylabel('Salinity [PSU]', fontsize=fsize)
        
  plt.tick_params(axis='both', labelsize=fsize)
  plt.tight_layout()
  plt.xlim(lon_tsg.min()-0.5, lon_tsg.max()+0.5)      
  plt.xticks([-72, -70, -68, -66], 
             ['72$^{\circ}$W', '70$^{\circ}$W', '68$^{\circ}$W', '66$^{\circ}$W'])
             
  plt.savefig(dirfig_ppr + 'sim_Oleander_back_alt_BD_int_' + date_to_plot + '.png', dpi=200)
          

                          