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
Open file with fronts detected on step 4 and make figures for paper:
    plot some examples of transects with fronts
    
    FIGURE FOR PAPER!! (Figure 3)

written by Bàrbara Barceló-Llull on 16-April-2020 at Mallorca

"""

if __name__ == '__main__':
    
  plt.close('all')

  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/'+ \
                'fig_Oleander_BiggerRegion/fronts/'
                
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all_Bigger_Domain/post-processing/'
  
  
  
  ''' File with the detected fronts '''  
  
  #file_dic = 'Oleander_BD_front_detection.pkl'
   
  file_dic = 'Oleander_BD_front_detection_5km.pkl'
  
  f = open(dir_dic + file_dic, 'rb')
  dict_fronts = pickle.load(f)
  f.close() 

  dates_strings = list(dict_fronts.keys())


 
  
  ''' 5 examples of transects '''  

  dates_to_plot = ['20160924', '20170212', '20180114', '20180602', '20190106']
  smin, smax = 29.5, 37
   
  fsize = 10
  lw = 1.5
  fig = plt.figure(figsize=(8,10))
  
  n=0
  
  enum = ['(a) ', '(b) ','(c) ','(d) ','(e) ',]  
  
  for ie, date_to_plot in enumerate(dates_to_plot):
      
    dict_date = dict_fronts[date_to_plot]
    print('>> DATE...', date_to_plot)
    
 
    
    dist_int = dict_date['dist_int']
    lon_int  = dict_date['lon_int']
    lat_int  = dict_date['lat_int']

    sal_tsg          = dict_date['sal_tsgi']
    sal_tsg_pfronts  = dict_date['sal_tsgi_pfronts']
    sal_tsg_nfronts  = dict_date['sal_tsgi_nfronts']

    sal_1way         = dict_date['sal_1wayi']
    sal_1way_pfronts = dict_date['sal_1wayi_pfronts']
    sal_1way_nfronts = dict_date['sal_1wayi_nfronts']
    
    sal_sat          = dict_date['sal_sati']
    sal_sat_pfronts  = dict_date['sal_sati_pfronts']
    sal_sat_nfronts  = dict_date['sal_sati_nfronts']
    
    # east of -70.7
    
    lonmin_sarg1 = -70.7
    print('mean tsg east of -70.7... ', 
          np.nanmean(sal_tsg[lon_int>lonmin_sarg1]))
    print('std tsg east of -70.7... ', 
          np.nanstd(sal_tsg[lon_int>lonmin_sarg1]))  
    print('')

    print('mean adv east of -70.7... ', 
          np.nanmean(sal_1way[lon_int>lonmin_sarg1]))
    print('std adv east of -70.7... ', 
          np.nanstd(sal_1way[lon_int>lonmin_sarg1]))  
    print('')
    
    print('mean smap east of -70.7... ', 
          np.nanmean(sal_sat[lon_int>lonmin_sarg1]))
    print('std smap east of -70.7... ', 
          np.nanstd(sal_sat[lon_int>lonmin_sarg1]))  
    print('')
    
    print('')
    
    lonmin_sarg2 = -69
    print('mean tsg east of -69... ', 
          np.nanmean(sal_tsg[lon_int>lonmin_sarg2]))  
    print('std tsg east of -69... ', 
          np.nanstd(sal_tsg[lon_int>lonmin_sarg2]))      
    print('')
    print('mean adv east of -69... ', 
          np.nanmean(sal_1way[lon_int>lonmin_sarg2]))  
    print('std adv east of -69... ', 
          np.nanstd(sal_1way[lon_int>lonmin_sarg2]))      
    print('')
    print('mean sat east of -69... ', 
          np.nanmean(sal_sat[lon_int>lonmin_sarg2]))  
    print('std sat east of -69... ', 
          np.nanstd(sal_sat[lon_int>lonmin_sarg2]))      
    print('')
    print('')
    
    lonmin_sarg3 = -70.5
    print('mean tsg east of -70.5... ', 
          np.nanmean(sal_tsg[lon_int>lonmin_sarg3]))
    print('std tsg east of -70.5... ', 
          np.nanstd(sal_tsg[lon_int>lonmin_sarg3]))
    print('')
    print('mean adv east of -70.5... ', 
          np.nanmean(sal_1way[lon_int>lonmin_sarg3]))
    print('std adv east of -70.5... ', 
          np.nanstd(sal_1way[lon_int>lonmin_sarg3]))
    print('')
    print('')
      
    if date_to_plot == '20180602':
        
       # interpolate gap
       
       f_1way_int = interpolate.interp1d(lon_int[~np.isnan(sal_1way)], 
                                          sal_1way[~np.isnan(sal_1way)], 
                                         kind='linear', fill_value="extrapolate")
       # range of tsg data with adv data
       rmin = lon_int[~np.isnan(sal_1way)].min()
       rmax = lon_int[~np.isnan(sal_1way)].max()
       
       cond = np.logical_and(lon_int>=rmin, lon_int<=rmax)
       
       
       sal_1way_int = f_1way_int(lon_int[cond])
       
#       plt.figure()
#       plt.plot(lon_ori, sal_1wayo, 'ob')
#       plt.plot(lon_ori[cond], sal_1wayo_int, 'xr')
       
       sal_1way = np.ones(lon_int.shape) * np.nan
       sal_1way[cond] = sal_1way_int
       
       
    n = n+1
    t = len(dates_to_plot)
    
    ax = plt.subplot(t,1,n)
    plt.plot(lon_int, sal_tsg, linestyle='-', color='k' ,
                   linewidth=lw-0.5,label=r'S$_{TSG}$', zorder=1)

    plt.plot(lon_int, sal_1way, linestyle='--', color='darkviolet', 
             linewidth=lw, 
             label=r'S$_{adv}$', zorder=3)
  
    plt.plot(lon_int, sal_sat, linestyle='dashdot', color='orange' ,
             linewidth=lw, 
             label=r'S$_{SMAP}$', zorder=2)  

     
    #plt.xlabel('Longitude', fontsize=fsize)
    plt.ylabel('Salinity [PSU]', fontsize=fsize)
        
    plt.tick_params(axis='both', labelsize=fsize)
    plt.tight_layout()
    plt.xlim(lon_int.min()-0.5, lon_int.max()+0.5)   
    
    date_datetime = datetime.strptime(date_to_plot, '%Y%m%d') 
    plt.text(0.01, 0.895, enum[ie] + date_datetime.strftime("%d %b %Y"),
             transform=ax.transAxes)
    
    if n==t:  
      plt.legend(loc=4,prop={'size': fsize})
      plt.xticks([-72, -70, -68, -66], 
             ['72$^{\circ}$W', '70$^{\circ}$W', '68$^{\circ}$W', '66$^{\circ}$W'])
    else:    
      plt.xticks([])
      # Hide the right and top spines
      ax.spines['bottom'].set_visible(False)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim(smin, smax)
  plt.savefig(dirfig + 'figppr_Oleander_5transects.png', dpi=200)
          


  ''' 3 examples and fronts '''

#  dates_to_plot = ['20170212', '20180114', '20190106']
#
#
#  fsize = 10
#  lw = 1.5
#  fig = plt.figure(figsize=(8,10))
#  
#  n=0
#  
#    
#  for date_to_plot in dates_to_plot:
#      
#    dict_date = dict_fronts[date_to_plot]
#  
#
#    dist_int = dict_date['dist_int']
#    lon_int  = dict_date['lon_int']
#    lat_int  = dict_date['lat_int']
#
#    sal_tsg          = dict_date['sal_tsgi']
#    sal_tsg_pfronts  = dict_date['sal_tsgi_pfronts']
#    sal_tsg_nfronts  = dict_date['sal_tsgi_nfronts']
#
#    sal_1way         = dict_date['sal_1wayi']
#    sal_1way_pfronts = dict_date['sal_1wayi_pfronts']
#    sal_1way_nfronts = dict_date['sal_1wayi_nfronts']
#    
#    sal_sat          = dict_date['sal_sati']
#    sal_sat_pfronts  = dict_date['sal_sati_pfronts']
#    sal_sat_nfronts  = dict_date['sal_sati_nfronts']
#    
#  
#    n = n+1
#    t = len(dates_to_plot)
#    
#    ax = plt.subplot(t,1,n)
#    plt.plot(lon_int, sal_tsg, linestyle='-', color='0.8' ,
#                   linewidth=lw-0.5,label=r'S$_{TSG}$', zorder=1)
#    plt.plot(lon_int, sal_tsg_pfronts, linestyle='-', color='k' ,
#                   linewidth=lw+1, zorder=1)  
#    plt.plot(lon_int, sal_tsg_nfronts, linestyle='-', color='k' ,
#                   linewidth=lw+1, zorder=1) 
#    
#    plt.plot(lon_int, sal_1way, linestyle='-', color='violet', 
#             linewidth=lw, 
#             label=r'S$_{adv}$', zorder=3)
#    plt.plot(lon_int, sal_1way_pfronts, linestyle='-', color='indigo' ,
#                   linewidth=lw+1, zorder=4)  
#    plt.plot(lon_int, sal_1way_nfronts, linestyle='-', color='indigo' ,
#                   linewidth=lw+1, zorder=4) 
#  
#    plt.plot(lon_int, sal_sat, linestyle='-', color='gold' ,
#             linewidth=lw, 
#             label=r'S$_{SMAP}$', zorder=2)  
#    plt.plot(lon_int, sal_sat_pfronts, linestyle='-', color='darkgoldenrod' ,
#                   linewidth=lw+1, zorder=3)  
#    plt.plot(lon_int, sal_sat_nfronts, linestyle='-', color='darkgoldenrod' ,
#                   linewidth=lw+1, zorder=3) 
#    
#    #plt.xlabel('Longitude', fontsize=fsize)
#    plt.ylabel('Salinity [PSU]', fontsize=fsize)
#        
#    plt.tick_params(axis='both', labelsize=fsize)
#    plt.tight_layout()
#    plt.xlim(lon_int.min()-0.5, lon_int.max()+0.5)   
#
#    plt.text(0.01, 0.895, date_to_plot,transform=ax.transAxes)
#
#    
#    if n==t:  
#      plt.legend(loc=4,prop={'size': fsize})
#      plt.xticks([-72, -70, -68, -66], 
#             ['72$^{\circ}$W', '70$^{\circ}$W', '68$^{\circ}$W', '66$^{\circ}$W'])
#    else:    
#      plt.xticks([])
#  plt.savefig(dirfig + 'figppr_Oleander_3transects_and_fronts.png', dpi=200)
#
