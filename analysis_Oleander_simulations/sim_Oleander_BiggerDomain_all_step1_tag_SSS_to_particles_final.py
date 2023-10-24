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

'''
Tag SSS from SMAP to each particle at T0 of the backward simulation.

Here we open the files from the Lagrangian simulation, 
tag the salinity from SMAP at T0 position, 
and save the salinity data into the netcdf file to make 
comparisons in other codes.

Backward simulations for the Oleander data (Gulf Stream).

written by Bàrbara Barceló-Llull on 5-November-2019 at Mallorca

changed for the new version of SMAP (v4.0)

> on 21-Nov-2019:
    
changed for the massive simulations


> updated on 13 April 2020 for the simulations at the Bigger Domain
  FINAL CODE FOR PAPER!
  
  All simulations with Delayed-time currents as Oleander transects 
  are within the period in which we have dt altimetry currents: 
  [20160924 ,..., 20190227]. 
  Dt currents from 2015-03-31 to 2019-05-13.
'''


    
if __name__ == '__main__':
  
  plt.close('all')
    
  ''' Directories backward simulations '''

  dirsave     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_Oleander_all_Bigger_Domain/' #save outputs directory
  diroutputs  = dirsave + 'sim_outputs/'
  dirsave_sss = dirsave + 'sss_tagged/'
  
  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_Oleander_BiggerRegion/'
    
  
  ''' Set parameters for the code '''

  # running days
  rds = [7]

  # save figures?
  save = False #False#True

  # FINAL BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46

  # first day of the simulation in the GS 
  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0)) 
  # date_alt_fin = mdates.date2num(datetime(2019, 5, 13, 0, 0, 0))  
  
  date_ori = date_ori_dt  
  
  ''' File with the simulation dates '''  
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'
  #Bigger Domain and quality controlled.
  file_dic = 'Oleander_TSG_per_transect_QCed_BiggerDomain.pkl' 

  f = open(dir_dic + file_dic, 'rb')
  dict_tsg_all = pickle.load(f)
  f.close() 

  dates_strings_all = list(dict_tsg_all.keys())
  print(dates_strings_all)

  #after 2019-01-12 we don't have delayed-time altimetry data
  dates_strings = dates_strings_all.copy()
  dates_strings.sort(key=int)
  print(dates_strings)
  
  
  ''' Loop for each simulation '''
  
  for idt, date in enumerate(dates_strings):

    date_time_obj = datetime.strptime(date, '%Y%m%d')
    print('')
    print('Date of simulation...', date_time_obj)
  
    
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))
         

    for rd in rds:

      filename = 'sim_Oleander_back_alt_BD_'+mdates.num2date(date_fin).strftime("%Y%m%d")+ \
             '_' + '%02i'% rd + 'days.nc'
    
      print('')
      print('Openning... ' + filename)
      

      ''' Open simulation data '''
    
      ind_tf = 0  #first observation, beginning of simulation at TF
      ind_ti = -1 #last observation, end of simulation at T0
    
      nc    = netcdf.Dataset(diroutputs + filename, 'r')
    
      # read initial position at TF
      traj_tf  = nc.variables['trajectory'][:, ind_tf].data  # (traj, obs)
      lat_tf   = nc.variables['lat'][:, ind_tf].data
      lon_tf   = nc.variables['lon'][:, ind_tf].data   
      time_tf  = nc.variables['time'][:, ind_tf].data   #"seconds since 2017-01-01T00:00:00.000000000"

      # read final position at T0
      traj_ti  = nc.variables['trajectory'][:, ind_ti].data  # (traj, obs)
      lat_ti   = nc.variables['lat'][:, ind_ti].data
      lon_ti   = nc.variables['lon'][:, ind_ti].data
      time_ti  = nc.variables['time'][:, ind_ti].data
      nc.close()

      # time_ti in python format
      time_ini_days  = time_ti[~np.isnan(time_ti)][0]/(3600*24) #from seconds to days
      time_ini_py    = date_ori + time_ini_days
      
#      time_fin_py = date_ori + time_tf[~np.isnan(time_tf)][0]/(3600*24)
#      
#      print('fin date...',mdates.num2date(time_fin_py).strftime("%Y-%m-%d"))
#      
      print('')
      print('ini date...',mdates.num2date(time_ini_py).strftime("%Y-%m-%d"))
      print('')
      
      ''' Search SMAP data for the initial day (T0) of the simulation '''

      #dir_SSS  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_70km/'
      
      dir_SSS  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
      
      lon_smap_180, lat_smap, sss_smap, time_smap = \
              sto.sim_open_smap_v4(time_ini_py, dir_SSS, lonmin, lonmax, latmin, latmax,)

      lon_smap = lon_smap_180 + 360
    
      

      '''
      Interpolate the corresponding SSS_SMAP data to each particle initial position
      '''
  
      # Griddata to interpolate the SMAP data onto the particle positions

      sss_int_nonan = interpolate.griddata((lon_smap.ravel(), lat_smap.ravel()), sss_smap.ravel(),
                            (lon_ti, lat_ti), method='linear')  

      # mask values at the coast
    
      sss_int = np.copy(sss_int_nonan)
      sss_int[sss_int<30] = np.nan
    
      smin = np.nanmin(sss_int)
      smax = np.nanmax(sss_int)
      

      ''' figure original SMAP data '''
      
#      plt.figure()
#      pc = plt.pcolor(lon_smap[0,:], lat_smap[:,0],sss_smap, 
#                    vmin=smin, vmax=smax, 
#                    cmap=plt.cm.jet)
#      plt.axis('image')
#      plt.xlim(lonmin+ 360 - 2, lonmax+ 360+4)
#      plt.ylim(latmin-3, latmax+4)
#      plt.title('SSS SMAP on ' + mdates.num2date(time_smap).strftime("%Y-%m-%d"))  
#      #plt.tight_layout()
#      plt.colorbar()
#      if save == True:
#        plt.savefig(dirfig + 'simOl_BF_alt_' + filename[-18:-3] + \
#                '_SSS_SMAP.png', dpi=200) 



      '''
      plot the advected field at T0 
      (SMAP interpolated at particle position)
      '''
      
#      plt.figure()
#      sc=plt.scatter(lon_ti, lat_ti, c=sss_int, 
#                  vmin=smin, vmax=smax, 
#                cmap=plt.cm.jet)
#      plt.scatter(lon_ti[sss_int<25], lat_ti[sss_int<25], c='k', marker='+')
#      plt.axis('image')
#      plt.xlim(lonmin+ 360 - 2, lonmax+ 360+4)
#      plt.ylim(latmin-3, latmax+4)
#      plt.title('SSS adv at Ti running time=' + filename[-9:-3] )  
#      #plt.tight_layout()
#      plt.colorbar(sc)
##      if save == True:
##        plt.savefig(dirfig + 'simOl_B_alt_' + filename[-18:-3] + 
##                '_L'+ '%02i'% scale +'km_Sadv_Tf.png', dpi=200)
        
      
      ''' plot the advected field at TF '''
#      plt.figure()
#      plt.scatter(lon_tf, lat_tf, c=sss_int, 
#                  vmin=smin, vmax=smax, 
#                cmap=plt.cm.jet)
#      plt.axis('image')
#      plt.xlim(lonmin+ 360 - 2, lonmax+ 360+4)
#      plt.ylim(latmin-3, latmax+4)
#      plt.title('SSS adv at TF running time=' + filename[-9:-3] )  
#      #plt.tight_layout()
#      plt.colorbar()
#      if save == True:
#        plt.savefig(dirfig + 'simOl_B_alt_' + filename[-18:-3] + '_Sadv_Tf.png', dpi=200)
#     
      
      '''
      Save SSS from SMAP to the .nc file
      '''

      # If a new file doesn't exist, create the file to save simulation data + sss_int
      name_file = dirsave_sss + 'SimOl_B_alt_BD_'+filename[-18:-3]+'_sss.nc'

      if path.exists(name_file) == False:
  
        shutil.copyfile(diroutputs + filename, name_file)
        print('')
        print('file for sss created,...', name_file)

        # Create variable
        nc = netcdf.Dataset(name_file, 'a', format='NETCDF3_CLASSIC')

        # create variable...
        nc.createVariable('sss_adv',  'f8', ('traj'))
  
        # Write in variable attributes...
        nc.variables['sss_adv'].long_name       = 'SSS advected'
        nc.variables['sss_adv'].units           = '1e-3'

        nc.close() 
    
        
      else:
        print('ATENTION: Remove previous file first!')     
        print('DATA NOT SAVED!!')

      # Save data...    
      nc = netcdf.Dataset(name_file, 'a', format='NETCDF3_CLASSIC')
      nc.variables['sss_adv'][:]   = sss_int
      nc.close() 
      