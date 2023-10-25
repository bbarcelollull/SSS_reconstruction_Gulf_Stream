#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Toolbox for the salinity Lagrangian simulation analysis 

written by Bàrbara Barceló Llull on 26-07-2019 at APL-UW
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
import pickle
from matplotlib              import dates as mdates
from datetime                import datetime, timedelta
from scipy                   import interpolate
import sss_toolbox           as to

def rmsd(x, y):
    
  '''
  Infer the root mean square difference between 2 variables x and y
  '''
    
  return np.sqrt(np.nanmean(np.square(x - y)))

def std(x):
    
  '''
  Infer the standard deviation of the variable x.
  x can be the difference of two variables.
  '''
    
  std = np.sqrt(np.nanmean(abs(x - np.nanmean(x))**2))
  return std

def open_back_sim_file(dirsave, filename):
  
  '''
  Open backward simulation file and read initial and final data.
  '''  
    
  ind_tf = 0  #first observation, beginning of simulation at TF
  ind_ti = -1 #last observation, end of simulation at T0
    
  nc    = netcdf.Dataset(dirsave + filename, 'r')
    
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

    
  # SSS tagged to each particle
  sss       = nc.variables['sss_adv'][:].data
  nc.close()
  
  # from 360 format to -180... 180
  lon_ti[lon_ti>180] = lon_ti[lon_ti>180] - 360
  lon_tf[lon_tf>180] = lon_tf[lon_tf>180] - 360
  
  return traj_ti, lat_ti, lon_ti, time_ti, traj_tf, lat_tf, lon_tf, time_tf, sss


def open_sim_file(dirsave, filename):
  
  '''
  Open simulation file and read initial and final data.
  '''  
    
  ind_ini = 0  #first observation, beginning of simulation
  ind_fin = -1 #last observation, end of simulation
  
  nc    = netcdf.Dataset(dirsave + filename, 'r')
  # read initial position
  traj_ini  = nc.variables['trajectory'][:, ind_ini].data  # (traj, obs)
  lat_ini   = nc.variables['lat'][:, ind_ini].data
  lon_ini   = nc.variables['lon'][:, ind_ini].data   
  time_ini  = nc.variables['time'][:, ind_ini].data   #"seconds since 2017-01-01T00:00:00.000000000"

  # read final position
  traj_fin  = nc.variables['trajectory'][:, ind_fin].data  # (traj, obs)
  lat_fin   = nc.variables['lat'][:, ind_fin].data
  lon_fin   = nc.variables['lon'][:, ind_fin].data
  time_fin  = nc.variables['time'][:, ind_fin].data
  # SSS tagged to each particle
  sss       = nc.variables['sss_int'][:].data
  nc.close()
  
  # from 360 format to -180... 180
  lon_ini[lon_ini>180] = lon_ini[lon_ini>180] - 360
  lon_fin[lon_fin>180] = lon_fin[lon_fin>180] - 360
  
  return traj_ini, lat_ini, lon_ini, time_ini, traj_fin, lat_fin, lon_fin, time_fin, sss


def open_tsg(dir_savedic, file_tsg):
    
  '''
  Open TSG data from Kyla data set.
  '''
  
  f = open(dir_savedic + file_tsg, 'rb')
  dict_tsg_all = pickle.load(f)
  f.close() 

  time_tsg_all = dict_tsg_all['time']
  lon_tsg_all  = dict_tsg_all['lon']
  lat_tsg_all  = dict_tsg_all['lat']
  sal_tsg_all  = dict_tsg_all['sal'] # practical salinity
  
  return time_tsg_all, lon_tsg_all, lat_tsg_all, sal_tsg_all


def simtime2datetime(time_ini, time_fin, year):
    
  '''
  Convert time_ini and time_fin to python format
  
  (Now they are in seconds since 1-1-2017 (day 0),
  that is the date of the first altimetry currents that we have used
  to do the simulations for the year 2017) 
  
  We can use the same function for other years of simulations
  as long as the altimetry currents start on 1-1-year.
  '''
  
  # num days of the beggining of the simulation from 1-1-2017 (day 0)
  time_ini_days = time_ini[0]/(3600*24)
  
  # num days of the end of the simulation from from 1-1-2017 (day 0)
  time_fin_days = time_fin[0]/(3600*24) 

  # first date
  d_ini_or    = mdates.date2num(datetime(year, 1, 1, 0, 0, 0))
  
  # Sum these days to d_ini_or
  timepi = d_ini_or + time_ini_days  
  timepf = d_ini_or + time_fin_days  
  
  print('Initial simulation time...', mdates.num2date(timepi))
  print('Final simulation time...', mdates.num2date(timepf))  
  
  return timepi, timepf

def simtime2datetime_allsims(date_ori, time_sim):   
      
      ''' 
      Convert simulation time to python format.
      Input: 
          date_ori: initial date of the currents fields in the simulations
          time_sim: time to convert (in seconds since date_ori)
      '''    
      
      time_days  = np.nanmin(time_sim)/(3600*24) #from seconds to days
      time_py    = date_ori + time_days
      
      return time_py
  

def sim_delimit_tsg(timepf, dt, lonmin, lonmax, latmin, latmax, \
                 time_tsg_all, lon_tsg_all, lat_tsg_all, sal_tsg_all):
  
  '''
  Search the TSG data that correpond to the final date of the simulation 
  +- dt [timepf-dt, timepf+dt] and within the simulation domain.
  '''  
  
  indt = np.logical_and(time_tsg_all>=timepf-dt, time_tsg_all<=timepf+dt) # all day period
  
  time_tsga = time_tsg_all[indt]
  lon_tsga  = lon_tsg_all[indt]
  lat_tsga  = lat_tsg_all[indt]
  sal_tsga  = sal_tsg_all[indt]  
  
  #only within the region of study  
  indlon = np.logical_and(lon_tsga>=lonmin, lon_tsga<=lonmax)

  time_tsg_lon = time_tsga[indlon]
  lon_tsg_lon  = lon_tsga[indlon]
  lat_tsg_lon  = lat_tsga[indlon]
  sal_tsg_lon  = sal_tsga[indlon]  

  indlat = np.logical_and(lat_tsg_lon>=latmin, lat_tsg_lon<=latmax)

  time_tsg = time_tsg_lon[indlat]
  lon_tsg  = lon_tsg_lon[indlat]
  lat_tsg  = lat_tsg_lon[indlat]
  sal_tsg  = sal_tsg_lon[indlat]    
  
  return time_tsg, lon_tsg, lat_tsg, sal_tsg


def sim_file_smap(time_fin):
    
  '''
  Function to find the file of the smap data that corresponds to
  time_fin; this time is in the simulation format.
  
  Only valid for 2017!
  '''
  
  time_fin_days      = time_fin[0]/(3600*24) #from seconds to days
  time_fin_days_smap = time_fin_days + 1     #smap counts the 1-1-2017 as day 1
  day_fins_folder_smap    = time_fin_days_smap - 4   #time of the folder of smap (4 days before the desired data)

  file_SSS = '/2017/' + '%03i'% day_fins_folder_smap + \
           '/RSS_smap_SSS_L3_8day_running_70km_2017_'+ '%03i'% time_fin_days_smap +'_FNL_v03.0.nc'
  
  return file_SSS

def sim_file_smap_all_years(time_py):
    
  '''
  Function to find the file of the smap data that corresponds to
  time_py; this time is in python format.
  
  '''
  
  # Search the year of time_py
  yy = mdates.num2date(time_py).year
  
  # Search the day of the year, smap counts the 1-1-year as day 1
  day_of_year = time_py - mdates.date2num(datetime(yy, 1, 1)) + 1 
  
  #time of the folder of smap (4 days before the desired data)
  day_folder    = day_of_year - 4   

  file_SSS = '/'+np.str(yy)+'/' + '%03i'% day_folder + \
           '/RSS_smap_SSS_L3_8day_running_70km_'+np.str(yy)+'_'+ \
           '%03i'% day_of_year +'_FNL_v03.0.nc'
  
  return file_SSS

def sim_file_smap_all_years_v4(time_py):
    
  '''
  Function to find the file of the smap data that corresponds to
  time_py; this time is in python format.
  
  '''
  
  # Search the year of time_py
  yy = mdates.num2date(time_py).year
  
  # Search the day of the year, smap counts the 1-1-year as day 1
  day_of_year = time_py - mdates.date2num(datetime(yy, 1, 1)) + 1 
  
  if day_of_year > 4:
    #time of the folder of smap (4 days before the desired data)
    day_folder    = day_of_year - 4   
    
    file_SSS = '/'+np.str(yy)+'/' + '%03i'% day_folder + \
           '/RSS_smap_SSS_L3_8day_running_'+np.str(yy)+'_'+ \
           '%03i'% day_of_year +'_FNL_v04.0.nc'  
           
  elif day_of_year <=4:
      
    yy = yy-1
      
    days_year = mdates.date2num(datetime(yy, 12, 31)) - \
                  mdates.date2num(datetime(yy-1, 12, 31))
                  
    day_folder = day_of_year - 4 + days_year
      
    file_SSS = '/'+np.str(yy)+'/' + '%03i'% day_folder + \
           '/RSS_smap_SSS_L3_8day_running_'+np.str(yy+1)+'_'+ \
           '%03i'% day_of_year +'_FNL_v04.0.nc'        

  
  return file_SSS

def sim_delimit_smap_region(lonmin, lonmax, latmin, latmax, \
                      lon_smap_glo, lat_smap_glo, sss_smap_glo):
    
    ''' 
    Function that delimits the SMAP data inside the simulation domain.
    '''
    
    # delimit SMAP data to the simulation area
    lonmin = lonmin -10 
    lonmax = lonmax +10
    latmin = latmin -10
    latmax = latmax +10  
    
    ilon = np.logical_and(lon_smap_glo[0,:]>=lonmin, lon_smap_glo[0,:]<=lonmax)
    ilat = np.logical_and(lat_smap_glo[:,0]>=latmin, lat_smap_glo[:,0]<=latmax)
    
    lon_smap = lon_smap_glo[ilat,:][:,ilon]
    lat_smap = lat_smap_glo[ilat,:][:,ilon]
    sss_smap = sss_smap_glo[ilat,:][:,ilon]
    
    return lon_smap, lat_smap, sss_smap    

  
def sim_open_smap(time_py, dir_SSS, lonmin, lonmax, latmin, latmax):
    
  '''
  Search SMAP file on time_py, open it, delimit SMAP to the simulation domain.
  time_py in python format
  '''  
    
  # SMAP file 
  file_SSS = sim_file_smap_all_years(time_py)

  # Open data
  lon_smap_glo, lat_smap_glo, time_smap, sss_smap_glo = to.read_SMAP(dir_SSS + file_SSS)

  # delimit SMAP domain 
  lon_smap, lat_smap, sss_smap = \
                      sim_delimit_smap_region(lonmin, lonmax, latmin, latmax, \
                      lon_smap_glo, lat_smap_glo, sss_smap_glo)
  
  return lon_smap, lat_smap, sss_smap, time_smap

def sim_open_smap_v4(time_py, dir_SSS, lonmin, lonmax, latmin, latmax):
    
  '''
  Search SMAP file on time_py, open it, delimit SMAP to the simulation domain.
  time_py in python format
  '''  
    
  # SMAP file 
  file_SSS = sim_file_smap_all_years_v4(time_py)

  # Open data
  lon_smap_glo, lat_smap_glo, time_smap, sss_smap_glo = \
                                   to.read_SMAP(dir_SSS + file_SSS)

  # delimit SMAP domain 
  lon_smap, lat_smap, sss_smap = \
                      sim_delimit_smap_region(lonmin, lonmax, latmin, latmax, \
                      lon_smap_glo, lat_smap_glo, sss_smap_glo)
  
  return lon_smap, lat_smap, sss_smap, time_smap


def search_domain(date_fin):
      
  # Domain for the simulation

  if date_fin == mdates.date2num(datetime(2017, 2, 10, 0, 0, 0)): #offshore
      
    lonmin, lonmax = -138, -130
    latmin, latmax =  28, 33
    
  elif date_fin == mdates.date2num(datetime(2017, 2, 11, 0, 0, 0)):
  
    lonmin, lonmax = -135, -125  
    latmin, latmax =  29, 34
  
    
  elif date_fin == mdates.date2num(datetime(2017, 2, 12, 0, 0, 0)): #nearshore
      
    lonmin, lonmax = -127, -122
    latmin, latmax =  32, 37

  elif date_fin == mdates.date2num(datetime(2017, 4, 2, 0, 0, 0)): #nearshore
      
    lonmin, lonmax = -126, -121 #-125, -121
    latmin, latmax =  28, 33 #29, 32

  elif date_fin == mdates.date2num(datetime(2017, 4, 12, 0, 0, 0)): #nearshore
      
    lonmin, lonmax = -127, -121#-125.5, -122
    latmin, latmax =  30, 36 #32, 34.5

  return lonmin, lonmax, latmin, latmax


def search_domain_sim_chl(date_fin, satellite):
      
  # Domain for the simulation

  if date_fin == mdates.date2num(datetime(2017, 2, 15, 0, 0, 0)): 
      
    if satellite == 'modis':
        
      lonmin, lonmax = -151, -146
      latmin, latmax = 31, 33
      
    elif satellite == 'viirs':

      lonmin, lonmax = -135, -131.5
      latmin, latmax = 25, 28 

  return lonmin, lonmax, latmin, latmax

def file_l3_smap_from_date(date_fin):
  
  '''
  Function to find the file of the L3 smap data that corresponds to
  a desired date in datetime format.
  '''

  year        = date_fin.year
  day_of_year = date_fin.timetuple().tm_yday
  
  day_folder_smap    = day_of_year - 4   #time of the folder of smap (4 days before the desired data)

  file_SSS = '/'+np.str(year)+'/' + '%03i'% day_folder_smap + \
           '/RSS_smap_SSS_L3_8day_running_70km_2017_'+ '%03i'% day_of_year +'_FNL_v03.0.nc'
  
  return file_SSS


def read_l3_smap_from_date(date_fin, dir_SSS, lonmin, lonmax, latmin, latmax):
    
  '''
  Search SMAP file, open it, delimit SMAP to the simulation domain.
  '''  
    
  # SMAP file 
  file_SSS = file_l3_smap_from_date(date_fin)

  # Open data
  lon_smap_glo, lat_smap_glo, time_smap, sss_smap_glo = to.read_SMAP(dir_SSS + file_SSS)

  # delimit SMAP domain 
  lon_smap, lat_smap, sss_smap = \
                      sim_delimit_smap_region(lonmin, lonmax, latmin, latmax, \
                      lon_smap_glo, lat_smap_glo, sss_smap_glo)
  
  return lon_smap, lat_smap, sss_smap, time_smap

def open_sim_backward(file, date_ori):  
    
      ind_tf = 0  #first observation, beginning of simulation at TF
      ind_ti = -1 #last observation, end of simulation at T0
    
      nc    = netcdf.Dataset(file, 'r')
    
      # read initial position at TF
      traj_tf  = nc.variables['trajectory'][:, ind_tf].data  # (traj, obs)
      lat_tf   = nc.variables['lat'][:, ind_tf].data
      lon_tf   = nc.variables['lon'][:, ind_tf].data   
      time_tf  = nc.variables['time'][:, ind_tf].data   #seconds since date_ori

      # read final position at T0
      traj_ti  = nc.variables['trajectory'][:, ind_ti].data  # (traj, obs)
      lat_ti   = nc.variables['lat'][:, ind_ti].data
      lon_ti   = nc.variables['lon'][:, ind_ti].data
      time_ti  = nc.variables['time'][:, ind_ti].data  #seconds since date_ori

      # SSS tagged to each particle
      sss       = nc.variables['sss_adv'][:].data
      nc.close()

      time_ini_py    = simtime2datetime_allsims(date_ori, time_ti)
      time_fin_py    = simtime2datetime_allsims(date_ori, time_tf)   
  
      # from 360 format to -180... 180
      lon_tf[lon_tf>180] = lon_tf[lon_tf>180] - 360
      lon_ti[lon_ti>180] = lon_ti[lon_ti>180] - 360
      
      # mask where sss is -9999 (salinity coming from the land)
      sss[sss<-100] = np.nan

        
      print('')
      print('ini date...',mdates.num2date(time_ini_py).strftime("%Y-%m-%d"))
      print('')
      print('fin date...',mdates.num2date(time_fin_py).strftime("%Y-%m-%d"))
      print('')
      
      return lat_tf, lon_tf, time_fin_py, lat_ti, lon_ti, time_ini_py, sss