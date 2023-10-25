#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
from netCDF4                 import Dataset
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
from scipy.ndimage           import gaussian_filter
from scipy.interpolate import griddata
from mpl_toolkits.basemap    import Basemap, shiftgrid
import matplotlib.gridspec   as gridspec
import rt_anatools_3         as rt

'''
Plot the gradient of SMAP data for each simulation day and then
average for each season: winter vs. summer. More gradient in summer?

written by Bàrbara Barceló-Llull on 4-Jan-2020 at Mallorca

adapted on 31-March-2020 for the Bigger Domain, paper figure.

I will include in the figure of the gradient, the mean SMAP SSS. 

FIGURES FOR PAPER!! Figures 7 and 9.

'''


def download_bathymetry():
    
  # read in etopo5 topography/bathymetry.
  url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
  etopodata = Dataset(url)

  topoin = etopodata.variables['ROSE'][:]
  lons = etopodata.variables['ETOPO05_X'][:]
  lats = etopodata.variables['ETOPO05_Y'][:]
  # shift data so lons go from -180 to 180 instead of 20 to 380.
  topoin,lons = shiftgrid(180.,topoin,lons,start=False) 

  return topoin, lons, lats

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
      
#      print('')
#      print('SSS - center of the product time interval...') 
#      print(mdates.num2date(time).strftime("%Y-%m-%d %H:%M"))  
      
      return time

def sim_delimit_smap_region_1degmore(lonmin, lonmax, latmin, latmax, \
                      lon_smap_glo, lat_smap_glo, sss_smap_glo):
    
    ''' 
    Function that delimits the SMAP data inside the simulation domain.
    '''
    
    # delimit SMAP data to the simulation area
    lonmin = lonmin -1
    lonmax = lonmax +1
    latmin = latmin -1
    latmax = latmax +1
    
    ilon = np.logical_and(lon_smap_glo[0,:]>=lonmin, lon_smap_glo[0,:]<=lonmax)
    ilat = np.logical_and(lat_smap_glo[:,0]>=latmin, lat_smap_glo[:,0]<=latmax)
    
    lon_smap = lon_smap_glo[ilat,:][:,ilon]
    lat_smap = lat_smap_glo[ilat,:][:,ilon]
    sss_smap = sss_smap_glo[ilat,:][:,ilon]
    
    return lon_smap, lat_smap, sss_smap   


def plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize):
  
      
      # decor
      bm.drawcoastlines()
      bm.drawrivers(zorder=6)
      bm.fillcontinents(color='0.8', lake_color='0.8', zorder=5)
      
      
      parallels = np.arange(latmin,latmax,5)
      bm.drawparallels(parallels,labels=[1, 0, 0, 0],
                             fontsize=fsize-1, linewidth=0.1, zorder=8)
      meridians = np.arange(-80, lonmax, 5)
      #meridians = [-72, -70, -68]
      bm.drawmeridians(meridians,labels=[0, 0, 0, 1],
                             fontsize=fsize-1, linewidth=0.1, zorder=9)
      

def extract_data_for_season_new(season_months, var, time_all):  
  
  '''
  In the new version of this function, var is a list 
  with len(var) = len(time_all)
  
  Also, var_season will be a list.
  
  '''  
  print('extrating data for the months... ', season_months)  
  
  # append data for the desired season
  var_season  = [] #np.zeros(shape=(var.shape[0],var.shape[1]))
  time_season = [] 
  
  for mm in season_months:
    for it, tt in enumerate(time_all):
        
        time_object = mdates.num2date(tt)
        
        if np.logical_or(
                np.logical_or(time_object.month == season_months[0],
                              time_object.month == season_months[1]),
                              time_object.month == season_months[2]):
                    
  
            var_season.append(var[it])
            time_season.append(tt)
            
  return var_season, time_season  

def extract_data_for_season_3D(season_months, var, time_all):  
  
  '''
  In the new version of this function, var is a list 
  with len(var) = len(time_all)
  
  Also, var_season will be a list.
  
  '''  
  print('extrating data for the months... ', season_months)  
  
  # append data for the desired season
  var_season  = [] #np.zeros(shape=(var.shape[0],var.shape[1]))
  time_season = [] 
  
  for mm in season_months:
    for it, tt in enumerate(time_all):
        
        time_object = mdates.num2date(tt)
        
        if np.logical_or(
                np.logical_or(time_object.month == season_months[0],
                              time_object.month == season_months[1]),
                              time_object.month == season_months[2]):
                    
  
            var_season.append(var[:,:, it])
            time_season.append(tt)
            
  return var_season, time_season  
      
if __name__ == '__main__':
  
  plt.close('all')
    
  dir_SSS  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  dirfig   = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_weekly_sims_BiggerDomain/SMAP_maps_gradient/' 

  

  # BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46
  
  lat_mean = np.mean([latmin, latmax])
  
  ''' convertion factors: degrees to meters '''
  
  lat_mean_rad = np.nanmean(lat_mean) * np.pi/180 # lat mean in radians
  factor_diflon_deg2m, factor_diflat_deg2m = to.length_lon_lat_degs(lat_mean_rad)
       

  ''' File with the simulation dates '''
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_weekly_all/'
  file_dic = 'sim_weekly_relase_dates_corrected_gaps.pkl'

  f = open(dir_dic + file_dic, 'rb')
  dict_release_dates = pickle.load(f)
  f.close() 

  dates_release_all = dict_release_dates['date_release'][:-9] #- 0.5 #starting at 00:00 instead than at 12:00
  
  '''
  Create the basemap
  '''
  bm = Basemap(projection = 'merc',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')      

  mag_all  = []
  date_all = []
  sss_all  = []

  ''' Download bathymetry '''
  
  topo_global, lonbat, latbat = download_bathymetry()
  

  #limit data within domain
  
  inds_lat = np.logical_and(latbat >= 27.5, latbat <= latmax)
  inds_lon = np.logical_and(lonbat >= lonmin, lonbat <= lonmax)
  
  lonbat2d, latbat2d = np.meshgrid(lonbat[inds_lon], latbat[inds_lat])
  
  topo = topo_global[inds_lat,:][:, inds_lon]
  
  ''' Open MDT '''
  
  dir_mld  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/MDT_CNES-CLS18_fromAviso/' 
  file_mld = 'dataset-mdt-cnes-cls18-global_1586181917100.nc'

  nc    = netcdf.Dataset(dir_mld + file_mld, 'r')
  
  lon_mdt  = nc.variables['longitude'][:]
  lat_mdt  = nc.variables['latitude'][:]
  u_mdt    = nc.variables['u'][:]
  v_mdt    = nc.variables['v'][:]
  mdt      = nc.variables['mdt'][:]
  
  nc.close()

  ''' Open altimetry file '''
  
  dir_netcdf    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SLA_GS_big_region/' 
  # Data from 2015-03-31 to 2019-05-13
  file_alt_dt = 'dataset-duacs-rep-global-merged-allsat-phy-l4_1583502357526.nc'
  # Data from 2019-05-14 to 2020-01-15
  file_alt_nrt = 'dataset-duacs-nrt-global-merged-allsat-phy-l4_1583920119159.nc'
  
  time_dt, lat_dt, lon_dt, adt_dt, ugos_dt, \
        vgos_dt, sla_dt, ugosa_dt, vgosa_dt = to.read_altimetry(dir_netcdf + file_alt_dt)

  time_nrt, lat_nrt, lon_nrt, adt_nrt, ugos_nrt, \
        vgos_nrt, sla_nrt, ugosa_nrt, vgosa_nrt = to.read_altimetry(dir_netcdf + file_alt_nrt) 
   
  lon_alt = lon_dt
  lat_alt = lat_dt

  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0))  
  date_ori_nrt = mdates.date2num(datetime(2019, 5, 14, 0, 0, 0))  

  
  ugeo_all = np.zeros(shape=lon_alt.shape)
  vgeo_all = np.zeros(shape=lon_alt.shape)
  sla_all  = np.zeros(shape=lon_alt.shape)
  adt_all  = np.zeros(shape=lon_alt.shape)
  
  for date in dates_release_all:
      
    date_time_obj = mdates.num2date(date)
    
    print('')
    print('Date of simulation...', date_time_obj)


    ''' Open and save altimetry data '''
   
    if date < date_ori_nrt:
       
          date_ori = date_ori_dt
          
          # altimetry data for this date
          ind_alt = np.where(time_dt == int(date))
          
          sla_date = sla_dt[ind_alt].squeeze()
          adt_date = adt_dt[ind_alt].squeeze()
          ug_date  = ugos_dt[ind_alt].squeeze()
          vg_date  = vgos_dt[ind_alt].squeeze()
          adt_date = adt_dt[ind_alt].squeeze()
          
          # remove mask
          
          ug_date_nomask = np.copy(ug_date.data)
          ug_date_nomask[ug_date.mask==True]=np.nan

          vg_date_nomask = np.copy(vg_date.data)
          vg_date_nomask[vg_date.mask==True]=np.nan 
          
          sla_date_nomask = np.copy(sla_date.data)
          sla_date_nomask[sla_date.mask==True]=np.nan   

          adt_date_nomask = np.copy(adt_date.data)
          adt_date_nomask[adt_date.mask==True]=np.nan   
          
          # save into a big array
          ugeo_all = np.dstack((ugeo_all, ug_date_nomask))
          vgeo_all = np.dstack((vgeo_all, vg_date_nomask))
          sla_all  = np.dstack((sla_all,  sla_date_nomask))
          adt_all  = np.dstack((adt_all,  adt_date_nomask))
          
    elif  date > date_ori_nrt:
          date_ori = date_ori_nrt
          
          ind_alt = np.where(time_nrt == int(date))
          
          sla_date = sla_nrt[ind_alt].squeeze()
          adt_date = adt_nrt[ind_alt].squeeze()
          ug_date  = ugos_nrt[ind_alt].squeeze()
          vg_date  = vgos_nrt[ind_alt].squeeze()
          adt_date = adt_nrt[ind_alt].squeeze()

          # remove mask
          
          ug_date_nomask = np.copy(ug_date.data)
          ug_date_nomask[ug_date.mask==True]=np.nan

          vg_date_nomask = np.copy(vg_date.data)
          vg_date_nomask[vg_date.mask==True]=np.nan 
          
          sla_date_nomask = np.copy(sla_date.data)
          sla_date_nomask[sla_date.mask==True]=np.nan   
          
          adt_date_nomask = np.copy(adt_date.data)
          adt_date_nomask[adt_date.mask==True]=np.nan  
          
         # save into a big array
          ugeo_all = np.dstack((ugeo_all, ug_date_nomask))
          vgeo_all = np.dstack((vgeo_all, vg_date_nomask))
          sla_all  = np.dstack((sla_all,  sla_date_nomask))
          adt_all  = np.dstack((adt_all,  adt_date_nomask))


    ''' Search SMAP data for this date '''


    dir_SSS  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
      
    lon_smap_180, lat_smap, sss_smap, time_smap = \
              sto.sim_open_smap_v4(date, dir_SSS, lonmin, lonmax, latmin, latmax,)

    lon_smap = lon_smap_180 #lon_smap_180 + 360
                     
    ''' Magnitude of the gradient '''
          
    # grid spacing in km
    dx_deg = 0.25 # in deg
    dy_deg = 0.25 # in deg
          
    dx_km = dx_deg * factor_diflon_deg2m /1000 #in km
    dy_km = dy_deg * factor_diflat_deg2m /1000 #in km
          
   
    vgrad = np.gradient(sss_smap, dy_km, dx_km) 
    # vgrad2 = np.gradient(sss_smap) = vgrad2[0]/dy_m vgrad2[1]/dx_m
    mag = np.sqrt(vgrad[0]**2 + vgrad[1]**2)  
          
#          sss_prova= np.array([np.arange(0,5), np.arange(0,5)])
#          grad_prova = np.gradient(sss_prova) #grad[0] is dlat and grad[1] is dlon
#          grad_prova2 = np.gradient(sss_prova, 2, 4)
          
          # figure to check
#          x_smap, y_smap = bm(lon_smap, lat_smap)
#          
#          
#          plt.figure()
#          plt.pcolor(x_smap, y_smap, mag, cmap=plt.cm.jet)
#          plt.colorbar()
#          #plt.streamplot(x_smap, y_smap, vgrad[0], vgrad[1], color='w')
#          plt.quiver(x_smap, y_smap, vgrad[0], vgrad[1], color='w')
#          
#          plot_decor_map(bm,lonmin, lonmax, latmin, latmax)
      
        
    ''' Save gradient magnitude for each simulation date '''
    mag_all.append(mag)
    date_all.append(date)
    sss_all.append(sss_smap)
          


  ''' Temporal average '''

  # convert list of 2d array into 3d array
  mag_smap_all_3D = np.dstack(mag_all)
  sss_smap_all_3D = np.dstack(sss_all)
  
  mag_smap_all_mean  = np.nanmean(mag_smap_all_3D, axis=2)
  sss_smap_all_mean  = np.nanmean(sss_smap_all_3D, axis=2)


  ''' Figure temporal mean GRADIENT AND MEAN for PAPER'''
  x_smap, y_smap = bm(lon_smap, lat_smap)
  lon_mdt_2d, lat_mdt_2d = np.meshgrid(lon_mdt, lat_mdt)
  x_mdt, y_mdt   = bm(lon_mdt_2d, lat_mdt_2d)
  xbat, ybat     = bm(lonbat2d, latbat2d)
  
  
  fsize = 12
  
  # find coordinates of the desired 200 m isobath
  plt.figure()
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  # isobath
  cs1 = plt.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=0.9,
              linestyles='-')  
  seg = 0
  isob = cs1.allsegs[0][seg]

  # plot isobath
  plt.plot(isob[:,0], isob[:,1], 'r',linewidth=0.9,
              linestyle='-')  
  plt.close()
  

  # FIGURE!!
  
  fig = plt.figure(figsize=(8,8/2))
  
  wp  = 0.44
  hp  = 0.9 - 0.01
  ws  = wp + 0.08
  wcb = 0.015
  hcb = 0.85  - 0.01
  
  dd = 0.05
  
  ax0 = plt.axes([0, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  axcb0 = plt.axes([wp-0.015 , (1-hcb)/2, wcb, hcb])


  ax1 = plt.axes([0.5-0.01, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  axcb1 = plt.axes([0.5-0.01 + wp-0.015 , (1-hcb)/2, wcb, hcb])
#  cs1 = ax1.contour(xbat, ybat, topo, levels=[-1000],
#              colors='c',linewidths=0.9,
#              linestyles='-')    
  # plot 200 m isobath
  lw=0.9
  ax1.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')

  # plot MDT = 0.5 m (Gulf Stream mean position)
  
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed')   
  
  pc0=ax1.pcolor(x_smap, y_smap, mag_smap_all_mean, vmin=0.003, vmax=0.03, 
                 cmap=plt.cm.magma_r)
  
  # qv0 = ax1.quiver(x_mdt[::5,::5], y_mdt[::5,::5], 
  #            u_mdt.squeeze()[::5,::5], v_mdt.squeeze()[::5,::5],
  #            color='b', scale=9)
  # ax1.quiverkey(qv0, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)
  
  ax1.set_title('S$_{SMAP}$ grad. mag.', fontsize=fsize)     


  #ax1.clabel(cs1, cs1.levels,  fmt='%1.0f', fontsize=fsize-4)#, inline=True, fontsize=fsize-2)
 
  pc1 = ax0.pcolor(x_smap, y_smap, sss_smap_all_mean, vmin=32, vmax=37,
             cmap=plt.cm.Spectral_r) # Spectral_r

  qv1 = ax0.quiver(x_mdt[::5,::5], y_mdt[::5,::5], 
             u_mdt.squeeze()[::5,::5], v_mdt.squeeze()[::5,::5],
             color='b', scale=9)
  ax0.quiverkey(qv1, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)

  ax0.set_title('S$_{SMAP}$', fontsize=fsize)  
  
  
  cb0 = plt.colorbar(pc0, cax=axcb1)
  cb0.ax.set_title('[PSU/km]', fontsize=fsize-1, pad=10)
  cb0.ax.tick_params(labelsize=fsize-1) 

  cb1 = plt.colorbar(pc1, cax=axcb0)
  cb1.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb1.ax.tick_params(labelsize=fsize-1)   
  

  fig.text(0.01, 0.96, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.52+0.01, 0.96, '(b)' , fontsize = fsize+1, zorder=200)
             
  #plt.tight_layout()
  fig.savefig(dirfig + 'SMAP_temp_avg_grad_and_sss_BiggerDomain.png', dpi=200)
  
  
       
  ''' late winter vs. summer gradient averages  '''
  
  winter_months = [2, 3, 4]
  summer_months = [6, 7, 8]
  
  # smap grad winter
  mag_smap_winter, time_smap_winter  = extract_data_for_season_new( \
                                        winter_months,  \
                                        mag_all, \
                                        date_all)
  
  # smap grad summer
  mag_smap_summer, time_smap_summer  = extract_data_for_season_new( \
                                        summer_months,  \
                                        mag_all, \
                                        date_all)             
  
  # convert list of 2d array into 3d array
  mag_smap_winter_3D = np.dstack(mag_smap_winter)
  mag_smap_summer_3D = np.dstack(mag_smap_summer)
  
  # Plot grad means
  mag_smap_winter_mean = np.nanmean(mag_smap_winter_3D, axis=2)
  mag_smap_summer_mean = np.nanmean(mag_smap_summer_3D, axis=2)


  ''' late winter vs. summer SSS averages  '''
  
  # smap grad winter
  sss_smap_winter, time_smap_winter  = extract_data_for_season_new( \
                                        winter_months,  \
                                        sss_all, \
                                        date_all)
  
  # smap grad summer
  sss_smap_summer, time_smap_summer  = extract_data_for_season_new( \
                                        summer_months,  \
                                        sss_all, \
                                        date_all)             
  
  # convert list of 2d array into 3d array
  sss_smap_winter_3D = np.dstack(sss_smap_winter)
  sss_smap_summer_3D = np.dstack(sss_smap_summer)
  
  # Plot grad means
  sss_smap_winter_mean = np.nanmean(sss_smap_winter_3D, axis=2)
  sss_smap_summer_mean = np.nanmean(sss_smap_summer_3D, axis=2)  


  ''' late winter vs. summer (ugeo, vgeo)  '''
  
  ugeo_winter, time_ugeo_winter  = extract_data_for_season_3D( \
                                        winter_months,  \
                                        ugeo_all, \
                                        date_all)
  
  vgeo_winter, time_vgeo_winter  = extract_data_for_season_3D( \
                                        winter_months,  \
                                        vgeo_all, \
                                        date_all)             

  ugeo_summer, time_ugeo_summer  = extract_data_for_season_3D( \
                                        summer_months,  \
                                        ugeo_all, \
                                        date_all)
  
  vgeo_summer, time_vgeo_summer  = extract_data_for_season_3D( \
                                        summer_months,  \
                                        vgeo_all, \
                                        date_all)      
  # convert list of 2d array into 3d array
  ugeo_winter_3D = np.dstack(ugeo_winter)
  vgeo_winter_3D = np.dstack(vgeo_winter)
  ugeo_summer_3D = np.dstack(ugeo_summer)
  vgeo_summer_3D = np.dstack(vgeo_summer)
  
  # Plot grad means
  ugeo_winter_mean = np.nanmean(ugeo_winter_3D, axis=2)
  vgeo_winter_mean = np.nanmean(vgeo_winter_3D, axis=2)
  
  ugeo_summer_mean = np.nanmean(ugeo_summer_3D, axis=2) 
  vgeo_summer_mean = np.nanmean(vgeo_summer_3D, axis=2) 
  

  ''' late winter vs. summer ADT '''
  
  adt_winter, time_adt_winter  = extract_data_for_season_3D( \
                                        winter_months,  \
                                        adt_all, \
                                        date_all)
       

  adt_summer, time_adt_summer  = extract_data_for_season_3D( \
                                        summer_months,  \
                                        adt_all, \
                                        date_all)
  
 
  # convert list of 2d array into 3d array
  adt_winter_3D = np.dstack(adt_winter)
  adt_summer_3D = np.dstack(adt_summer)

  # Plot seasonal average of ADT
  adt_winter_mean = np.nanmean(adt_winter_3D, axis=2)
  
  adt_summer_mean = np.nanmean(adt_summer_3D, axis=2) 
  
  
  plt.figure()
  plt.contour(lon_alt, lat_alt, adt_winter_mean)
  plt.contour(lon_alt, lat_alt, adt_winter_mean, levels=[0.5],
             colors= 'r', linewidths=lw, linestyles='dashed')   
  plt.title('winter')
  
  plt.figure()
  plt.contour(lon_alt, lat_alt, adt_summer_mean)  
  plt.contour(lon_alt, lat_alt, adt_summer_mean, levels=[0.5],
             colors= 'r', linewidths=lw, linestyles='dashed')     
  plt.title('summer')
  
  plt.figure()
  plt.contour(lon_alt, lat_alt, adt_winter_mean, levels=[0.5],
             colors= 'r', linewidths=lw, linestyles='dashed')  
  plt.contour(lon_alt, lat_alt, adt_summer_mean, levels=[0.5],
             colors= 'b', linewidths=lw, linestyles='dashed')     
  
  
 

  ''' Figure winter vs. summer GRADIENT AND MEAN for PAPER'''
  x_alt, y_alt = bm(lon_alt, lat_alt)
 
      
      
  fsize = 10
  fig = plt.figure(figsize=(7,8))
   
  gs = gridspec.GridSpec(2, 3, width_ratios=[20, 20, 1])
  
  # GRADIENT 
  
  ax0 = plt.subplot(gs[1,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  ax1 = plt.subplot(gs[1, 1]) 
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  
  axcb = plt.subplot(gs[1, 2])

  # >> smap grad winter <<
  
  pc1=ax0.pcolor(x_smap, y_smap, mag_smap_winter_mean, vmin=0.003, vmax=0.03, 
                 cmap=plt.cm.magma_r)
  ax0.set_title('S$_{SMAP}$ grad. mag. - Winter', fontsize=fsize)     

  # isobath
  cs0 = ax0.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  
  # adt winter
  ax0.contour(x_alt, y_alt, adt_winter_mean, levels=[0.5],
             colors= 'r', linewidths=lw, linestyles='dashed')    
  
  # >> smap grad summer <<
  
  ax1.pcolor(x_smap, y_smap, mag_smap_summer_mean, vmin=0.003, vmax=0.03,
             cmap=plt.cm.magma_r) # Spectral_r
  ax1.set_title('S$_{SMAP}$ grad. mag. - Summer', fontsize=fsize)  
  
  # isobath
  cs1 = ax1.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-') 
  
  # adt summer
  ax1.contour(x_alt, y_alt, adt_summer_mean, levels=[0.5],
             colors= 'b', linewidths=lw, linestyles='dashed')   
  
  cb = plt.colorbar(pc1, cax=axcb)
  cb.ax.set_title('[PSU/km]', fontsize=fsize-1, pad=10)
  cb.ax.tick_params(labelsize=fsize-1) 
  
  
  # >> smap seasonal average <<
  
  ax2 = plt.subplot(gs[0,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  ax3 = plt.subplot(gs[0, 1]) 
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  
  axcb2 = plt.subplot(gs[0, 2])

  
  pc2=ax2.pcolor(x_smap, y_smap, sss_smap_winter_mean, vmin=32, vmax=37, 
                 cmap=plt.cm.Spectral_r)
  ax2.set_title('S$_{SMAP}$ - Winter', fontsize=fsize)  
   
  qv0 = ax2.quiver(x_alt[::2,::2], y_alt[::2,::2], 
             ugeo_winter_mean[::2,::2], vgeo_winter_mean[::2,::2],
             color='b',  scale=9)
  ax2.quiverkey(qv0, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)

  
  ax3.pcolor(x_smap, y_smap, sss_smap_summer_mean, vmin=32, vmax=37,
             cmap=plt.cm.Spectral_r) # Spectral_r
  ax3.set_title('S$_{SMAP}$ - Summer', fontsize=fsize)  
  
  qv1 = ax3.quiver(x_alt[::2,::2], y_alt[::2,::2], 
             ugeo_summer_mean[::2,::2], vgeo_summer_mean[::2,::2],
             color='b', scale=9)
  ax3.quiverkey(qv1, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)
  
  
  cb = plt.colorbar(pc2, cax=axcb2)
  cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb.ax.tick_params(labelsize=fsize-1)   

  fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.43, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
  fig.text(0.01, 0.48, '(c)' , fontsize = fsize+1, zorder=200)
  fig.text(0.43, 0.48, '(d)' , fontsize = fsize+1, zorder=200)
             
  plt.tight_layout()
  fig.savefig(dirfig + 'SMAP_seasonal_gradient_and_mean_BiggerDomain.png', dpi=200)
  
  