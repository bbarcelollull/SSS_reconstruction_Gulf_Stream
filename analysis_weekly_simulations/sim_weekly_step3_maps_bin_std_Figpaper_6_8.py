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
import matplotlib.gridspec   as gridspec
import glob
import sss_toolbox           as to
from scipy.interpolate       import griddata
from scipy                   import stats
from mpl_toolkits.basemap    import Basemap, shiftgrid
import rt_anatools_3         as rt
from matplotlib.ticker       import FormatStrFormatter
from netCDF4                 import Dataset

"""
Bin analysis: Horizontal STD in bins of different sizes.


Figures for PAPER!! Figures 6 and 8.

written by Bàrbara Barceló-Llull on 2-April-2020 in Mallorca

"""

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

def infer_binned_std(px, py, val, binx, biny):    

    
    cond_nans = np.isnan(val)
    
    px_nonans  = px[~cond_nans]
    py_nonans  = py[~cond_nans]    
    val_nonans = val[~cond_nans]

    
    stats_variable = stats.binned_statistic_2d(px_nonans, py_nonans,  
                                         values=val_nonans, 
                                         statistic='std',
                                         bins=[binx, biny])
    std_variable = stats_variable.statistic.T
    
    return std_variable

def infer_binned_count(px, py, val, binx, biny):    

    
    cond_nans = np.isnan(val)
    
    px_nonans  = px[~cond_nans]
    py_nonans  = py[~cond_nans]    
    val_nonans = val[~cond_nans]

    
    stats_variable = stats.binned_statistic_2d(px_nonans, py_nonans,  
                                         values=px_nonans, 
                                         statistic='count',
                                         bins=[binx, biny])
    count_variable = stats_variable.statistic.T
    
    return count_variable
  
def plot_temp_evol_std(loni, lati):
    
  # center of the bins
  
  binx_center = binx[:-1] + dx/2
  biny_center = biny[:-1] + dy/2
 
  binx_center_s = binx_s[:-1] + dxsmall/2
  biny_center_s = biny_s[:-1] + dysmall/2
  
  
  indx = np.argmin(np.abs(binx_center-loni))
  indy = np.argmin(np.abs(biny_center-lati))

  indx_s = np.argmin(np.abs(binx_center_s-loni))
  indy_s = np.argmin(np.abs(biny_center_s-lati))

  print('')
  print('Plotting data at...')
  print('')
  print('Big bins... ')
  print('Lat=', binx_center[indx])
  print('Lon=', biny_center[indy])
  print('')
  print('Small bins... ')
  print('Lat=', binx_center_s[indx_s])
  print('Lon=', biny_center_s[indy_s])
  
  fig = plt.figure(figsize=(10,8))
  ax1=plt.subplot(211)
  plt.plot(time_adv_all, std_adv_all[indy, indx, :], label='SSS$_{adv}$ dx=' + np.str(dx))
  plt.plot(time_smap_all, std_smap_all[indy, indx, :], label='SSS$_{SMAP}$ dx=' + np.str(dx))
  plt.legend()
  plt.ylabel('2D STD salinity')
  plt.title('lon = ' + np.str(loni) + '  lat = ' + np.str(lati))
  # Format the y axis
  ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
  ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))  #%b month: Jan
 
  ax2=plt.subplot(212)   
  plt.plot(time_adv_all, std_adv_s_all[indy_s, indx_s, :], label='SSS$_{adv}$ dx=' + np.str(dxsmall))
  plt.legend()
  plt.ylabel('2D STD salinity')
  ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
  ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
  
  fig.savefig(dirfig + 'maps_SSS_comp_2D_STD_temp_evol'+ 
                       '_loni_' + np.str(loni) +
                       '_lati_' + np.str(lati) + '.png', dpi=200) 
    
def nan_where_std_0(std_adv_s):       
   # nan where std is 0 (because all nans in bin)
   std_adv_s_masked = np.ma.copy(std_adv_s)
   
   std_adv_s_masked =  np.ma.masked_where(std_adv_s_masked==0, std_adv_s_masked)
   
   std_adv_s_masked.data[std_adv_s_masked.mask==True] = np.nan
   
   return std_adv_s_masked.data


def extract_data_for_season(season_months, var, time_all):  
    
  print('extrating data for the months... ', season_months)  
  
  # append data for the desired season
  var_season  = np.zeros(shape=(var.shape[0],var.shape[1]))
  time_season = [] 
  
  for mm in season_months:
    for it, tt in enumerate(time_all):
        
        time_object = mdates.num2date(tt)
        
        if np.logical_or(
                np.logical_or(time_object.month == season_months[0],
                              time_object.month == season_months[1]),
                              time_object.month == season_months[2]):
                    
            #print(time_object.month)        
            # save data       
            var_season = np.dstack((var_season, var[:,:,it]))
            time_season = np.append(time_season, tt)
  
  return var_season, time_season    
 
def infer_2D_std_for_each_date(binx, biny, time_all, 
                            lon, lat, var):
  
  if  var.shape[2] == time_all.shape[0]:
      print('Correct shape')
      
  else:
      print('Variable and time have different sizes!')
      
  std_2dmap_all = np.zeros(shape=(biny.shape[0]-1, binx.shape[0]-1))

  for ii in np.arange(time_all.shape[0]):

    # sss smap 
    std_2dmap = infer_binned_std(lon.flatten(), 
                                lat.flatten(), 
                                var[:,:,ii].flatten(),
                                binx, biny)

    # save into a big array
    std_2dmap_all = np.dstack((std_2dmap_all, 
                              nan_where_std_0(std_2dmap)))

  return std_2dmap_all[:,:,1:]


def infer_2D_count_for_each_date(binx, biny, time_all, 
                            lon, lat, var):
  
  if  var.shape[2] == time_all.shape[0]:
      print('Correct shape')
      
  else:
      print('Variable and time have different sizes!')
      
  count_2dmap_all = np.zeros(shape=(biny.shape[0]-1, binx.shape[0]-1))

  for ii in np.arange(time_all.shape[0]):

    # sss smap 
    count_2dmap = infer_binned_count(lon.flatten(), 
                                lat.flatten(), 
                                var[:,:,ii].flatten(),
                                binx, biny)

    # save into a big array
    count_2dmap_all = np.dstack((count_2dmap_all,count_2dmap))

  return count_2dmap_all[:,:,1:]


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
      

def figure_season_averages(std_summer_mean, std_winter_mean, 
                           binx, biny, smin, smax, name_var, save_name):  

  # fig parameters
  binx_2d, biny_2d = np.meshgrid(binx, biny)
  xx_big, yy_big = bm(binx_2d, biny_2d)
  
  
  fsize = 14
 
  # start figure
  fig = plt.figure(figsize=(12,4))
  gs = gridspec.GridSpec(1, 3, width_ratios=[20, 20, 1])
  
  #define axis and plot decor
  ax0 = plt.subplot(gs[0,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  ax1 = plt.subplot(gs[0, 1]) 
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  axcb = plt.subplot(gs[0, 2])    
    
    
  pc1=ax0.pcolor(xx_big, yy_big, std_winter_mean, vmin=smin, vmax=smax, 
                 cmap=plt.cm.magma_r)

  ax0.set_title(name_var +  ' STD (feb-mar-abr)') 
    
  ax1.pcolor(xx_big, yy_big, std_summer_mean, vmin=smin, vmax=smax, 
                 cmap=plt.cm.magma_r)
 
  ax1.set_title(name_var + ' STD (jun-jul-aug)')  
  
  cb = plt.colorbar(pc1, cax=axcb)
  cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb.ax.tick_params(labelsize=fsize-2)   
  
  fig.savefig(dirfig + 'maps_SSS_STD_seasons_' + save_name + '.png', dpi=200)
  
  
def mask_no_enough_data(count_adv_s_all, std_adv_s_all_mean, count_threshold):
  
  '''
  Mask array where counts are less than threshold times the maximum
  number of data points inside bins.
  
  count_adv_s_all    = counts
  std_adv_s_all_mean = array to mask
  count_threshold    = count threshold in %
  
  '''  
  #counts 
  count_adv_s_sum    = np.nansum(count_adv_s_all, axis=2)
  max_sum_adv_s      = np.nanmax(count_adv_s_sum)
  
  count_adv_s_sum_perc = count_adv_s_sum*100/max_sum_adv_s

  
  mask_counts_adv_s = np.where(count_adv_s_sum_perc<count_threshold)
  
  count_adv_s_sum_perc_mk = np.copy(count_adv_s_sum_perc)
  count_adv_s_sum_perc_mk[mask_counts_adv_s] = np.nan

  # mask std data where bins have less than 80% of the maximum number of data inside bins
  std_adv_s_all_mean_mk = np.copy(std_adv_s_all_mean)
  std_adv_s_all_mean_mk[mask_counts_adv_s] = np.nan
  
  return std_adv_s_all_mean_mk
    
if __name__ == '__main__':
    
  plt.close('all')
    
  ''' Directories backward simulations '''

  dirsave     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_weekly_all_Bigger_Domain/' #save outputs directory
  diroutputs  = dirsave + 'sim_outputs/'
  dirsave_sss = dirsave + 'sss_tagged/'
  
  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_weekly_sims_BiggerDomain/bin_maps/'
    
  dir_SSS  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'

  files    = sorted(glob.glob(dirsave_sss + '*.nc'))
  
  ''' first day of the simulation in the GS '''
  
  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0))  
  date_ori_nrt = mdates.date2num(datetime(2019, 5, 14, 0, 0, 0))   
  
  # DOMAIN FOR weekly and Oleander simulations BIGGER DOMAIN
  lonmin, lonmax = -82, -63 #-74, -66 #-74, -68
  latmin, latmax =  25, 46 #35, 40 #37, 40 
  
  # Grid of the advected particles
  dx = 0.04

  lon_adv, lat_adv = np.meshgrid(np.arange(lonmin, lonmax, dx), np.arange(latmin, latmax, dx))


  '''
  Create the basemap
  '''
  bm = Basemap(projection = 'merc',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')  


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
  
  
  sss_adv_all   = np.zeros(shape=lon_adv.shape) #remove first index 3rd axis
  time_adv_all  = []
  
  sss_smap_all  = np.zeros(shape=(164, 156)) #remove first index 3rd axis
  time_smap_all = []

  date_all = []
  
  for file in files: 
      
      print(file[80:])
      

      date_file = datetime.strptime(file[-22:-14], '%Y%m%d')
      
      time_smap_all = np.append(time_smap_all, 
                                mdates.date2num(date_file))
      
      time_alt_all = np.copy(time_smap_all)
      
      ''' Open and save altimetry data '''
      date = mdates.date2num(date_file)
      date_all.append(date)
      
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

          
      ''' Open simulation data '''
      
      lat_tf, lon_tf, time_fin_py, lat_ti, lon_ti, time_ini_py, sss = \
                                                      sto.open_sim_backward(file, date_ori)

      sss_adv = griddata((lon_tf,lat_tf ), sss, (lon_adv, lat_adv), method='linear')
      
      
      sss_adv_all = np.dstack((sss_adv_all, sss_adv))
      
      time_adv_all = np.append(time_adv_all, time_fin_py)
      
#      plt.figure()
#      plt.pcolor(lon_adv, lat_adv, sss_adv, vmin=35.5, vmax=37)
#      plt.colorbar()
#
#      plt.figure()
#      plt.scatter(lon_tf, lat_tf, c=sss, vmin=35.5, vmax=37)
#      plt.colorbar()      
      
      
      ''' Search SMAP data for the initial and final day of the simulation '''

      lon_smapi, lat_smapi, sss_smapi, time_smapi = \
              sto.sim_open_smap_v4(time_ini_py, dir_SSS, lonmin, lonmax, latmin, latmax,)

      
      lon_smapf, lat_smapf, sss_smapf, time_smapf = \
              sto.sim_open_smap_v4(time_fin_py, dir_SSS, lonmin, lonmax, latmin, latmax,)

      # save final day smap data
      sss_smap_all = np.dstack((sss_smap_all, sss_smapf))
 
  
  ''' Infer 2D horizontal STD in bins with different sizes '''
  
  # 0.6º 
  dx = 0.75 #0.6
  dy = dx
  
#  binx = np.arange(lon_adv.min(), lon_adv.max(), dx)
#  biny = np.arange(lat_adv.min(), lat_adv.max(), dy)
  
  binx = np.arange(lonmin, lonmax+dx, dx)
  biny = np.arange(latmin, latmax+dy, dy) #corners of the boxes
  
  std_smap_all = infer_2D_std_for_each_date(binx, biny, time_smap_all, 
                            lon_smapf, lat_smapf, sss_smap_all[:,:, 1:]) #remove first data!!

#  count_smap_all = infer_2D_count_for_each_date(binx, biny, time_smap_all, 
#                            lon_smapf, lat_smapf, sss_smap_all[:,:, 1:])

  std_adv_all = infer_2D_std_for_each_date(binx, biny, time_adv_all, 
                            lon_adv, lat_adv, sss_adv_all[:,:, 1:])

  count_adv_all = infer_2D_count_for_each_date(binx, biny, time_adv_all, 
                            lon_adv, lat_adv, sss_adv_all[:,:, 1:]) 
  
  
  # 0.25º
  dxsmall = 0.04*6 #to have the same number of data inside each bin 
  dysmall = dxsmall
  
  binx_s = np.arange(lonmin, lonmax+dxsmall, dxsmall)
  biny_s = np.arange(latmin, latmax+dysmall, dysmall)

 
  std_adv_s_all = infer_2D_std_for_each_date(binx_s, biny_s, time_adv_all, 
                            lon_adv, lat_adv, sss_adv_all[:,:, 1:])

  count_adv_s_all = infer_2D_count_for_each_date(binx_s, biny_s, time_adv_all, 
                            lon_adv, lat_adv, sss_adv_all[:,:, 1:])  
  

#  dxsmall = 0.04*6
#  dysmall = dxsmall
#  
#  binx_s = np.arange(lonmin, lonmax+dxsmall, dxsmall)
#  biny_s = np.arange(latmin, latmax+dysmall, dysmall)  
#  scountpr = infer_binned_count(lon_adv, lat_adv, sss_adv_all[:,:, 1:][:,:,0], 
#                                binx_s, biny_s)   
#  
#  plt.figure()
#  plt.pcolor(binx_s, biny_s, scountpr)
#  plt.colorbar()    
  
  ''' Temporal average '''
  
  std_smap_all_mean  = np.nanmean(std_smap_all, axis=2)
  std_adv_all_mean   = np.nanmean(std_adv_all, axis=2)  
  std_adv_s_all_mean = np.nanmean(std_adv_s_all, axis=2) 

#plt.figure()
#plt.pcolor(binx_s, biny_s,std_adv_s_all[:,:,0])
#plt.colorbar()  
#
#plt.figure()
#plt.pcolor(binx_s, biny_s, count_adv_s_all[:,:,0])
#plt.colorbar()  
#
#
#plt.figure()
#plt.pcolor(lon_adv, lat_adv,sss_adv_all[:,:, 1:][:,:,0])
#plt.colorbar()
  
  ''' Mask where no enough data '''

  count_threshold = 80 #%80 of the max data inside bins

  # big bin size
  std_adv_all_mean_mk = mask_no_enough_data(count_adv_all, 
                                 std_adv_all_mean, count_threshold)
  
  # small bin size
  std_adv_s_all_mean_mk = mask_no_enough_data(count_adv_s_all, 
                                 std_adv_s_all_mean, count_threshold)

  # mask also the SMAP data using the mask for Sadv to compare better

  std_smap_all_mean_mk = np.copy(std_smap_all_mean)
  std_smap_all_mean_mk[np.isnan(std_adv_all_mean_mk)] = np.nan
  
 
  
 
  ''' Plot MDT to see what contour defines the GS position ''' 
  
  lon_mdt_2d, lat_mdt_2d = np.meshgrid(lon_mdt, lat_mdt)
  x_mdt, y_mdt = bm(lon_mdt_2d, lat_mdt_2d)    
  
  plt.figure()
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, 16)
  
  cs= plt.contour(x_mdt, y_mdt, mdt.squeeze())
  plt.clabel(cs, cs.levels, inline=True, fontsize=10)
  qv0 = plt.quiver(x_mdt[::5,::5], y_mdt[::5,::5], 
             u_mdt.squeeze()[::5,::5], v_mdt.squeeze()[::5,::5],
             color='b', scale=9)
  plt.quiverkey(qv0, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200) 
  
  
  
  
  ''' Download bathymetry '''
  
  topo_global, lonbat, latbat = download_bathymetry()
  

  #limit data within domain
  
  inds_lat = np.logical_and(latbat >= 27.5, latbat <= latmax)
  inds_lon = np.logical_and(lonbat >= lonmin, lonbat <= lonmax)
  
  lonbat2d, latbat2d = np.meshgrid(lonbat[inds_lon], latbat[inds_lat])
  
  topo = topo_global[inds_lat,:][:, inds_lon]


  ''' Find coordinates of the desired 200 m isobath '''
  
  xbat, ybat = bm(lonbat2d, latbat2d)


  plt.figure()
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, 16)

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
  
  
  ''' Figure temporal average STD Smap 0.75, Sadv 0.75, Sadv 0.25 '''

  
  # fig parameters
  binx_2d, biny_2d = np.meshgrid(binx, biny)
  xx_big, yy_big = bm(binx_2d, biny_2d)
  
  binx_s_2d, biny_s_2d = np.meshgrid(binx_s, biny_s)
  xx_s, yy_s = bm(binx_s_2d, biny_s_2d)  
  
  lon_mdt_2d, lat_mdt_2d = np.meshgrid(lon_mdt, lat_mdt)
  x_mdt, y_mdt = bm(lon_mdt_2d, lat_mdt_2d)

  
  
  sminb, smaxb= 0.1, 0.5
  smins, smaxs= 0.1, 0.3
  
  fsize = 16
 
  # start figure
  fig = plt.figure(figsize=(14,5))
  
  wp  = 0.24
  hp  = 0.9 - 0.01
  ws  = wp + 0.08
  wcb = 0.015
  hcb = 0.85  - 0.01
  
  dd = 0.01
  
  ax0 = plt.axes([0.05-dd, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  ax1 = plt.axes([0.01 + ws-dd, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  axcb1 = plt.axes([0.01 + ws + wp + 0.02-dd, (1-hcb)/2, wcb, hcb])
  
  
  ax2 = plt.axes([0.01 + ws*2 + 0.03, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  axcb2 = plt.axes([0.01 + ws*2 + 0.03 + wp + 0.02, (1-hcb)/2, wcb, hcb])
  
  
  # plot 200 m isobath in all subplots
  lw = 1.2
  ax0.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  ax1.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  ax2.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  
  # plot MDT = 0.5 m (Gulf Stream mean position)

  ax0.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 
  ax2.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 
  
  #plt.colorbar(cax=cax)    
  
  # plot data big bins
  pc1=ax0.pcolor(xx_big, yy_big, std_smap_all_mean_mk, vmin=sminb, vmax=smaxb, 
                 cmap=plt.cm.magma_r)
  ax0.set_title('S$_{SMAP}$ ' + np.str(dx) + '$^\circ$x ' \
                + np.str(dx) + '$^\circ$', fontsize=fsize ) 
  
  ax1.pcolor(xx_big, yy_big, std_adv_all_mean_mk, vmin=sminb, vmax=smaxb, 
                 cmap=plt.cm.magma_r)
  ax1.set_title('S$_{adv}$ ' + np.str(dx) + '$^\circ$x ' \
                + np.str(dx) + '$^\circ$', fontsize=fsize )  
  

  cb1 = plt.colorbar(pc1, cax=axcb1)
  cb1.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb1.ax.tick_params(labelsize=fsize-1)   
  cb1.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
  # plot data small bins

  pc2 = ax2.pcolor(xx_s, yy_s, std_adv_s_all_mean_mk, vmin=smins, vmax=smaxs, 
                 cmap=plt.cm.magma_r)
  ax2.set_title('S$_{adv}$ ' + np.str(dxsmall) + '$^\circ$x ' \
                + np.str(dxsmall) + '$^\circ$', fontsize=fsize ) 
  #ax2.contour(x_mdt, y_mdt, mdt.squeeze(),10)  
#  ax2.quiver(x_mdt[::5,::5], y_mdt[::5,::5], 
#             u_mdt.squeeze()[::5,::5], v_mdt.squeeze()[::5,::5],
#             color='b')
  
  cb2 = plt.colorbar(pc2, cax=axcb2)
  cb2.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb2.ax.tick_params(labelsize=fsize-1)  
  cb2.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  fig.text(0.015, 0.97 - 0.01, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.295, 0.97 - 0.01, '(b)' , fontsize = fsize+1, zorder=200)
  fig.text(0.645, 0.97 - 0.01, '(c)' , fontsize = fsize+1, zorder=200)

#  pos1 = ax.get_position() # get the original position 
#  pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width / 2.0, pos1.height / 2.0] 
#  ax.set_position(pos2)
  
  #plt.tight_layout()#w_pad = -2)
  fig.savefig(dirfig + 'maps_SSS_STD_temporal_average_paper.png', dpi=200)
  
  



  ''' Extract data for each season '''
  
  winter_months = [2, 3, 4] #[12, 1, 2]
  summer_months = [6, 7, 8]
  
  # smap winter
  sss_smap_winter, time_smap_winter  = extract_data_for_season( \
                                        winter_months,  \
                                        sss_smap_all[:, :, 1:], \
                                        time_smap_all)

  # smap summer
  sss_smap_summer, time_smap_summer  = extract_data_for_season( \
                                        summer_months,  \
                                        sss_smap_all[:, :, 1:], \
                                        time_smap_all)  
  
  # adv winter
  sss_adv_winter, time_adv_winter  = extract_data_for_season( \
                                        winter_months,  \
                                        sss_adv_all[:, :, 1:], \
                                        time_adv_all)

  # adv summer
  sss_adv_summer, time_adv_summer  = extract_data_for_season( \
                                        summer_months,  \
                                        sss_adv_all[:, :, 1:], \
                                        time_adv_all)    
  
  
  ''' Infer std for each of these variables and seasons '''
  
  # Big bin,  Sadv
  std_adv_summer = infer_2D_std_for_each_date(binx, biny, time_adv_summer, 
                            lon_adv, lat_adv, sss_adv_summer[:,:,1:])
  
  std_adv_winter = infer_2D_std_for_each_date(binx, biny, time_adv_winter, 
                            lon_adv, lat_adv, sss_adv_winter[:,:,1:])

  # Big bin, smap
  std_smap_summer = infer_2D_std_for_each_date(binx, biny, time_smap_summer, 
                            lon_smapf, lat_smapf, sss_smap_summer[:,:,1:])

  std_smap_winter = infer_2D_std_for_each_date(binx, biny, time_smap_winter, 
                            lon_smapf, lat_smapf, sss_smap_winter[:,:,1:])
  
  
  # Small bin, adv
  std_adv_summer_s = infer_2D_std_for_each_date(binx_s, biny_s, time_adv_summer, 
                            lon_adv, lat_adv, sss_adv_summer[:,:,1:])
  
  std_adv_winter_s = infer_2D_std_for_each_date(binx_s, biny_s, time_adv_winter, 
                            lon_adv, lat_adv, sss_adv_winter[:,:,1:])

  count_adv_summer_s = infer_2D_count_for_each_date(binx_s, biny_s, time_adv_summer, 
                            lon_adv, lat_adv, sss_adv_summer[:,:,1:])
  
  count_adv_winter_s = infer_2D_count_for_each_date(binx_s, biny_s, time_adv_winter, 
                            lon_adv, lat_adv, sss_adv_winter[:,:,1:]) 
 
    
  ''' average over the seasons '''
  
  # Big bin, adv and smap
  std_adv_summer_mean    = np.nanmean(std_adv_summer, axis=2)
  std_adv_winter_mean    = np.nanmean(std_adv_winter, axis=2)

  std_smap_summer_mean   = np.nanmean(std_smap_summer, axis=2)
  std_smap_winter_mean   = np.nanmean(std_smap_winter, axis=2) 

  # Small bin, adv
  std_adv_summer_s_mean  = np.nanmean(std_adv_summer_s, axis=2)
  std_adv_winter_s_mean  = np.nanmean(std_adv_winter_s, axis=2)
  
 
    
  ''' Mask where no enough data '''
  

  count_threshold = 80 #%80 of the max data inside bins

  
  # small bin size
  std_adv_summer_s_mean_mk = mask_no_enough_data(count_adv_summer_s, 
                                 std_adv_summer_s_mean, count_threshold)

  std_adv_winter_s_mean_mk = mask_no_enough_data(count_adv_winter_s, 
                                 std_adv_winter_s_mean, count_threshold)
 


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
  
  
  # plt.figure()
  # plt.contour(lon_alt, lat_alt, adt_winter_mean)
  # plt.contour(lon_alt, lat_alt, adt_winter_mean, levels=[0.5],
  #            colors= 'r', linewidths=lw, linestyles='dashed')   
  # plt.title('winter')
  
  # plt.figure()
  # plt.contour(lon_alt, lat_alt, adt_summer_mean)  
  # plt.contour(lon_alt, lat_alt, adt_summer_mean, levels=[0.5],
  #            colors= 'r', linewidths=lw, linestyles='dashed')     
  # plt.title('summer')
  
  # plt.figure()
  # plt.contour(lon_alt, lat_alt, adt_winter_mean, levels=[0.5],
  #            colors= 'r', linewidths=lw, linestyles='dashed')  
  # plt.contour(lon_alt, lat_alt, adt_summer_mean, levels=[0.5],
  #            colors= 'b', linewidths=lw, linestyles='dashed')     
  
  

  ''' Figure winter vs. summer STD (only Sadv)'''

 # fig parameters
  binx_2d, biny_2d = np.meshgrid(binx, biny)
  xx_big, yy_big = bm(binx_2d, biny_2d)
  
  binx_s_2d, biny_s_2d = np.meshgrid(binx_s, biny_s)
  xx_s, yy_s = bm(binx_s_2d, biny_s_2d)  
  
  x_alt, y_alt = bm(lon_alt, lat_alt)
  
  sminb, smaxb= 0.1, 0.5
  smins, smaxs= 0.1, 0.3
  
  
  fsize = 10
  fig = plt.figure(figsize=(7,4))
   
  gs = gridspec.GridSpec(1, 3, width_ratios=[20, 20, 1])

  # Sadv 
  
  ax2 = plt.subplot(gs[0,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  ax3 = plt.subplot(gs[0, 1]) 
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  
  axcb2 = plt.subplot(gs[0, 2])
  
  # plot 200 m isobath in all subplots
  lw=0.9
  ax2.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  ax3.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
 
  # adt winter
  ax2.contour(x_alt, y_alt, adt_winter_mean, levels=[0.5],
             colors= 'r', linewidths=lw, linestyles='dashed') 

  # adt summer
  ax3.contour(x_alt, y_alt, adt_summer_mean, levels=[0.5],
             colors= 'b', linewidths=lw, linestyles='dashed')   
   
  # ax2.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
  #            colors= 'c', linewidths=lw, linestyles='dashed') 
  # ax3.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
  #            colors= 'c', linewidths=lw, linestyles='dashed')   
  
  pc2 = ax2.pcolor(xx_s, yy_s, std_adv_winter_s_mean_mk, vmin=smins, vmax=smaxs, 
                 cmap=plt.cm.magma_r)
  ax2.set_title('S$_{adv}$ ' + np.str(dxsmall) + '$^\circ$x ' \
                + np.str(dxsmall) + '$^\circ$  - Winter', fontsize=fsize ) 
  
  ax3.pcolor(xx_s, yy_s, std_adv_summer_s_mean_mk, vmin=smins, vmax=smaxs, 
                 cmap=plt.cm.magma_r)
  ax3.set_title('S$_{adv}$ ' + np.str(dxsmall) + '$^\circ$x ' \
                + np.str(dxsmall) + '$^\circ$  - Summer', fontsize=fsize ) 
  
  
  
  cb = plt.colorbar(pc2, cax=axcb2)
  cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb.ax.tick_params(labelsize=fsize-1)   
  cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
  
  fig.text(0.01, 0.95, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.43, 0.95, '(b)' , fontsize = fsize+1, zorder=200)

             
  plt.tight_layout()
  fig.savefig(dirfig + 'maps_SSS_STD_seasonal_average_paper.png', dpi=200)
  
    
  
  
