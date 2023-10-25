#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
from matplotlib              import dates as mdates 
from datetime                import datetime, timedelta
import pickle
from mpl_toolkits.basemap    import Basemap, shiftgrid
import sss_toolbox           as to
from scipy                   import stats
from matplotlib.ticker       import FormatStrFormatter
from netCDF4                 import Dataset
import matplotlib.gridspec   as gridspec


"""
Open CMEMS Mercator model data, read dates in which there are
weekly simulations.


Plot mean and gradient magnitude
Plot STD 0.75º x 0.75º and 0.24º x 0.24

Figure 10 in the paper.

written by Bàrbara Barceló-Llull on 18 August 2020 in Mallorca
"""


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

def read_cmems_model_surf(dir_model, file_model):
    
  nc      = netcdf.Dataset(dir_model+file_model, 'r') 
   
  depm         = nc.variables['depth'][:] # only upper layer = 0.494025 m
  time_m_orig  = nc.variables['time'][:]  # hours since 1950-01-01 00:00:00
  latm1d       = nc.variables['latitude'][:]  
  lonm1d       = nc.variables['longitude'][:]   
  um           = nc.variables['uo'][:]  #eastward_sea_water_velocity
  vm           = nc.variables['vo'][:]  #northward_sea_water_velocity
  ptm          = nc.variables['thetao'][:]  #sea_water_potential_temperature
  psm          = nc.variables['so'][:]  #sea_water_salinity,Practical Salinity Unit
   
  nc.close()  

  # lon and lat 2d array
  lonm, latm = np.meshgrid(lonm1d, latm1d)   


  # time in python format
  # first date
  d_ini_or    = mdates.date2num(datetime(1950, 1, 1, 0, 0, 0))  
  

  # Original time from hours to days 
  time_days = np.zeros(time_m_orig.shape)
  
  for ind, tt in enumerate(time_m_orig):
      time_days[ind] = mdates.hours(tt) # Return seconds as days.
    
    
  # Sum these days to d_ini_or
  timem = d_ini_or + time_days

  print(' ')
  print('ini time model...', mdates.num2date(timem.min()))
  print('fin time model...', mdates.num2date(timem.max()))
  print(' ')
  
  return timem, lonm, latm, depm, psm, ptm, um, vm

def infer_grad_magnitude(sal_model, dx_deg, dy_deg, lat_mean):
  ''' 
  
  Function to infer the magnitude of the 2D gradient.  
  
  dx_deg and dy_deg grid: spacing in degrees.
  lat_mean: mean latitude in degrees.
  
  sal_model: 2D variable to infer the magnitude of the gradient.
  '''
    
  # convertion factors: degrees to meters
  
  lat_mean_rad = np.nanmean(lat_mean) * np.pi/180 # lat mean in radians
  
  factor_diflon_deg2m, factor_diflat_deg2m = to.length_lon_lat_degs(lat_mean_rad)
       
     
  # grid spacing from degrees to km 
          
  dx_km = dx_deg * factor_diflon_deg2m /1000 #in km
  dy_km = dy_deg * factor_diflat_deg2m /1000 #in km
          
    
  # Magnitude of the gradient
  
  vgrad = np.gradient(sal_model, dy_km, dx_km) 
  # vgrad2 = np.gradient(sss_smap) = vgrad2[0]/dy_m vgrad2[1]/dx_m
  mag = np.sqrt(vgrad[0]**2 + vgrad[1]**2)   

  return mag   

def infer_2D_std_for_each_date(binx, biny, time_all, 
                            lon, lat, var):
  
  if  var.shape[2] == time_all.shape[0]:
      print('Data have the correct shape.')
      
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

def nan_where_std_0(std_adv_s):       
   # nan where std is 0 (because all nans in bin)
   std_adv_s_masked = np.ma.copy(std_adv_s)
   
   std_adv_s_masked =  np.ma.masked_where(std_adv_s_masked==0, std_adv_s_masked)
   
   std_adv_s_masked.data[std_adv_s_masked.mask==True] = np.nan
   
   return std_adv_s_masked.data

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

if __name__ == '__main__':
    
  plt.close('all')
  
  dirfig   = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_weekly_sims_BiggerDomain/model/' 
  
  
  ''' Open CMEMS Model data'''
  
  dir_model    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/CMEMS_global_reanalysis/' 
  # reanalysis from 2015-04-07 to 2018-12-25
  file_model_re  = 'global-reanalysis-phy-001-030-daily_1597743581842.nc'

  timem_re, lonm_re, latm_re, depm_re, psm_re, ptm_re, um_re, vm_re = \
                          read_cmems_model_surf(dir_model, file_model_re)
 
    
  # analysis from 2018-12-26 to 2020-01-15
  file_model_an  = 'global-analysis-forecast-phy-001-024_1597746265540.nc'

  timem_an, lonm_an, latm_an, depm_an, psm_an, ptm_an, um_an, vm_an = \
                          read_cmems_model_surf(dir_model, file_model_an)
 

  date_start_an = mdates.date2num(datetime(2018, 12, 26, 0, 0, 0))  

  # Parameters of the model to infer the gradient
  dx_deg = 0.08332825 # in deg
  dy_deg = 0.08333588 # in deg
   
   
  ''' File with the simulation dates '''
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_weekly_all/'
  file_dic = 'sim_weekly_relase_dates_corrected_gaps.pkl'

  f = open(dir_dic + file_dic, 'rb')
  dict_release_dates = pickle.load(f)
  f.close() 

  dates_release_all = dict_release_dates['date_release'][:-9] #- 0.5 #starting at 00:00 instead than at 12:00

  print('')    
  print('Date of simulation...', mdates.num2date(dates_release_all.min()))  
  print('Date of simulation...', mdates.num2date(dates_release_all.max()))       
  print('')
    
  '''
  Create the basemap
  '''

  # BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46  
  lat_mean = np.mean([latmin, latmax]) # to infer the gradient

  bm = Basemap(projection = 'merc',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')      
  
  ''' Open and save model data for each simulation end date '''

  sal_model_all = []
  u_model_all   = []
  v_model_all   = []   
  mag_model_all = []
  
  for date in dates_release_all:

   
    date_time_obj = mdates.num2date(date)
    
    print('')
    print('Date of simulation...', date_time_obj)
   

    ''' Open and save model reanalysis or analysis data '''


    if date < date_start_an:
        
          # use reanalysis data
       
          ind_mod = np.where(timem_re == date)
    
          print(' ')
          print('reanalysis for day... ', mdates.num2date(timem_re[ind_mod]))
          print(' ')
     
          sal_model = psm_re[ind_mod].squeeze()
          u_model   = um_re[ind_mod].squeeze()
          v_model   = vm_re[ind_mod].squeeze()
           

          
    elif  date >= date_start_an:
        
          # use analysis data 
        
          ind_mod = np.where(timem_an == date)
    
          print(' ')
          print('analysis for day... ', mdates.num2date(timem_an[ind_mod]))
          print(' ')
     
          sal_model = psm_an[ind_mod].squeeze()
          u_model   = um_an[ind_mod].squeeze()
          v_model   = vm_an[ind_mod].squeeze()
 
    # remove mask
    sal_model_nomask = np.copy(sal_model.data)
    sal_model_nomask[sal_model.mask==True]=np.nan

    u_model_nomask = np.copy(u_model.data)
    u_model_nomask[u_model.mask==True]=np.nan 
          
    v_model_nomask = np.copy(v_model.data)
    v_model_nomask[v_model.mask==True]=np.nan   

    # save into a big array
    sal_model_all.append(sal_model_nomask) 
    u_model_all.append(u_model_nomask)
    v_model_all.append(v_model_nomask)

                          
    ''' Infer and save gradient magnitude for each date '''

    mag_model = infer_grad_magnitude(sal_model_nomask, dx_deg, dy_deg, lat_mean)
   
    mag_model_all.append(mag_model)

          
  # arrays with all data from the model: sal_model_all, mag_model_all
          
          
  ''' Temporal average '''

  # convert list of 2d array into 3d array
  mag_model_all_3D = np.dstack(mag_model_all)
  sal_model_all_3D = np.dstack(sal_model_all)
  
  mag_model_all_mean  = np.nanmean(mag_model_all_3D, axis=2)
  sal_model_all_mean  = np.nanmean(sal_model_all_3D, axis=2)          
          
     
  ''' Bin STD 0.75º x 0.75º'''
  
  dx_big = 0.75 
  dy_big = dx_big
  
  binx_big = np.arange(lonmin, lonmax+dx_big, dx_big)
  biny_big = np.arange(latmin, latmax+dy_big, dy_big) #corners of the boxes
  
  std_model_big_all = infer_2D_std_for_each_date(binx_big, biny_big, dates_release_all, 
                            lonm_re, latm_re, sal_model_all_3D)

  std_model_big_all_mean  = np.nanmean(std_model_big_all, axis=2)


  ''' Bin STD 0.24º x 0.24º'''
  
  dxsmall = 0.04*6 
  dysmall = dxsmall
  
  binx_s = np.arange(lonmin, lonmax+dxsmall, dxsmall)
  biny_s = np.arange(latmin, latmax+dysmall, dysmall)
  
  
  std_model_small_all = infer_2D_std_for_each_date(binx_s, biny_s, dates_release_all, 
                            lonm_re, latm_re, sal_model_all_3D)
  
  std_model_small_all_mean  = np.nanmean(std_model_small_all, axis=2)
  
  
  
  ''' Open data for the example '''
  
  select_dates = ['20170212']
  
  
  for idt, date in enumerate(select_dates): 

      
    date_time_obj = datetime.strptime(date, '%Y%m%d')
    
    print('')
    print ('date... ', date_time_obj)
    
      
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))

    ''' Open model data for this date '''
    
    ind_mod = np.where(timem_re == date_fin+0.5)  
  
    print(' ')
    print('model data for day... ', mdates.num2date(timem_re[ind_mod]))
    print(' ')
     
    sal_model = psm_re[ind_mod].squeeze()
    u_model   = um_re[ind_mod].squeeze()
    v_model   = vm_re[ind_mod].squeeze()
    
    
    

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
  
  ''' Open MDT '''
  
  dir_mld  = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/MDT_CNES-CLS18_fromAviso/' 
  file_mld = 'dataset-mdt-cnes-cls18-global_1586181917100.nc'

  nc    = netcdf.Dataset(dir_mld + file_mld, 'r')
  
  lon_mdt  = nc.variables['longitude'][:]
  lat_mdt  = nc.variables['latitude'][:]
  u_mdt    = nc.variables['u'][:]
  v_mdt    = nc.variables['v'][:]
  mdt      = nc.variables['mdt'][:]  
  
  
  ''' Figure bin STD '''
  lon_mdt_2d, lat_mdt_2d = np.meshgrid(lon_mdt, lat_mdt)
  x_mdt, y_mdt = bm(lon_mdt_2d, lat_mdt_2d)   
  
  binx_2d, biny_2d = np.meshgrid(binx_big, biny_big)
  xx_big, yy_big = bm(binx_2d, biny_2d)
  
  binx_s_2d, biny_s_2d = np.meshgrid(binx_s, biny_s)
  xx_s, yy_s = bm(binx_s_2d, biny_s_2d)  
  
  fsize = 12
  sminb, smaxb= 0.1, 0.5
  smins, smaxs= 0.1, 0.3  

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

  axcb0 = plt.axes([wp-0.015-0.01 , (1-hcb)/2, wcb, hcb])

  ax1 = plt.axes([0.5-0.01, 0.05, wp, hp])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  axcb1 = plt.axes([0.5-0.01-0.01 + wp-0.015 , (1-hcb)/2, wcb, hcb])
  
  
  # plot data big bins
  pc0=ax0.pcolor(xx_big, yy_big, std_model_big_all_mean, vmin=sminb, vmax=smaxb, 
                 cmap=plt.cm.magma_r)
  ax0.set_title('S$_{model}$ ' + np.str(dx_big) + '$^\circ$x ' \
                + np.str(dy_big) + '$^\circ$', fontsize=fsize ) 
  
  # plot data small bins    
  pc1=ax1.pcolor(xx_s, yy_s, std_model_small_all_mean, vmin=smins, vmax=smaxs, 
                 cmap=plt.cm.magma_r)
  ax1.set_title('S$_{model}$ ' + np.str(dxsmall) + '$^\circ$x ' \
                + np.str(dysmall) + '$^\circ$', fontsize=fsize ) 
 
  # plot 200 m isobath in all subplots
  lw = 1.2
  ax0.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
  ax1.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')

  # plot MDT = 0.5 m (Gulf Stream mean position)

  ax0.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 
    
  # cb big bins
  cb0 = plt.colorbar(pc0, cax=axcb0)
  cb0.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb0.ax.tick_params(labelsize=fsize-1)  
  cb0.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  # cb small bins
  cb1 = plt.colorbar(pc1, cax=axcb1)
  cb1.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb1.ax.tick_params(labelsize=fsize-1)   
  cb1.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  fig.text(0.015, 0.97 - 0.01, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.5, 0.97 - 0.01, '(b)' , fontsize = fsize+1, zorder=200)

 

  fig.savefig(dirfig + 'S_model_bins_STD.png', dpi=200)
  
  
     
  ''' Figure temporal average Salinity model and  gradient for revision'''
 
  x_model, y_model = bm(lonm_re, latm_re)

  
  fsize = 12
  

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


  
  pc0=ax1.pcolor(x_model, y_model, mag_model_all_mean, vmin=0.003, vmax=0.03, 
                 cmap=plt.cm.magma_r)
  
  # qv0 = ax1.quiver(x_mdt[::5,::5], y_mdt[::5,::5], 
  #            u_mdt.squeeze()[::5,::5], v_mdt.squeeze()[::5,::5],
  #            color='b', scale=9)
  # ax1.quiverkey(qv0, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)

  # plot 200 m isobath in all subplots
  lw = 1.2
  ax1.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')

  # plot MDT = 0.5 m (Gulf Stream mean position)
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 

  
  ax1.set_title('S$_{model}$ grad. mag.', fontsize=fsize)     


 
  pc1 = ax0.pcolor(x_model, y_model, sal_model_all_mean, vmin=32, vmax=37,
             cmap=plt.cm.Spectral_r) # Spectral_r


  ax0.set_title('S$_{model}$', fontsize=fsize)  
  
  
  cb0 = plt.colorbar(pc0, cax=axcb1)
  cb0.ax.set_title('[PSU/km]', fontsize=fsize-1, pad=10)
  cb0.ax.tick_params(labelsize=fsize-1) 

  cb1 = plt.colorbar(pc1, cax=axcb0)
  cb1.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb1.ax.tick_params(labelsize=fsize-1)   
  

  fig.text(0.01, 0.96, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.5, 0.96, '(b)' , fontsize = fsize+1, zorder=200)
             
  #plt.tight_layout()
  fig.savefig(dirfig + 'S_model_temp_avg_grad_and_sss_revision.png', dpi=200)
  



  ''' Figure for PAPER'''
 
  #x_model, y_model = bm(lonm, latm)  
  x_model, y_model = bm(lonm_re, latm_re)
  
  
  fsize = 13
  fig = plt.figure(figsize=(12,8))
   
  gs = gridspec.GridSpec(2, 3)



  # ----- S EXAMPLE -----
  ax1 = fig.add_subplot(gs[0, 0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
  smin = 32#np.nanmin(sss)
  smax = 37#np.nanmax(sss)       
  pc1 = ax1.pcolor(x_model[0,:], y_model[:,0],sal_model, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
  # ax1.scatter(xtsg, ytsg, c=sal_tsg, 
  #           vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0,
  #           edgecolors='k',zorder=3)  
        
  ax1.set_title('S$_{model}$ on 12 Feb. 2017', fontsize=fsize)        

  cb1 = plt.colorbar(pc1)  
  cb1.ax.tick_params(labelsize=fsize-1) 
  cb1.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)


        
  # ----- S mean -----
    
  ax2 = fig.add_subplot(gs[0, 1])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
   
  pc2 = ax2.pcolor(x_model, y_model, sal_model_all_mean, vmin=32, vmax=37,
             cmap=plt.cm.Spectral_r) # Spectral_r

  ax2.set_title('S$_{model}$ mean', fontsize=fsize)     
   
 
  cb2 = plt.colorbar(pc2)
  cb2.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb2.ax.tick_params(labelsize=fsize-1)           

  # ----- S gradient -----
        
  ax3 = fig.add_subplot(gs[0, 2])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  pc3=ax3.pcolor(x_model, y_model, mag_model_all_mean, vmin=0.003, vmax=0.03, 
                 cmap=plt.cm.magma_r)
  
  # plot 200 m isobath in all subplots
  lw = 1.2
  ax3.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')

  # plot MDT = 0.5 m (Gulf Stream mean position)
  ax3.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 

  
  ax3.set_title('S$_{model}$ grad. mag. mean', fontsize=fsize)    

  
  cb3 = plt.colorbar(pc3)
  cb3.ax.set_title('[PSU/km]', fontsize=fsize-1, pad=10)
  cb3.ax.tick_params(labelsize=fsize-1)   

  # ----- S 0.75ºx0.75º -----
  ax4 = fig.add_subplot(gs[1, 1])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
   
         
  # plot data big bins
  pc4=ax4.pcolor(xx_big, yy_big, std_model_big_all_mean, vmin=sminb, vmax=smaxb, 
                 cmap=plt.cm.magma_r)
  ax4.set_title('S$_{model}$ ' + np.str(dx_big) + '$^\circ$x ' \
                + np.str(dy_big) + '$^\circ$', fontsize=fsize ) 

  # plot 200 m isobath in all subplots
  lw = 1.2
  ax4.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')
 
  # plot MDT = 0.5 m (Gulf Stream mean position)
  ax4.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed')     
  
  

  cb4 = plt.colorbar(pc4)
  cb4.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb4.ax.tick_params(labelsize=fsize-1)  
  cb4.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

  
  # ----- S 0.24ºx0.24º -----
  ax5 = fig.add_subplot(gs[1, 2])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  
  # plot data small bins    
  pc5=ax5.pcolor(xx_s, yy_s, std_model_small_all_mean, vmin=smins, vmax=smaxs, 
                 cmap=plt.cm.magma_r)
  ax5.set_title('S$_{model}$ ' + np.str(dxsmall) + '$^\circ$x ' \
                + np.str(dysmall) + '$^\circ$', fontsize=fsize ) 
 

  ax5.plot(isob[:,0], isob[:,1], 'c',linewidth=lw,
              linestyle='-')


  ax5.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed') 



  cb5 = plt.colorbar(pc5)
  cb5.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
  cb5.ax.tick_params(labelsize=fsize-1)  
  cb5.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))  
  

  fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
  fig.text(0.33, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
  fig.text(0.66, 0.97, '(c)' , fontsize = fsize+1, zorder=200)
  fig.text(0.33, 0.48, '(d)' , fontsize = fsize+1, zorder=200)
  fig.text(0.66, 0.48, '(e)' , fontsize = fsize+1, zorder=200)  

             
  plt.tight_layout()
  fig.savefig(dirfig + 'S_model_all_for_paper.png', dpi=200)
        