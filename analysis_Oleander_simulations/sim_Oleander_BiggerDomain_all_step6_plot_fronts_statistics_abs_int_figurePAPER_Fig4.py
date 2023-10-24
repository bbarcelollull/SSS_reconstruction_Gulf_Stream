#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
from netCDF4                 import Dataset
import pickle
from matplotlib              import dates as mdates
from datetime                import datetime, timedelta
from scipy                   import interpolate
import sim_toolbox           as sto
import geopy.distance
import matplotlib.gridspec   as gridspec
from scipy                   import stats
from matplotlib.ticker       import PercentFormatter
from mpl_toolkits.basemap    import Basemap, shiftgrid
from scipy                   import interpolate
"""
Open file with fronts detected on step 4 and make figures for paper:
    number of fronts and absolute intensity of fronts
    
    FIGURE FOR PAPER (Figure 4)

written by Bàrbara Barceló-Llull on 14-April-2020 at Mallorca
modified on 4 May 2020 to plot fronts on interpolated data.

"""


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

def plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize):
  
      
      # decor
      bm.drawcoastlines()
      bm.drawrivers(zorder=6)
      bm.fillcontinents(color='0.8', lake_color='0.8', zorder=5)
      
      
      parallels = np.arange(latmin,latmax,5)
      bm.drawparallels(parallels,labels=[1, 0, 0, 0],
                             fontsize=fsize-1, linewidth=0.1, zorder=8)
      meridians = np.arange(lonmin-3, lonmax, 5)
      #meridians = [-72, -70, -68]
      bm.drawmeridians(meridians,labels=[0, 0, 0, 1],
                             fontsize=fsize-1, linewidth=0.1, zorder=9)


def detect_chunks(a):
    
    return [a[s] for s in np.ma.clump_unmasked(np.ma.masked_invalid(a))]

def intensity_fronts(stsg_fronts, dist, stsg_fronts_list):
 
    # distance of the fronts    
    cond_tsg                 = np.isnan(stsg_fronts)
    dist_front_tsg           = np.copy(dist[1:]) 
    dist_front_tsg[cond_tsg] = np.nan
    
    list_dist_front = detect_chunks(dist_front_tsg)
       
    int_abs_all   =[]
    int_grad_all  = []
    mean_dist_all = []
    for il, s_front in enumerate(stsg_fronts_list):
        
        int_abs     = np.max(s_front) - np.min(s_front)
        int_abs_all = np.append(int_abs_all, int_abs)
        
        dist_front   = list_dist_front[il]
        
        distance_over_front = np.max(dist_front) - np.min(dist_front)
        
        int_grad     = int_abs/dist_front
        int_grad_all = np.append(int_grad_all, int_grad)
    
        mean_dist_all = np.append(mean_dist_all, np.mean(dist_front))
        
    return int_abs_all, int_grad_all, mean_dist_all

def front_in_grid(var): 
    varcopy = np.copy(var)
    
    varcopy[np.isnan(var)] = 0
    varcopy[~np.isnan(var)] = 1
    
    return varcopy



def map_stats_new_intensity(lon_ori, lat_ori, sal_fronts):
    #lon_ori, lat_ori, sal_1wayo_nfronts
    ''' Number of fronts per transect and intensity of fronts
    inferred as the difference between maximum and minimum salinity
    divided by length of fronts '''
    
#    plt.figure()
#    plt.scatter(lon_ori,lat_ori, c=sal_fronts)
#    plt.colorbar()
    
    if len(sal_fronts[~np.isnan(sal_fronts)]) >0:
        
      ''' 
      map of the number of fronts in each bin and mean intensity 
      also return the mean lon, mean lat and absolute intensity of each front
      '''
    
      lon_fronts_tsg_p = np.copy(lon_ori)
      lat_fronts_tsg_p = np.copy(lat_ori)
    
      lon_fronts_tsg_p[np.isnan(sal_fronts)] = np.nan
      lat_fronts_tsg_p[np.isnan(sal_fronts)] = np.nan
    
      lon_fronts_tsg_p_list = detect_chunks(lon_fronts_tsg_p)
      lat_fronts_tsg_p_list = detect_chunks(lat_fronts_tsg_p)
      sal_tsgo_pfronts_list = detect_chunks(sal_fronts)
    
      lon_front_mean = []
      lat_front_mean = []
      sal_front_int  = []
      for ifront in np.arange(len(lon_fronts_tsg_p_list)):
        
        lon_front = lon_fronts_tsg_p_list[ifront]
        lat_front = lat_fronts_tsg_p_list[ifront]
        sal_front = sal_tsgo_pfronts_list[ifront]
        
        intensity_front_abs = np.max(sal_front) - np.min(sal_front)
        #print('abs int...', intensity_front_abs )
        # distance between (latmax, lonmin) and (latmin, lonmax)

        coords_1 = (lon_front.min(), lat_front.max())
        coords_2 = (lon_front.max(), lat_front.min())
        
        len_front  = geopy.distance.geodesic(coords_1, coords_2).km 
        #print('length...', len_front )
        
        #plt.figure()
        #plt.subplot(211)
        #plt.scatter(lon_front, lat_front, c=sal_front)
        #plt.colorbar()
        #plt.subplot(212)
        #plt.plot(lon_front, sal_front)
        
        # absolute front intensity/front length
        intensity_front_rel = intensity_front_abs/len_front #PSU/km
        #print('rel int...', intensity_front_rel )
        #print('')
        
        
        lon_front_mean = np.append(lon_front_mean, np.mean(lon_front))
        lat_front_mean = np.append(lat_front_mean, np.mean(lat_front))
        sal_front_int  = np.append(sal_front_int, intensity_front_rel)
        
      # 2D statistics for this transect  
      retnum = stats.binned_statistic_2d(lon_front_mean, lat_front_mean, 
                                          None, 'count', bins=[lon_grid,lat_grid])
      map_num = retnum.statistic.T


      retint = stats.binned_statistic_2d(lon_front_mean, lat_front_mean, 
                                          sal_front_int, 'median', bins=[lon_grid,lat_grid])
      map_int = retint.statistic.T
      
    else:
       lon_front_mean = [] 
       lat_front_mean = []
       sal_front_int = []
       map_num = np.zeros(lon_grid_2d[:-1,:-1].shape) 
       map_int = np.ones(lon_grid_2d[:-1,:-1].shape) * np.nan
       
        
    return lon_front_mean, lat_front_mean, sal_front_int, map_num, map_int


def map_stats(lon_ori, lat_ori, sal_fronts):
    
    if len(sal_fronts[~np.isnan(sal_fronts)]) >0:
        
      ''' 
      map of the number of fronts in each bin and mean intensity 
      also return the mean lon, mean lat and absolute intensity of each front
      '''
    
      lon_fronts_tsg_p = np.copy(lon_ori)
      lat_fronts_tsg_p = np.copy(lat_ori)
    
      lon_fronts_tsg_p[np.isnan(sal_fronts)] = np.nan
      lat_fronts_tsg_p[np.isnan(sal_fronts)] = np.nan
    
      lon_fronts_tsg_p_list = detect_chunks(lon_fronts_tsg_p)
      lat_fronts_tsg_p_list = detect_chunks(lat_fronts_tsg_p)
      sal_tsgo_pfronts_list = detect_chunks(sal_fronts)
    
      lon_front_mean = []
      lat_front_mean = []
      sal_front_int  = []
      for ifront in np.arange(len(lon_fronts_tsg_p_list)):
        
        lon_front = lon_fronts_tsg_p_list[ifront]
        lat_front = lat_fronts_tsg_p_list[ifront]
        sal_front = sal_tsgo_pfronts_list[ifront]
        
        intensity_front = np.max(sal_front) - np.min(sal_front)
        
        lon_front_mean = np.append(lon_front_mean, np.mean(lon_front))
        lat_front_mean = np.append(lat_front_mean, np.mean(lat_front))
        sal_front_int  = np.append(sal_front_int, intensity_front)
        
      # 2D statistics for this transect  
      retnum = stats.binned_statistic_2d(lon_front_mean, lat_front_mean, 
                                          lon_front_mean, 'count', bins=[lon_grid,lat_grid])
      map_num = retnum.statistic.T


      retint = stats.binned_statistic_2d(lon_front_mean, lat_front_mean, 
                                          sal_front_int, 'median', bins=[lon_grid,lat_grid])
      map_int = retint.statistic.T
      
    else:
       lon_front_mean = [] 
       lat_front_mean = []
       sal_front_int = []
       map_num = np.zeros(lon_grid_2d[:-1,:-1].shape) 
       map_int = np.ones(lon_grid_2d[:-1,:-1].shape) * np.nan
       
        
    return lon_front_mean, lat_front_mean, sal_front_int, map_num, map_int

  
def map_all(lon_fronts, lat_fronts, sal_fronts, lon_grid, lat_grid):  
    
  if len(sal_fronts) >0:
    
    retnum = stats.binned_statistic_2d(lon_fronts, lat_fronts, 
                                          lon_fronts, 'count', bins=[lon_grid,lat_grid])
    map_num = retnum.statistic.T


    #retint = stats.binned_statistic_2d(lon_fronts, lat_fronts, 
    #                                      sal_fronts, 'median', bins=[lon_grid,lat_grid])
    #map_median_int = retint.statistic.T    
  
    # Better to plot mean instead of median
    retint = stats.binned_statistic_2d(lon_fronts, lat_fronts, 
                                          sal_fronts, 'mean', bins=[lon_grid,lat_grid])
    map_mean_int = retint.statistic.T      
    
  else: #for empty list
      
    print('')  
    print('list is empty')
    map_num = np.zeros(lon_grid_2d[:-1,:-1].shape) 
    map_mean_int = np.ones(lon_grid_2d[:-1,:-1].shape) * np.nan 
       
  return map_num, map_mean_int

def merge_all_fronts(lon_p, lon_n, lat_p, lat_n, sal_p, sal_n):  
    
  lon_total = np.append(lon_p, lon_n)
  lat_total = np.append(lat_p, lat_n)
  sal_total = np.append(sal_p, sal_n)
  
  return lon_total, lat_total, sal_total

def mask_not_enough_trans_in_bins(map_num_tsg_p, cond_trans):
    
  map_num_tsg_p_mask = np.copy(map_num_tsg_p)
  map_num_tsg_p_mask[cond_trans==True] = np.nan
  
  return map_num_tsg_p_mask

if __name__ == '__main__':
    
  plt.close('all')

  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/'+ \
                'fig_Oleander_BiggerRegion/fronts/'
                
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all_Bigger_Domain/post-processing/'
  
  plot_fig = False #True
  
  
  ''' File with the detected fronts '''  
  

  file_dic = 'Oleander_BD_front_detection_5km.pkl'


  f = open(dir_dic + file_dic, 'rb')
  dict_fronts = pickle.load(f)
  f.close() 

  dates_strings = list(dict_fronts.keys())
  
  ''' bin grid '''
  
  # DOMAIN FOR BINS
  lonmin, lonmax = -75, -63 #-82, -63 
  latmin, latmax =  32, 41 #25, 46
  
  lon_grid = np.arange(lonmin, lonmax+0.01, 0.25)
  lat_grid = np.arange(latmin, latmax+0.01, 0.25)
    
  lon_grid_2d, lat_grid_2d = np.meshgrid(lon_grid, lat_grid)
  
  map_number_p = np.zeros(lon_grid_2d.shape)
  map_number_n = np.zeros(lon_grid_2d.shape)
  map_intens = np.zeros(lon_grid_2d.shape)

  '''
  Create the basemap
  '''
  bm = Basemap(projection = 'merc',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')      
  
  ''' Loop for each date '''
  
  dates_py = [] 

  lon_front_mean_tsg_p_all = []
  lat_front_mean_tsg_p_all = []
  sal_front_int_tsg_p_all  = []

  lon_front_mean_1w_p_all = []
  lat_front_mean_1w_p_all = []
  sal_front_int_1w_p_all  = []

  lon_front_mean_sat_p_all = []
  lat_front_mean_sat_p_all = []
  sal_front_int_sat_p_all  = []
    
  lon_front_mean_tsg_n_all = []
  lat_front_mean_tsg_n_all = []
  sal_front_int_tsg_n_all  = []
        
  lon_front_mean_1w_n_all = []
  lat_front_mean_1w_n_all = []
  sal_front_int_1w_n_all  = []

            
  lon_front_mean_sat_n_all = []
  lat_front_mean_sat_n_all = []
  sal_front_int_sat_n_all  = []
        
  map_tran_in_bins_all = np.ones((len(dates_strings), 
                                  lat_grid.shape[0]-1, 
                                  lon_grid.shape[0]-1))*np.nan
  
  for idt, date in enumerate(dates_strings):

    date_time_obj = datetime.strptime(date, '%Y%m%d')
    
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))  
    
    print('')    
    print(date)
    
    ''' Open and save data for this date '''

    # save dates
    dates_py = np.append(dates_py, date_fin)

    # open dictionary for this date    
    dict_date = dict_fronts[date]

 
#    dist_ori = dict_date['dist_ori']
#    lon_ori  = dict_date['lon_ori']
#    lat_ori  = dict_date['lat_ori']
    
    dist_int = dict_date['dist_int']
    lon_int  = dict_date['lon_int']
    lat_int  = dict_date['lat_int']

    sal_tsgi          = dict_date['sal_tsgi']
    sal_tsgi_pfronts  = dict_date['sal_tsgi_pfronts']
    sal_tsgi_nfronts  = dict_date['sal_tsgi_nfronts']

    sal_1wayi         = dict_date['sal_1wayi']
    sal_1wayi_pfronts = dict_date['sal_1wayi_pfronts']
    sal_1wayi_nfronts = dict_date['sal_1wayi_nfronts']
    
    sal_sati          = dict_date['sal_sati']
    sal_sati_pfronts  = dict_date['sal_sati_pfronts']
    sal_sati_nfronts  = dict_date['sal_sati_nfronts']
    
    ''' Interpolate lon_ori, lat_ori onto dist_int '''
    # not necessary with new fronts file, it's provisional
    
#    flon = interpolate.interp1d(dist_ori, lon_ori)
#    lon_int = flon(dist_int)
#    
#    flat = interpolate.interp1d(dist_ori, lat_ori)
#    lat_int = flat(dist_int)
#    
#    plt.figure()
#    plt.scatter(lon_ori, lat_ori, s=5, marker = 'o', color='k')
#    plt.scatter(lon_int, lat_int, s=5, marker = 'x', color='r')
#    
#    
    ''' 
    In which bin do we have data? 1 = data, 0 = no data within each bin,
    then sum number of transects inside bins.
    '''
    
    # number of data points inside each bin
    
    nonan = ~np.isnan(lon_int)
    stats_data_in_bins = stats.binned_statistic_2d(lon_int[nonan], lat_int[nonan], 
                                          lon_int[nonan], 'count', bins=[lon_grid,lat_grid])
    map_data_in_bins = stats_data_in_bins.statistic.T   
    
    # 0 = no data, 1 = data
    
    map_tran_in_bins = np.copy(map_data_in_bins)
    map_tran_in_bins[map_data_in_bins>0] = 1
    
    map_tran_in_bins_all[idt] = map_tran_in_bins
    
#    plt.figure()
#    plt.plot(lon_ori, lat_ori)
#    plt.pcolor(lon_grid, lat_grid, map_tran_in_bins)
#    plt.colorbar()
   
    
    ''' 
    mean position and median absolute intensity of each front 
    in each transect (date)
    sal_front_int_tsg_p is the intensity of each front for this date!!
    '''
    
    # TSG positive
    lon_front_mean_tsg_p, lat_front_mean_tsg_p, sal_front_int_tsg_p, \
    map_num_tsg_p_tran, map_int_tsg_p_tran = \
                              map_stats(lon_int, lat_int, sal_tsgi_pfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_tsg_p_all = np.append(lon_front_mean_tsg_p_all, lon_front_mean_tsg_p )
    lat_front_mean_tsg_p_all = np.append(lat_front_mean_tsg_p_all, lat_front_mean_tsg_p )
    sal_front_int_tsg_p_all  = np.append(sal_front_int_tsg_p_all, sal_front_int_tsg_p )
    
    # 1 way positive
    
    lon_front_mean_1w_p, lat_front_mean_1w_p, sal_front_int_1w_p, \
    map_num_1w_p_tran, map_int_1w_p_tran = \
                              map_stats(lon_int, lat_int, sal_1wayi_pfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_1w_p_all = np.append(lon_front_mean_1w_p_all, lon_front_mean_1w_p )
    lat_front_mean_1w_p_all = np.append(lat_front_mean_1w_p_all, lat_front_mean_1w_p )
    sal_front_int_1w_p_all  = np.append(sal_front_int_1w_p_all, sal_front_int_1w_p )

#    plt.figure(figsize=(16,8))
#    plt.plot(sal_front_int_1w_p, 'ob')
#    plt.title(date)
#    
#    if len(sal_front_int_1w_p)>0:
#      if np.nanmin(sal_front_int_1w_p) <0.02:
#        aaa
      
    # sat positive
    lon_front_mean_sat_p, lat_front_mean_sat_p, sal_front_int_sat_p, \
    map_num_sat_p_tran, map_int_sat_p_tran = \
                              map_stats(lon_int, lat_int, sal_sati_pfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_sat_p_all = np.append(lon_front_mean_sat_p_all, lon_front_mean_sat_p )
    lat_front_mean_sat_p_all = np.append(lat_front_mean_sat_p_all, lat_front_mean_sat_p )
    sal_front_int_sat_p_all  = np.append(sal_front_int_sat_p_all, sal_front_int_sat_p )
    
    
    
    # TSG negative
    
    lon_front_mean_tsg_n, lat_front_mean_tsg_n, sal_front_int_tsg_n, \
    map_num_tsg_n_tran, map_int_tsg_n_tran = \
                              map_stats(lon_int, lat_int, sal_tsgi_nfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_tsg_n_all = np.append(lon_front_mean_tsg_n_all, lon_front_mean_tsg_n )
    lat_front_mean_tsg_n_all = np.append(lat_front_mean_tsg_n_all, lat_front_mean_tsg_n )
    sal_front_int_tsg_n_all  = np.append(sal_front_int_tsg_n_all, sal_front_int_tsg_n )
    
    # 1 way negative
    
    lon_front_mean_1w_n, lat_front_mean_1w_n, sal_front_int_1w_n, \
    map_num_1w_n_tran, map_int_1w_n_tran = \
                              map_stats(lon_int, lat_int, sal_1wayi_nfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_1w_n_all = np.append(lon_front_mean_1w_n_all, lon_front_mean_1w_n )
    lat_front_mean_1w_n_all = np.append(lat_front_mean_1w_n_all, lat_front_mean_1w_n )
    sal_front_int_1w_n_all  = np.append(sal_front_int_1w_n_all, sal_front_int_1w_n )

#    plt.figure(figsize=(16,8))
#    plt.plot(sal_front_int_1w_n, 'ob')
#    plt.title(date)
#    
#    if np.nanmin(sal_front_int_1w_n) <0.02:
#      aaa
      
    # sat negative
    lon_front_mean_sat_n, lat_front_mean_sat_n, sal_front_int_sat_n, \
    map_num_sat_n_tran, map_int_sat_n_tran = \
                              map_stats(lon_int, lat_int, sal_sati_nfronts)
    
    # list of mean lon, mean lat and intensity of each front detected
    lon_front_mean_sat_n_all = np.append(lon_front_mean_sat_n_all, lon_front_mean_sat_n )
    lat_front_mean_sat_n_all = np.append(lat_front_mean_sat_n_all, lat_front_mean_sat_n )
    sal_front_int_sat_n_all  = np.append(sal_front_int_sat_n_all, sal_front_int_sat_n )


    lon_front_mean_sat_n, lat_front_mean_sat_n, sal_front_int_sat_n, \
    map_num_sat_n_tran, map_int_sat_n_tran = \
                              map_stats(lon_int, lat_int, sal_sati_nfronts)

    

    
    ''' Plot figures for each date'''

#    if plot_fig == True:
#
#      ''' Figure fronts on the smoothed data '''
#      
#      cmin = np.nanmin(sal_tsgo)-0.5 # 0.5
#      cmax = np.nanmax(sal_tsgo)+1.2
#    
#      plt.figure(figsize=(10,10))
#      plt.subplot(311)
#      plt.plot(lon_int, sal_tsgi, '-', color='0.7', linewidth =2,  
#             label='TSG')
#    
#      plt.plot(lon_int, sal_tsgi_pfronts, '-', color='r', linewidth =4, 
#             label='TSG')    
#      plt.plot(lon_int, sal_tsgi_nfronts, '-', color='b', linewidth =4, 
#             label='TSG')        
#      plt.ylabel('TSG')
#      #plt.xlim(0,500) 
#      plt.xlim(lon_int.min()-0.1,lon_int.max()+0.1) 
#      #plt.xlim(dist_ori.min(), dist_ori.max()) 
#      plt.ylim(cmin, cmax)
#      #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
#      plt.title('Front detection - ' + mdates.num2date(date_fin).strftime("%Y%m%d"))
#
#      plt.subplot(312)
#      plt.plot(lon_int, sal_1wayi,  '-', color='0.7', linewidth =2,  label='1 way')
#      plt.plot(lon_int, sal_1wayi_pfronts, '-', color='r', linewidth =4,label='1 way')
#      plt.plot(lon_int, sal_1wayi_nfronts, '-', color='b', linewidth =4,label='1 way')
#    
#      plt.ylabel('1-way')
#      #plt.xlim(0,500) 
#      plt.xlim(lon_int.min()-0.1,lon_int.max()+0.1) 
#      #plt.xlim(dist_ori.min(), dist_ori.max()) 
#      plt.ylim(cmin, cmax)
#      #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
#
#    
#      plt.subplot(313)
#      plt.plot(lon_int, sal_sati,  '-', color='0.7', linewidth =2,  label='smap')
#      plt.plot(lon_int, sal_sati_pfronts,  '-', '-', color='r', linewidth =4, label='smap')
#      plt.plot(lon_int, sal_sati_nfronts,  '-', '-', color='b', linewidth =4, label='smap')
#
#      plt.xlabel('Longitude')
#      plt.ylabel('SMAP')
#      #plt.xlim(0,500) 
#      plt.xlim(lon_int.min()-0.1,lon_int.max()+0.1) 
#      #plt.xlim(dist_ori.min(), dist_ori.max()) 
#      plt.ylim(cmin, cmax)
#      #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
#    
#      plt.savefig(dirfig + 'Oleander_fronts_original_sign'  + 
#                 '_' + mdates.num2date(date_fin).strftime("%Y%m%d") + '_5km.png', dpi=200)
#       
#      plt.close('all')    


  ''' Statistics of the total number of fronts (for all trans) 
  and mean intensity of all the fronts detected within each bin '''
  
  # tsg positive
  map_num_tsg_p, map_int_tsg_p = map_all(lon_front_mean_tsg_p_all, 
                                         lat_front_mean_tsg_p_all, \
                           sal_front_int_tsg_p_all, lon_grid, lat_grid)
     
  # 1 way positive 
  map_num_1w_p, map_int_1w_p = map_all(lon_front_mean_1w_p_all, 
                                       lat_front_mean_1w_p_all, \
                           sal_front_int_1w_p_all, lon_grid, lat_grid)


  # sat positive 
  map_num_sat_p, map_int_sat_p = map_all(lon_front_mean_sat_p_all, 
                                         lat_front_mean_sat_p_all, \
                           sal_front_int_sat_p_all, lon_grid, lat_grid)
      
  # tsg negative
  map_num_tsg_n, map_int_tsg_n = map_all(lon_front_mean_tsg_n_all, lat_front_mean_tsg_n_all, \
                           sal_front_int_tsg_n_all, lon_grid, lat_grid)
     
  # 1 way negative 
  map_num_1w_n, map_int_1w_n = map_all(lon_front_mean_1w_n_all, lat_front_mean_1w_n_all, \
                           sal_front_int_1w_n_all, lon_grid, lat_grid)

 
  # sat negative 
  map_num_sat_n, map_int_sat_n = map_all(lon_front_mean_sat_n_all, lat_front_mean_sat_n_all, \
                           sal_front_int_sat_n_all, lon_grid, lat_grid)
  
  
  # tsg total, merge positive and negative fronts, and statistics
  lon_total_tsg, lat_total_tsg, sal_total_tsg = \
          merge_all_fronts(lon_front_mean_tsg_p_all, lon_front_mean_tsg_n_all, 
                           lat_front_mean_tsg_p_all, lat_front_mean_tsg_n_all, 
                           sal_front_int_tsg_p_all,  sal_front_int_tsg_n_all)

  map_num_tsg_t, map_int_tsg_t = map_all(lon_total_tsg, lat_total_tsg, \
                           sal_total_tsg, lon_grid, lat_grid)
  
  
  # 1 way total
  lon_total_1w, lat_total_1w, sal_total_1w = \
          merge_all_fronts(lon_front_mean_1w_p_all, lon_front_mean_1w_n_all, 
                           lat_front_mean_1w_p_all, lat_front_mean_1w_n_all, 
                           sal_front_int_1w_p_all,  sal_front_int_1w_n_all)          

  map_num_1w_t, map_int_1w_t = map_all(lon_total_1w, lat_total_1w, \
                           sal_total_1w, lon_grid, lat_grid)

#  plt.figure(figsize=(16,8))
#  plt.plot(sal_front_int_1w_p_all, 'ob')
#  plt.plot(sal_front_int_1w_n_all, 'or')
  
  # sat total
  lon_total_sat, lat_total_sat, sal_total_sat = \
          merge_all_fronts(lon_front_mean_sat_p_all, lon_front_mean_sat_n_all, 
                           lat_front_mean_sat_p_all, lat_front_mean_sat_n_all, 
                           sal_front_int_sat_p_all,  sal_front_int_sat_n_all) 
          
  map_num_sat_t, map_int_sat_t = map_all(lon_total_sat, lat_total_sat,  \
                           sal_total_sat, lon_grid, lat_grid)
  

  ''' Download bathymetry '''
  
  topo_global, lonbat, latbat = download_bathymetry()
  

  #limit data within domain
  
  inds_lat = np.logical_and(latbat >= latmin, latbat <= latmax)
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
  
  
  
  
  ''' Number of transects per bin '''
  
  num_trans_in_bins = np.sum(map_tran_in_bins_all, axis=0)
  
  
  ''' Don't plot bins with at least 5% of transects crossing it '''
  th = 15/100 # At least 80 % of max # of transects 
  
  max_trans_in_bins = num_trans_in_bins.max()
  
  num_trans_rel = num_trans_in_bins/max_trans_in_bins
  
  cond_trans = num_trans_rel<th

  map_num_tsg_t_mask = mask_not_enough_trans_in_bins(map_num_tsg_t, cond_trans)
  map_int_tsg_t_mask = mask_not_enough_trans_in_bins(map_int_tsg_t, cond_trans)

  map_num_1w_t_mask = mask_not_enough_trans_in_bins(map_num_1w_t, cond_trans)
  map_int_1w_t_mask = mask_not_enough_trans_in_bins(map_int_1w_t, cond_trans)

  map_num_sat_t_mask = mask_not_enough_trans_in_bins(map_num_sat_t, cond_trans)
  map_int_sat_t_mask = mask_not_enough_trans_in_bins(map_int_sat_t, cond_trans)

  ''' plot intensity only if number of fronts/number of transects > 0.15 '''
  
  thnum = 0.15
  
  map_int_tsg_t_mkmk = np.copy(map_int_tsg_t_mask)
  map_int_1w_t_mkmk  = np.copy(map_int_1w_t_mask)
  map_int_sat_t_mkmk = np.copy(map_int_sat_t_mask)
  
  map_int_tsg_t_mkmk[map_num_tsg_t_mask/num_trans_in_bins<thnum] = np.nan
  map_int_1w_t_mkmk[map_num_1w_t_mask/num_trans_in_bins<thnum]   = np.nan
  map_int_sat_t_mkmk[map_num_sat_t_mask/num_trans_in_bins<thnum] = np.nan

  
  ''' Figure for paper: number and intensity of fronts - no mask '''


  x_grid_2d, y_grid_2d= bm(lon_grid_2d, lat_grid_2d)
  xbat, ybat = bm(lonbat2d, latbat2d)


  lon_mdt_2d, lat_mdt_2d = np.meshgrid(lon_mdt, lat_mdt)
  x_mdt, y_mdt = bm(lon_mdt_2d, lat_mdt_2d)    
  
  
  fsize = 14


  
  # --- start figure ---
  
  fig = plt.figure(figsize=(10,10))
  
  gs = gridspec.GridSpec(4, 3,  height_ratios=[1,20, 1,20])


  ''' number total - FOR PAPER '''
  
  nmin, nmax = 0,1
  
  ax1 = plt.subplot(gs[1,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  
  #isobath
  lw = 0.8
  cs1 = ax1.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
  #ax1.clabel(cs1, cs1.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)

  
  # GS mean position
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1) 
  
  
  
  pc=plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                    map_num_tsg_t_mask/num_trans_in_bins, #map_num_tsg_t.sum(),
                    vmin=nmin, vmax=nmax, cmap=plt.cm.magma_r) #cmap=plt.cm.magma_r)  
  plt.title('S$_{TSG}$', fontsize = fsize)

#  plt.axis('image')
#  plt.ylabel('Latitude')
  
  ax2 = plt.subplot(gs[1,1])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
#  cs2 = ax2.contour(xbat, ybat, topo, levels=[-1000],
#              colors='k',linewidths=0.2,
#              linestyles='solid')
  
  # isobath
  cs2 = ax2.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
 
 # ax2.clabel(cs2, cs2.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)

  # GS mean position
  ax2.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1) 
  
  
  plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                 map_num_1w_t_mask/num_trans_in_bins, #map_num_1w_t.sum(),
                 vmin=nmin, vmax=nmax, cmap=plt.cm.magma_r) #cmap=plt.cm.magma_r)      
  plt.title('S$_{adv}$', fontsize = fsize)
#  plt.axis('image')
  
  
  ax4 = plt.subplot(gs[1,2])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  # isobath
  cs4 = ax4.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
 
  #ax4.clabel(cs4, cs4.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)
  
  # GS mean position
  ax4.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1)   
  
  plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                 map_num_sat_t_mask/num_trans_in_bins, #map_num_sat_t.sum(),
                 vmin=nmin, vmax=nmax, cmap=plt.cm.magma_r) #cmap=plt.cm.magma_r)      
  plt.title('S$_{SMAP}$', fontsize = fsize)
#  plt.axis('image')
#  plt.xlabel('Longitude')  
  
  axcb = plt.subplot(gs[0,:])
  cbar = plt.colorbar(pc,  cax=axcb, orientation='horizontal')   
  #cbar.ax.xaxis.set_major_formatter(PercentFormatter(1, 0))    
  #cbar.ax.set_title('Number of fronts [%]', fontsize = fsize)
  cbar.ax.set_title('Occurrence of salinity fronts', fontsize = fsize)
  cbar.ax.tick_params(labelsize=fsize-1) 
  
  xt, yt = 0.03, 0.05
  ax1.text(xt, yt, '(a)',
           transform=ax1.transAxes, fontsize=fsize)
  ax2.text(xt, yt, '(b)',
           transform=ax2.transAxes, fontsize=fsize)
  ax4.text(xt, yt, '(c)',
           transform=ax4.transAxes, fontsize=fsize)



  ''' mean absolute intentisity of fronts  - FOR PAPER '''
  

  #x_grid_2d, y_grid_2d= bm(lon_grid_2d, lat_grid_2d)
  #fsize = 14
  smin, smax = 0.5,1.4 #0, 1.5
  #fig = plt.figure(figsize=(8,8))
  
  #gs = gridspec.GridSpec(3, 2,  height_ratios=[1,20,20])

  ax1 = plt.subplot(gs[3,0])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  #isobath
  cs1 = ax1.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
 
  #ax1.clabel(cs1, cs1.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)
  
  # GS mean position
  ax1.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1) 
  
  pc=plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                    map_int_tsg_t_mkmk,
                    vmin=smin, vmax=smax, cmap=plt.cm.magma_r )#plt.cm.gnuplot2_r) #cmap=plt.cm.magma_r)  
  plt.title('S$_{TSG}$', fontsize = fsize)
#  plt.axis('image')
#  plt.ylabel('Latitude')
  
  
  ax2 = plt.subplot(gs[3,1])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)

  # isobath
  cs2 = ax2.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
 
  #ax2.clabel(cs2, cs2.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)

  # GS mean position
  ax2.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1)   
  
  plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                 map_int_1w_t_mkmk,
                 vmin=smin, vmax=smax, cmap=plt.cm.magma_r) #plt.cm.gnuplot2_r) #cmap=plt.cm.magma_r)      
  plt.title('S$_{adv}$', fontsize = fsize)
#  plt.axis('image')
  

  
  ax4 = plt.subplot(gs[3,2])
  plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
  
  # isobath
  cs4 = ax4.contour(xbat, ybat, topo, levels=[-200],
              colors='c',linewidths=lw,
              linestyles='solid', zorder=1)
 
  #ax4.clabel(cs4, cs4.levels,  fmt='%1.0f', fontsize=fsize-5)#, inline=True, fontsize=fsize-2)
  
  # GS mean position
  ax4.contour(x_mdt, y_mdt, mdt.squeeze(), levels=[0.5],
             colors= 'c', linewidths=lw, linestyles='dashed', zorder=1) 
  
  plt.pcolormesh(x_grid_2d[:-1,:-1], y_grid_2d[:-1,:-1], 
                 map_int_sat_t_mkmk,
                 vmin=smin, vmax=smax, cmap=plt.cm.magma_r) #cmap=plt.cm.magma_r)      
  plt.title('S$_{SMAP}$', fontsize = fsize)
#  plt.axis('image')
#  plt.xlabel('Longitude')

  xt, yt = 0.03, 0.05
  ax1.text(xt, yt, '(d)',
           transform=ax1.transAxes, fontsize=fsize)
  ax2.text(xt, yt, '(e)',
           transform=ax2.transAxes, fontsize=fsize)
  ax4.text(xt, yt, '(f)',
           transform=ax4.transAxes, fontsize=fsize)
  
  axcb = plt.subplot(gs[2,:])
  cbar = plt.colorbar(pc,  cax=axcb, orientation='horizontal')     
  cbar.ax.set_title('Mean intensity of salinity fronts [PSU/km]', 
                    fontsize = fsize)
  cbar.ax.tick_params(labelsize=fsize-1) 
  plt.tight_layout()
  
  fig.savefig(dirfig + 'Oleander_BD_fronts_number_and_intensity_abs_5km_mean_intensity_int.png', dpi=200)

  
  print('Fronts TSG...', bm(452589, 568563,inverse=True))
  print('Fronts SMAP...', bm(418647, 613706,inverse=True))
  print('Intensity...', bm(411550, 628260,inverse=True))
