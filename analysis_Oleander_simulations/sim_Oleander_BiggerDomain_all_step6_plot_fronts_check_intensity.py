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

"""
Open file with fronts detected on step 4 and check that
the relative intensity of fronts on interpolated data
is correct

written by Bàrbara Barceló-Llull on 22-April-2020 at Mallorca

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
    
    plt.figure()
    plt.scatter(lon_ori,lat_ori, c=sal_fronts)
    plt.colorbar()
    
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
        print('abs int...', intensity_front_abs )
        # distance between (latmax, lonmin) and (latmin, lonmax)

        coords_1 = (lon_front.min(), lat_front.max())
        coords_2 = (lon_front.max(), lat_front.min())
        
        len_front  = geopy.distance.geodesic(coords_1, coords_2).km 
        print('length...', len_front )
        
        plt.figure()
        plt.subplot(211)
        plt.scatter(lon_front, lat_front, c=sal_front)
        plt.colorbar()
        plt.subplot(212)
        plt.plot(lon_front, sal_front)
        
        # absolute front intensity/front length
        intensity_front_rel = intensity_front_abs/len_front #PSU/km
        print('rel int...', intensity_front_rel )
        print('')
        
        
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

  
def map_all(lon_fronts, lat_fronts, sal_fronts, lon_grid, lat_grid):  
    
  if len(sal_fronts) >0:
    
    retnum = stats.binned_statistic_2d(lon_fronts, lat_fronts, 
                                          None, 'count', bins=[lon_grid,lat_grid])
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

def mask_not_enough_trans_in_bins(map_num_tsg_p):
    
  map_num_tsg_p_mask = np.copy(map_num_tsg_p)
  map_num_tsg_p_mask[cond_trans==True] = np.nan
  
  return map_num_tsg_p_mask

def check_intensity_fronts_int(dist_int, sal_fronts):
    
    #lon_ori, lat_ori, sal_1wayo_nfronts
    ''' Number of fronts per transect and intensity of fronts
    inferred as the difference between maximum and minimum salinity
    divided by length of fronts '''
    
    plt.figure()
    plt.plot(dist_int, sal_fronts)
    
    if len(sal_fronts[~np.isnan(sal_fronts)]) >0:
        
      ''' 
      map of the number of fronts in each bin and mean intensity 
      also return the mean lon, mean lat and absolute intensity of each front
      '''
    
      dist_fronts = np.copy(dist_int)
    
      dist_fronts[np.isnan(sal_fronts)] = np.nan
    
      dist_fronts_list = detect_chunks(dist_fronts)
      sal_fronts_list  = detect_chunks(sal_fronts)
    

      
      for ifront in np.arange(len(dist_fronts_list)):
        
        dist_front = dist_fronts_list[ifront]
        sal_front  = sal_fronts_list[ifront]
        
        intensity_front_abs = np.max(sal_front) - np.min(sal_front)
        print('abs int...', intensity_front_abs )
        
        len_front  = dist_front.max() - dist_front.min()
        print('length...', len_front )
        
        plt.figure()
        plt.plot(dist_front, sal_front)
        
        # absolute front intensity/front length
        intensity_front_rel = intensity_front_abs/len_front #PSU/km
        print('rel int...', intensity_front_rel )
        print('')
        
        
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
  
  
  for idt, date in enumerate(['20160928']):

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

 
    dist_int = dict_date['dist_int']
    lon_ori  = dict_date['lon_ori']
    lat_ori  = dict_date['lat_ori']

#    sal_tsgi          = dict_date['sal_tsgi']
#    sal_tsgi_pfronts  = dict_date['sal_tsgi_pfronts']
#    sal_tsgi_nfronts  = dict_date['sal_tsgi_nfronts']

    sal_1wayi         = dict_date['sal_1wayi']
    sal_1wayi_pfronts = dict_date['sal_1wayi_pfronts']
    sal_1wayi_nfronts = dict_date['sal_1wayi_nfronts']

    sal_1wayo         = dict_date['sal_1wayo']
    sal_1wayo_pfronts = dict_date['sal_1wayo_pfronts']
    sal_1wayo_nfronts = dict_date['sal_1wayo_nfronts']
    
    plt.figure()
    plt.subplot(211)
    plt.plot(dist_int, sal_1wayi, '-k')
    plt.plot(dist_int[1:], sal_1wayi_pfronts, '-r')
    plt.plot(dist_int[1:], sal_1wayi_nfronts, '-b')
    plt.subplot(212)
    plt.plot(lon_ori, sal_1wayo, '-k')
    plt.plot(lon_ori, sal_1wayo_pfronts, '-r')
    plt.plot(lon_ori, sal_1wayo_nfronts, '-b')
    
#    sal_sati          = dict_date['sal_sati']
#    sal_sati_pfronts  = dict_date['sal_sati_pfronts']
#    sal_sati_nfronts  = dict_date['sal_sati_nfronts']
    
    check_intensity_fronts_int(dist_int[1:], sal_1wayi_nfronts) 
    
    check_intensity_fronts_int(dist_int[1:], sal_1wayi_pfronts)
        
#
#def check_intensity_fronts_orig(dist_int, sal_fronts):
#    
#    #lon_ori, lat_ori, sal_1wayo_nfronts
#    ''' Number of fronts per transect and intensity of fronts
#    inferred as the difference between maximum and minimum salinity
#    divided by length of fronts '''
#    
#    plt.figure()
#    plt.plot(dist_int, sal_fronts)
#    
#    if len(sal_fronts[~np.isnan(sal_fronts)]) >0:
#        
#      ''' 
#      map of the number of fronts in each bin and mean intensity 
#      also return the mean lon, mean lat and absolute intensity of each front
#      '''
#    
#      dist_fronts = np.copy(dist_int)
#    
#      dist_fronts[np.isnan(sal_fronts)] = np.nan
#    
#      dist_fronts_list = detect_chunks(dist_fronts)
#      sal_fronts_list  = detect_chunks(sal_fronts)
#    
#
#      
#      for ifront in np.arange(len(dist_fronts_list)):
#        
#        dist_front = dist_fronts_list[ifront]
#        sal_front  = sal_fronts_list[ifront]
#        
#        intensity_front_abs = np.max(sal_front) - np.min(sal_front)
#        print('abs int...', intensity_front_abs )
#        
#        len_front  = dist_front.max() - dist_front.min()
#        print('length...', len_front )
#        
#        plt.figure()
#        plt.plot(dist_front, sal_front)
#        
#        # absolute front intensity/front length
#        intensity_front_rel = intensity_front_abs/len_front #PSU/km
#        print('rel int...', intensity_front_rel )
#        print('')
