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
Open file with fronts detected on step 4 and make figures for the revised paper:
    scatter plots with the fronts detected on each Salinity field
    and fronts detected in TSG/Adv TSG/SMAP 
    
    FIGURE FOR PAPER (Figure 5)

written by Bàrbara Barceló-Llull on 14-April-2020 at Mallorca
modified on 14 August 2020: scatter plot for revision


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

def mask_to_identify_fronts(sal_pfronts_all, sal_nfronts_all):
    
  # merge positive and negative fronts and 
  # identify with 1 where there is a front in a grid point
  # and nan where no front
    
  cond_front_p = ~np.isnan(sal_pfronts_all)
  cond_front_n = ~np.isnan(sal_nfronts_all)
  
  
  cond_front = np.logical_or(cond_front_p==True, 
                             cond_front_n==True)
  
  
  mask_front = np.ones(sal_pfronts_all.shape) * np.nan
  mask_front[cond_front] = 1
  
  
  return mask_front

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

  size_dates = len(dates_strings)
  size_data  = 5001
    
  sal_tsgi_all          = np.ones((size_dates, size_data))*np.nan
  sal_tsgi_pfronts_all  = np.ones((size_dates, size_data))*np.nan
  sal_tsgi_nfronts_all  = np.ones((size_dates, size_data))*np.nan

  sal_1wayi_all         = np.ones((size_dates, size_data))*np.nan
  sal_1wayi_pfronts_all = np.ones((size_dates, size_data))*np.nan
  sal_1wayi_nfronts_all = np.ones((size_dates, size_data))*np.nan
    
  sal_sati_all          = np.ones((size_dates, size_data))*np.nan
  sal_sati_pfronts_all  = np.ones((size_dates, size_data))*np.nan
  sal_sati_nfronts_all  = np.ones((size_dates, size_data))*np.nan
    

  
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


    
    dist_int = dict_date['dist_int']
    lon_int  = dict_date['lon_int']
    lat_int  = dict_date['lat_int']

    sal_tsgi_all[idt,:]          = dict_date['sal_tsgi']
    sal_tsgi_pfronts_all[idt,:]  = dict_date['sal_tsgi_pfronts']
    sal_tsgi_nfronts_all[idt,:]  = dict_date['sal_tsgi_nfronts']

    sal_1wayi_all[idt,:]         = dict_date['sal_1wayi']
    sal_1wayi_pfronts_all[idt,:] = dict_date['sal_1wayi_pfronts']
    sal_1wayi_nfronts_all[idt,:] = dict_date['sal_1wayi_nfronts']
    
    sal_sati_all[idt,:]          = dict_date['sal_sati']
    sal_sati_pfronts_all[idt,:]  = dict_date['sal_sati_pfronts']
    sal_sati_nfronts_all[idt,:]  = dict_date['sal_sati_nfronts']
    
  #xx, yy  = np.meshgrid(dates_py, lat_int)
  xx, yy  = np.meshgrid(dates_py, lon_int)
  smin = 30
  smax = 37

  fsize = 12
  
  fig = plt.figure(figsize=(12,5))
  
  
  gs = gridspec.GridSpec(1, 3)

  
  ax1 = plt.subplot(gs[0,0])  

  ax1.scatter(xx, yy, c=sal_tsgi_pfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax1.scatter(xx, yy, c=sal_tsgi_nfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax1.set_title('S$_{TSG}$', fontsize = fsize)
  
  
  ax2 = plt.subplot(gs[0,1])  

  ax2.scatter(xx, yy, c=sal_1wayi_pfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax2.scatter(xx, yy, c=sal_1wayi_nfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax2.set_title('S$_{adv}$', fontsize = fsize)  
  
  
  ax3 = plt.subplot(gs[0,2]) 
  ax3.scatter(xx, yy, c=sal_sati_pfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax3.scatter(xx, yy, c=sal_sati_nfronts_all.T, s=7, vmin=smin, vmax=smax)
  ax3.set_title('S$_{SMAP}$', fontsize = fsize)  
  
  plt.tight_layout()
  
  
  
  ''' 1 = front in that point '''

  mask_front_tsg = mask_to_identify_fronts(sal_tsgi_pfronts_all, 
                                           sal_tsgi_nfronts_all)
  
  mask_front_adv = mask_to_identify_fronts(sal_1wayi_pfronts_all, 
                                           sal_1wayi_nfronts_all)

  mask_front_sat = mask_to_identify_fronts(sal_sati_pfronts_all, 
                                           sal_sati_nfronts_all)

  
  ''' 1 = front in that point in both S_tsg and S_adv '''

  cond_fronts_tsg = ~np.isnan(mask_front_tsg)
  cond_fronts_adv = ~np.isnan(mask_front_adv)
  
  cond_fronts_merge = np.logical_and(cond_fronts_tsg==True, 
                                     cond_fronts_adv==True)

  mask_merge = np.ones(mask_front_tsg.shape) * np.nan
  mask_merge[cond_fronts_merge] = 1

  ''' 1 = front in that point in both S_tsg and S_SMAP '''

  cond_fronts_tsg = ~np.isnan(mask_front_tsg)
  cond_fronts_sat = ~np.isnan(mask_front_sat)
  
  cond_fronts_merge_tsg_smap = np.logical_and(cond_fronts_tsg==True, 
                                     cond_fronts_sat==True)

  mask_merge_tsg_smap = np.ones(mask_front_tsg.shape) * np.nan
  mask_merge_tsg_smap[cond_fronts_merge_tsg_smap] = 1  
  
  ''' Figure '''

  fsize = 16
  marker_size = 1
  
  fig = plt.figure(figsize=(12,4))
  
  gs = gridspec.GridSpec(1, 3)

  ax1 = plt.subplot(gs[0,0])  

  ax1.scatter(xx, yy, c=mask_front_tsg.T, s=marker_size, marker='s', 
              cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  ax1.set_title('Fronts in S$_{TSG}$', fontsize = fsize)
  #ax1.set_ylabel('Latitude [$^{\circ}$N]')
  ax1.set_ylabel('Longitude [$^{\circ}$W]', fontsize=fsize-2)
  date_form = mdates.DateFormatter("%b-%Y")
  ax1.xaxis.set_major_formatter(date_form)
  ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  #ax1.set_ylim(32., 40.5)
  ax1.set_ylim(-74., -64)
  ax1.set_yticklabels((abs(ax1.get_yticks().astype(int))))
  ax1.invert_yaxis()
  ax1.tick_params(labelsize=fsize-2) 
  
  ax2 = plt.subplot(gs[0,1])  

  ax2.scatter(xx, yy, c=mask_front_adv.T, s=marker_size, marker='s', 
              cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  ax2.set_title('Fronts in S$_{adv}$', fontsize = fsize)  
  ax2.xaxis.set_major_formatter(date_form)
  ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  #ax2.set_ylim(32., 40.5)
  ax2.set_ylim(-74., -64)
  ax2.set_yticklabels((abs(ax1.get_yticks().astype(int))))
  ax2.invert_yaxis()
  ax2.tick_params(labelsize=fsize-2) 
  # fronts in TSG and Sadv
  
  
  ax2.scatter(xx, yy, c=mask_merge.T, s=marker_size, marker='s', 
              cmap=plt.cm.Wistia, vmin=0, vmax=1)
  # ax4.set_title('Fronts in S$_{TSG}$ and S$_{adv}$', fontsize = fsize)  
  # ax4.set_ylabel('Latitude [$^{\circ}$N]')  
  # ax4.xaxis.set_major_formatter(date_form)
  # ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  # ax4.set_ylim(32., 40.5)

    
  
  ax3 = plt.subplot(gs[0,2]) 
  ax3.scatter(xx, yy, c=mask_front_sat.T, s=marker_size, marker='s', 
              cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  ax3.set_title('Fronts in S$_{SMAP}$', fontsize = fsize)  
  ax3.xaxis.set_major_formatter(date_form)
  ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  #ax3.set_ylim(32., 40.5)
  ax3.set_ylim(-74., -64)
  ax3.set_yticklabels((abs(ax1.get_yticks().astype(int)))) 
  ax3.invert_yaxis()
  ax3.tick_params(labelsize=fsize-2) 
  
  ax3.scatter(xx, yy, c=mask_merge_tsg_smap.T, s=marker_size, marker='s', 
              cmap=plt.cm.Wistia, vmin=0, vmax=1)

    
  xt, yt = -0.1, 1.1 #0.03, 0.05
  ax1.text(xt, yt, '(a)',
           transform=ax1.transAxes, fontsize=fsize)
  ax2.text(xt, yt, '(b)',
           transform=ax2.transAxes, fontsize=fsize)
  ax3.text(xt, yt, '(c)',
           transform=ax3.transAxes, fontsize=fsize)


  #fig.autofmt_xdate()
  plt.tight_layout()
  
  fig.savefig(dirfig + 'PAPERrev_Oleander_BD_fronts_scatter_comparison.png', dpi=500)


  ''' Figure with 5 subplots '''
  
  # fsize = 12
  
  # fig = plt.figure(figsize=(12,8))
  
  # gs = gridspec.GridSpec(2, 3)

  # ax1 = plt.subplot(gs[0,0])  

  # ax1.scatter(xx, yy, c=mask_front_tsg.T, s=1, marker='s', 
  #             cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  # ax1.set_title('Fronts in S$_{TSG}$', fontsize = fsize)
  # ax1.set_ylabel('Latitude [$^{\circ}$N]')
  # date_form = mdates.DateFormatter("%b-%Y")
  # ax1.xaxis.set_major_formatter(date_form)
  # ax1.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  # ax1.set_ylim(32., 40.5)
  
  # ax2 = plt.subplot(gs[0,1])  

  # ax2.scatter(xx, yy, c=mask_front_adv.T, s=1, marker='s', 
  #             cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  # ax2.set_title('Fronts in S$_{adv}$', fontsize = fsize)  
  # ax2.xaxis.set_major_formatter(date_form)
  # ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  # ax2.set_ylim(32., 40.5)
    
  
  # ax3 = plt.subplot(gs[0,2]) 
  # ax3.scatter(xx, yy, c=mask_front_sat.T, s=1, marker='s', 
  #             cmap=plt.cm.YlGnBu, vmin=0, vmax=1)
  # ax3.set_title('Fronts in S$_{SMAP}$', fontsize = fsize)  
  # ax3.xaxis.set_major_formatter(date_form)
  # ax3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  # ax3.set_ylim(32., 40.5)
    
  
  
  # ax4 = plt.subplot(gs[1,1])  
  
  # ax4.scatter(xx, yy, c=mask_merge.T, s=1, marker='s', 
  #             cmap=plt.cm.Wistia, vmin=0, vmax=1)
  # ax4.set_title('Fronts in S$_{TSG}$ and S$_{adv}$', fontsize = fsize)  
  # ax4.set_ylabel('Latitude [$^{\circ}$N]')  
  # ax4.xaxis.set_major_formatter(date_form)
  # ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))

  # #ax4.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,6)))
  # ax4.set_ylim(32., 40.5)

  
  # ax5 = plt.subplot(gs[1,2])  
  # ax5.scatter(xx, yy, c=mask_merge_tsg_smap.T, s=1, marker='s', 
  #             cmap=plt.cm.Wistia, vmin=0, vmax=1)
  # ax5.set_title('Fronts in S$_{TSG}$ and S$_{SMAP}$', fontsize = fsize)  
  # ax5.xaxis.set_major_formatter(date_form)
  # ax5.xaxis.set_major_locator(mdates.MonthLocator(bymonth=range(1,12,12)))
  # ax5.set_ylim(32., 40.5)
    
  
  # xt, yt = -0.1, 1.05 #0.03, 0.05
  # ax1.text(xt, yt, '(a)',
  #          transform=ax1.transAxes, fontsize=fsize)
  # ax2.text(xt, yt, '(b)',
  #          transform=ax2.transAxes, fontsize=fsize)
  # ax3.text(xt, yt, '(c)',
  #          transform=ax3.transAxes, fontsize=fsize)
  # ax4.text(xt, yt, '(d)',
  #          transform=ax4.transAxes, fontsize=fsize) 
  # ax5.text(xt, yt, '(e)',
  #          transform=ax5.transAxes, fontsize=fsize)  
  
  
  # #fig.autofmt_xdate()
  # plt.tight_layout()
  
  # fig.savefig(dirfig + 'PAPERrev_Oleander_BD_fronts_scatter_comparison.png', dpi=200)

  