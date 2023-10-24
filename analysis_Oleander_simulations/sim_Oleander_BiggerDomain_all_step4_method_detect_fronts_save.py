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
Code to detect and quantify fronts. 

Remove the large scale trend (gradient) and analyze the small scale 
variability, gradients X are fronts. (Desprès et al., 2011)

Separate the type of front (positive or negative).

written by Bàrbara Barceló-Llull on 11-Nov-2019 at Mallorca

adapted for the long simulation on 26-Nov-2019

> adapted for the Bigger Domain and finals simulations on 14-Apr-2020
"""

def distance2apoint(lon, lat, lonp, latp):
  dist = np.zeros(lon.shape)
  for k in np.arange(lon.shape[0]):
    coords_1 = (lon[k], lat[k])
    coords_2 = (lonp, latp)
    #dist[k]  = geopy.distance.vincenty(coords_1, coords_2).km  
    dist[k]  = geopy.distance.geodesic(coords_1, coords_2).km  
    
  return dist

def segments_TSG_trans_mean_std(dist, sal_tsg, lon_tsg, lat_tsg,time_tsg):    
    # Divide data in segments 20 km long
    
    dchunks = np.arange(0, dist.max(), 20)
    
    mean_all = []
    std_all  = []
    
    lon_all  = []
    lat_all  = []
    time_all = []
    dist_all = []
    
    for ii, dchunk in enumerate(dchunks[:-1]):
        
        ind_chunk = np.logical_and(dist>=dchunks[ii], dist<=dchunks[ii+1])
        print('')
        print('Data within distmin=', dchunks[ii])
        print('Data within distmax=', dchunks[ii+1])
        
        
        dist_seg = dist[ind_chunk]
        sal_seg  = sal_tsg[ind_chunk]
        lon_seg  = lon_tsg[ind_chunk]
        lat_seg  = lat_tsg[ind_chunk]
        time_seg = time_tsg[ind_chunk]
        
        
        print('size of data within this chunk...', dist_seg.shape)
        print('')
        
        #plt.figure()
        #plt.plot(dist_seg, sal_seg,'o-')
        
        mean_seg = np.nanmean(sal_seg)
        std_seg  = np.nanstd(sal_seg)
        
        mean_all = np.append(mean_all, mean_seg)
        std_all  = np.append(std_all, std_seg)
        
        # save mid position
        lon_all = np.append(lon_all, np.nanmean(lon_seg))
        lat_all = np.append(lat_all, np.nanmean(lat_seg))
        dist_all = np.append(dist_all, dchunks[ii+1]-10)
        
        #save mid time
        time_all = np.append(time_all, np.nanmean(time_seg))
        
    return  mean_all, std_all, lon_all, lat_all, dist_all, time_all


def edit_variable(sal_tsg, dist, scale):
    
    ''' For TSG even a scale of 2 km removes values that are good 
        Skip this step and assume that the TSG data is good'''
    
    # Repeate editing to keep removing bad values --> Iterative editing
    stsg_edited = sal_tsg.copy() # initialize the variable
    N_out     = 10                  # to start the first editing (it could be any number)
    
    # ITERATIVE EDITING
    # keep editing until N_out = 0
    while N_out>0: 
        
        # Editing: remove 3 times the variance of the data
        # to remove the points that go too far from the other values
        # 1)smooth 2)look at the residuals
        stsg_editsmooth=rt_anatools.loess_smooth_handmade(stsg_edited,1/scale,\
                                                        t=dist)
#        plt.figure()
#        plt.plot(dist, sal_tsg,'o-')
#        plt.plot(dist, stsg_editsmooth,'o-r')
        
        stsg_editresiduals = stsg_edited-stsg_editsmooth
        
#        plt.figure()
#        plt.plot(dist, stsg_editresiduals,'o-')
        
        id_outlier = np.where(np.abs(stsg_editresiduals)>3*np.nanstd(stsg_editresiduals))
        
        # Remove outliers
        stsg_edited[id_outlier[0]] = np.nan
        N_out = len(id_outlier[0])
        
#    plt.figure()
#    plt.plot(dist, sal_tsg,'o-')
#    plt.plot(dist, stsg_edited,'x-')
    
    return stsg_edited

#def detect_fronts_gradient(stsg_residuals, gradient_limit, lim_points):
#    
#    dif_stsg    = np.diff(stsg_residuals)
#    ind_front   = np.where(np.abs(dif_stsg)>=gradient_limit)
#    ind_nofront = np.where(np.abs(dif_stsg)<gradient_limit)
#    
##    plt.figure(figsize=(10,6))
##    plt.plot(dist, stsg_residuals, 'o-')
##    #plt.plot(dist[1:][ind_front], stsg_residuals[1:][ind_front], 'x')
##    plt.plot(dist[1:], stsg_fronts_raw, 'x')
#
#    stsg_fronts = front_condition_consecutive_points(stsg_residuals, ind_nofront, lim_points) 
#    
#    return stsg_fronts

#def detect_fronts_gradient_old(stsg_residuals, gradient_limit, lim_points):
#    
#    dif_stsg      = np.diff(stsg_residuals)
#    
#    ind_front_p   = np.where(dif_stsg >= gradient_limit)
#    ind_nofront_p = np.where(dif_stsg < gradient_limit)
#    
#    ind_front_n   = np.where(dif_stsg <= -gradient_limit)
#    ind_nofront_n = np.where(dif_stsg >  -gradient_limit)
#
#    stsg_fronts_p = front_condition_consecutive_points(stsg_residuals, \
#                                                       ind_nofront_p, lim_points)
#    
#    stsg_fronts_n = front_condition_consecutive_points(stsg_residuals, \
#                                                       ind_nofront_n, lim_points) 
#    
#    return stsg_fronts_p, stsg_fronts_n

def detect_fronts_gradient(stsg_residuals, gradient_limit, lim_points):
    ''' 
    Detect fronts following a gradient-based method
    
    adapted to save both limits of a front on 4 May 2020
    '''
    
    dif_stsg      = np.diff(stsg_residuals)
    
    ind_front_p   = np.where(dif_stsg >= gradient_limit)
    ind_nofront_p = np.where(dif_stsg < gradient_limit)
    
    ind_front_n   = np.where(dif_stsg <= -gradient_limit)
    ind_nofront_n = np.where(dif_stsg >  -gradient_limit)
     
    # when dif>= gradient, the limits of this gradient are good points 
    # (both delimit the front)
    # if in dif i=1 is a front, in sal i=1 and i=2 are good
    
    stsg_res_fronts_p_raw = np.ones(stsg_residuals.shape) * np.nan
    stsg_res_fronts_n_raw = np.ones(stsg_residuals.shape) * np.nan
    
    #pre-limit of the p front
    stsg_res_fronts_p_raw[:-1][ind_front_p] = \
                                   stsg_residuals[:-1][ind_front_p]    
    
    # post-limit of the p front
    stsg_res_fronts_p_raw[1:][ind_front_p] = \
                                   stsg_residuals[1:][ind_front_p]

    #pre-limit of the n front
    stsg_res_fronts_n_raw[:-1][ind_front_n] = \
                                   stsg_residuals[:-1][ind_front_n]    
    
    # post-limit of the n front
    stsg_res_fronts_n_raw[1:][ind_front_n] = \
                                   stsg_residuals[1:][ind_front_n]
                                   
                                   
    stsg_fronts_p = front_condition_consecutive_points(stsg_res_fronts_p_raw, \
                                                       lim_points)
    
    stsg_fronts_n = front_condition_consecutive_points(stsg_res_fronts_n_raw, \
                                                       lim_points) 
    
    return stsg_fronts_p, stsg_fronts_n

def front_condition_consecutive_points(stsg_fronts_raw, lim_points):    

    ''' 
    check if for each point, within the  neighbours there are lim_points 
    consecutive detected points
     
    adapted to save both limits of a front on 4 May 2020
    '''
    
    
    stsg_fronts = np.copy(stsg_fronts_raw)
    
    for ii, st in enumerate(stsg_fronts_raw):
        if np.isnan(st) == False:
            
            # check if a desired number of neighbours are no nan
            
            stsg_neighbours = stsg_fronts_raw[ii-lim_points+1:ii+lim_points]
            num_nonan = len(stsg_neighbours[~np.isnan(stsg_neighbours)])
            
            if num_nonan < lim_points:
                
                # reject this value, it may be noise
                stsg_fronts[ii] = np.nan
                
            
            elif num_nonan >= lim_points:
                
                # check all the combinations of consecutive values
                len_chunks = []
                for jj in np.arange(lim_points):
                    
                    # check each combination
                    s_chunk = stsg_fronts_raw[ii-jj : ii+lim_points-jj]
                    
                    # save the length of each combination
                    len_chunks = np.append(len_chunks, len(s_chunk[~np.isnan(s_chunk)]))
                    
                    
                # if no combination have more tan lim_points values (no nan), reject
                # if one combination is valid, then keep this value
                if np.any(len_chunks==lim_points) == False:
                     
                    # reject this value, it may be noise
                    stsg_fronts[ii] = np.nan 
                     
            
#    plt.figure(figsize=(10,6))
#    plt.plot(dist, stsg_residuals, 'o-b')
#    plt.plot(dist[1:], stsg_fronts_raw, '*r')   
#    plt.plot(dist[1:], stsg_fronts, 'xg')   
    
    return stsg_fronts

#def front_condition_consecutive_points_old(stsg_residuals, ind_nofront, lim_points):    
#    # over at least 3 successive points
#    
#    stsg_fronts_raw = np.copy(stsg_residuals[1:])
#    stsg_fronts_raw[ind_nofront] = np.nan
#    
#    # check if for each point, within the 10 neighbours there are 5 
#    # consecutive detected points
#    
#    
#    stsg_fronts = np.copy(stsg_fronts_raw)
#    
#    for ii, st in enumerate(stsg_fronts_raw):
#        if np.isnan(st) == False:
#            
#            # check if a desired number of neighbours are no nan
#            
#            stsg_neighbours = stsg_fronts_raw[ii-lim_points+1:ii+lim_points]
#            num_nonan = len(stsg_neighbours[~np.isnan(stsg_neighbours)])
#            
#            if num_nonan < lim_points:
#                
#                # reject this value, it may be noise
#                stsg_fronts[ii] = np.nan
#                
#            
#            elif num_nonan >= lim_points:
#                
#                # check all the combinations of consecutive values
#                len_chunks = []
#                for jj in np.arange(lim_points):
#                    
#                    # check each combination
#                    s_chunk = stsg_fronts_raw[ii-jj : ii+lim_points-jj]
#                    
#                    # save the length of each combination
#                    len_chunks = np.append(len_chunks, len(s_chunk[~np.isnan(s_chunk)]))
#                    
#                    
#                # if no combination have more tan lim_points values (no nan), reject
#                # if one combination is valid, then keep this value
#                if np.any(len_chunks==lim_points) == False:
#                     
#                    # reject this value, it may be noise
#                    stsg_fronts[ii] = np.nan 
#                     
#            
##    plt.figure(figsize=(10,6))
##    plt.plot(dist, stsg_residuals, 'o-b')
##    plt.plot(dist[1:], stsg_fronts_raw, '*r')   
##    plt.plot(dist[1:], stsg_fronts, 'xg')   
#    
#    return stsg_fronts

def mask_where_nans(dist_ori, s1way_adv_ori, s1way_adv_raw, dist):
    
    mask_nans_1way = np.ones(dist_ori.shape)
    mask_nans_1way[np.isnan(s1way_adv_ori)] = 0
    
    f_mask_1way = interpolate.interp1d(dist_ori, mask_nans_1way, 
                                       kind='nearest',fill_value="extrapolate")
    maski_1way  = f_mask_1way(dist)  
    
    maski_1way_nan = np.copy(maski_1way)
    maski_1way_nan[maski_1way_nan<1] = np.nan
    
    s1way_adv = s1way_adv_raw * maski_1way_nan

    return s1way_adv


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

def mask_gaps(dist_ori, dist, dx, sal_smooth_raw, sal_ori):
    

  # gaps in dist and also nans in orignal data 
  dist_nans = dist_ori[~np.isnan(sal_ori)]
  
  # search gaps (3dx?)
  dif_dist       = np.diff(dist_nans)
  inds_gaps_pre  = np.where(np.abs(dif_dist) > 3*dx)[0]
  inds_gaps_post = inds_gaps_pre + 1       
          
  dist_pre  = dist_nans[inds_gaps_pre]
  dist_post = dist_nans[inds_gaps_post.astype(int)]

  # mask where there are no data, the gaps
        
  sal_masked = np.copy(sal_smooth_raw)
           
  for ii, dpre in enumerate(dist_pre):
    cond = np.logical_and(dist>min(dpre, dist_post[ii]),
                                  dist<max(dpre, dist_post[ii]))
            
    sal_masked[cond] = np.nan
            
  return sal_masked

def mask_fronts_original(dist, dist_ori, cond_fronts):
    
      ''' mask = nan when no front, mask =1 when there is a front '''

      dist_int_fronts_all = np.copy(dist)#[1:])
      dist_int_fronts_all[cond_fronts] = np.nan
      dist_fronts_list = detect_chunks(dist_int_fronts_all)
      
      mask_ori_fronts = np.ones(dist_ori.shape)*np.nan
      
      for ifr in np.arange(len(dist_fronts_list)):
          
          distif = dist_fronts_list[ifr]
          
          ind_ori_front = np.where(np.logical_and(dist_ori >= distif.min(),
                                  dist_ori <= distif.max()))
            
          # mask = 1 when front, mask=nan when no front  
          mask_ori_fronts[ind_ori_front] = 1
          
      return  mask_ori_fronts  
    
if __name__ == '__main__':
    
  plt.close('all')
  
  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/'+ \
                'fig_Oleander_BiggerRegion/fronts/'
  # v4 SMAP data              
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all_Bigger_Domain/post-processing/'
  
  ''' Parameters for the detection '''
  
  dx     = 0.2 # km
  dist   = np.arange(0, 1000+dx, dx)# km
  
  L_small = 10 #10 #km
  L_large = 200#140#70 #km SMAP smooth scale

  #gradient_limit = (0.2/20) * dx # this limit is based on original data, not valid for residuals
  gradient_limit = (0.02) * dx#(0.04) * dx # gradient on residuals 0.02 PSU/km

  lim_points = round(5/dx) #3 2.5km minimum length of fronts
  
  plot_fig = True
  
  ''' File with Backward advected field interpolated onto TSG position '''
  
  file_dic = 'SimOl_B_alt_BD_sal_int_to_TSG_final.pkl'

  f = open(dir_dic + file_dic, 'rb')
  dict_1way = pickle.load(f)
  f.close()  
  
  
  ''' List of simulated dates after doing the comparison '''
  
  dates_strings = list(dict_1way.keys())
  
  
  ''' Loop for each simulation '''

  # FINAL BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46
  
#  lon_min = -73.9
#  lat_max = 40.0
  
  dic_all_dates = {}
  
  select_dates_to_check = ['20160928']
  
  for date in dates_strings:

    date_time_obj = datetime.strptime(date, '%Y%m%d')

    print('')
    print ('date... ', date_time_obj)
    
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))  
    

    
    ''' Open Backward data for this date and 7 days of running time '''   
    
    dict_1way_date = dict_1way[date]['07']
              
    s1way_adv_ori = dict_1way_date['s_advi']
    s_sat_ori     = dict_1way_date['s_smapi']
    
    # TSG data
    sal_tsg_ori  = dict_1way_date['sal_tsg']
    
    time_tsg = dict_1way_date['time_tsg']
    lon_tsg  = dict_1way_date['lon_tsg']
    lat_tsg  = dict_1way_date['lat_tsg']
    
    # convert (lon, lat) to distance [km], different for each transect!!
    #dist_ori = distance2apoint(lon_tsg, lat_tsg, lon_min, lat_max)
    
    # we don't want to work on distance, so here is only to infer the smooth of the fields
    # (lon_tsg.min(), lat_tsg.max()) will be the value closer to the coast
    dist_ori = distance2apoint(lon_tsg, lat_tsg, lon_tsg.min(), lat_tsg.max())

    
    
    ''' Interpolate data to a regular distance axis (intelligent interpolation)'''

    print('Interpolating data for this date...')
    print('') 
      
    # Smooth to have represented the same scales in all fields
    sal_tsg_raw  = rt_anatools.loess_smooth_handmade(sal_tsg_ori,1/L_small,\
                                                 t_final=dist,t=dist_ori)
    
    s_sat_raw  = rt_anatools.loess_smooth_handmade(s_sat_ori,1/L_small,\
                                                 t_final=dist,t=dist_ori)
    
    s1way_adv_raw  = rt_anatools.loess_smooth_handmade(s1way_adv_ori,1/L_small,\
                                                 t_final=dist,t=dist_ori)
    

    # mask where nans in original data
    
    sal_tsg   = mask_gaps(dist_ori, dist, dx, sal_tsg_raw,   sal_tsg_ori)
    s_sat     = mask_gaps(dist_ori, dist, dx, s_sat_raw,     s_sat_ori)
    s1way_adv = mask_gaps(dist_ori, dist, dx, s1way_adv_raw, s1way_adv_ori)

    print('')
    print('Data interpolated!')
    print('') 
    
#    plt.figure(figsize=(8,10))
#    plt.subplot(311)
#    plt.plot(dist_ori, sal_tsg_ori, '-k', linewidth=3)
#    #plt.plot(dist, sal_tsg_raw, '-r')
#    plt.plot(dist, sal_tsg, '--y', linewidth=1.5)
#    plt.title('TSG ' + date)
#
#    plt.subplot(312)
#    plt.plot(dist_ori, s1way_adv_ori, '-k', linewidth=3)
#    plt.plot(dist, s1way_adv, '--y', linewidth=1.5)
#    plt.title('adv')
#    
#    plt.subplot(313)
#    plt.plot(dist_ori, s_sat_ori, '-k', linewidth=3)
#    plt.plot(dist, s_sat, '--y', linewidth=1.5)
#    plt.title('sat')
#    
#    plt.tight_layout()
#    plt.savefig(dirfig + 'Front_detection_step1_interp_and_filter_'+ date +'.png', dpi=200)
#

    
    
    ''' Only continue if there are data for each field '''

    if np.logical_and(np.logical_and(len(sal_tsg[~np.isnan(sal_tsg)]) > 0,
                                   len(s1way_adv[~np.isnan(s1way_adv)]) > 0 ),
                                   len(s_sat[~np.isnan(s_sat)]) > 0 ):
                         
      print('Inferring residuals for this date...')
      print('')     
                         
      ''' Infer the large scale trend for each salinity field'''
 

      # smooth for the desired scale (10 km)
      stsg_large_scale  = rt_anatools.loess_smooth_handmade(sal_tsg,1/L_large,\
                                                         t_final=dist,t=dist)
      s1way_large_scale_ori = rt_anatools.loess_smooth_handmade(s1way_adv,1/L_large,\
                                                         t_final=dist,t=dist)
      ssat_large_scale = rt_anatools.loess_smooth_handmade(s_sat,1/L_large,\
                                                         t_final=dist,t=dist)
    
      # nan where no original data
      s1way_large_scale  = np.copy(s1way_large_scale_ori)
      s1way_large_scale[np.isnan(s1way_adv)] = np.nan
    
    
#      plt.figure(figsize=(8,10))
#      plt.subplot(311)
#      plt.plot(dist, sal_tsg, '-k', linewidth=2)
#      plt.plot(dist, stsg_large_scale, '-r', linewidth=2)
#      plt.title('TSG')
#    
#      plt.subplot(312)
#      plt.plot(dist, s1way_adv, '-k', linewidth=2)
#      plt.plot(dist, s1way_large_scale, '-r', linewidth=2)
#      plt.title('ADV')    
#
#      plt.subplot(313)
#      plt.plot(dist, s_sat, '-k', linewidth=2)
#      plt.plot(dist, ssat_large_scale, '-r', linewidth=2)
#      plt.title('SMAP')
#      
#      plt.tight_layout()
#      plt.savefig(dirfig + 'Front_detection_step2_large_scale_'+ date +'.png', dpi=200)
#      
      
      
      
      ''' Infer the residuals for each field '''
  
      stsg_residuals       = sal_tsg - stsg_large_scale
      s1way_residuals_noqc = s1way_adv - s1way_large_scale
      ssat_residuals       = s_sat - ssat_large_scale
    
    # Apply QC to the residuals before writting the QC for the advected fields
#    max_residual = 0.5
#    
#    s1way_residuals = np.copy(s1way_residuals_noqc)
#    s2way_residuals = np.copy(s2way_residuals_noqc)
#    
#    s1way_residuals[np.abs(s1way_residuals_noqc)>max_residual] = np.nan
#    s2way_residuals[np.abs(s2way_residuals_noqc)>max_residual] = np.nan
    
      s1way_residuals = s1way_residuals_noqc
      
      print('Residuals inferred!')
      print('')  
      
#    plt.figure()
#    plt.plot(dist, s1way_residuals_noqc, 'ob')
#    plt.plot(dist, s1way_residuals, 'xr')
      
      ''' Plot the residuals '''
    
      if plot_fig == True:
        
        plt.figure(figsize=(10,10))
        plt.subplot(311)
        plt.plot(dist, sal_tsg, color='k', linewidth =2, label='TSG') 
        plt.plot(dist, s_sat, color='y', linewidth =2, linestyle =  '--', label='smap')    
        plt.plot(dist, s1way_adv, color='r', linewidth =1.5 ,label='1 way') 
 
        plt.legend()
        #plt.xlabel('Distance [km]')
        plt.ylabel('Original Salinity')
        plt.xlim(dist.min(), dist.max())
        plt.title('Front detection - ' + mdates.num2date(date_fin).strftime("%Y%m%d"))
    
    
        plt.subplot(313)
        plt.plot(dist, stsg_residuals, '-', color='k', linewidth =2, markersize=4, 
             label='TSG')
        plt.plot(dist, ssat_residuals, color='y', linewidth =2, linestyle =  '--', 
               label='smap')
        plt.plot(dist, s1way_residuals, color='r', linewidth =1.5 ,label='1 way')
    
        #plt.xlabel('Distance [km]')
        plt.ylabel('Residual Salinity')
        plt.xlim(dist.min(), dist.max())
    
        plt.subplot(312)
        plt.plot(dist, stsg_large_scale, '-', color='k', linewidth =2, markersize=4, 
             label='TSG')
        plt.plot(dist, ssat_large_scale, color='y', linewidth =2, linestyle =  '--', 
              label='smap')
        plt.plot(dist, s1way_large_scale, color='r', linewidth =1.5 ,label='1 way')
      
        plt.xlabel('Distance [km]')
        plt.ylabel('Large-scale Salinity')
        plt.xlim(dist.min(), dist.max())    
 
        plt.tight_layout()
        plt.savefig(dirfig + 'Front_detection_step3_residuals_' + date + '.png', dpi=200)
        

      ''' Detect front '''
      
      print('Detecting fronts for this date...')
      print('')
    
      stsg_fronts_p, stsg_fronts_n  = detect_fronts_gradient(stsg_residuals,  gradient_limit, lim_points)
      ssat_fronts_p, ssat_fronts_n  = detect_fronts_gradient(ssat_residuals,  gradient_limit, lim_points)
      s1way_fronts_p, s1way_fronts_n = detect_fronts_gradient(s1way_residuals, gradient_limit, lim_points)

      print('Fronts detected!')
      print('')
      
      # plot
      cmin = 0.5
      res_lim = 1.5
    
      plt.figure(figsize=(10,10))
      plt.subplot(311)
      plt.plot(dist, stsg_residuals, '-', color='0.7', linewidth =2,  
             label='TSG')
    
      plt.plot(dist, stsg_fronts_p, '-', color='r', linewidth =4, 
             label='TSG')    
      plt.plot(dist, stsg_fronts_n, '-', color='b', linewidth =4, 
             label='TSG')      
      plt.ylabel('TSG')
      plt.xlim(dist.min(), dist.max()) 
      plt.ylim(-cmin, cmin)
      plt.yticks([-res_lim,0,res_lim], [-res_lim,0,res_lim])
      plt.title('Front detection Method1 - ' + mdates.num2date(date_fin).strftime("%Y%m%d"))

      plt.subplot(312)
      plt.plot(dist, s1way_residuals,  '-', color='0.7', linewidth =2,  label='1 way')
      plt.plot(dist, s1way_fronts_p, '-', color='r', linewidth =4,label='1 way')
      plt.plot(dist, s1way_fronts_n, '-', color='b', linewidth =4,label='1 way')
    
      plt.ylabel('ADV')
      plt.xlim(dist.min(), dist.max()) 
      plt.ylim(-cmin, cmin)
      plt.yticks([-res_lim,0,res_lim], [-res_lim,0,res_lim])
    
      plt.subplot(313)
      plt.plot(dist, ssat_residuals,  '-', color='0.7', linewidth =2,  label='smap')
      plt.plot(dist, ssat_fronts_p,  '-', '-', color='r', linewidth =4, label='smap')
      plt.plot(dist, ssat_fronts_n,  '-', '-', color='b', linewidth =4, label='smap')

      plt.xlabel('Distance [km]')
      plt.ylabel('SMAP')
      plt.xlim(dist.min(), dist.max())   
      plt.ylim(-cmin, cmin)
      plt.yticks([-res_lim,0,res_lim], [-res_lim,0,res_lim])
    
      plt.savefig(dirfig + 'Front_detection_step4_detect_fronts_' + date + '.png', dpi=200)
      
      
    
      ''' Fronts on the interpolated and smoothed ("original") data ''' 
      
      # positive fronts
    
      cond_tsg_p = np.isnan(stsg_fronts_p)
      cond_sat_p = np.isnan(ssat_fronts_p)
      cond_1w_p  = np.isnan(s1way_fronts_p)
    
      stsg_fronts_int_p = np.copy(sal_tsg)
      ssat_fronts_int_p = np.copy(s_sat)
      s1w_fronts_int_p  = np.copy(s1way_adv)
    
      stsg_fronts_int_p[cond_tsg_p] = np.nan
      ssat_fronts_int_p[cond_sat_p] = np.nan
      s1w_fronts_int_p[cond_1w_p]   = np.nan

      # negative fronts
    
      cond_tsg_n = np.isnan(stsg_fronts_n)
      cond_sat_n = np.isnan(ssat_fronts_n)
      cond_1w_n  = np.isnan(s1way_fronts_n)
    
      stsg_fronts_int_n = np.copy(sal_tsg)
      ssat_fronts_int_n = np.copy(s_sat)
      s1w_fronts_int_n  = np.copy(s1way_adv)
    
      stsg_fronts_int_n[cond_tsg_n] = np.nan
      ssat_fronts_int_n[cond_sat_n] = np.nan
      s1w_fronts_int_n[cond_1w_n]   = np.nan


      ''' Fronts on the original data ''' 
      
      mask_ori_tsg_pf = mask_fronts_original(dist, dist_ori, cond_tsg_p)
      mask_ori_sat_pf = mask_fronts_original(dist, dist_ori, cond_sat_p)
      mask_ori_s1w_pf = mask_fronts_original(dist, dist_ori, cond_1w_p)
       
      mask_ori_tsg_nf = mask_fronts_original(dist, dist_ori, cond_tsg_n)
      mask_ori_sat_nf = mask_fronts_original(dist, dist_ori, cond_sat_n)
      mask_ori_s1w_nf = mask_fronts_original(dist, dist_ori, cond_1w_n)
       
      sal_tsg_ori_pf = sal_tsg_ori   * mask_ori_tsg_pf
      sal_sat_ori_pf = s_sat_ori     * mask_ori_sat_pf
      sal_1w_ori_pf  = s1way_adv_ori * mask_ori_s1w_pf

      sal_tsg_ori_nf = sal_tsg_ori   * mask_ori_tsg_nf
      sal_sat_ori_nf = s_sat_ori     * mask_ori_sat_nf
      sal_1w_ori_nf  = s1way_adv_ori * mask_ori_s1w_nf
      
       
#      plt.figure()
#      plt.subplot(211)
#      plt.plot(dist_ori, sal_tsg_ori, 'k-')
#      plt.plot(dist_ori, sal_tsg_ori_pf, 'r-')
#      plt.plot(dist_ori, sal_tsg_ori_nf, 'b-')
#      plt.subplot(212)
#      plt.plot(dist, sal_tsg, 'k-')
#      plt.plot(dist[1:], stsg_fronts_int_p, 'r-')  
#      plt.plot(dist[1:], stsg_fronts_int_n, 'b-')  
#      
#      plt.figure()
#      plt.subplot(211)
#      plt.plot(dist_ori, s1way_adv_ori, 'k-')
#      plt.plot(dist_ori, sal_1w_ori_pf, 'r-')
#      plt.plot(dist_ori, sal_1w_ori_nf, 'b-')
#      plt.subplot(212)
#      plt.plot(dist, s1way_adv, 'k-')
#      plt.plot(dist[1:], s1w_fronts_int_p, 'r-')  
#      plt.plot(dist[1:], s1w_fronts_int_n, 'b-')        
      
      '''
      Interpolate lon_ori, lat_ori onto dist_int 
      to have lon_int, lat_int
      '''
    
      flon = interpolate.interp1d(dist_ori, lon_tsg, 
                                  bounds_error=False,
                                  fill_value= (np.nan, np.nan))
      lon_int = flon(dist)
    
      flat = interpolate.interp1d(dist_ori, lat_tsg, 
                                  bounds_error=False,
                                  fill_value= (np.nan, np.nan))
      lat_int = flat(dist)
    
#      plt.figure()
#      plt.scatter(lon_tsg, lat_tsg, s=5, marker = 'o', color='k')
#      plt.scatter(lon_int, lat_int, s=5, marker = 'x', color='r')
             
    
      # Figure
      
      if plot_fig == True:    
        
        
        cmin = 30 #np.nanmin(sal_tsg)-0.1 # 0.5
        cmax = 37 #np.nanmax(sal_tsg)+1
    
        plt.figure(figsize=(10,10))
        plt.subplot(311)
        plt.plot(dist, sal_tsg, '-', color='0.7', linewidth =2,  
             label='TSG')
    
        plt.plot(dist, stsg_fronts_int_p, '-', color='r', linewidth =4, 
             label='TSG')    
        plt.plot(dist, stsg_fronts_int_n, '-', color='b', linewidth =4, 
             label='TSG')        
        plt.ylabel('TSG')
        plt.xlim(dist.min(), dist.max()) 
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
        plt.title('Front detection - ' + date)

        plt.subplot(312)
        plt.plot(dist, s1way_adv,  '-', color='0.7', linewidth =2,  label='1 way')
        plt.plot(dist, s1w_fronts_int_p, '-', color='r', linewidth =4,label='1 way')
        plt.plot(dist, s1w_fronts_int_n, '-', color='b', linewidth =4,label='1 way')
    
        plt.ylabel('ADV')
        plt.xlim(dist.min(), dist.max()) 
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])

    
        plt.subplot(313)
        plt.plot(dist, s_sat,  '-', color='0.7', linewidth =2,  label='smap')
        plt.plot(dist, ssat_fronts_int_p,  '-', '-', color='r', linewidth =4, label='smap')
        plt.plot(dist, ssat_fronts_int_n,  '-', '-', color='b', linewidth =4, label='smap')

        plt.xlabel('Distance [km]')
        plt.ylabel('SMAP')
        plt.xlim(dist.min(), dist.max())   
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
        
        plt.tight_layout()
        plt.savefig(dirfig + 'Front_detection_step5_plot_fronts_int_' + date + '.png', dpi=200)

#
#        # Fronts on original data 
#        
        plt.figure(figsize=(10,10))
        plt.subplot(311)
        plt.plot(dist_ori, sal_tsg_ori, '-', color='0.7', linewidth =2,  
             label='TSG')
    
        plt.plot(dist_ori, sal_tsg_ori_pf, '-', color='r', linewidth =4, 
             label='TSG')    
        plt.plot(dist_ori, sal_tsg_ori_nf, '-', color='b', linewidth =4, 
             label='TSG')        
        plt.ylabel('TSG')
        plt.xlim(dist.min(), dist.max()) 
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
        plt.title('Front detection - ' + date)

        plt.subplot(312)
        plt.plot(dist_ori, s1way_adv_ori,  '-', color='0.7', linewidth =2,  label='1 way')
        plt.plot(dist_ori, sal_1w_ori_pf, '-', color='r', linewidth =4,label='1 way')
        plt.plot(dist_ori, sal_1w_ori_nf, '-', color='b', linewidth =4,label='1 way')
    
        plt.ylabel('1-way')
        plt.xlim(dist.min(), dist.max()) 
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])

    
        plt.subplot(313)
        plt.plot(dist_ori, s_sat_ori,  '-', color='0.7', linewidth =2,  label='smap')
        plt.plot(dist_ori, sal_sat_ori_pf,  '-', '-', color='r', linewidth =4, label='smap')
        plt.plot(dist_ori, sal_sat_ori_nf,  '-', '-', color='b', linewidth =4, label='smap')

        plt.xlabel('Distance [km]')
        plt.ylabel('SMAP')
        plt.xlim(dist.min(), dist.max())   
        plt.ylim(cmin, cmax)
        #plt.yticks([-0.25,0,0.25], [-0.25,0,0.25])
        
        plt.tight_layout()
        plt.savefig(dirfig + 'Front_detection_step6_plot_fronts_orign_' + date + '.png', dpi=200)
        
        plt.close('all')


      ''' Save data for each date '''
      # we neet to save also lon_int and lat_int!
      
      dic_tran = {'dist_int'        : dist,
                  'lon_int'         : lon_int,
                  'lat_int'         : lat_int,
                  
                  'sal_tsgi'          : sal_tsg,
                  'sal_tsgi_pfronts'  : stsg_fronts_int_p,
                  'sal_tsgi_nfronts'  : stsg_fronts_int_n,
                  
                  'sal_1wayi'         : s1way_adv,
                  'sal_1wayi_pfronts' : s1w_fronts_int_p, 
                  'sal_1wayi_nfronts' : s1w_fronts_int_n,

                  'sal_sati'          : s_sat,
                  'sal_sati_pfronts'  : ssat_fronts_int_p,
                  'sal_sati_nfronts'  : ssat_fronts_int_n,
                  
                  'dist_ori'          : dist_ori,
                  'lon_ori'           : lon_tsg,
                  'lat_ori'           : lat_tsg,
                  
                  'sal_tsgo'          : sal_tsg_ori,
                  'sal_tsgo_pfronts'  : sal_tsg_ori_pf,
                  'sal_tsgo_nfronts'  : sal_tsg_ori_nf,
                  
                  'sal_1wayo'         : s1way_adv_ori,
                  'sal_1wayo_pfronts' : sal_1w_ori_pf, 
                  'sal_1wayo_nfronts' : sal_1w_ori_nf,
                  
                  'sal_sato'          : s_sat_ori,
                  'sal_sato_pfronts'  : sal_sat_ori_pf,
                  'sal_sato_nfronts'  : sal_sat_ori_nf,                  
                  # also fronts on the original (raw) data
                  }    
    
      dic_all_dates.update({date: dic_tran})          
      
    
  ''' Save fronts detected data '''
  
  f = open(dir_dic + 'Oleander_BD_front_detection_5km.pkl','wb')
  pickle.dump(dic_all_dates,f)
  f.close()     
