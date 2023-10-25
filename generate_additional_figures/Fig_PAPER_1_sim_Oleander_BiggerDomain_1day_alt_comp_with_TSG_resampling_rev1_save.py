#usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
import pickle
from matplotlib              import dates as mdates
from datetime                import datetime
import sim_toolbox           as sto
from scipy.interpolate       import griddata
import Tools_OI              as toi 


"""  

> written by Bàrbara Barceló-Llull on 14 August 2020:
    
    To create the final Figure 1 for the revised manuscript: 
    coarse grain the reconstructed field to the original resolution
    
    Here we resample and smooth S_adv and save data in a .pkl file
    named S_adv_resampled_smoothed_20170212.pkl 
    
    The final Figure 1 is created in the code: 
    Fig_PAPER_1_sim_Oleander_BiggerDomain_1day_alt_comp_with_TSG_resampling_rev1_plot.py
        

"""



if __name__ == '__main__':
    
  plt.close('all')
  
  ''' Directories '''
  
  dirsave     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all/sss_tagged/'

  # v4 SMAP data        
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'

  
  ''' Set parameters for the code '''

  # running days
  rds = [7]

  # TRYING A BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46
  
      
  # first day of the simulation in the GS 
  # when using delayed time altimetry data for the massive simulations
  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0))  
  date_ori_nrt = mdates.date2num(datetime(2019, 5, 14, 0, 0, 0))   

  
  i_change     = 132 # '20190110' last date of dt 
  
  
  ''' File with the simulation dates '''  
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'
  file_dic   = 'Oleander_TSG_per_transect_all_BiggerDomain.pkl'


  f = open(dir_dic + file_dic, 'rb')
  dict_tsg_all = pickle.load(f)
  f.close() 

  
  
  ''' For the simulation on the selected date... '''
  
  dic_all_dates    = {}
  
  select_dates = ['20170212']
  
  n_short = 0
  
  for idt, date in enumerate(select_dates): #dates_strings):#dates_strings:

      
    date_time_obj = datetime.strptime(date, '%Y%m%d')
    
    print('')
    print ('date... ', date_time_obj)
    
    if idt <= i_change:
      date_ori = date_ori_dt  
      #print('using dt ini date')
      
    elif idt > i_change:
      date_ori = date_ori_nrt
      #print('using nrt ini date')
      
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))


    for rd in rds:

             
        filename = 'sim_Oleander_back_alt_BD_'+ mdates.num2date(date_fin).strftime("%Y%m%d")+ \
                  '_' + '%02i'% rd + 'days'+'_sss.nc'
    
        print('')
        print('Openning...' + filename)
      
      
        ''' Open simulation data '''
      

        ind_tf = 0  #first observation, beginning of simulation at TF
        ind_ti = -1 #last observation, end of simulation at T0
    
        nc    = netcdf.Dataset(dirsave + filename, 'r')
    
        # read initial position at TF
        #traj_tf  = nc.variables['trajectory'][:, ind_tf].data  # (traj, obs)
        lat_tf   = nc.variables['lat'][:, ind_tf].data
        lon_tf   = nc.variables['lon'][:, ind_tf].data   
        time_tf  = nc.variables['time'][:, ind_tf].data   #"seconds since 2017-01-01T00:00:00.000000000"

        # read final position at T0
        #traj_ti  = nc.variables['trajectory'][:, ind_ti].data  # (traj, obs)
        #lat_ti   = nc.variables['lat'][:, ind_ti].data
        #lon_ti   = nc.variables['lon'][:, ind_ti].data
        #time_ti  = nc.variables['time'][:, ind_ti].data

    
        # SSS tagged to each particle
        sss       = nc.variables['sss_adv'][:].data
        nc.close()
   
      
        # from 360 format to -180... 180
        lon_tf[lon_tf>180] = lon_tf[lon_tf>180] - 360
  
        # rename variables
        #time_ini = time_ti
        time_fin = time_tf

        lon_fin = lon_tf
        lat_fin = lat_tf
      
      
        # mask where sss is -9999 (salinity coming from the land)
        sss[sss<-100] = np.nan
  
        # convert time of the simulation to python format
      
        timepf = sto.simtime2datetime_allsims(date_ori, time_fin)
        #timepi = sto.simtime2datetime_allsims(date_ori, time_ini)

        #print('ini date...',mdates.num2date(timepi).strftime("%Y-%m-%d"))
        print('end date...',mdates.num2date(timepf).strftime("%Y-%m-%d"))


        ''' Search SMAP data for the final day of the simulation '''

        lon_smap, lat_smap, sss_smap, time_smap = \
              sto.sim_open_smap_v4(timepf, dir_SSS, lonmin, lonmax, latmin, latmax,)
      
      
    
        ''' COARSE GRAIN advected field using griddata (revision 1) '''
       
        # reduce grid smap to the region of the reconstruction
       
        ii = np.where(np.logical_and(lon_smap[0,:]>=lon_tf.min(), 
                     lon_smap[0,:]<=lon_tf.max()))
       
        jj =  np.where(np.logical_and(lat_smap[:,0]>=lat_tf.min(), 
                     lat_smap[:,0]<=lat_tf.max()))
       
        lon_smap_red = lon_smap[jj[0], :][:, ii[0]]
        lat_smap_red = lat_smap[jj[0], :][:, ii[0]]
       
        # resampling
        sss_rsp = griddata((lon_tf, lat_tf), sss, 
                          (lon_smap_red, lat_smap_red), method='nearest')
       

        ''' optimal interpolation to smooth the field '''
        
        lx  = 70 #km
        ly  = 70 #km
        ang    = 0
        eps    = 0.03
       
        print('start optimal interpolation to smooth the resampled field... ')
        print('')


        lonikm2d, latikm2d, sss_rsp_s, errint = toi.compute_OI_2d(
                  lon_smap_red.flatten(), lat_smap_red.flatten(), 
                  sss_rsp.flatten(),
                  lon_smap_red[0,:], lat_smap_red[:,0], 
                  lx, ly, ang, eps,'plane')
       
        print('optimal interpolation done... ')
        print('')
        
        

        '''
          Save S advected, resampled and smoothed to an .pkl file
        '''

        dict_smooth = {'lon_smap_red'  : lon_smap_red,
                        'lat_smap_red'  : lat_smap_red,
                        'sss_rsp_s'     : sss_rsp_s}
    
    
        # save dict_smooth 

        f = open(dirsave + 'S_adv_resampled_smoothed_20170212.pkl','wb')
        pickle.dump(dict_smooth,f)
        f.close()      

