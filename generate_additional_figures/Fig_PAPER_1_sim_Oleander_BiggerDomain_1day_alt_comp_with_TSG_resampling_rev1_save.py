#usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
import pickle
from matplotlib              import dates as mdates
from datetime                import datetime, timedelta
from scipy                   import interpolate
import sim_toolbox           as sto
from mpl_toolkits.basemap    import Basemap
import sss_toolbox           as to
from scipy.interpolate       import griddata
import Tools_OI              as toi 
from scipy.ndimage           import gaussian_filter
import pickle


"""
Code to open the simulation file with SSS tagged, 
(file written in "*_tag_SSS_to_particles.py")
and compare with TSG data. 


Infer the rmsd for Sadv vs. Stsg,
and Ssmap vs. Stsg

adapted to Oleander simulations backward simulation


written by Bàrbara Barceló-Llull on 6-Nov-2019 at Mallorca

> on 21-Nov-2019:
    
changed for the massive simulations, only the simulations where
the delayed-time altimetry currents where used

simulations from: sim_Oleander_backward_altimetry_all

QC:
    advected fields: check interpolation where no data
    TSG data: only long transects and remove short transects inside
    long transects as in 20180323 and 20181122
    
> on 9-3-2020:

Backward simulations for the simulations at the Gulf Stream
in a Bigger domain and for 1 day ('20170212') using 
currents from altimetry.
  

> 14 August 2020:
    
    revision 1: coarse grain the reconstructed field to 
    the original resolution
    
    resample and smooth S_adv and save data in a file:
        
        S_adv_resampled_smoothed_20170212.pkl
        
        in /Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/sim_Oleander_all/sss_tagged/
        


"""
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


if __name__ == '__main__':
    
  plt.close('all')
  
  dirsave     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all/sss_tagged/'

  dirfig      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/' + \
                'fig_Oleander_BiggerRegion/'
  dirfigppr     = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/fig_paper/'


  # v4 SMAP data        
      
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  
  dir_savedic = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all/post-processing/'

  
  ''' Set parameters for the code '''

  # running days
  rds = [7]

  # TRYING A BIGGER DOMAIN
  lonmin, lonmax = -82, -63 
  latmin, latmax =  25, 46

  
  #Create the basemap
  
  bm = Basemap(projection = 'merc',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')   
      
  # first day of the simulation in the GS 
  # when using delayed time altimetry data for the massive simulations
  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0))  
  date_ori_nrt = mdates.date2num(datetime(2019, 5, 14, 0, 0, 0))   

  
  i_change     = 132 # '20190110' last date of dt 
  
  ''' File with the simulation dates '''  
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'
  #file_dic = 'Oleander_TSG_per_transect_all.pkl'
  file_dic = 'Oleander_TSG_per_transect_all_BiggerDomain.pkl'


  f = open(dir_dic + file_dic, 'rb')
  dict_tsg_all = pickle.load(f)
  f.close() 
#
#  dates_strings_all = list(dict_tsg_all.keys())
#  print(dates_strings_all)
#
#  #after 2019-01-12 we don't have delayed-time altimetry data
#  dates_strings = dates_strings_all.copy()
#  dates_strings.sort(key=int)
#  print(dates_strings)
  

  ''' Open altimetry file - for figure for paper'''
  
  dir_netcdf    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SLA_GS_big_region/' 
  # from 2015-03-31 to 2019-05-13
  file_alt_dt  = 'dataset-duacs-rep-global-merged-allsat-phy-l4_1583502357526.nc'
  # from 2019-05-14 to 2020-01-15
  file_alt_nrt = 'dataset-duacs-nrt-global-merged-allsat-phy-l4_1583920119159.nc'
  
  
  time_dt, lat_dt, lon_dt, adt_dt, ugos_dt, \
        vgos_dt, sla_dt, ugosa_dt, vgosa_dt = to.read_altimetry(dir_netcdf + file_alt_dt)

  time_nrt, lat_nrt, lon_nrt, adt_nrt, ugos_nrt, \
        vgos_nrt, sla_nrt, ugosa_nrt, vgosa_nrt = to.read_altimetry(dir_netcdf + file_alt_nrt) 
   
  lon_alt = lon_dt
  lat_alt = lat_dt
  
  
  ''' Loop for each simulation '''
  
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

    

    ''' Open TSG data for this date'''
      
    dict_tsg_date = dict_tsg_all[date]
      
    time_tsg = dict_tsg_date['time']
    lon_tsg  = dict_tsg_date['lon']
    lat_tsg  = dict_tsg_date['lat']
    sal_tsg  = dict_tsg_date['sal']
      
    if len(sal_tsg[~np.isnan(sal_tsg)]) == 0:
        
        print('empty transect')
        n_short = n_short + 1  
        
        
    elif len(sal_tsg[~np.isnan(sal_tsg)]) > 0:
        
     tot_dist = lon_tsg[~np.isnan(sal_tsg)].max()-lon_tsg[~np.isnan(sal_tsg)].min()

     if tot_dist < 1.5: 
          
        print('profile too short!')
        print('')
        
        n_short = n_short + 1  
    
     else:
        
      # particular QC for some transects:
      if date == '20161001':

            time_tsg = time_tsg[1:]
            lon_tsg  = lon_tsg[1:]
            lat_tsg  = lat_tsg[1:]
            sal_tsg  = sal_tsg[1:]


      elif date == '20160928':
            time_tsg = time_tsg[1:]
            lon_tsg  = lon_tsg[1:]
            lat_tsg  = lat_tsg[1:]
            sal_tsg  = sal_tsg[1:]            

      elif date == '20161005':            
            bad_ind = 2880
            
            time_tsg = np.append(time_tsg[:bad_ind], time_tsg[bad_ind+1:])
            lon_tsg  = np.append(lon_tsg[:bad_ind], lon_tsg[bad_ind+1:])
            lat_tsg  = np.append(lat_tsg[:bad_ind], lat_tsg[bad_ind+1:])
            sal_tsg  = np.append(sal_tsg[:bad_ind], sal_tsg[bad_ind+1:])
            
      elif date == '20180323':
            bad_ind = 40
            
            time_tsg = time_tsg[bad_ind+1:]
            lon_tsg  = lon_tsg[bad_ind+1:]
            lat_tsg  = lat_tsg[bad_ind+1:]
            sal_tsg  = sal_tsg[bad_ind+1:]
  
            
      elif date == '20181122':   
            bad_ind = 1085

            time_tsg = time_tsg[bad_ind+1:]
            lon_tsg  = lon_tsg[bad_ind+1:]
            lat_tsg  = lat_tsg[bad_ind+1:]
            sal_tsg  = sal_tsg[bad_ind+1:]

      elif date == '20180125':   

            time_tsg = time_tsg[:-1]
            lon_tsg  = lon_tsg[:-1]
            lat_tsg  = lat_tsg[:-1]
            sal_tsg  = sal_tsg[:-1]
            
      dic_all_runtime  = {}

      for rd in rds:

             
        filename = 'sim_Oleander_back_alt_BD_'+ mdates.num2date(date_fin).strftime("%Y%m%d")+ \
                  '_' + '%02i'% rd + 'days'+'_sss.nc'
    
        print('')
        print('Openning...' + filename)
      
      
        ''' Open simulation data '''
      
#      traj_ini, lat_ini, lon_ini, time_ini, \
#      traj_fin, lat_fin, lon_fin, time_fin, sss = sto.open_back_sim_file(dirsave, filename)


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

  
        ''' altimetry data for this date (dt) '''
        
        ind_alt_fin = np.where(time_dt == timepf)
        #ind_alt_ini = np.where(time_dt == timepi)

        #sla_date_ini = sla_dt[ind_alt_ini].squeeze()
        #adt_date_ini = adt_dt[ind_alt_ini].squeeze()
        #ug_date_ini  = ugos_dt[ind_alt_ini].squeeze()
        #vg_date_ini  = vgos_dt[ind_alt_ini].squeeze()
                
        sla_date_fin = sla_dt[ind_alt_fin].squeeze()
        adt_date_fin = adt_dt[ind_alt_fin].squeeze()
        ug_date_fin  = ugos_dt[ind_alt_fin].squeeze()
        vg_date_fin  = vgos_dt[ind_alt_fin].squeeze()


  
        ''' Interpolate Sadv to Stsg position '''

  
        cond_nans = np.isnan(lon_fin)
        sadvi = interpolate.griddata((lon_fin[~cond_nans], 
                                lat_fin[~cond_nans]), 
                                sss[~cond_nans].ravel(),
                            (lon_tsg, lat_tsg), method='linear')
    

  

 
        ''' Search SMAP data for the final day of the simulation '''

        lon_smap, lat_smap, sss_smap, time_smap = \
              sto.sim_open_smap_v4(timepf, dir_SSS, lonmin, lonmax, latmin, latmax,)
      
       

  
        ''' Interpolate SSS from SMAP to TSG positions '''
  
        ssmapfi_nonan = interpolate.griddata((lon_smap.ravel(), lat_smap.ravel()), sss_smap.ravel(),
                            (lon_tsg, lat_tsg), method='linear')
  
        ssmapfi = np.copy(ssmapfi_nonan)
        ssmapfi[ssmapfi<20] = np.nan
      
  
  
        ''' ------- Comparison S_adv with TSG and SMAP salinity data ------- '''
  
        '''
        Infer statistics between (sadvi, sal_tsg) and (sadvi, ssmapfi)
        correlation coefficient, rms, std
        '''

        corr_adv  = np.corrcoef(sadvi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])[1,0] 
        corr_smap = np.corrcoef(ssmapfi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])[1,0] 

        rms_adv  = sto.rmsd(sadvi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])
        rms_smap = sto.rmsd(ssmapfi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])
 
        std_adv  = sto.std(sadvi[~np.isnan(sadvi)] - sal_tsg[~np.isnan(sadvi)])
        std_smap = sto.std(ssmapfi[~np.isnan(sadvi)]- sal_tsg[~np.isnan(sadvi)])    
  
        print('')
        print('Correlation Sadv and Stsg...', corr_adv)
        print('Correlation Ssmap and Stsg...', corr_smap)
        print('')

        print('rms Sadv and Stsg...', rms_adv)
        print('rms Ssmap and Stsg...', rms_smap)
        print('')  

        print('std Sadv and Stsg...', std_adv)
        print('std Ssmap and Stsg...', std_smap)
        print('')  
        
        ''' limits of salinity for plots '''
        smin = 32#np.nanmin(sss)
        smax = 37#np.nanmax(sss)
        
        
 

    
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
       
       
        ''' Figure S_SMAP and S_adv_resampled '''
       
        fig = plt.figure(figsize=(10,4)) #(figsize=(12,12))
        gs  = fig.add_gridspec(1, 4, width_ratios=[20,20,20,1])
      
        
        ax1 = fig.add_subplot(gs[0, 0])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
        
        ax2 = fig.add_subplot(gs[0, 1])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)        
        
        ax3 = fig.add_subplot(gs[0, 2])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)          
         
    
        # reconstruction
        
        sc=ax1.scatter(xfin, yfin, c=sss[:].flatten(), 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0)#, edgecolors='m')

        ax1.set_title('S$_{adv}$', fontsize=fsize) 


        # reconstruction resampled
        xrsp, yrsp = bm(lon_smap_red, lat_smap_red)
        
        pc = ax2.pcolor(xrsp[0,:], yrsp[:,0],sss_rsp, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
        
        ax2.set_title('S$_{adv}$ resampled', fontsize=fsize) 
        
        # smap
        pc = ax3.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)               

        ax3.set_title('S$_{SMAP}$', fontsize=fsize)
 
        axcb = fig.add_subplot(gs[0, 3])
        cb = plt.colorbar(pc, cax=axcb)  
        cb.ax.tick_params(labelsize=fsize-1) 
        cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
        
        fig.text(0.01, 0.95, '(a)' , fontsize = fsize+1, zorder=200)
        fig.text(0.35, 0.95, '(b)' , fontsize = fsize+1, zorder=200)
        fig.text(0.65, 0.95, '(c)' , fontsize = fsize+1, zorder=200)
       
        plt.tight_layout()
        
        plt.savefig(dirfigppr + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_maps_resampling.png', dpi=200)
               
        
        
       
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
          Save S advected, resampled and smoothed to an .nc file
        '''

        dict_smooth = {'lon_smap_red'  : lon_smap_red,
                        'lat_smap_red'  : lat_smap_red,
                        'sss_rsp_s'     : sss_rsp_s}
    
    
        # save dict_smooth 

        f = open(dirsave + 'S_adv_resampled_smoothed_20170212.pkl','wb')
        pickle.dump(dict_smooth,f)
        f.close()      


        
        # ''' try gaussian smoother '''
        
        
        # sss_rsp_s2 = gaussian_filter(sss_rsp, sigma=0.5)
        
        
        # ''' Figure S_SMAP and S_adv_resampled and smoothed 2 '''
       
        # fig = plt.figure(figsize=(7,8)) #(figsize=(12,12))
        # gs  = fig.add_gridspec(2, 3, width_ratios=[20,20,1])
      
        
        # ax1 = fig.add_subplot(gs[0, 0])
        # plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
        
        # ax2 = fig.add_subplot(gs[0, 1])
        # plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)        
        
        # ax3 = fig.add_subplot(gs[1, 0])
        # plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)          
         
        # ax4 = fig.add_subplot(gs[1, 1])
        # plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)   
        
        # axcb = fig.add_subplot(gs[:, 2])
        
        # # reconstruction
        
        # sc=ax1.scatter(xfin, yfin, c=sss[:].flatten(), 
        #     vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0)#, edgecolors='m')

        # ax1.set_title('S$_{adv}$', fontsize=fsize) 


        # # reconstruction resampled
        # xrsp, yrsp = bm(lon_smap_red, lat_smap_red)
        
        # pc = ax2.pcolor(xrsp[0,:], yrsp[:,0],sss_rsp, 
        #             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
        
        # ax2.set_title('S$_{adv}$ resampled', fontsize=fsize) 
        
        
        # # reconstruction resampled and smoothed

        # pc = ax3.pcolor(xrsp[0,:], yrsp[:,0],sss_rsp_s2, 
        #             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
        
        
        # ax3.set_title('S$_{adv}$ resampled and smoothed', fontsize=fsize) 
        
        # # smap
        # pc = ax4.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
        #             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)               

        # ax4.set_title('S$_{SMAP}$', fontsize=fsize)
 
    
        
        # cb = plt.colorbar(pc, cax=axcb)  
        # cb.ax.tick_params(labelsize=fsize-1) 
        # cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
        
        # fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
        # fig.text(0.45, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
        # fig.text(0.01, 0.49, '(c)' , fontsize = fsize+1, zorder=200)
        # fig.text(0.45, 0.49, '(d)' , fontsize = fsize+1, zorder=200)
       
        # plt.tight_layout()
        
        # plt.savefig(dirfigppr + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_maps_resamp_smoothed2.png', dpi=200)
                    
