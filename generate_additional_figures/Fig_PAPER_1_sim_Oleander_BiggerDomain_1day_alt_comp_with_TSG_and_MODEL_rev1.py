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

"""
Code to open the simulation file with SSS tagged, 
(file written in "*_tag_SSS_to_particles.py")
and compare with TSG data. 

simulations from: sim_Oleander_backward_altimetry_all

QC:
    advected fields: check interpolation where no data
    TSG data: only long transects and remove short transects inside
    long transects as in 20180323 and 20181122
    
> on 9-3-2020:

Backward simulations for the simulations at the Gulf Stream
in a Bigger domain and for 1 day ('20170212') using 
currents from altimetry.
  

> on 18 August 2020:
    
    Add comparison with CMEMS (Mercator) model for revision 1

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
  

  ''' Open CMEMS Model data'''
  
  dir_model    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/CMEMS_global_reanalysis/' 
  # reanalysis from 2015-04-07 to 2018-12-25
  file_model  = 'global-reanalysis-phy-001-030-daily_1597743581842.nc'


  timem, lonm, latm, depm, psm, ptm, um, vm = \
                read_cmems_model_surf(dir_model, file_model)
                
  
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

    ''' Open model data for this date '''
    
    ind_mod = np.where(timem == date_fin+0.5)
    
    print(' ')
    print('model data for day... ', mdates.num2date(timem[ind_mod]))
    print(' ')
     
    sal_model = psm[ind_mod].squeeze()
    u_model   = um[ind_mod].squeeze()
    v_model   = vm[ind_mod].squeeze()
    
    # lonm, latm

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

                
        sla_date_fin = sla_dt[ind_alt_fin].squeeze()
        adt_date_fin = adt_dt[ind_alt_fin].squeeze()
        ug_date_fin  = ugos_dt[ind_alt_fin].squeeze()
        vg_date_fin  = vgos_dt[ind_alt_fin].squeeze()
        
        ''' Make a map for the comparison Sadv vs. Stsg '''
      
        smin = 32#np.nanmin(sss)
        smax = 37#np.nanmax(sss)
    
    
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

        
        ''' Only maps, comparison S_adv with S_model '''
        
        
        fsize = 13
        ss = 1
        
        xfin, yfin = bm(lon_fin[:].flatten(), lat_fin[:].flatten())
        xtsg, ytsg = bm(lon_tsg, lat_tsg)
        xsmap, ysmap = bm(lon_smap, lat_smap)
        

        xx_alt, yy_alt   = bm(lon_alt, lat_alt)
        x_model, y_model = bm(lonm, latm)
        
        xgm, ygm         = bm(-73.1, 44)
        xmab, ymab       = bm(-77.8, 40.2)
        xch, ych         = bm(-79.2, 35.7)
        
        
        # Start figure
        fig = plt.figure(figsize=(15,6)) #(figsize=(12,12))
        gs  = fig.add_gridspec(1, 4, width_ratios=[20, 20,20,1])

        # S_smap
        ax1 = fig.add_subplot(gs[0, 0])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
        pc = ax1.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
        ax1.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0,
            edgecolors='k',zorder=3)  
        # qv1=ax1.quiver(xx_alt[::ss, ::ss], yy_alt[::ss, ::ss], 
        #                     ug_date_fin[::ss, ::ss], vg_date_fin[::ss, ::ss], 
        #                     color='steelblue', scale=15,alpha=0.5, zorder=2)
        # ax1.quiverkey(qv1, 0.07, 0.93, 1, '1 m/s', 
        #               coordinates='figure', color='k', alpha=1)                  
            

        ax1.set_title('S$_{SMAP}$', fontsize=fsize)
   
        # region labels 
        bbox_props = dict(boxstyle="rarrow", fc='w', ec="w", lw=0.5, alpha=0.8)
        ax1.text(xgm, ygm, 'Gulf of Maine', ha="center", va="center", 
            bbox=bbox_props, fontsize=fsize-4, color='b',zorder=100)    
        ax1.text(xmab, ymab, 'Middle Atl. Bight', ha="center", va="center", 
            bbox=bbox_props, fontsize=fsize-4, color='b',zorder=100)    
        ax1.text(xch, ych, 'Cape Hatteras', ha="center", va="center", 
            bbox=bbox_props, fontsize=fsize-4, color='b',zorder=100) 

        # S_adv subplot 
        
        ax2 = fig.add_subplot(gs[0, 1])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
        sc=ax2.scatter(xfin, yfin, c=sss[:].flatten(), 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0)#, edgecolors='m')
        ax2.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0, edgecolors='k')  

        ax2.set_title('S$_{adv}$', fontsize=fsize) 

        # S_model subplot
        
        ax3 = fig.add_subplot(gs[0, 2])
        plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
        
        pc = ax3.pcolor(x_model[0,:], y_model[:,0],sal_model, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
        ax3.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0,
            edgecolors='k',zorder=3)  
        # qv1=ax1.quiver(xx_alt[::ss, ::ss], yy_alt[::ss, ::ss], 
        #                    ug_date_fin[::ss, ::ss], vg_date_fin[::ss, ::ss], 
        #                    color='steelblue', scale=15,alpha=0.5, zorder=2)
            
        # #ax1.quiverkey(qv1, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)                
        # ax1.quiverkey(qv1, 0.07, 0.945, 1, '1 m/s', 
        #               coordinates='figure', color='k', alpha=1)                  
            

        ax3.set_title('S$_{model}$', fontsize=fsize)
   

        axcb = fig.add_subplot(gs[0, 3])
        cb = plt.colorbar(pc, cax=axcb)  
        cb.ax.tick_params(labelsize=fsize-1) 
        cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
        
        fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
        fig.text(0.315, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
        fig.text(0.62, 0.97, '(c)' , fontsize = fsize+1, zorder=200)
        
        plt.tight_layout()
        
        plt.savefig(dirfigppr + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_MODEL.png', dpi=200)
        
    
    
    

#         ''' Only maps, adapted for the PAPER!!! '''
        
        
#         fsize = 11
#         ss = 1
#         xx_alt, yy_alt   = bm(lon_alt, lat_alt)
#         xgm, ygm         = bm(-73.1, 44)
#         xmab, ymab       = bm(-77.8, 40.2)
#         xch, ych         = bm(-79.2, 35.7)
        
#         fig = plt.figure(figsize=(10,6)) #(figsize=(12,12))
#         gs  = fig.add_gridspec(1, 3, width_ratios=[20,20,1])
      
        
#         ax1 = fig.add_subplot(gs[0, 0])
#         plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
#         pc = plt.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
#                     vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
#         plt.scatter(xtsg, ytsg, c=sal_tsg, 
#             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0,
#             edgecolors='k',zorder=3)  
#         qv1=ax1.quiver(xx_alt[::ss, ::ss], yy_alt[::ss, ::ss], 
#                            ug_date_fin[::ss, ::ss], vg_date_fin[::ss, ::ss], 
#                            color='steelblue', scale=15,alpha=0.5, zorder=2)
            
#         #ax1.quiverkey(qv1, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)                
#         ax1.quiverkey(qv1, 0.07, 0.945, 1, '1 m/s', 
#                       coordinates='figure', color='k', alpha=1)                  
            

#         plt.title('S$_{SMAP}$', fontsize=fsize)
   
  
#         # region labels 
#         bbox_props = dict(boxstyle="rarrow", fc='w', ec="w", lw=0.5, alpha=0.8)
#         ax1.text(xgm, ygm, 'Gulf of Maine', ha="center", va="center", 
#             bbox=bbox_props, fontsize=fsize-2, color='b',zorder=100)    
#         ax1.text(xmab, ymab, 'Middle Atl. Bight', ha="center", va="center", 
#             bbox=bbox_props, fontsize=fsize-2, color='b',zorder=100)    
#         ax1.text(xch, ych, 'Cape Hatteras', ha="center", va="center", 
#             bbox=bbox_props, fontsize=fsize-2, color='b',zorder=100)  
        
#         ax2 = fig.add_subplot(gs[0, 1])
#         plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
#         sc=plt.scatter(xfin, yfin, c=sss[:].flatten(), 
#             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0)#, edgecolors='m')
#         plt.scatter(xtsg, ytsg, c=sal_tsg, 
#             vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0, edgecolors='k')  

# #        qv2=ax2.quiver(xx_alt, yy_alt, 
# #                           ug_date_fin, vg_date_fin, scale=15, color='b')
# #            
# #        ax2.quiverkey(qv2, 0.05, 0.935, 1, '1 m/s', coordinates='figure', color='k')                  
# #            


#         plt.title('S$_{adv}$', fontsize=fsize) 

#         axcb = fig.add_subplot(gs[0, 2])
#         cb = plt.colorbar(pc, cax=axcb)  
#         cb.ax.tick_params(labelsize=fsize-1) 
#         cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
        
#         fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
#         fig.text(0.45, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
        
#         plt.tight_layout()
        
#         plt.savefig(dirfigppr + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_maps.png', dpi=200)
        
    