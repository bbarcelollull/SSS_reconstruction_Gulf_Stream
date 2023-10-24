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
from mpl_toolkits.basemap    import Basemap
import sss_toolbox           as to

"""
Code to open the simulation file with SSS tagged, 
(file written in "*_tag_SSS_to_particles.py")
and compare with TSG data. 


Backward simulations for the simulations at the Gulf Stream
in a Bigger domain and using 
currents from altimetry and TSG quality controlled.


> Save advected and SMAP data interpolated onto 
TSG data points for each simulation to detect fronts.


> old FIGURE 1 for paper: maps of SMAP vs. Adv 
with TSG and altimetry currents ('20170212') 


adapted by Bàrbara Barceló-Llull on 13-April-2020 in Mallorca  

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
                'sim_Oleander_all_Bigger_Domain/sss_tagged/'

  dirfig_ppr      = '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/' + \
                'fig_Oleander_BiggerRegion/forpaper/'
  
  dirfig_all      =  '/Users/bbarcelo/HOME_SCIENCE/Figures/2019_SSS_Ladvection/' + \
                'fig_Oleander_BiggerRegion/comparison/'          


  # v4 SMAP data        
      
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  
  dir_savedic = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all_Bigger_Domain/post-processing/'


  ''' Do you want to save a figure for each simulation? '''
  
  save_comp_all = True

  ''' Do you want to save interpolated data? '''
  
  save_interp_data = True  
  
  ''' Set parameters for the code '''

  # running days
  rds = [7]

  # FINAL BIGGER DOMAIN
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
  date_ori_dt  = mdates.date2num(datetime(2015, 3, 31, 0, 0, 0))  
  #date_ori_nrt = mdates.date2num(datetime(2019, 5, 14, 0, 0, 0))   

  date_ori = date_ori_dt  

  
  ''' File with the simulation dates '''  
  
  dir_dic    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'
  #Bigger Domain and quality controlled.
  file_dic = 'Oleander_TSG_per_transect_QCed_BiggerDomain.pkl' 


  f = open(dir_dic + file_dic, 'rb')
  dict_tsg_all = pickle.load(f)
  f.close() 

  dates_strings_all = list(dict_tsg_all.keys())
  print(dates_strings_all)

  dates_strings = dates_strings_all.copy()
  dates_strings.sort(key=int)
  print(dates_strings)
  
  
  ''' Open altimetry file - for figure for paper'''
  
  dir_netcdf    = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SLA_GS_big_region/' 
  # from 2015-03-31 to 2019-05-13
  file_alt_dt  = 'dataset-duacs-rep-global-merged-allsat-phy-l4_1583502357526.nc'

  
  time_dt, lat_dt, lon_dt, adt_dt, ugos_dt, \
        vgos_dt, sla_dt, ugosa_dt, vgosa_dt = to.read_altimetry(dir_netcdf + file_alt_dt)


  lon_alt = lon_dt
  lat_alt = lat_dt
  
  
  ''' Loop for each simulation '''
  
  dic_all_dates    = {}

  date_to_check = ['20170212']
  
  for date in dates_strings:

    date_time_obj = datetime.strptime(date, '%Y%m%d')
    
    print('')
    print ('date... ', date_time_obj)
    
      
    year  = date_time_obj.year
    month = date_time_obj.month
    day   = date_time_obj.day
    
    date_fin = mdates.date2num(datetime(year, month, day, 0, 0, 0))

    ''' Open TSG data for this date (already QC-ed) '''
      
    dict_tsg_date = dict_tsg_all[date]
      
    time_tsg = dict_tsg_date['time']
    lon_tsg  = dict_tsg_date['lon']
    lat_tsg  = dict_tsg_date['lat']
    sal_tsg  = dict_tsg_date['sal']
     
    ''' Open data from simulation '''
         
    dic_all_runtime  = {}

    for rd in rds:

        filename = 'SimOl_B_alt_BD_'+ mdates.num2date(date_fin).strftime("%Y%m%d")+ \
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

        # SSS tagged to each particle
        sss       = nc.variables['sss_adv'][:].data
        nc.close()
      
#      smin = np.nanmin(sss)
#      smax = np.nanmax(sss)
#      # plot the advected field at TF
#      plt.figure()
#      plt.scatter(lon_tf, lat_tf, c=sss, 
#                  vmin=smin, vmax=smax, 
#                cmap=plt.cm.jet)
#      plt.axis('image')
#      plt.xlim(lonmin+ 360 - 2, lonmax+ 360+4)
#      plt.ylim(latmin-3, latmax+4)
#      plt.title('SSS adv at TF running time=' + filename[-9:-3] )  
#      #plt.tight_layout()
#      plt.colorbar()      
#      
      
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

  
        ''' altimetry data for this date (dt) FOR FIGURE 1 PAPER '''
        
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
        
        ''' Make a map for the comparison Sadv vs. Stsg '''
      
        smin = 32#np.nanmin(sss)
        smax = 37#np.nanmax(sss)
    
#    
#      plt.figure()#plt.figure(figsize=(12,6))
#      sc=plt.scatter(lon_fin[:].flatten(), lat_fin[:].flatten(), c=sss[:].flatten(), 
#            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0)#, edgecolors='m')
#      plt.scatter(lon_tsg, lat_tsg, c=sal_tsg, 
#            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0, edgecolors='k')  
#      plt.axis('image')
#      plt.xlim(lonmin, lonmax)
#      plt.ylim(latmin, latmax)
#      plt.title('SSS advected after '+ np.str(rd)+ 
#            ' days of simulation ending on '+
#              mdates.num2date(timepf).strftime("%Y-%m-%d")+
#              ' + TSG data')  
#      plt.colorbar(sc)
#      plt.xlim(lonmin - 2, lonmax+4)
#      plt.ylim(latmin-3, latmax+4)     
#      plt.savefig(dirfig + 'SimOl_BF_' + filename[-28:-12] + filename[-8:-3]+ \
#                '_Sadv_Tf_and_TSG.png', dpi=200)

  
        ''' Interpolate Sadv to Stsg position '''
       
        

        cond_nans = np.isnan(lon_fin)
        sadvi = interpolate.griddata((lon_fin[~cond_nans], 
                                lat_fin[~cond_nans]), 
                                sss[~cond_nans].ravel(),
                            (lon_tsg, lat_tsg), method='linear')
       
      # Map advected field and interpolated at TSG positions
  
#        plt.figure(figsize=(12,6))
#        plt.scatter(lon_fin[~cond_nans], lat_fin[~cond_nans], s=2, c='m')
#        plt.scatter(lon_fin[~cond_nans], lat_fin[~cond_nans],
#                    c=sss[~cond_nans], s=1, cmap=plt.cm.jet)
#        plt.scatter(lon_tsg, lat_tsg, c=sadvi, vmin=32.3, vmax=37, cmap=plt.cm.jet, 
#              linewidths=0.05, edgecolors='w')        
#        
#        plt.figure(figsize=(12,6))
#        plt.scatter(lon_fin, lat_fin, c=sss, vmin=32.3, vmax=37, cmap=plt.cm.jet)  
#        plt.scatter(lon_tsg, lat_tsg, marker='x')
#        plt.scatter(lon_tsg, lat_tsg, c=sadvi, vmin=32.3, vmax=37, cmap=plt.cm.jet, 
#              linewidths=0.05, edgecolors='w')
#        
#        plt.title('SSS advected at TSG position')
#        plt.colorbar()
#      
#        plt.figure(figsize=(12,6))
#        plt.plot(lon_tsg, sadvi, 'ob')
#        plt.plot(lon_tsg, sal_tsg, 'xr')
#        aaa
      
 
        ''' Search SMAP data for the final day of the simulation '''

        lon_smap, lat_smap, sss_smap, time_smap = \
              sto.sim_open_smap_v4(timepf, dir_SSS, lonmin, lonmax, latmin, latmax,)
      
       
#      plt.figure()
#      pc = plt.pcolor(lon_smap[0,:], lat_smap[:,0],sss_smap, 
#                      vmin=smin, vmax=smax, cmap=plt.cm.jet)
#      plt.axis('image')
#      plt.xlim(lonmin - 2, lonmax+4)
#      plt.ylim(latmin-3, latmax+4)
#      plt.colorbar()
#      plt.title('SSS SMAP on ' + mdates.num2date(time_smap).strftime("%Y-%m-%d"))  
#      
  
        ''' Interpolate SSS from SMAP to TSG positions '''
  
        ssmapfi_nonan = interpolate.griddata((lon_smap.ravel(), lat_smap.ravel()), sss_smap.ravel(),
                            (lon_tsg, lat_tsg), method='linear')
  
        ssmapfi = np.copy(ssmapfi_nonan)
        ssmapfi[ssmapfi<20] = np.nan
      
      # check if still out
#      plt.figure()
#      plt.plot(lon_tsg, ssmapfi_nonan, 'xr')
#      plt.plot(lon_tsg, ssmapfi, 'o-y', 
#           markersize=4,linewidth=1, 
#           label=r'S$_{SMAPf}$ ' + '(rmsd = %0.2f)'% rms_smap) 
      
      
      # Map SSS SMAP field and interpolated
  
#      plt.figure()
#      pc = plt.pcolor(lon_smap[0,:], lat_smap[:,0],sss_smap, 
#                    vmin=smin, vmax=smax, cmap=plt.cm.jet)  
#      plt.scatter(lon_tsg, lat_tsg, c=sal_tsg, 
#            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0, edgecolors='k')  
#      plt.title('SSS from SMAP on '+ 
#            mdates.num2date(time_smap).strftime("%Y-%m-%d")+
#            ' + TSG data')
#      plt.axis('image')   
#      plt.xlim(lonmin - 2, lonmax+4)
#      plt.ylim(latmin-3, latmax+4)
#    
#      plt.colorbar(pc)
#      plt.savefig(dirfig + 'SimOl_BF_' + filename[-28:-12] + filename[-8:-3]+
#                '_SSMAP_Tf_and_TSG.png', dpi=200)  
#      
  
        ''' ------- Comparison S_adv with TSG and SMAP salinity data ------- '''
  
        '''
        Infer statistics between (sadvi, sal_tsg) and (sadvi, ssmapfi)
        correlation coefficient, rms, std
        '''

#        corr_adv  = np.corrcoef(sadvi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])[1,0] 
#        corr_smap = np.corrcoef(ssmapfi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])[1,0] 
#
#        rms_adv  = sto.rmsd(sadvi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])
#        rms_smap = sto.rmsd(ssmapfi[~np.isnan(sadvi)], sal_tsg[~np.isnan(sadvi)])
# 
#        std_adv  = sto.std(sadvi[~np.isnan(sadvi)] - sal_tsg[~np.isnan(sadvi)])
#        std_smap = sto.std(ssmapfi[~np.isnan(sadvi)]- sal_tsg[~np.isnan(sadvi)])    
#  
#        print('')
#        print('Correlation Sadv and Stsg...', corr_adv)
#        print('Correlation Ssmap and Stsg...', corr_smap)
#        print('')
#
#        print('rms Sadv and Stsg...', rms_adv)
#        print('rms Ssmap and Stsg...', rms_smap)
#        print('')  
#
#        print('std Sadv and Stsg...', std_adv)
#        print('std Ssmap and Stsg...', std_smap)
#        print('')  
        
        if save_comp_all == True:
            
          '''
          Figure: maps S_adv vs. S_SMAP (end of simulation) 
          and plot to compare fields at TSG position.
          '''
        
          xfin, yfin = bm(lon_fin[:].flatten(), lat_fin[:].flatten())
          xtsg, ytsg = bm(lon_tsg, lat_tsg)
          xsmap, ysmap = bm(lon_smap, lat_smap)
        
          fsize = 14
  
          fig = plt.figure(figsize=(10,10)) #(figsize=(12,12))
          gs  = fig.add_gridspec(2, 2)
      
          ax1 = fig.add_subplot(gs[0, 0])
          plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
          sc=plt.scatter(xfin, yfin, c=sss[:].flatten(), 
            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0)#, edgecolors='m')
          plt.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0, edgecolors='k')  
          #plt.axis('image')
          #plt.xlim(lonmin, lonmax)
          #plt.ylim(latmin, latmax)
          plt.title('S$_{adv}$ '+ mdates.num2date(timepf).strftime("%Y-%m-%d")+ ' '\
                + np.str(rd)+  ' days '+ '+ TSG data')  
          plt.colorbar(sc)
          #plt.xlim(lonmin - 1, lonmax+1)
          #plt.ylim(latmin-1, latmax+1)

          ax2 = fig.add_subplot(gs[0, 1])
          plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
          pc = plt.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
                    vmin=smin, vmax=smax, cmap=plt.cm.jet)  
          plt.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.jet, linewidths=0, edgecolors='k')  
          plt.title('S$_{SMAP}$ SMAP on '+ 
            mdates.num2date(time_smap).strftime("%Y-%m-%d")+
            ' + TSG data')
        #plt.axis('image')   
        #plt.xlim(lonmin - 1, lonmax+1)
        #plt.ylim(latmin-1, latmax+1)
    
          plt.colorbar(pc)      
      
          ax3 = fig.add_subplot(gs[1, :])
          plt.plot(lon_tsg, sal_tsg, 'o-k', 
                   markersize=4,linewidth=1,label=r'S$_{tsg}$')


          plt.plot(lon_tsg, ssmapfi, 'o-y', 
             markersize=4,linewidth=1, 
             label=r'S$_{SMAPf}$')  
      
          plt.plot(lon_tsg, sadvi, 'o-', color='steelblue', 
             markersize=4,linewidth=1, 
             label=r'S$_{adv}$')
           
#        plt.plot(lon_tsg, ssmapfi, 'o-y', 
#           markersize=4,linewidth=1, 
#           label=r'S$_{SMAPf}$ ' + '(rmsd = %0.2f)'% rms_smap)  
#      
#        plt.plot(lon_tsg, sadvi, 'o-', color='steelblue', 
#           markersize=4,linewidth=1, 
#           label=r'S$_{adv}$ ' + '(rmsd = %0.2f)'% rms_adv)
      
          plt.legend(loc=4,prop={'size': fsize})
          plt.title('Does the simulation improve the SSS field?' + 
              ' (simulation: ' + filename[-28:-3] + ')', fontsize=fsize)
          plt.xlabel('Longitude', fontsize=fsize)
          plt.ylabel('SSS', fontsize=fsize)
        
          plt.tick_params(axis='both', labelsize=fsize)
          plt.tight_layout()
          plt.xlim(lon_tsg.min()-0.5, lon_tsg.max()+0.5)

        
          plt.savefig(dirfig_all + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_comp.png', dpi=200)
        
          
          
        if date == '20170212': 
            
          ''' Figure 1 paper: maps on '20170212' '''
        
        
          fsize = 12
          ss = 1
          xx_alt, yy_alt   = bm(lon_alt, lat_alt)
        
          fig = plt.figure(figsize=(10,6)) #(figsize=(12,12))
          gs  = fig.add_gridspec(1, 3, width_ratios=[20,20,1])
      
        
          ax1 = fig.add_subplot(gs[0, 0])
          plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
          pc = plt.pcolor(xsmap[0,:], ysmap[:,0],sss_smap, 
                    vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r,zorder=1)  
          plt.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0,
            edgecolors='k',zorder=3)  
          qv1=ax1.quiver(xx_alt[::ss, ::ss], yy_alt[::ss, ::ss], 
                           ug_date_fin[::ss, ::ss], vg_date_fin[::ss, ::ss], 
                           color='steelblue', scale=15,alpha=0.5, zorder=2)
            
          #ax1.quiverkey(qv1, 0.89, 0.03, 1, '1 m/s', color='k', zorder=200)                
          ax1.quiverkey(qv1, 0.07, 0.945, 1, '1 m/s', 
                      coordinates='figure', color='k', alpha=1)                  
            

          plt.title('S$_{SMAP}$ on '+ 
            mdates.num2date(time_smap).strftime("%d-%b-%Y"), fontsize=fsize)
   
        
          ax2 = fig.add_subplot(gs[0, 1])
          plot_decor_map(bm,lonmin, lonmax, latmin, latmax, fsize)
          sc=plt.scatter(xfin, yfin, c=sss[:].flatten(), 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0)#, edgecolors='m')
          plt.scatter(xtsg, ytsg, c=sal_tsg, 
            vmin=smin, vmax=smax, cmap=plt.cm.Spectral_r, linewidths=0, edgecolors='k')  
    

          plt.title('S$_{adv}$ on '+ \
                  mdates.num2date(timepf).strftime("%d-%b-%Y"),fontsize=fsize) #"%Y-%m-%d") )  

          axcb = fig.add_subplot(gs[0, 2])
          cb = plt.colorbar(pc, cax=axcb)  
          cb.ax.tick_params(labelsize=fsize-1) 
          cb.ax.set_title('[PSU]', fontsize=fsize-1, pad=10)
        
          fig.text(0.01, 0.97, '(a)' , fontsize = fsize+1, zorder=200)
          fig.text(0.45, 0.97, '(b)' , fontsize = fsize+1, zorder=200)
        
          plt.tight_layout()
        
          plt.savefig(dirfig_ppr + 'sim_Oleander_back_alt_BD_' + filename[-22:-14] + filename[-14:-7] + '_maps.png', dpi=200)
        
          

      
        ''' Save data for this transect in a dictionary '''
        # dictionary of the data of this transect:
        dic_tran = {'lon_tsg'  : lon_tsg,
                  'lat_tsg'  : lat_tsg,
                  'sal_tsg'  : sal_tsg,
                  'time_tsg' : time_tsg,
                  's_advi'   : sadvi,
                  's_smapi'  : ssmapfi}
      
        # add to dictionary for all running times
        dic_all_runtime.update({'%02i'% rd: dic_tran})
      
      
        # add to dictionary with all TSG transects -> [transect]['07'][var]
        dic_all_dates.update({date: dic_all_runtime}) 
      
        plt.close('all') 
    
  
  if save_interp_data == True:
      
    ''' Save Sadv and Ssmap on TSG position to detect fronts '''
  
    f = open(dir_savedic + 'SimOl_B_alt_BD_sal_int_to_TSG_final.pkl','wb')
    pickle.dump(dic_all_dates,f)
    f.close()  