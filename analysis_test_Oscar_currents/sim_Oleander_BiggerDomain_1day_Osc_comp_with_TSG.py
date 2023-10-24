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


"""
Code to open the simulation file with SSS tagged, 
(file written in "*_tag_SSS_to_particles.py")
and compare with TSG data. 


Infer the rmsd for Sadv vs. Stsg,
and Ssmap vs. Stsg

adapted to Oleander simulations backward simulation


written by Bàrbara Barceló-Llull on 6-Nov-2019 at Mallorca


QC:
    advected fields: check interpolation where no data
    TSG data: only long transects and remove short transects inside
    long transects as in 20180323 and 20181122
    
> on 10-3-2020:

Backward simulations for the simulations at the Gulf Stream
in a Bigger domain and for 1 day ('20170212') using 
currents from Oscar.
  

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

  # v4 SMAP data        
      
  dir_SSS     = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS/SMAP_L3_RSS_V4_8day_SCI/'
  
  dir_savedic = '/Users/bbarcelo/HOME_SCIENCE/Data/2019_SSS_Ladvection/'+ \
                'sim_Oleander_all/post-processing/'

  name_simulation = 'sim_Oleander_back_Oscar_BD_'
  
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
  date_ori_dt = mdates.date2num(datetime(2017, 1, 1, 0, 0, 0))     

  
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

             
        filename = name_simulation + mdates.num2date(date_fin).strftime("%Y%m%d")+ \
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

  

          
        ''' Make a map for the comparison Sadv vs. Stsg '''
      
        smin = 32#np.nanmin(sss)
        smax = 37#np.nanmax(sss)
    
    
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

  
#      sadvi = interpolate.griddata((lon_fin[~np.isnan(sss)], 
#                                lat_fin[~np.isnan(sss)]), 
#                                sss[~np.isnan(sss)].ravel(),
#                            (lon_tsg, lat_tsg), method='linear')
        #cond_nans = np.logical_or(np.isnan(lon_fin), np.isnan(sss))
        cond_nans = np.isnan(lon_fin)
        sadvi = interpolate.griddata((lon_fin[~cond_nans], 
                                lat_fin[~cond_nans]), 
                                sss[~cond_nans].ravel(),
                            (lon_tsg, lat_tsg), method='linear')
    
      # Map advected field and interpolated at TSG positions
  
#        plt.figure(figsize=(12,6))
#        plt.scatter(lon_fin, lat_fin, c=sss, vmin=32.3, vmax=37, cmap=plt.cm.jet)  
#        plt.scatter(lon_tsg, lat_tsg)
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
  
        '''
        Plot to compare S_adv and S_SMAP (end of simulation) at the TSG position
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
        plt.title('SSS from SMAP on '+ 
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

        
        plt.savefig(dirfig + name_simulation + filename[-22:-14] + filename[-14:-7] + '_comp.png', dpi=200)
        
        ''' Only maps '''
        
        fig = plt.figure(figsize=(12,6)) #(figsize=(12,12))
        gs  = fig.add_gridspec(1, 2)
      
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
        plt.title('SSS from SMAP on '+ 
            mdates.num2date(time_smap).strftime("%Y-%m-%d")+
            ' + TSG data')
        #plt.axis('image')   
        #plt.xlim(lonmin - 1, lonmax+1)
        #plt.ylim(latmin-1, latmax+1)
    
        plt.colorbar(pc)      
      

        
        plt.savefig(dirfig + name_simulation + filename[-22:-14] + filename[-14:-7] + '_maps.png', dpi=200)
        
 