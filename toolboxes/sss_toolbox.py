#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Toolbox for the SSS analysis 

written by Bàrbara Barceló Llull on 06-05-2019 at APL-UW
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import netCDF4               as netcdf
from datetime                import datetime
from matplotlib              import dates as mdates
from mpl_toolkits.basemap    import Basemap
from matplotlib.colors       import LinearSegmentedColormap


def read_nc_sss(dirnc, file):
   
   nc    = netcdf.Dataset(dirnc + file, 'r')
   time  = nc.variables['time'][:]  # (time)   
   lat   = nc.variables['lat'][:]   # (lat, lon)
   lon   = nc.variables['lon'][:]   # (lat, lon)
   sss   = nc.variables['sss'][:]   #(time, lat, lon)
   nc.close()

   # nan invalid values
   cond_invalid = np.logical_or(sss<28, sss>39)
   sss[cond_invalid] = np.nan
   
   return time, lat, lon, sss

def read_nc_OISST(dirnc, file):
   
   nc        = netcdf.Dataset(dirnc + file, 'r')
   time_orig = nc.variables['time'][:]        # (time) seconds since 1970-01-01T00:00:00Z
   lat1d     = nc.variables['latitude'][:]    # (lat)
   lon360    = nc.variables['longitude'][:]   # (lon) from 0.125 to 359.875 by 0.25
   sst       = nc.variables['sst'][:].squeeze() # (time, depth, latitude, longitude)
   err       = nc.variables['err'][:].squeeze()
   anom      = nc.variables['anom'][:].squeeze()
   nc.close()
   
   
   # lon from "0 to 360" format to "-180 to 180" format
   lon1d = np.copy(lon360)
   lon1d[lon1d>180] = lon1d[lon1d>180] - 360
   
   
   # lon and lat 2d array
   lon, lat = np.meshgrid(lon1d, lat1d)
   
   
   # time in python format
   # first date
   d_ini_or    = mdates.date2num(datetime(1970, 1, 1, 0, 0, 0))  
  
   # Return seconds as days
   time_days = np.zeros(time_orig.shape)
  
   for ind, tt in enumerate(time_orig):
      time_days[ind] = mdates.seconds(tt) # Return seconds as days.
    
   # Sum these days to d_ini_or
   time = d_ini_or + time_days
      
#   time = time[0]
#   print('file time...', mdates.num2date(time).strftime("%y/%m/%d %H:%M"))   
   
   
   return time, lat, lon, sst, err, anom

def read_nc_ssh(files_ssh):
    
  time_ssh = []
  
  for ii, file in enumerate(files_ssh):
      
      nc        = netcdf.Dataset(file, 'r')     
      time_old  = nc.variables['time'][:]     # "days since 1950-01-01 00:00:00"
      lat1d     = nc.variables['latitude'][:]   
      lon360    = nc.variables['longitude'][:].squeeze() # (lon) from 0.125 to 359.875 by 0.25
      adt       = nc.variables['adt'][:]        
      ugos      = nc.variables['ugos'][:].squeeze() #"Absolute geostrophic velocity: zonal component"
      vgos      = nc.variables['vgos'][:].squeeze()
      sla       = nc.variables['sla'][:].squeeze()
      ugosa     = nc.variables['ugosa'][:].squeeze()
      vgosa     = nc.variables['vgosa'][:].squeeze()
      nc.close()      

      # time in python format
      # first date
      d_ini_or    = mdates.date2num(datetime(1950, 1, 1, 0, 0, 0))  
  
      # Sum time to original date
      time = d_ini_or + time_old
      
#      print(file)
#      print('time min... ', mdates.num2date(time.min()) )
#      print('time max... ', mdates.num2date(time.max()) )

      if ii==0:
          # lon from "0 to 360" format to "-180 to 180" format
          lon1d = np.copy(lon360)
          lon1d[lon1d>180] = lon1d[lon1d>180] - 360

          # lon and lat 2d array
          lon, lat = np.meshgrid(lon1d, lat1d) 
          
          adt_all   = adt
          ugos_all  = ugos
          vgos_all  = vgos
          sla_all   = sla
          ugosa_all = ugosa
          vgosa_all = vgosa
          time_ssh = np.append(time_ssh, time)
          
      elif ii>0:
          
          ind_time  = time>time_min #avoid overlap
          
          adt_all   = np.append(adt_all,   adt[ind_time,:,:],   axis=0)
          ugos_all  = np.append(ugos_all,  ugos[ind_time,:,:],  axis=0)
          vgos_all  = np.append(vgos_all,  vgos[ind_time,:,:],  axis=0)
          sla_all   = np.append(sla_all,   sla[ind_time,:,:],   axis=0)
          ugosa_all = np.append(ugosa_all, ugosa[ind_time,:,:], axis=0)
          vgosa_all = np.append(vgosa_all, vgosa[ind_time,:,:], axis=0)        
          time_ssh  = np.append(time_ssh,  time[ind_time])
          
#          print('original time shape', time.shape)
#          print('after removing overlapping',time[ind_time].shape )
          
      # to evoid overlap with the next file
      time_min = time.max() 
      

  return time_ssh, lon, lat, adt_all, ugos_all, vgos_all
  


def read_altimetry(file_alt):
    
   nc      = netcdf.Dataset(file_alt, 'r') 
   time_orig  = nc.variables['time'][:]  # days since 1950-01-01 00:00:00
   lat1d   = nc.variables['latitude'][:]  
   lon360  = nc.variables['longitude'][:]   # 0-360
   adt     = nc.variables['adt'][:]  #(time, lat, lon)
   ugos    = nc.variables['ugos'][:]  #(time, lat, lon) Absolute geostrophic velocity
   vgos    = nc.variables['vgos'][:]  #(time, lat, lon) Absolute geostrophic velocity
   sla     = nc.variables['sla'][:]  #(time, lat, lon)
   ugosa   = nc.variables['ugosa'][:]  #(time, lat, lon)
   vgosa   = nc.variables['vgosa'][:]  #(time, lat, lon)
   
   nc.close()    

   # lon from "0 to 360" format to "-180 to 180" format
   lon1d = np.ma.copy(lon360)
   lon1d[lon1d>180] = lon1d[lon1d>180] - 360

   # lon and lat 2d array
   lon, lat = np.meshgrid(lon1d, lat1d)

   # time in python format
   # first date
   d_ini_or    = mdates.date2num(datetime(1950, 1, 1, 0, 0, 0))  
  
   # Time is in days already
   time_days = time_orig
    
   # Sum these days to d_ini_or
   time = d_ini_or + time_days

    
   return time, lat, lon, adt, ugos, vgos, sla, ugosa, vgosa

def remove_mask_alt(time, lat, lon, adt, ugos, vgos, sla, ugosa, vgosa):
    
   adt.data[adt.mask==True]     = np.nan 
   ugos.data[ugos.mask==True]   = np.nan
   vgos.data[vgos.mask==True]   = np.nan 
   sla.data[sla.mask==True]     = np.nan
   ugosa.data[ugosa.mask==True] = np.nan
   vgosa.data[vgosa.mask==True] = np.nan

   return time.data, lat.data, lon.data, adt.data, ugos.data, vgos.data, sla.data, ugosa.data, vgosa.data

#def read_altimetry_nomask(file_alt):
#    
#   nc      = netcdf.Dataset(file_alt, 'r') 
#   time_orig  = nc.variables['time'][:]  # days since 1950-01-01 00:00:00
#   lat1d   = nc.variables['latitude'][:]  
#   lon360  = nc.variables['longitude'][:]   # 0-360
#   adt     = nc.variables['adt'][:]  #(time, lat, lon)
#   ugos    = nc.variables['ugos'][:]  #(time, lat, lon) Absolute geostrophic velocity
#   vgos    = nc.variables['vgos'][:]  #(time, lat, lon) Absolute geostrophic velocity
#   sla     = nc.variables['sla'][:]  #(time, lat, lon)
#   ugosa   = nc.variables['ugosa'][:]  #(time, lat, lon)
#   vgosa   = nc.variables['vgosa'][:]  #(time, lat, lon)
#   
#   nc.close()    
#
#   # lon from "0 to 360" format to "-180 to 180" format
#   lon1d = np.ma.copy(lon360)
#   lon1d[lon1d>180] = lon1d[lon1d>180] - 360
#
#   # lon and lat 2d array
#   lon, lat = np.meshgrid(lon1d, lat1d)
#
#   # time in python format
#   # first date
#   d_ini_or    = mdates.date2num(datetime(1950, 1, 1, 0, 0, 0))  
#  
#   # Time is in days already
#   time_days = time_orig
#    
#   # Sum these days to d_ini_or
#   time = d_ini_or + time_days
#
#    
#   return time.data, lat.data, lon.data, adt.data, ugos.data, vgos.data, sla.data, ugosa.data, vgosa.data

def create_basemap():
    
#  lonmin = -124 #124
#  lonmax = -122.5
#  latmin = 37.2
#  latmax = 38.4
  
  lonmin = -135 #124
  lonmax = -109
  latmin = 22
  latmax = 50
  
  fsize = 12
  # Create Basemap
  print ('starting Basemap...')
  m     = Basemap(llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,
              urcrnrlon=lonmax,resolution='h',projection='merc') 
    
  return m


def plot_coast(m, fsize):
    
    m.drawcoastlines()
    m.drawrivers(linewidth=0.5, color='k', zorder=6)
    m.fillcontinents(color='0.8', lake_color='0.8', zorder=5)
    parallels = np.arange(20.,50.,5)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],fontsize=fsize-1, zorder=8)
    meridians = np.arange(-135,-100,5)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],fontsize=fsize-1, zorder=9)

def create_basemap_sim_GS():
  
  ''' Basemap for the Lagrangian simulation region '''
  
  lonmin = -90 #-135
  lonmax = -50
  latmin = 15#20
  latmax = 50   

  print ('starting Basemap...')
  m     = Basemap(llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,
              urcrnrlon=lonmax,resolution='i',projection='merc') 
  return m


def create_basemap_sim():
  
  ''' Basemap for the Lagrangian simulation region '''
  
  lonmin = -150 
  lonmax = -109
  latmin = 22
  latmax = 50  

  print ('starting Basemap...')
  m     = Basemap(llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,
              urcrnrlon=lonmax,resolution='h',projection='merc') 
  return m

def plot_coast_sim(m, fsize):
    
    ''' Plot the coast for the Lagrangian simulation region '''
    
    m.drawcoastlines()
    #m.drawrivers(linewidth=0.5, color='k', zorder=6)
    m.fillcontinents(color='0.8', lake_color='0.8', zorder=5)
    parallels = np.arange(20.,50.,10)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],fontsize=fsize-1, zorder=8)
    meridians = np.arange(-150,-100,10)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],fontsize=fsize-1, zorder=9)

def plot_coast_Oleander(m, fsize):
    
    ''' Plot the coast for the Lagrangian simulation region '''
    
    m.drawcoastlines()
    #m.drawrivers(linewidth=0.5, color='k', zorder=6)
    m.fillcontinents(color='0.8', lake_color='0.8', zorder=5)
    parallels = np.arange(30.,50.,2)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],fontsize=fsize-1, zorder=8)
    meridians = np.arange(-80,-60,2)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],fontsize=fsize-1, zorder=9)



def cmap_div():
   # from: https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale

   my_cmap_div = LinearSegmentedColormap.from_list("", ['darkblue' ,
                                                 'steelblue',
                                                 'mediumseagreen',
                                                 'white',
                                                 'gold',
                                                 'salmon',
                                                 'darkred',])

   return my_cmap_div


def memade_cmap():
   # from: https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale

#   my_cmap_div = LinearSegmentedColormap.from_list("", ['darkblue' ,
#                                                 'steelblue',
#                                                 'limegreen',
#                                                 'white',
#                                                 'gold',
#                                                 'salmon',
#                                                 'darkred',])

    my_cmap = LinearSegmentedColormap.from_list("", ['darkblue' ,
                                                 'steelblue',
                                                 'mediumseagreen',#'limegreen',
                                                 'palegoldenrod',#'lightgoldenrodyellow',#'palegoldenrod', #'lemonchiffon',#'beige',
                                                 'gold',
                                                 'salmon',
                                                 'darkred',])      
    # check cmap 
#    x, y = np.meshgrid(np.linspace(-1,1,50), np.linspace(-1,1,50))
#    g = np.tile(np.linspace(0,1,50),(50,1))
#    
#    plt.figure()
#    plt.pcolormesh(x,y,g, cmap=my_cmap)
#    plt.colorbar()
#    plt.show()
    
    return my_cmap    

def check_year(time, year_checked) :
    
    ''' Return the indices of the desired year '''
    
    ind_yy = []
    for ind, tt in enumerate(time):
        yy = mdates.num2date(tt).year
 
        if yy == year_checked:
            ind_yy = np.append(ind_yy, int(ind))
            
    return ind_yy

def check_month(time, month_checked) :
    
    ''' Return the indices of the desired month '''
    
    ind_mm = []
    for ind, tt in enumerate(time):
        mm = mdates.num2date(tt).month
 
        if mm == month_checked:
            ind_mm = np.append(ind_mm, int(ind))
            
    return ind_mm    

def plot_map(m, lon, lat, sss, cmin, cmax, title_plot, dir_fig, filename, fsize):  
    
    x, y = m(lon, lat)
    
    fig = plt.figure(figsize=(8,8))
    ax  = plt.subplot(111)  
    
    plot_coast(m, fsize)  
    
    #plt.pcolormesh(x,y,sss,cmap=plt.cm.jet)
    plt.pcolormesh(x,y,sss,cmap=memade_cmap())
    #plt.pcolormesh(x,y,sss,cmap=sns.diverging_palette(128, 240, as_cmap=True))# cmap=plt.cm.jet)
    plt.clim(vmin=cmin, vmax=cmax)
    
    plt.title(title_plot)
    
    plt.colorbar()
    plt.tight_layout()
    
    fig.savefig(dir_fig + filename, dpi=200)   
    
    
def length_lon_lat_degs(lat_mean_rad):
      
  ''' 
  Function to infer the length of a degree of longitude
  and the length of a degree of latitude
  at a specific latitude position. 
  Assuming that the Earth is an ellipsoid.
  
  input:  latitude in RADIANS!!!
  output: length_deg_lon, length_deg_lat
  
  from:
      https://en.wikipedia.org/wiki/Longitude#Length_of_a_degree_of_longitude
      https://en.wikipedia.org/wiki/Latitude#Length_of_a_degree_of_latitude
  '''  
  
  ''' Earth parameters '''
  
  a = 6378137.0                # m (equatorial radius)
  b = 6356752.3142             # m  (polar radius)
  
  ecc_2 = (a**2 - b**2) / a**2 # eccentricity squared


  ''' The length of a degree of longitude is... '''
  
  divident_lon = (a*np.pi/180) * np.cos(lat_mean_rad)
  divisor_lon  = np.sqrt(1 - (ecc_2*np.sin(lat_mean_rad)*np.sin(lat_mean_rad)))
  
  length_deg_lon = divident_lon / divisor_lon 
  
  
  ''' The length of a degree of latitude is... '''
  
  divident_lat = (a*np.pi/180) * (1 - ecc_2)
  divisor_lat  = (1 - (ecc_2 * np.sin(lat_mean_rad) * np.sin(lat_mean_rad)))**(3/2)
  
  length_deg_lat = divident_lat / divisor_lat
  
  
  return length_deg_lon, length_deg_lat    



def save_currents_simulation(savedirnc, filename, time_new, lon, lat, u_int, v_int):
    
  ''' 
  Save a Netcdf file for the currents to do the Lagrangian simulation with Parcels
  '''  
  nc = netcdf.Dataset(savedirnc + filename, 'w', format='NETCDF3_CLASSIC')
    
  # Create the dimensions...
  nc.createDimension('lat',  lat.shape[0]) 
  nc.createDimension('lon',  lon.shape[0])    
  nc.createDimension('time', time_new.shape[0]) 
  
  # Create the variables...
 # 'f8' (64-bit floating point)
#  nc.createVariable('lon',     'f4', ('dlat', 'dlon'))
#  nc.createVariable('lat',     'f4', ('dlat', 'dlon'))
  nc.createVariable('lon',     'f4', ('lon'))
  nc.createVariable('lat',     'f4', ('lat'))
  nc.createVariable('time',    'f8', ('time'))
  
  nc.createVariable('ug_int',  'f4', ('time', 'lat', 'lon'))#,fill_value=tem.get_fill_value())
  nc.createVariable('vg_int',  'f4', ('time', 'lat', 'lon'))
  

  # Write in variable attributes...
  nc.variables['time'].long_name       = 'time count'
  nc.variables['time'].units           = 'seconds'
  
  nc.variables['lon'].long_name        = 'longitude'
  nc.variables['lon'].units            = 'degrees_east'

  nc.variables['lat'].long_name        = 'latitude'
  nc.variables['lat'].units            = 'degrees_north'
  
  nc.variables['ug_int'].long_name     = 'Zonal geostrophic velocity interpolated'
  nc.variables['ug_int'].units         = 'm/s'

  nc.variables['vg_int'].long_name     = 'Meridional geostrophic velocity interpolated'
  nc.variables['vg_int'].units         = 'm/s'
  
  
  nc.close()   

  # ------------------------------------------
  # Save data
  # ------------------------------------------
     
  nc = netcdf.Dataset(savedirnc + filename, 'a', format='NETCDF3_CLASSIC')
  nc.variables['time'][:]   = time_new
  nc.variables['lon'][:]    = lon # parcels need 1D
  nc.variables['lat'][:]    = lat #      parcels need 1D
  nc.variables['ug_int'][:] = u_int
  nc.variables['vg_int'][:] = v_int

  nc.close() 
  
  
def read_SMAP(file):
        
      nc   = netcdf.Dataset(file, 'r')
      lon360  = nc.variables['lon'][:].data # center longitude of grid cell 
      lat1d  = nc.variables['lat'][:].data # center latitude of grid cell
      # reference time of analyzed variable field corresponding 
      # to center of the product time interval
      time_orig = nc.variables['time'][:].data # seconds since 2000-01-01T00:00:00Z
      nobs      = nc.variables['nobs'][:].data # Number of observations for L3 average
      sss_masked       = nc.variables['sss_smap'][:]
       # Reference sea surface salinity from HYCOM: check what is this!
      sss_ref   = nc.variables['sss_ref'][:]
      # Ancillary sea surface temperature (from CMC): check what is this!
      #surtep    = nc.variables['sss_ref'][:] 
      nc.close() 
    
      lon1d = np.copy(lon360)
      lon1d[lon1d>180] = lon1d[lon1d>180] - 360
    
      lon, lat = np.meshgrid(lon1d, lat1d)
      
      # remove mask in sss_masked
      sss = np.copy(sss_masked.data)
      sss[sss_masked.mask==True] = np.nan
      
      ''' convert time to python format: '''
    
      # first date
      d_ini_or    = mdates.date2num(datetime(2000, 1, 1, 0, 0, 0))  
  
      # Return seconds as days
      time_days = np.zeros(time_orig.shape)
  
      for ind, tt in enumerate(time_orig):
         time_days[ind] = mdates.seconds(tt) # Return seconds as days.
    
      # Sum these days to d_ini_or
      time = d_ini_or + time_days
      
      time = time[0]
      print('')
      print('SSS - center of the product time interval...') 
      print(mdates.num2date(time).strftime("%Y-%m-%d %H:%M"))
        
      
      return lon, lat, time, sss    

def delimit_region(lon_smap_glo, lat_smap_glo, sss_smap_glo):
    
    # delimit SMAP data to the simulation area
    lonmin = -150 #-135
    lonmax = -109
    latmin = 22
    latmax = 50    
    
    ilon = np.logical_and(lon_smap_glo[0,:]>=lonmin, lon_smap_glo[0,:]<=lonmax)
    ilat = np.logical_and(lat_smap_glo[:,0]>=latmin, lat_smap_glo[:,0]<=latmax)
    
    lon_smap = lon_smap_glo[ilat,:][:,ilon]
    lat_smap = lat_smap_glo[ilat,:][:,ilon]
    sss_smap = sss_smap_glo[ilat,:][:,ilon]
    
    return lon_smap, lat_smap, sss_smap     