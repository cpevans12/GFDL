#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 14:35:15 2025

@author: cpe28
"""

import netCDF4 as nc
import numpy as np
import glob as glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import get_cmap
from scipy.stats import theilslopes
from matplotlib.colors import BoundaryNorm

image_dir = '/Volumes/Elements/GFDL_project/images/new_images/'
# =============================================================================
# 
# =============================================================================
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
nclimgrid_prcptot = np.array(nclimgrid_data.variables['prcptot'])
mask = np.isclose(nclimgrid_prcptot,1e+20)
nclimgrid_prcptot[mask] = np.nan
lat = np.array(nclimgrid_data.variables['lat'])
lon = np.array(nclimgrid_data.variables['lon'])

slopes = np.full((nclimgrid_prcptot.shape[1], nclimgrid_prcptot.shape[2]), np.nan)
years = np.arange(nclimgrid_prcptot.shape[0])
for i in range(nclimgrid_prcptot.shape[1]):
    for j in range(nclimgrid_prcptot.shape[2]):
        time_series = nclimgrid_prcptot[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
nclimgrid_slopes = slopes



livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/ClimdexPhard_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
lat = np.array(livneh_data.variables['lat'])
lon = np.array(livneh_data.variables['lon'])
livneh_prcptot = np.array(livneh_data.variables['prcptot'])
mask = np.isclose(livneh_prcptot,1e+20)
livneh_prcptot[mask] = np.nan
mask = np.isnan(nclimgrid_prcptot)
livneh_prcptot[mask] = np.nan

slopes = np.full((livneh_prcptot.shape[1], livneh_prcptot.shape[2]), np.nan)
intercepts = np.full((livneh_prcptot.shape[1], livneh_prcptot.shape[2]), np.nan)
years = np.arange(livneh_prcptot.shape[0])

for i in range(livneh_prcptot.shape[1]):
    for j in range(livneh_prcptot.shape[2]):
        time_series = livneh_prcptot[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
        intercepts[i,j] = intercept

livneh_slopes = slopes

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,nclimgrid_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),extend='both')
plt.savefig(image_dir+'nclimgrid_prcptot_slopes_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()



fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,livneh_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),extend='both')
plt.savefig(image_dir+'livneh_prcptot_slopes_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

from matplotlib.colorbar import ColorbarBase
levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(bottom=0.5)

# Create the colorbar
cb = ColorbarBase(ax, cmap=get_cmap("BrBG"), norm=norm, orientation='horizontal', extend='both')

cb.set_label(r'Theil Slope (mm/year)')
plt.savefig(image_dir+"prcptot_slope_nclimgrid_colorbar.png",dpi=300,bbox_inches='tight')
plt.show()
plt.close()

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']

star_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/ClimdexPhard_STAR-*CONUS16thD2021_1985-2014_ann.nc")
star_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['star_'+str(models[i])+'_prcptot'] = np.array(star_data.variables['prcptot'])
    mask = np.isclose(globals()['star_'+str(models[i])+'_prcptot'],1e+20)
    globals()['star_'+str(models[i])+'_prcptot'][mask] = np.nan
    if i in {15}:
        lat = np.array(star_data.variables['lat'])
        lon = np.array(star_data.variables['lon'])
    
    slopes = np.full((globals()['star_'+str(models[i])+'_prcptot'].shape[1], globals()['star_'+str(models[i])+'_prcptot'].shape[2]), np.nan)
    years = np.arange(globals()['star_'+str(models[i])+'_prcptot'].shape[0])

    for j in range(globals()['star_'+str(models[i])+'_prcptot'].shape[1]):
        for k in range(globals()['star_'+str(models[i])+'_prcptot'].shape[2]):
            time_series = globals()['star_'+str(models[i])+'_prcptot'][:,j,k]
            if np.all(np.isnan(time_series)):
                continue
            slope, intercept, _,_ = theilslopes(time_series,years)
            slopes[j,k] = slope
    
    globals()['star_slope_'+str(models[i])+'_prcptot'] = slopes


all_star = np.empty([16,30,444,922])
for i in range(16):
    all_star[i,:,:,:] = globals()['star_'+str(models[i])+'_prcptot']

star_ensemble_prcptot_avg = np.nanmean(all_star,axis=0)

slopes = np.full((star_ensemble_prcptot_avg.shape[1], star_ensemble_prcptot_avg.shape[2]), np.nan)
years = np.arange(star_ensemble_prcptot_avg.shape[0])
for i in range(star_ensemble_prcptot_avg.shape[1]):
    for j in range(star_ensemble_prcptot_avg.shape[2]):
        time_series = star_ensemble_prcptot_avg[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
star_ensemble_prcptot_slopes = slopes        


levels = np.arange(-5,5.1,0.5)
norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,star_ensemble_prcptot_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),extend='both')
plt.savefig(image_dir+'star_prcptot_slope_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

levels = np.arange(-5,5.1,0.5)
norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(bottom=0.5)

# Create the colorbar
cb = ColorbarBase(ax, cmap=get_cmap("BrBG"), norm=norm, orientation='horizontal', extend='both')

cb.set_label(r'Theil Slope (mm/year)')
plt.savefig(image_dir+"prcptot_slope_star_colorbar.png",dpi=300,bbox_inches='tight')
plt.show()
plt.close()

# =============================================================================
# 
# =============================================================================
star_ensemble_slopes_matrix = np.empty_like(nclimgrid_slopes)
for i in range(444):
    for j in range(922):
        if np.isnan(nclimgrid_slopes[i,j]):
            star_ensemble_slopes_matrix[i,j] = np.nan
        else:
            model = star_ensemble_prcptot_slopes[i,j]
            nclim = nclimgrid_slopes[i,j]
            if model > 0. and nclim > 0.:
                star_ensemble_slopes_matrix[i,j] = 2
            elif model > 0. and nclim < 0.:
                star_ensemble_slopes_matrix[i,j] = 1
            elif model < 0. and nclim > 0.:
                star_ensemble_slopes_matrix[i,j] = -1
            elif model < 0. and nclim < 0.:
                star_ensemble_slopes_matrix[i,j] = -2


contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
cb = plt.contourf(lon, lat, star_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
cbar.set_label(r'$\frac{STAR}{nClimGrid}$',fontsize=12)
cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
# plt.savefig(image_dir+"star_ensemble_prcptot_slopes_directions_comparison.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

from matplotlib.colors import ListedColormap
contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = BoundaryNorm(contours,ncolors=4)
colors = ['#884c1b','#efd8ae','#a7d3d1','#005b5b']

# 2. Make a ListedColormap
cmap = ListedColormap(colors)

fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(bottom=0.5)

# Create the colorbar
cb = ColorbarBase(ax, cmap=cmap,norm=norm, orientation='horizontal',ticks=[-1.5, -0.5, 0.5, 1.5])
cb.set_label(r'$\frac{Model}{Observations}$',fontsize=12)
cb.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
plt.savefig(image_dir+"prcptot_slope_comparison_colorbar.png",dpi=300,bbox_inches='tight')
plt.show()
plt.close()



from scipy.stats import kendalltau

nclimgrid_tau = np.empty_like(nclimgrid_slopes)
nclimgrid_pval = np.empty_like(nclimgrid_slopes)
mask = np.isnan(nclimgrid_slopes)
nclimgrid_tau[mask] = np.nan
nclimgrid_pval[mask] = np.nan
years = np.arange(1985,2015)
for i in range(nclimgrid_tau.shape[0]):
    for j in range(nclimgrid_tau.shape[1]):
        if np.isnan(nclimgrid_slopes[i,j]):
            nclimgrid_tau[i,j] = np.nan
            nclimgrid_pval[i,j] = np.nan
        else:
            time_series = nclimgrid_prcptot[:,i,j]
            tau, p_value = kendalltau(years, time_series)
            nclimgrid_tau[i,j] = tau
            nclimgrid_pval[i,j] = p_value        
nclimgrid_masked_slopes = np.where(nclimgrid_pval<0.1,nclimgrid_slopes,np.nan)


levels = np.arange(-15,15.1,0.5)
sig_norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

positive_mask = nclimgrid_masked_slopes > 0.5                     # Where slope is significantly positive
negative_mask = nclimgrid_masked_slopes < -0.5 




fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
cb = plt.contourf(lon, lat, star_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cb2 = plt.contourf(lon,lat,nclimgrid_masked_slopes,levels=levels,colors='none',hatches=['++++'],transform=ccrs.PlateCarree(),norm=sig_norm,extend='both')
plt.savefig(image_dir+"star_ensemble_prcptot_slopes_directions_comparison_significance_version1.png",dpi=300,bbox_inches='tight',pad_inches=0,transparent=True)
plt.show()
plt.close()


import numpy as np
import matplotlib.pyplot as plt

# Define the custom colormap
colors = [
    (0.0, 'red'),    # Far negative
    (0.5, 'white'),  # Zero
    (1.0, 'blue')    # Far positive
]

custom_cmap = LinearSegmentedColormap.from_list('red_white_blue', colors)


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
cb = plt.contourf(lon, lat, star_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cb2 = plt.contourf(lon,lat,nclimgrid_masked_slopes,levels=levels,colors='none',hatches=['++++'],transform=ccrs.PlateCarree(),norm=sig_norm,extend='both')
overlay = plt.pcolormesh(lon,lat,nclimgrid_masked_slopes,cmap=custom_cmap,alpha=0.4,transform=ccrs.PlateCarree())
plt.savefig(image_dir+"star_ensemble_prcptot_slopes_directions_comparison_significance_version2.png",dpi=300,bbox_inches='tight',pad_inches=0,transparent=True)
plt.show()
plt.close()

# =============================================================================
# 
# =============================================================================

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']

letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

loca_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/ClimdexPhard_LOCA2-*_historical_r1i1p1f1_gn_CONUS16thD2021_1985-2014_ann.nc")
loca_files.sort()

all_loca = np.empty([16,30,444,922])
for i in range(16):
    loca_data = nc.Dataset(loca_files[i])
    globals()['LOCA2_'+str(models[i])+'_prcptot'] = np.array(loca_data.variables['prcptot'])
    mask = np.isclose(globals()['LOCA2_'+str(models[i])+'_prcptot'],1e+20)
    globals()['LOCA2_'+str(models[i])+'_prcptot'][mask] = np.nan
    mask = np.isnan(nclimgrid_prcptot)
    globals()['LOCA2_'+str(models[i])+'_prcptot'][mask]= np.nan
    if i in {15}:
        lat = np.array(loca_data.variables['lat'])
        lon = np.array(loca_data.variables['lon'])
    all_loca[i,:,:,:] = globals()['LOCA2_'+str(models[i])+'_prcptot']
    slopes = np.full((globals()['LOCA2_'+str(models[i])+'_prcptot'].shape[1], globals()['LOCA2_'+str(models[i])+'_prcptot'].shape[2]), np.nan)
    years = np.arange(globals()['LOCA2_'+str(models[i])+'_prcptot'].shape[0])

    for j in range(globals()['LOCA2_'+str(models[i])+'_prcptot'].shape[1]):
        for k in range(globals()['LOCA2_'+str(models[i])+'_prcptot'].shape[2]):
            time_series = globals()['LOCA2_'+str(models[i])+'_prcptot'][:,j,k]
            if np.all(np.isnan(time_series)):
                continue
            slope, intercept, _,_ = theilslopes(time_series,years)
            slopes[j,k] = slope
    
    globals()['loca2_slope_'+str(models[i])+'_prcptot'] = slopes

loca_ens_precip = np.nanmin(all_loca,axis=0)

slopes = np.full((loca_ens_precip.shape[1], loca_ens_precip.shape[2]), np.nan)
years = np.arange(loca_ens_precip.shape[0])
for i in range(loca_ens_precip.shape[1]):
    for j in range(loca_ens_precip.shape[2]):
        time_series = loca_ens_precip[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
loca_ens_precip_slopes = slopes


levels = np.arange(-5,5.1,0.5)
norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,loca_ens_precip_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),extend='both')
plt.savefig(image_dir+'loca_prcptot_slope_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

loca_ensemble_slopes_matrix = np.empty_like(livneh_slopes)
for i in range(444):
    for j in range(922):
        if np.isnan(livneh_slopes[i,j]):
            loca_ensemble_slopes_matrix[i,j] = np.nan
        else:
            model = loca_ens_precip_slopes[i,j]
            liv = livneh_slopes[i,j]
            if model > 0. and liv > 0.:
                loca_ens_precip_slopes[i,j] = 2
            elif model > 0. and liv < 0.:
                loca_ens_precip_slopes[i,j] = 1
            elif model < 0. and liv > 0.:
                loca_ens_precip_slopes[i,j] = -1
            elif model < 0. and liv < 0.:
                loca_ens_precip_slopes[i,j] = -2

from scipy.stats import kendalltau

livneh_tau = np.empty_like(livneh_slopes)
livneh_pval = np.empty_like(livneh_slopes)
mask = np.isnan(livneh_slopes)
livneh_tau[mask] = np.nan
livneh_pval[mask] = np.nan
years = np.arange(1985,2015)
for i in range(livneh_tau.shape[0]):
    for j in range(livneh_tau.shape[1]):
        if np.isnan(livneh_slopes[i,j]):
            livneh_tau[i,j] = np.nan
            livneh_pval[i,j] = np.nan
        else:
            time_series = livneh_prcptot[:,i,j]
            tau, p_value = kendalltau(years, time_series)
            livneh_tau[i,j] = tau
            livneh_pval[i,j] = p_value        
livneh_masked_slopes = np.where(livneh_pval<0.1,livneh_slopes,np.nan)

levels = np.arange(-15,15.1,0.5)
sig_norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
cb = plt.contourf(lon, lat, loca_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cb2 = plt.contourf(lon,lat,livneh_masked_slopes,levels=levels,colors='none',hatches=['++++'],transform=ccrs.PlateCarree(),norm=sig_norm,extend='both')
plt.savefig(image_dir+"loca_ensemble_prcptot_slopes_directions_comparison_significance_version1.png",dpi=300,bbox_inches='tight',pad_inches=0,transparent=True)
plt.show()
plt.close()

# =============================================================================
# 
# =============================================================================

