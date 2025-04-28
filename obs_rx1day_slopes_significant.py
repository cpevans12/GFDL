##this looks at the slope of rx1day in the observational products and then finds the significance and only plots where the slope is significant

import netCDF4 as nc
import numpy as np
import glob as glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import get_cmap
from scipy.stats import theilslopes
from matplotlib.colors import BoundaryNorm

image_dir = '/Volumes/Elements/GFDL_project/images/'
# =============================================================================
# =============================================================================
# # 
# =============================================================================
# =============================================================================
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
nclimgrid_rx1day = np.array(nclimgrid_data.variables['rx1day'])
mask = np.isclose(nclimgrid_rx1day,1e+20)
nclimgrid_rx1day[mask] = np.nan

slopes = np.full((nclimgrid_rx1day.shape[1], nclimgrid_rx1day.shape[2]), np.nan)
intercepts = np.full((nclimgrid_rx1day.shape[1], nclimgrid_rx1day.shape[2]), np.nan)
years = np.arange(nclimgrid_rx1day.shape[0])

for i in range(nclimgrid_rx1day.shape[1]):
    for j in range(nclimgrid_rx1day.shape[2]):
        time_series = nclimgrid_rx1day[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
        intercepts[i,j] = intercept
nclimgrid_slopes = slopes
nclimgrid_intercepts = intercepts
###
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/ClimdexPhard_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
lat = np.array(livneh_data.variables['lat'])
lon = np.array(livneh_data.variables['lon'])
livneh_rx1day = np.array(livneh_data.variables['rx1day'])
mask = np.isclose(livneh_rx1day,1e+20)
livneh_rx1day[mask] = np.nan
mask = np.isnan(nclimgrid_rx1day)
livneh_rx1day[mask] = np.nan

slopes = np.full((livneh_rx1day.shape[1], livneh_rx1day.shape[2]), np.nan)
intercepts = np.full((livneh_rx1day.shape[1], livneh_rx1day.shape[2]), np.nan)
years = np.arange(livneh_rx1day.shape[0])

for i in range(livneh_rx1day.shape[1]):
    for j in range(livneh_rx1day.shape[2]):
        time_series = livneh_rx1day[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
        intercepts[i,j] = intercept

livneh_slopes = slopes
livneh_intercepts = intercepts
###
prism_rx1day = np.load("/Volumes/Elements/GFDL_project/data/prism_yrly_rx1day_regridded.npy")
mask = np.isnan(nclimgrid_rx1day)
prism_rx1day[mask] = np.nan

slopes = np.full((prism_rx1day.shape[1], prism_rx1day.shape[2]), np.nan)
intercepts = np.full((prism_rx1day.shape[1], prism_rx1day.shape[2]), np.nan)
years = np.arange(prism_rx1day.shape[0])

for i in range(prism_rx1day.shape[1]):
    for j in range(prism_rx1day.shape[2]):
        time_series = prism_rx1day[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
        intercepts[i,j] = intercept

prism_slopes = slopes
prism_intercepts = intercepts
# =============================================================================
# significance
# =============================================================================

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
            time_series = nclimgrid_rx1day[:,i,j]
            tau, p_value = kendalltau(years, time_series)
            nclimgrid_tau[i,j] = tau
            nclimgrid_pval[i,j] = p_value        
nclimgrid_masked_slopes = np.where(nclimgrid_pval<0.1,nclimgrid_slopes,np.nan)


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
            time_series = livneh_rx1day[:,i,j]
            tau, p_value = kendalltau(years, time_series)
            livneh_tau[i,j] = tau
            livneh_pval[i,j] = p_value        
livneh_masked_slopes = np.where(livneh_pval<0.1,livneh_slopes,np.nan)


prism_tau = np.empty_like(prism_slopes)
prism_pval = np.empty_like(prism_slopes)
mask = np.isnan(prism_slopes)
prism_tau[mask] = np.nan
prism_pval[mask] = np.nan
years = np.arange(1985,2015)
for i in range(prism_tau.shape[0]):
    for j in range(prism_tau.shape[1]):
        if np.isnan(prism_slopes[i,j]):
            prism_tau[i,j] = np.nan
            prism_pval[i,j] = np.nan
        else:
            time_series = prism_rx1day[:,i,j]
            tau, p_value = kendalltau(years, time_series)
            prism_tau[i,j] = tau
            prism_pval[i,j] = p_value        
prism_masked_slopes = np.where(prism_pval<0.1,prism_slopes,np.nan)

total_cells = np.count_nonzero(~np.isnan(livneh_slopes))
count_livneh = np.count_nonzero(~np.isnan(livneh_masked_slopes))
percent_livneh = (count_livneh/total_cells)*100.

total_cells = np.count_nonzero(~np.isnan(nclimgrid_slopes))
count_nclimgrid = np.count_nonzero(~np.isnan(nclimgrid_masked_slopes))
percent_nclimgrid = (count_nclimgrid/total_cells)*100.

total_cells = np.count_nonzero(~np.isnan(prism_slopes))
count_prism = np.count_nonzero(~np.isnan(prism_masked_slopes))
percent_prism = (count_prism/total_cells)*100.


variables = [livneh_masked_slopes,nclimgrid_masked_slopes,prism_masked_slopes]
percents = [percent_livneh,percent_nclimgrid,percent_prism]
letters = ["(a)","(b)","(c)"]
titles = ["Livneh","nClimGrid","PRISM"]

levels = np.arange(-3,3.1,0.5)
norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig = plt.figure(figsize=(14,6))
for i in range(3):
    ax = fig.add_subplot(1,3,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-126,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='50m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    ax.set_title(titles[i])
    cb = plt.contourf(lon,lat,variables[i],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9,pad=0.1)
    cbar.set_label(r"$m_{\mathrm{TS}}$ rx1day (p<0.1)",fontsize=10)
    cbar.ax.tick_params(labelsize=7)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    ax.text(0.8, 0.27, "Percent\nSignificant:\n"+f"{percents[i]:.2f}%", transform=ax.transAxes, fontsize=7,verticalalignment='top', bbox=props)
plt.savefig(image_dir+"livneh_nclimgrid_prism_rx1day_slopes_significant_p01.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
