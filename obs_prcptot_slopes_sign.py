#this code compares the direction of the slope between the observational products, so if they all agree, one value is assigned, if only 2 agree and one is different, a different value is assigned, etc.
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
nclimgrid_prcptot = np.array(nclimgrid_data.variables['prcptot'])
mask = np.isclose(nclimgrid_prcptot,1e+20)
nclimgrid_prcptot[mask] = np.nan

slopes = np.full((nclimgrid_prcptot.shape[1], nclimgrid_prcptot.shape[2]), np.nan)
intercepts = np.full((nclimgrid_prcptot.shape[1], nclimgrid_prcptot.shape[2]), np.nan)
years = np.arange(nclimgrid_prcptot.shape[0])

for i in range(nclimgrid_prcptot.shape[1]):
    for j in range(nclimgrid_prcptot.shape[2]):
        time_series = nclimgrid_prcptot[:,i,j]
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
livneh_intercepts = intercepts
###
prism_prcptot = np.load("/Volumes/Elements/GFDL_project/data/prism_yrly_prcptot_regridded.npy")
mask = np.isnan(nclimgrid_prcptot)
prism_prcptot[mask] = np.nan

slopes = np.full((prism_prcptot.shape[1], prism_prcptot.shape[2]), np.nan)
intercepts = np.full((prism_prcptot.shape[1], prism_prcptot.shape[2]), np.nan)
years = np.arange(prism_prcptot.shape[0])

for i in range(prism_prcptot.shape[1]):
    for j in range(prism_prcptot.shape[2]):
        time_series = prism_prcptot[:,i,j]
        if np.all(np.isnan(time_series)):
            continue
        slope, intercept, _,_ = theilslopes(time_series,years)
        slopes[i,j] = slope
        intercepts[i,j] = intercept

prism_slopes = slopes
prism_intercepts = intercepts
# =============================================================================
# =============================================================================
# # 
# =============================================================================
# =============================================================================
variables = [livneh_slopes,nclimgrid_slopes,prism_slopes]

letters = ["(a)","(b)","(c)"]
titles = ["Livneh","nClimGrid","PRISM"]

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)
states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')

fig = plt.figure(figsize=(14,6))
for i in range(3):
    ax = fig.add_subplot(1,3,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='50m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    ax.set_title(titles[i])
    cb = plt.contourf(lon,lat,variables[i],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9,pad=0.1)
    cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (mm/year)",fontsize=10)
    cbar.ax.tick_params(labelsize=7)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.savefig(image_dir+"livneh_nclimgrid_prism_prcptot_slopes.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
# =============================================================================
# 
# =============================================================================

slopes_matrix = np.empty_like(livneh_slopes)
for i in range(444):
    for j in range(922):
        # Handle missing data
        if np.isnan(livneh_slopes[i, j]):
            slopes_matrix[i, j] = np.nan
        else:
            # Extract the slope values
            liv = livneh_slopes[i, j]
            nclim = nclimgrid_slopes[i, j]
            prism = prism_slopes[i, j]

            # Assign specific codes based on the sign of slopes
            if liv > 0 and nclim > 0 and prism > 0:
                slopes_matrix[i, j] = 4  # All positive
            elif liv > 0 and nclim > 0 and prism < 0:
                slopes_matrix[i, j] = 3  # liv and nclim positive, prism negative
            elif liv > 0 and nclim < 0 and prism > 0:
                slopes_matrix[i, j] = 2  # liv and prism positive, nclim negative
            elif liv < 0 and nclim > 0 and prism > 0:
                slopes_matrix[i, j] = 1  # liv negative, nclim and prism positive
            elif liv > 0 and nclim < 0 and prism < 0:
                slopes_matrix[i, j] = -1  # liv positive, nclim and prism negative
            elif liv < 0 and nclim > 0 and prism < 0:
                slopes_matrix[i, j] = -2  # liv and prism negative, nclim positive
            elif liv < 0 and nclim < 0 and prism > 0:
                slopes_matrix[i, j] = -3  # liv and nclim negative, prism positive
            elif liv < 0 and nclim < 0 and prism < 0:
                slopes_matrix[i, j] = -4  # All negative


total_cells = np.count_nonzero(~np.isnan(slopes_matrix))
count1 = np.count_nonzero(slopes_matrix == 4.0)
count2 = np.count_nonzero(slopes_matrix == -4.0)
percent = ((count1+count2)/total_cells)*100.

count3 = np.count_nonzero(slopes_matrix == 3.0)
count4 = np.count_nonzero(slopes_matrix == 2.0)
count5 = np.count_nonzero(slopes_matrix == 1.0)
count6 = np.count_nonzero(slopes_matrix == -1.0)
count7 = np.count_nonzero(slopes_matrix == -2.0)
count8 = np.count_nonzero(slopes_matrix == -3.0)

percentages = [((count2/total_cells)*100.),((count8/total_cells)*100.),((count7/total_cells)*100.),((count6/total_cells)*100.),((count5/total_cells)*100.),((count4/total_cells)*100.),((count3/total_cells)*100.),((count1/total_cells)*100.)]

contours = [-4.1,-3.1,-2.1, -1.1,0, 1.1,2.1, 3.1,4.1]
norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-126, -65, 24, 50], crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
cb = plt.contourf(lon, lat, slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-3.5,-2.5,-1.5, -0.5, 0.5, 1.5,2.5,3.5],pad=0.1)
cbar.set_label("Livneh/nClimGrid/PRISM",fontsize=12)
cbar.ax.set_xticklabels(["-/-/-", "-/-/+","-/+/-","+/-/-","-/+/+","+/-/+", "+/+/-", "+/+/+"], fontsize=8)
for tick, pct in zip(cbar.ax.get_xticks(), percentages):cbar.ax.text(tick, 0.5, f"{pct:.1f}%", transform=cbar.ax.transData, ha='center', va='center', fontsize=9, color='k')
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
# ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.85, 0.15, "PA: "+f"{percent:.2f}", transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
plt.title("prcptot slope agreement",fontsize=12)
plt.savefig(image_dir+"livneh_nclimgrid_prism_prcptot_slopes_directions_comparison_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

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
            time_series = nclimgrid_prcptot[:,i,j]
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
            time_series = livneh_prcptot[:,i,j]
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
            time_series = prism_prcptot[:,i,j]
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

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)
states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='50m',facecolor='none')

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
    cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (p<0.1)",fontsize=10)
    cbar.ax.tick_params(labelsize=7)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    ax.text(0.8, 0.27, "Percent\nSignificant:\n"+f"{percents[i]:.2f}%", transform=ax.transAxes, fontsize=7,verticalalignment='top', bbox=props)
plt.savefig(image_dir+"livneh_nclimgrid_prism_prcptot_slopes_significant_p01.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
