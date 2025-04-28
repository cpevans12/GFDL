##finds the slope of prcptot for all the observations and models

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

nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
nclimgrid_prcptot = np.array(nclimgrid_data.variables['prcptot'])
mask = np.isclose(nclimgrid_prcptot,1e+20)
nclimgrid_prcptot[mask] = np.nan
lat = np.array(nclimgrid_data.variables['lat'])
lon = np.array(nclimgrid_data.variables['lon'])

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
################3
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

lat_index1 = find_nearest(lat, 40.5)
lat_index2 = find_nearest(lat, 47.5)
lon_index1 = find_nearest(lon,360-80.0)
lon_index2 = find_nearest(lon,360-66.0)


ne_nclimgrid_slopes = nclimgrid_slopes[lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]
northeast_slopes_cleaned_nclimgrid = ne_nclimgrid_slopes[~np.isnan(ne_nclimgrid_slopes)]

ne_livneh_slopes = livneh_slopes[lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]
northeast_slopes_cleaned_livneh = ne_livneh_slopes[~np.isnan(ne_livneh_slopes)]

ne_prism_slopes = prism_slopes[lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]
northeast_slopes_cleaned_prism = ne_prism_slopes[~np.isnan(ne_prism_slopes)]

all_trends = [northeast_slopes_cleaned_nclimgrid,northeast_slopes_cleaned_livneh,northeast_slopes_cleaned_prism]

plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
# Customize plot
plt.title("Distribution of prcptot Trends (Northeast)")
plt.ylabel("Trend (mm per year)")
plt.xticks([1, 2, 3], ["nClimGrid", "Livneh", "PRISM"])  # Adjust X-axis label to represent dataset
plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
plt.ylim(-1,1)
plt.show()
#
#
#
#
states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(121,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
ax.set_title("Livneh")
cb = plt.contourf(lon,lat,livneh_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (mm/year)",fontsize=10)
cbar.ax.tick_params(labelsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
ax = fig.add_subplot(122,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
ax.set_title("nClimGrid")
cb = plt.contourf(lon,lat,nclimgrid_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
cb.ax.tick_params(labelsize=8)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (mm/year)",fontsize=10)
cbar.ax.tick_params(labelsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(b)",transform=ax.transAxes)
plt.savefig(image_dir+"obs_prcptot_slopes.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


obs_percent_diff = ((livneh_slopes/nclimgrid_slopes)-1)*100.





variables = [livneh_slopes,nclimgrid_slopes,obs_percent_diff]
titles = ['Livneh','nClimGrid',r'$\frac{Livneh}{nClimGrid}$']
letters = ['(a)','(b)','(c)']
cbar_titles = [r"$m_{\mathrm{TS}}$ prcptot (mm/year)",r"$m_{\mathrm{TS}}$ prcptot (mm/year)","Percent Difference"]

fig = plt.figure(figsize=(6,14))
for i in range(3):
    ax = fig.add_subplot(3,1,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    if i in {0,1}:
        contour_levels = np.arange(-15,15.1,0.5)
        norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)
        cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
        cb.ax.tick_params(labelsize=8)
    else:
        norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
        contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]
        cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
        cb.ax.tick_params(labelsize=8)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(cbar_titles[i],fontsize=12)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(titles[i],fontsize=12)
plt.savefig(image_dir+"obs_prcptot_slopes_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()



slopes_matrix = np.empty_like(livneh_slopes)

for i in range(444):
    for j in range(922):
        if np.isnan(livneh_slopes[i,j]):
            slopes_matrix[i,j] = np.nan
        else:
            liv = livneh_slopes[i,j]
            nclim = nclimgrid_slopes[i,j]
            if liv > 0. and nclim > 0.:
                slopes_matrix[i,j] = 2
            elif liv > 0. and nclim < 0.:
                slopes_matrix[i,j] = 1
            elif liv < 0. and nclim > 0.:
                slopes_matrix[i,j] = -1
            elif liv < 0. and nclim < 0.:
                slopes_matrix[i,j] = -2



total_cells = np.count_nonzero(~np.isnan(slopes_matrix))
count1 = np.count_nonzero(slopes_matrix == 2.0)
count2 = np.count_nonzero(slopes_matrix == -2.0)
percent = ((count1+count2)/total_cells)*100.

contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125, -65, 24, 50], crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
cb = plt.contourf(lon, lat, slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
cbar.set_label(r'$\frac{Livneh}{nClimGrid}$',fontsize=12)
cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.85, 0.15, "PA: "+f"{percent:.2f}", transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
plt.title("prcptot slope direction",fontsize=12)
plt.savefig(image_dir+"obs_prcptot_slopes_directions_comparison_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
# =============================================================================
# =============================================================================
# # 
# =============================================================================
# =============================================================================

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']

letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

loca_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/ClimdexPhard_LOCA2-*_historical_r1i1p1f1_gn_CONUS16thD2021_1985-2014_ann.nc")
loca_files.sort()

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


all_loca2 = np.empty([16,30,444,922])
for i in range(16):
    all_loca2[i,:,:,:] = globals()['loca2_slope_'+str(models[i])+'_prcptot']

loca2_ensemble_slopes = np.nanmean(all_loca2,axis=(0,1))
loca2_ensemble_slopes_ne = loca2_ensemble_slopes[lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]



levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['slope_'+str(models[i])+'_prcptot'],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cb.ax.tick_params(labelsize=8)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
    cbar.set_label(r"$m_{\mathrm{TS}}$ (mm/year)",fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 Theil Slopes prcptot",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_prcptot_slopes.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


for i in range(16):
    globals()['LOCA2_'+str(models[i])+'_prcptot_slopes_percent_diff'] = ((globals()['slope_'+str(models[i])+'_prcptot']/livneh_slopes)-1)*100.

letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_prcptot_slopes_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("PiYG_r"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{Livneh}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 prcptot Slopes Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca2_all_models_slopes_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()



for h in range(16):
    globals()['LOCA2_'+str(models[h])+'_slopes_matrix'] = np.empty_like(livneh_slopes)
    for i in range(444):
        for j in range(922):
            if np.isnan(livneh_slopes[i,j]):
                globals()['LOCA2_'+str(models[h])+'_slopes_matrix'][i,j] = np.nan
            else:
                model = globals()['loca2_slope_'+str(models[h])+'_prcptot'][i,j]
                liv = livneh_slopes[i,j]
                if model > 0. and liv > 0.:
                    globals()['LOCA2_'+str(models[h])+'_slopes_matrix'][i,j] = 2
                elif model > 0. and liv < 0.:
                    globals()['LOCA2_'+str(models[h])+'_slopes_matrix'][i,j] = 1
                elif model < 0. and liv > 0.:
                    globals()['LOCA2_'+str(models[h])+'_slopes_matrix'][i,j] = -1
                elif model < 0. and liv < 0.:
                    globals()['LOCA2_'+str(models[h])+'_slopes_matrix'][i,j] = -2


percents = []
for i in range(16):
    total_area = np.count_nonzero(~np.isnan(globals()['LOCA2_'+str(models[i])+'_slopes_matrix']))
    count1 = np.count_nonzero(globals()['LOCA2_'+str(models[i])+'_slopes_matrix']==2.)
    count2 = np.count_nonzero(globals()['LOCA2_'+str(models[i])+'_slopes_matrix']== -2.)
    count3 = np.count_nonzero(globals()['LOCA2_'+str(models[i])+'_slopes_matrix']== -1.)
    count4 = np.count_nonzero(globals()['LOCA2_'+str(models[i])+'_slopes_matrix']== 1.)
    total = count1 + count2
    percent = (total/total_area)*100.
    percents.append(percent)

contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    if i in {8}:
        cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_slopes_matrix']+1e-10,contours,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm)
    else:
        cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_slopes_matrix'],contours,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
    cbar.set_label(r'$\frac{'+str(models[i])+'}{Livneh}$',fontsize=12)
    cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.82, 0.2, "PA: "+f"{percents[i]:.2f}", transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 prcptot Slopes Direction",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca2_all_models_slopes_directions.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


model_arrays = []
for i in range(16):
    model_arrays.append(globals()['loca2_slope_'+str(models[i])+'_prcptot'])
ensemble_average = np.nanmean(np.stack(model_arrays, axis=0), axis=0)

loca2_ensemble_slopes_matrix = np.empty_like(livneh_slopes)
for i in range(444):
    for j in range(922):
        if np.isnan(livneh_slopes[i,j]):
            loca2_ensemble_slopes_matrix[i,j] = np.nan
        else:
            model = ensemble_average[i,j]
            liv = livneh_slopes[i,j]
            if model > 0. and liv > 0.:
                loca2_ensemble_slopes_matrix[i,j] = 2
            elif model > 0. and liv < 0.:
                loca2_ensemble_slopes_matrix[i,j] = 1
            elif model < 0. and liv > 0.:
                loca2_ensemble_slopes_matrix[i,j] = -1
            elif model < 0. and liv < 0.:
                loca2_ensemble_slopes_matrix[i,j] = -2

total_area = np.count_nonzero(~np.isnan(loca2_ensemble_slopes_matrix))
count1 = np.count_nonzero(loca2_ensemble_slopes_matrix == 2.)
count2 = np.count_nonzero(loca2_ensemble_slopes_matrix == -2.)
loca2_ensemble_percent = ((count1+count2)/total_area)*100.


contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125, -65, 24, 50], crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
cb = plt.contourf(lon, lat, loca2_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
cbar.set_label(r'$\frac{LOCA2}{Livneh}$',fontsize=12)
cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.85, 0.15, "PA: "+f"{loca2_ensemble_percent:.2f}", transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
plt.title("LOCA2 Ensemble prcptot Slope Direction",fontsize=12)
plt.savefig(image_dir+"loca2_ensemble_prcptot_slopes_directions_comparison.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
# =============================================================================
# =============================================================================
# # 
# =============================================================================
# =============================================================================



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
    all_star[i,:,:,:] = globals()['star_slope_'+str(models[i])+'_prcptot']

for i in range(16):
    globals()['ne_star_slope_'+str(models[i])+'_prcptot'] = globals()['star_slope_'+str(models[i])+'_prcptot'][lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]
    globals()['ne_star_slope_'+str(models[i])+'_prcptot_cleaned']=globals()['ne_star_slope_'+str(models[i])+'_prcptot'][~np.isnan(globals()['ne_star_slope_'+str(models[i])+'_prcptot'])]
    
    globals()['ne_loca2_slope_'+str(models[i])+'_prcptot'] = globals()['loca2_slope_'+str(models[i])+'_prcptot'][lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]
    globals()['ne_loca2_slope_'+str(models[i])+'_prcptot_cleaned']=globals()['ne_loca2_slope_'+str(models[i])+'_prcptot'][~np.isnan(globals()['ne_loca2_slope_'+str(models[i])+'_prcptot'])]
    
    

star_ensemble_slopes = np.nanmean(all_star,axis=(0,1))
star_ensemble_slopes_ne = star_ensemble_slopes[lat_index1:(lat_index2+1),lon_index1:(lon_index2+1)]

northeast_slopes_cleaned_star = star_ensemble_slopes_ne[~np.isnan(star_ensemble_slopes_ne)]
northeast_slopes_cleaned_loca = loca2_ensemble_slopes_ne[~np.isnan(loca2_ensemble_slopes_ne)]

# =============================================================================
# 
# =============================================================================
all_trends = [northeast_slopes_cleaned_nclimgrid,northeast_slopes_cleaned_livneh,northeast_slopes_cleaned_prism,northeast_slopes_cleaned_loca,northeast_slopes_cleaned_star]

fig = plt.figure(figsize=(8,6))
plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
# Customize plot
plt.title("Distribution of Northeast prcptot Trends (Theil Slope)")
plt.ylabel(r"$m_{\mathrm{TS}}$ (mm/year)")
plt.xticks([1, 2, 3, 4, 5], ["nClimGrid", "Livneh", "PRISM","LOCA2", "STAR"])  # Adjust X-axis label to represent dataset
plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
plt.ylim(-1,3.5)
plt.savefig(image_dir+"boxplot_prcptot_slopes_obs_model_ensembles.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

all_trends = [northeast_slopes_cleaned_nclimgrid,northeast_slopes_cleaned_livneh,northeast_slopes_cleaned_prism]
titles = ["nClimGrid","Livneh","PRISM"]
for i in range(16):
    all_trends.append(globals()['ne_loca2_slope_'+str(models[i])+'_prcptot_cleaned'])
    titles.append(str(models[i]))

fig = plt.figure(figsize=(8,6))
plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
# Customize plot
plt.title("Distribution of Northeast prcptot Trends (Theil Slope)\nLOCA2 Models v. Obs.")
plt.ylabel(r"$m_{\mathrm{TS}}$ (mm/year)")
plt.xticks(np.arange(1, len(titles) + 1), titles, rotation=45, ha='right')  # Adjust X-axis label to represent dataset
plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
# plt.ylim(-1,3.5)
plt.savefig(image_dir+"boxplot_prcptot_slopes_obs_loca2_ind.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


all_trends = [northeast_slopes_cleaned_nclimgrid,northeast_slopes_cleaned_livneh,northeast_slopes_cleaned_prism]
titles = ["nClimGrid","Livneh","PRISM"]
for i in range(16):
    all_trends.append(globals()['ne_star_slope_'+str(models[i])+'_prcptot_cleaned'])
    titles.append(str(models[i]))

fig = plt.figure(figsize=(8,6))
plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
# Customize plot
plt.title("Distribution of Northeast prcptot Trends (Theil Slope)\nSTAR Models v. Obs.")
plt.ylabel(r"$m_{\mathrm{TS}}$ (mm/year)")
plt.xticks(np.arange(1, len(titles) + 1), titles, rotation=45, ha='right')  # Adjust X-axis label to represent dataset
plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
# plt.ylim(-1,3.5)
plt.savefig(image_dir+"boxplot_prcptot_slopes_obs_star_ind.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
# =============================================================================
# 
# =============================================================================

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['slope_'+str(models[i])+'_prcptot'],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cb.ax.tick_params(labelsize=8)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
    cbar.set_label(r"$m_{\mathrm{TS}}$ (mm/year)",fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR Theil Slopes prcptot",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_prcptot_slopes.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

for h in range(16):
    globals()['STAR_'+str(models[h])+'_slopes_matrix'] = np.empty_like(nclimgrid_slopes)
    for i in range(444):
        for j in range(922):
            if np.isnan(nclimgrid_slopes[i,j]):
                globals()['STAR_'+str(models[h])+'_slopes_matrix'][i,j] = np.nan
            else:
                model = globals()['star_slope_'+str(models[h])+'_prcptot'][i,j]
                nclim = nclimgrid_slopes[i,j]
                if model > 0. and nclim > 0.:
                    globals()['STAR_'+str(models[h])+'_slopes_matrix'][i,j] = 2
                elif model > 0. and nclim < 0.:
                    globals()['STAR_'+str(models[h])+'_slopes_matrix'][i,j] = 1
                elif model < 0. and nclim > 0.:
                    globals()['STAR_'+str(models[h])+'_slopes_matrix'][i,j] = -1
                elif model < 0. and nclim < 0.:
                    globals()['STAR_'+str(models[h])+'_slopes_matrix'][i,j] = -2



percents = []
for i in range(16):
    total_area = np.count_nonzero(~np.isnan(globals()['STAR_'+str(models[i])+'_slopes_matrix']))
    count1 = np.count_nonzero(globals()['STAR_'+str(models[i])+'_slopes_matrix']==2.)
    count2 = np.count_nonzero(globals()['STAR_'+str(models[i])+'_slopes_matrix']== -2.)
    total = count1 + count2
    percent = (total/total_area)*100.
    percents.append(percent)


contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['STAR_'+str(models[i])+'_slopes_matrix'],contours,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=12)
    cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.82, 0.2, "PA: "+f"{percents[i]:.2f}", transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR prcptot Slopes Direction",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_slopes_directions.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

model_arrays = []
for i in range(16):
    model_arrays.append(globals()['star_slope_'+str(models[i])+'_prcptot'])
ensemble_average_star = np.nanmean(np.stack(model_arrays, axis=0), axis=0)

star_ensemble_slopes_matrix = np.empty_like(nclimgrid_slopes)
for i in range(444):
    for j in range(922):
        if np.isnan(nclimgrid_slopes[i,j]):
            star_ensemble_slopes_matrix[i,j] = np.nan
        else:
            model = ensemble_average_star[i,j]
            nclim = nclimgrid_slopes[i,j]
            if model > 0. and nclim > 0.:
                star_ensemble_slopes_matrix[i,j] = 2
            elif model > 0. and nclim < 0.:
                star_ensemble_slopes_matrix[i,j] = 1
            elif model < 0. and nclim > 0.:
                star_ensemble_slopes_matrix[i,j] = -1
            elif model < 0. and nclim < 0.:
                star_ensemble_slopes_matrix[i,j] = -2

total_area = np.count_nonzero(~np.isnan(star_ensemble_slopes_matrix))
count1 = np.count_nonzero(star_ensemble_slopes_matrix == 2.)
count2 = np.count_nonzero(star_ensemble_slopes_matrix == -2.)
star_ensemble_percent = ((count1+count2)/total_area)*100.


contours = [-2.1, -1.1,0, 1.1, 2.1]
norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-125, -65, 24, 50], crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
cb = plt.contourf(lon, lat, star_ensemble_slopes_matrix, levels=contours,transform=ccrs.PlateCarree(), cmap=get_cmap("BrBG"), norm=norm)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
cbar.set_label(r'$\frac{STAR}{nClimGrid}$',fontsize=12)
cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.85, 0.15, "PA: "+f"{star_ensemble_percent:.2f}", transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
plt.title("STAR Ensemble prcptot Slope Direction",fontsize=12)
plt.savefig(image_dir+"star_ensemble_prcptot_slopes_directions_comparison.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()



for i in range(16):
    globals()[str(models[i])+"_ds_slope_matrix"] = np.empty_like(livneh_slopes)
    for j in range(444):
        for k in range(922):
            if np.isnan(livneh_slopes[j,k]):
                globals()[str(models[i])+"_ds_slope_matrix"][j,k] = np.nan
            else:
                loca = globals()['loca2_slope_'+str(models[i])+'_prcptot'][j,k]
                star = globals()['star_slope_'+str(models[i])+'_prcptot'][j,k]
                if loca > 0. and star > 0.:
                    globals()[str(models[i])+"_ds_slope_matrix"][j,k] = 2
                elif loca > 0. and star < 0.:
                    globals()[str(models[i])+"_ds_slope_matrix"][j,k] = 1
                elif loca < 0. and star > 0.:
                    globals()[str(models[i])+"_ds_slope_matrix"][j,k] = -1
                elif loca < 0. and star < 0.:
                    globals()[str(models[i])+"_ds_slope_matrix"][j,k] = -2


percent_ds_models = []
for i in range(16):
    total_area = np.count_nonzero(~np.isnan(globals()[str(models[i])+"_ds_slope_matrix"]))
    count1 = np.count_nonzero(globals()[str(models[i])+"_ds_slope_matrix"]==2.)
    count2 = np.count_nonzero(globals()[str(models[i])+"_ds_slope_matrix"]==-2.)
    percent = ((count1+count2)/total_area)*100.
    percent_ds_models.append(percent)


fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()[str(models[i])+"_ds_slope_matrix"],contours,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9, ticks=[-1.5, -0.5, 0.5, 1.5])
    cbar.set_label(r'$\frac{LOCA2}{STAR}$',fontsize=12)
    cbar.ax.set_xticklabels(["-/-", "-/+", "+/-", "+/+"], fontsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.82, 0.2, "PA: "+f"{percent_ds_models[i]:.2f}", transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2/STAR prcptot Slopes Direction",fontsize=20,y=0.93)
plt.savefig(image_dir+"all_models_slopes_directions_comparison.png",dpi=600,bbox_inches='tight')
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





states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(121,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
ax.set_title("Livneh")
cb = plt.contourf(lon,lat,livneh_masked_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (p<0.1)",fontsize=10)
cbar.ax.tick_params(labelsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(a)",transform=ax.transAxes)
ax = fig.add_subplot(122,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS, edgecolor='black')
ax.set_title("nClimGrid")
cb = plt.contourf(lon,lat,nclimgrid_masked_slopes,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
cb.ax.tick_params(labelsize=8)
cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
cbar.set_label(r"$m_{\mathrm{TS}}$ prcptot (p<0.1)",fontsize=10)
cbar.ax.tick_params(labelsize=8)
ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
ax.yaxis.tick_left()
ax.text(-0.06,0.95,"(b)",transform=ax.transAxes)
plt.savefig(image_dir+"obs_prcptot_slopes_significant_p01.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


# =============================================================================
# models significance
# =============================================================================

years = np.arange(1985,2015)
for i in range(16):
    globals()['loca2_'+str(models[i])+'_pvals'] = np.empty_like(globals()['loca2_slope_'+str(models[i])+'_prcptot'])
    globals()['loca2_'+str(models[i])+'_tau'] = np.empty_like(globals()['loca2_slope_'+str(models[i])+'_prcptot'])
    mask = np.isnan(globals()['loca2_slope_'+str(models[i])+'_prcptot'])
    globals()['loca2_'+str(models[i])+'_pvals'][mask] = np.nan
    globals()['loca2_'+str(models[i])+'_tau'][mask] = np.nan
    for j in range(globals()['loca2_slope_'+str(models[i])+'_prcptot'].shape[0]):
        for k in range(globals()['loca2_slope_'+str(models[i])+'_prcptot'].shape[1]):
            if np.isnan(globals()['loca2_slope_'+str(models[i])+'_prcptot'][j,k]):
                globals()['loca2_'+str(models[i])+'_pvals'][j,k] = np.nan
                globals()['loca2_'+str(models[i])+'_tau'][j,k] = np.nan
            else:
                time_series = globals()['LOCA2_'+str(models[i])+'_prcptot'][:,j,k]
                tau, p_value = kendalltau(years, time_series)
                globals()['loca2_'+str(models[i])+'_pvals'][j,k] = p_value
                globals()['loca2_'+str(models[i])+'_tau'][j,k] = tau
    globals()['loca2_'+str(models[i])+'_slopes_masked'] = np.where(globals()['loca2_'+str(models[i])+'_pvals']<0.1,globals()['loca2_slope_'+str(models[i])+'_prcptot'],np.nan)


years = np.arange(1985,2015)
for i in range(16):
    globals()['star_'+str(models[i])+'_pvals'] = np.empty_like(globals()['star_slope_'+str(models[i])+'_prcptot'])
    globals()['star_'+str(models[i])+'_tau'] = np.empty_like(globals()['star_slope_'+str(models[i])+'_prcptot'])
    mask = np.isnan(globals()['star_slope_'+str(models[i])+'_prcptot'])
    globals()['star_'+str(models[i])+'_pvals'][mask] = np.nan
    globals()['star_'+str(models[i])+'_tau'][mask] = np.nan
    for j in range(globals()['star_slope_'+str(models[i])+'_prcptot'].shape[0]):
        for k in range(globals()['star_slope_'+str(models[i])+'_prcptot'].shape[1]):
            if np.isnan(globals()['star_slope_'+str(models[i])+'_prcptot'][j,k]):
                globals()['star_'+str(models[i])+'_pvals'][j,k] = np.nan
                globals()['star_'+str(models[i])+'_tau'][j,k] = np.nan
            else:
                time_series = globals()['star_'+str(models[i])+'_prcptot'][:,j,k]
                tau, p_value = kendalltau(years, time_series)
                globals()['star_'+str(models[i])+'_pvals'][j,k] = p_value
                globals()['star_'+str(models[i])+'_tau'][j,k] = tau
    globals()['star_'+str(models[i])+'_slopes_masked'] = np.where(globals()['star_'+str(models[i])+'_pvals']<0.1,globals()['star_slope_'+str(models[i])+'_prcptot'],np.nan)
                

            


total_area = 444*922
number_of_non_nans = []
for i in range(16):
    number_of_non_nans.append(np.count_nonzero(~np.isnan(globals()['loca2_slope_'+str(models[i])+'_prcptot'])))

percent_significant = []
for i in range(16):
    count = np.count_nonzero((~np.isnan(globals()['star_'+str(models[i])+'_slopes_masked'])))
    percent = (count/number_of_non_nans[i])*100.
    percent_significant.append(percent)

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)
fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['loca2_'+str(models[i])+'_slopes_masked'],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cb.ax.tick_params(labelsize=8)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
    cbar.set_label(r"$m_{\mathrm{TS}}$ (p<0.1)",fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.84, 0.2, "p<0.1:\n"+f"{percent_significant[i]:.2f}"+"%", transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 Theil Slopes prcptot",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_prcptot_slopes_significant_p01.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

total_area = 444*922
number_of_non_nans = []
for i in range(16):
    number_of_non_nans.append(np.count_nonzero(~np.isnan(globals()['star_slope_'+str(models[i])+'_prcptot'])))

percent_significant = []
for i in range(16):
    count = np.count_nonzero((~np.isnan(globals()['star_'+str(models[i])+'_slopes_masked'])))
    percent = (count/number_of_non_nans[i])*100.
    percent_significant.append(percent)

levels = np.arange(-15,15.1,0.5)
norm = TwoSlopeNorm(vmin=-15, vcenter=0, vmax=15)
fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['star_'+str(models[i])+'_slopes_masked'],levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cb.ax.tick_params(labelsize=8)
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.9)
    cbar.set_label(r"$m_{\mathrm{TS}}$ (p<0.1)",fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.84, 0.2, "p<0.1:\n"+f"{percent_significant[i]:.2f}"+"%", transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR Theil Slopes prcptot",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_prcptot_slopes_significant_p01.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
