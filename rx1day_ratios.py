##this code looks at ratios (bias) between the observational products and the models

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import glob as glob

image_dir = '/Volumes/Elements/GFDL_project/images/'

states_provinces=cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']
# =============================================================================
# rx1dayAvg_annmax
# =============================================================================
star_data_rx1day = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/STAR/1985-2014/rx1day_EnsembleStats_ClimdexPhard-ann_STAR_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Max_ann.nc')
star_rx1dayMax_annmean= np.array(star_data_rx1day.variables['rx1dayMax_annmean'][0,:,:])
mask = np.isclose(star_rx1dayMax_annmean,1e+20)
star_rx1dayMax_annmean[mask] = np.nan
lon = np.array(star_data_rx1day.variables['lon'])
lat = np.array(star_data_rx1day.variables['lat'])
#
loca_data_rx1day = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/LOCA2/1985-2014/rx1day_EnsembleStats_ClimdexPhard-ann_LOCA2_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Max_ann.nc')
loca_rx1dayMax_annmean = np.array(loca_data_rx1day.variables['rx1dayMax_annmean'][0,:,:])
mask = np.isclose(loca_rx1dayMax_annmean,1e+20)
loca_rx1dayMax_annmean[mask] = np.nan
mask = np.isnan(star_rx1dayMax_annmean)
loca_rx1dayMax_annmean[mask] = np.nan
#
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
livneh_rx1dayMax_ann = np.array(livneh_data.variables['rx1dayMax_ann'][0,:,:])
mask = np.isclose(livneh_rx1dayMax_ann,1e+20)
livneh_rx1dayMax_ann[mask] = np.nan
mask = np.isnan(star_rx1dayMax_annmean)
livneh_rx1dayMax_ann[mask] = np.nan
#
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
nclimgrid_rx1dayMax_ann = np.array(nclimgrid_data.variables['rx1dayMax_ann'][0,:,:])
mask = np.isclose(nclimgrid_rx1dayMax_ann,1e+20)
nclimgrid_rx1dayMax_ann[mask] = np.nan
#

star_minus_nclimgrid = star_rx1dayMax_annmean/nclimgrid_rx1dayMax_ann
mask = np.isclose(star_minus_nclimgrid,1e+20)
star_minus_nclimgrid[mask] = np.nan
#
loca2_minus_livneh = loca_rx1dayMax_annmean/livneh_rx1dayMax_ann
mask = np.isclose(loca2_minus_livneh,1e+20)
loca2_minus_livneh[mask] = np.nan

obs_uncert = livneh_rx1dayMax_ann/nclimgrid_rx1dayMax_ann



# base_cmap = plt.get_cmap("bwr_r")

# # Create the custom colormap by combining parts
# colors = [
#     (0.0, base_cmap(0.0)),       # Start with the first color of reversed 'bwr' (blue)
#     (0.4, base_cmap(0.4)),       # First 40% is the first 40% of 'bwr_r'
#     (0.4, 'white'),              # Start of the white section
#     (0.6, 'white'),              # End of the white section (20%)
#     (0.6, base_cmap(0.6)),       # Start of the last 40% of 'bwr_r'
#     (1.0, base_cmap(1.0)),       # Last color of 'bwr_r' (red)
# ]

# # Create the colormap
# custom_cmap = LinearSegmentedColormap.from_list("custom_bwr_r", colors)

obs_uncert_percent_diff = (obs_uncert - 1) * 100
star_minus_nclimgrid_percent_diff = (star_minus_nclimgrid- 1) * 100
loca2_minus_livneh_percent_diff = (loca2_minus_livneh - 1) * 100

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
variables = [obs_uncert_percent_diff,star_minus_nclimgrid_percent_diff,loca2_minus_livneh_percent_diff]
titles = [r'$\frac{Livneh}{nClimGrid}$',r'$\frac{STAR}{nClimGrid}$',r'$\frac{LOCA2}{Livneh}$']
letters = ['(a)','(b)','(c)']

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]

fig = plt.figure(figsize=(6,14))
for i in range(3):
    ax = fig.add_subplot(3,1,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(titles[i],fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    if i in{0}:
        plt.title("rx1dayMax_annmean Percent Difference",fontsize=15)
plt.savefig(image_dir+"rx1dayMax_annmean_ratios_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


star_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-*CONUS16thD2021_1985-2014_30yr.nc')
star_files.sort()

loca_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_LOCA2-*CONUS16thD2021_1985-2014_30yr.nc')
loca_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['STAR_'+str(models[i])+'_rx1dayMax_ann'] = np.array(star_data.variables['rx1dayMax_ann'][0,:,:])
    mask = np.isclose(globals()['STAR_'+str(models[i])+'_rx1dayMax_ann'],1e+20)
    globals()['STAR_'+str(models[i])+'_rx1dayMax_ann'][mask] = np.nan
    #
    loca_data = nc.Dataset(loca_files[i])
    globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann'] = np.array(loca_data.variables['rx1dayMax_ann'][0,:,:])
    mask = np.isclose(globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann'],1e+20)
    globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann'][mask] = np.nan
    mask = np.isnan(globals()['STAR_'+str(models[i])+'_rx1dayMax_ann'])
    globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann'][mask] = np.nan


for i in range(16):
    globals()['STAR_'+str(models[i])+'_rx1dayMax_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_rx1dayMax_ann']/nclimgrid_rx1dayMax_ann)-1)*100.
    globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann_percent_diff'] = ((globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann']/livneh_rx1dayMax_ann)-1)*100.

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
    cb = plt.contourf(lon,lat,globals()['STAR_'+str(models[i])+'_rx1dayMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR rx1dayMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_rx1dayMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 rx1dayMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_rx1dayMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


for i in range(16):
    globals()['model_'+str(models[i])+'_rx1dayMax_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_rx1dayMax_ann'])/(globals()['LOCA2_'+str(models[i])+'_rx1dayMax_ann'])-1)*100.

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['model_'+str(models[i])+'_rx1dayMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{STAR}{LOCA2}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR/LOCA2 rx1dayMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"ds_all_models_rx1dayMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


# =============================================================================
# rx1dayAvg_annmean
# =============================================================================
star_data_rx1day = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/STAR/1985-2014/rx1day_EnsembleStats_ClimdexPhard-ann_STAR_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc')
star_rx1dayAvg_annmean= np.array(star_data_rx1day.variables['rx1dayAvg_annmean'][0,:,:])
mask = np.isclose(star_rx1dayAvg_annmean,1e+20)
star_rx1dayAvg_annmean[mask] = np.nan
lon = np.array(star_data_rx1day.variables['lon'])
lat = np.array(star_data_rx1day.variables['lat'])
#
loca_data_rx1day = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/LOCA2/1985-2014/rx1day_EnsembleStats_ClimdexPhard-ann_LOCA2_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc')
loca_rx1dayAvg_annmean = np.array(loca_data_rx1day.variables['rx1dayAvg_annmean'][0,:,:])
mask = np.isclose(loca_rx1dayAvg_annmean,1e+20)
loca_rx1dayAvg_annmean[mask] = np.nan
mask = np.isnan(star_rx1dayAvg_annmean)
loca_rx1dayAvg_annmean[mask] = np.nan
#
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
livneh_rx1dayAvg_ann = np.array(livneh_data.variables['rx1dayAvg_ann'][0,:,:])
mask = np.isclose(livneh_rx1dayAvg_ann,1e+20)
livneh_rx1dayAvg_ann[mask] = np.nan
mask = np.isnan(star_rx1dayAvg_annmean)
livneh_rx1dayAvg_ann[mask] = np.nan
#
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
nclimgrid_rx1dayAvg_ann = np.array(nclimgrid_data.variables['rx1dayAvg_ann'][0,:,:])
mask = np.isclose(nclimgrid_rx1dayAvg_ann,1e+20)
nclimgrid_rx1dayAvg_ann[mask] = np.nan
#

star_minus_nclimgrid = star_rx1dayAvg_annmean/nclimgrid_rx1dayAvg_ann
mask = np.isclose(star_minus_nclimgrid,1e+20)
star_minus_nclimgrid[mask] = np.nan
#
loca2_minus_livneh = loca_rx1dayAvg_annmean/livneh_rx1dayAvg_ann
mask = np.isclose(loca2_minus_livneh,1e+20)
loca2_minus_livneh[mask] = np.nan

obs_uncert = livneh_rx1dayAvg_ann/nclimgrid_rx1dayAvg_ann

obs_uncert_percent_diff = (obs_uncert - 1) * 100
star_minus_nclimgrid_percent_diff = (star_minus_nclimgrid- 1) * 100
loca2_minus_livneh_percent_diff = (loca2_minus_livneh - 1) * 100

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
variables = [obs_uncert_percent_diff,star_minus_nclimgrid_percent_diff,loca2_minus_livneh_percent_diff]
titles = [r'$\frac{Livneh}{nClimGrid}$',r'$\frac{STAR}{nClimGrid}$',r'$\frac{LOCA2}{Livneh}$']
letters = ['(a)','(b)','(c)']

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]

fig = plt.figure(figsize=(6,14))
for i in range(3):
    ax = fig.add_subplot(3,1,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(titles[i],fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    if i in{0}:
        plt.title("rx1dayAvg_annmean Percent Difference",fontsize=15)
plt.savefig(image_dir+"rx1dayAvg_annmean_ratios_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


star_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-*CONUS16thD2021_1985-2014_30yr.nc')
star_files.sort()

loca_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/rx1day_CalcTimeStats1_ClimdexPhard-ann_LOCA2-*CONUS16thD2021_1985-2014_30yr.nc')
loca_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann'] = np.array(star_data.variables['rx1dayAvg_ann'][0,:,:])
    mask = np.isclose(globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann'],1e+20)
    globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann'][mask] = np.nan
    #
    loca_data = nc.Dataset(loca_files[i])
    globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann'] = np.array(loca_data.variables['rx1dayAvg_ann'][0,:,:])
    mask = np.isclose(globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann'],1e+20)
    globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann'][mask] = np.nan
    mask = np.isnan(globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann'])
    globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann'][mask] = np.nan


for i in range(16):
    globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann']/nclimgrid_rx1dayAvg_ann)-1)*100.
    globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'] = ((globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann']/livneh_rx1dayAvg_ann)-1)*100.

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
    cb = plt.contourf(lon,lat,globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR rx1dayAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_rx1dayAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 rx1dayAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_rx1dayAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

for i in range(16):
    globals()['model_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_rx1dayAvg_ann'])/(globals()['LOCA2_'+str(models[i])+'_rx1dayAvg_ann'])-1)*100.

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['model_'+str(models[i])+'_rx1dayAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{STAR}{LOCA2}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR/LOCA2 rx1dayAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"ds_all_models_rx1dayAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

# =============================================================================
# prcptot_annmean
# =============================================================================
star_data_prcptot = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/STAR/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_STAR_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc')
star_prcptotAvg_annmean= np.array(star_data_prcptot.variables['prcptotAvg_annmean'][0,:,:])
mask = np.isclose(star_prcptotAvg_annmean,1e+20)
star_prcptotAvg_annmean[mask] = np.nan
#
loca_data_prcptot = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/LOCA2/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_LOCA2_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc')
loca_prcptotAvg_annmean = np.array(loca_data_prcptot.variables['prcptotAvg_annmean'][0,:,:])
mask = np.isclose(loca_prcptotAvg_annmean,1e+20)
loca_prcptotAvg_annmean[mask] = np.nan
mask = np.isnan(star_prcptotAvg_annmean)
loca_prcptotAvg_annmean[mask] = np.nan
#
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
livneh_prcptotAvg_ann = np.array(livneh_data.variables['prcptotAvg_ann'][0,:,:])
mask = np.isclose(livneh_prcptotAvg_ann,1e+20)
livneh_prcptotAvg_ann[mask] = np.nan
mask = np.isnan(star_prcptotAvg_annmean)
livneh_prcptotAvg_ann[mask] = np.nan
#
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
nclimgrid_prcptotAvg_ann = np.array(nclimgrid_data.variables['prcptotAvg_ann'][0,:,:])
mask = np.isclose(nclimgrid_prcptotAvg_ann,1e+20)
nclimgrid_prcptotAvg_ann[mask] = np.nan
#
star_minus_nclimgrid = star_prcptotAvg_annmean/nclimgrid_prcptotAvg_ann
mask = np.isclose(star_minus_nclimgrid,1e+20)
star_minus_nclimgrid[mask] = np.nan
#
loca2_minus_livneh = loca_prcptotAvg_annmean/livneh_prcptotAvg_ann
mask = np.isclose(loca2_minus_livneh,1e+20)
loca2_minus_livneh[mask] = np.nan

obs_uncert = livneh_prcptotAvg_ann/nclimgrid_prcptotAvg_ann

obs_uncert_percent_diff = (obs_uncert - 1) * 100
star_minus_nclimgrid_percent_diff = (star_minus_nclimgrid- 1) * 100
loca2_minus_livneh_percent_diff = (loca2_minus_livneh - 1) * 100

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
variables = [obs_uncert_percent_diff,star_minus_nclimgrid_percent_diff,loca2_minus_livneh_percent_diff]
titles = [r'$\frac{Livneh}{nClimGrid}$',r'$\frac{STAR}{nClimGrid}$',r'$\frac{LOCA2}{Livneh}$']
letters = ['(a)','(b)','(c)']

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]

fig = plt.figure(figsize=(6,14))
for i in range(3):
    ax = fig.add_subplot(3,1,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap('BrBG'),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(titles[i],fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    if i in{0}:
        plt.title("prcptotAvg_annmean Percent Difference",fontsize=15)
plt.savefig(image_dir+"prcptotAvg_annmean_ratios_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

star_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-*CONUS16thD2021_1985-2014_30yr.nc')
star_files.sort()

loca_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_LOCA2-*CONUS16thD2021_1985-2014_30yr.nc')
loca_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['STAR_'+str(models[i])+'_prcptotAvg_ann'] = np.array(star_data.variables['prcptotAvg_ann'][0,:,:])
    mask = np.isclose(globals()['STAR_'+str(models[i])+'_prcptotAvg_ann'],1e+20)
    globals()['STAR_'+str(models[i])+'_prcptotAvg_ann'][mask] = np.nan
    #
    loca_data = nc.Dataset(loca_files[i])
    globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann'] = np.array(loca_data.variables['prcptotAvg_ann'][0,:,:])
    mask = np.isclose(globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann'],1e+20)
    globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann'][mask] = np.nan
    mask = np.isnan(globals()['STAR_'+str(models[i])+'_prcptotAvg_ann'])
    globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann'][mask] = np.nan


for i in range(16):
    globals()['STAR_'+str(models[i])+'_prcptotAvg_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_prcptotAvg_ann']/nclimgrid_prcptotAvg_ann)-1)*100.
    globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann_percent_diff'] = ((globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann']/livneh_prcptotAvg_ann)-1)*100.

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
    cb = plt.contourf(lon,lat,globals()['STAR_'+str(models[i])+'_prcptotAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR prcptotAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_prcptotAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 prcptotAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_prcptotAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

for i in range(16):
    globals()['model_'+str(models[i])+'_prcptotAvg_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_prcptotAvg_ann'])/(globals()['LOCA2_'+str(models[i])+'_prcptotAvg_ann'])-1)*100.

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['model_'+str(models[i])+'_prcptotAvg_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{STAR}{LOCA2}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR/LOCA2 prcptotAvg_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"ds_all_models_prcptotAvg_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

# =============================================================================
# prcptot_annmax
# =============================================================================
star_data_prcptot = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/STAR/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_STAR_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Max_ann.nc')
star_prcptotMax_annmean= np.array(star_data_prcptot.variables['prcptotMax_annmean'][0,:,:])
mask = np.isclose(star_prcptotMax_annmean,1e+20)
star_prcptotMax_annmean[mask] = np.nan
#
loca_data_prcptot = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/LOCA2/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_LOCA2_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Max_ann.nc')
loca_prcptotMax_annmean = np.array(loca_data_prcptot.variables['prcptotMax_annmean'][0,:,:])
mask = np.isclose(loca_prcptotMax_annmean,1e+20)
loca_prcptotMax_annmean[mask] = np.nan
mask = np.isnan(star_prcptotMax_annmean)
loca_prcptotMax_annmean[mask] = np.nan
#
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
livneh_prcptotMax_ann = np.array(livneh_data.variables['prcptotMax_ann'][0,:,:])
mask = np.isclose(livneh_prcptotMax_ann,1e+20)
livneh_prcptotMax_ann[mask] = np.nan
mask = np.isnan(star_prcptotMax_annmean)
livneh_prcptotMax_ann[mask] = np.nan
#
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_30yr.nc')
nclimgrid_prcptotMax_ann = np.array(nclimgrid_data.variables['prcptotMax_ann'][0,:,:])
mask = np.isclose(nclimgrid_prcptotMax_ann,1e+20)
nclimgrid_prcptotMax_ann[mask] = np.nan
#
star_minus_nclimgrid = star_prcptotMax_annmean/nclimgrid_prcptotMax_ann
mask = np.isclose(star_minus_nclimgrid,1e+20)
star_minus_nclimgrid[mask] = np.nan
#
loca2_minus_livneh = loca_prcptotMax_annmean/livneh_prcptotMax_ann
mask = np.isclose(loca2_minus_livneh,1e+20)
loca2_minus_livneh[mask] = np.nan

obs_uncert = livneh_prcptotMax_ann/nclimgrid_prcptotMax_ann

obs_uncert_percent_diff = (obs_uncert - 1) * 100
star_minus_nclimgrid_percent_diff = (star_minus_nclimgrid- 1) * 100
loca2_minus_livneh_percent_diff = (loca2_minus_livneh - 1) * 100

variables = [obs_uncert_percent_diff,star_minus_nclimgrid_percent_diff,loca2_minus_livneh_percent_diff]
titles = [r'$\frac{Livneh}{nClimGrid}$',r'$\frac{STAR}{nClimGrid}$',r'$\frac{LOCA2}{Livneh}$']
letters = ['(a)','(b)','(c)']

# variables = [obs_uncert,star_minus_nclimgrid,loca2_minus_livneh]
norm = TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60)
contour_levels = [-60, -40, -20, -10, -5, 5, 10, 20, 40, 60]

fig = plt.figure(figsize=(6,14))
for i in range(3):
    ax = fig.add_subplot(3,1,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,variables[i],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(titles[i],fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    if i in{0}:
        plt.title("prcptotMax_annmean Percent Difference",fontsize=15)
plt.savefig(image_dir+"prcptotMax_annmean_ratios_v2.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

star_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-*CONUS16thD2021_1985-2014_30yr.nc')
star_files.sort()

loca_files = glob.glob('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/prcptot_CalcTimeStats1_ClimdexPhard-ann_LOCA2-*CONUS16thD2021_1985-2014_30yr.nc')
loca_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['STAR_'+str(models[i])+'_prcptotMax_ann'] = np.array(star_data.variables['prcptotMax_ann'][0,:,:])
    mask = np.isclose(globals()['STAR_'+str(models[i])+'_prcptotMax_ann'],1e+20)
    globals()['STAR_'+str(models[i])+'_prcptotMax_ann'][mask] = np.nan
    #
    loca_data = nc.Dataset(loca_files[i])
    globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann'] = np.array(loca_data.variables['prcptotMax_ann'][0,:,:])
    mask = np.isclose(globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann'],1e+20)
    globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann'][mask] = np.nan
    mask = np.isnan(globals()['STAR_'+str(models[i])+'_prcptotMax_ann'])
    globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann'][mask] = np.nan


for i in range(16):
    globals()['STAR_'+str(models[i])+'_prcptotMax_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_prcptotMax_ann']/nclimgrid_prcptotMax_ann)-1)*100.
    globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann_percent_diff'] = ((globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann']/livneh_prcptotMax_ann)-1)*100.

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
    cb = plt.contourf(lon,lat,globals()['STAR_'+str(models[i])+'_prcptotMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR prcptotMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"star_all_models_prcptotMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{'+str(models[i])+'}{nClimGrid}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("LOCA2 prcptotMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"loca_all_models_prcptotMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()


for i in range(16):
    globals()['model_'+str(models[i])+'_prcptotMax_ann_percent_diff'] = ((globals()['STAR_'+str(models[i])+'_prcptotMax_ann'])/(globals()['LOCA2_'+str(models[i])+'_prcptotMax_ann'])-1)*100.

fig = plt.figure(figsize=(20,16))
for i in range(16):
    ax = fig.add_subplot(4,4,(i+1),projection=ccrs.Mercator())
    ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, edgecolor='black')
    cb = plt.contourf(lon,lat,globals()['model_'+str(models[i])+'_prcptotMax_ann_percent_diff'],contour_levels,transform=ccrs.PlateCarree(),cmap=get_cmap("BrBG"),norm=norm,extend='both')
    cbar = fig.colorbar(cb, orientation='horizontal', shrink=0.7)
    cbar.set_label(r'$\frac{STAR}{LOCA2}$',fontsize=15)
    ax.set_xticks([-120.,-115., -110.,-105., -100.,-95.,-90.,-85.,-80.,-75.,-70.], crs=ccrs.PlateCarree())
    ax.set_xticklabels(["120°W","", "110°W","", "100°W","","90°W","","80°W","","70°W"],fontsize=10)
    ax.set_yticks([26.,28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.], crs=ccrs.PlateCarree())
    ax.set_yticklabels(["26°N","","30°N","","34°N","","38°N","","42°N","","46°N",""],fontsize=10)
    ax.yaxis.tick_left()
    ax.text(-0.06,0.95,letters[i],transform=ax.transAxes)
    plt.title(str(models[i]),fontsize=12)
plt.subplots_adjust(hspace=0.3)
plt.suptitle("STAR/LOCA2 prcptotMax_ann Percent Difference",fontsize=20,y=0.93)
plt.savefig(image_dir+"ds_all_models_prcptotMax_percent_diff.png",dpi=600,bbox_inches='tight')
plt.show()
plt.close()
