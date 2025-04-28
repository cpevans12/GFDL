#packages
import netCDF4 as nc
import numpy as np
import glob as glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import get_cmap
import geopandas as gpd

#this is where I've been saving images for the NPS site
image_dir = '/Volumes/Elements/GFDL_project/images/new_images/'
# =============================================================================
# 
# =============================================================================
#this is a shapefile from: https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
shape_file = "/Users/cpe28/Downloads/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp"
us_states = gpd.read_file(shape_file)

# Define the NCA Region states
northeast_states = [
    'Maine', 'New Hampshire', 'Vermont', 'Massachusetts',
    'Rhode Island', 'Connecticut', 'New York', 'Pennsylvania', 'New Jersey','Maryland','West Virginia','Delaware']
# Create a GeoDataFrame for U.S. states (from cartopy)
northeast_us_states = us_states[us_states['name'].isin(northeast_states)]
#
southeast_states = ['Virginia','Kentucky',
                  'North Carolina', 'South Carolina', 'Georgia', 'Florida', 
                  'Alabama', 'Mississippi', 'Tennessee','Louisiana','Arkansas']
southeast_us_states = us_states[us_states['name'].isin(southeast_states)]
#
midwest_states = ['Ohio','Michigan','Indiana','Illinois','Missouri','Iowa','Wisconsin','Minnesota']
midwest_us_states = us_states[us_states['name'].isin(midwest_states)]
#
south_plains_states = ['Texas','Oklahoma','Kansas']
south_plains_us_states = us_states[us_states['name'].isin(south_plains_states)]
#
north_plains_states = ['Montana','North Dakota','South Dakota','Nebraska','Wyoming']
north_plains_us_states = us_states[us_states['name'].isin(north_plains_states)]
#
southwest_states = ['Colorado','New Mexico','Arizona','Utah','Nevada','California']
southwest_us_states = us_states[us_states['name'].isin(southwest_states)]
#
northwest_states = ['Idaho','Oregon','Washington']
northwest_us_states = us_states[us_states['name'].isin(northwest_states)]

#using the Latitue and Longitude from ncligmrid to define the masks and regions
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
lat = np.array(nclimgrid_data.variables['lat'])
lon = np.array(nclimgrid_data.variables['lon'])

lon_adjust = lon-360.
lon_grid, lat_grid = np.meshgrid(lon_adjust, lat)

regions = {
    "Northeast": [
        'Maine', 'New Hampshire', 'Vermont', 'Massachusetts',
        'Rhode Island', 'Connecticut', 'New York', 'Pennsylvania', 
        'New Jersey', 'Maryland', 'West Virginia', 'Delaware'
    ],
    "Southeast": [
        'Virginia', 'Kentucky', 'North Carolina', 'South Carolina', 
        'Georgia', 'Florida', 'Alabama', 'Mississippi', 
        'Tennessee', 'Louisiana', 'Arkansas'
    ],
    "Midwest": [
        'Ohio', 'Michigan', 'Indiana', 'Illinois', 'Missouri', 
        'Iowa', 'Wisconsin', 'Minnesota'
    ],
    "South Great Plains": [
        'Texas', 'Oklahoma', 'Kansas'
    ],
    "North Great Plains": [
        'Montana', 'North Dakota', 'South Dakota', 'Nebraska', 
        'Wyoming'
    ],
    "Southwest": [
        'Colorado', 'New Mexico', 'Arizona', 'Utah', 'Nevada', 
        'California'
    ],
    "Northwest": [
        'Idaho', 'Oregon', 'Washington'
    ]
}

# Create a dictionary to hold the extents
extents = {}
for region, states in regions.items():
    region_states = us_states[us_states['name'].isin(states)]
    region_states = region_states[~region_states['type'].str.contains("County")]
    region_states = region_states[~region_states['type'].str.contains("Departamento")]
    region_states = region_states[~region_states['type'].str.contains("Oblast")]
    
    # Check if the region has any states
    if not region_states.empty:
        min_lon, min_lat, max_lon, max_lat = region_states.total_bounds
        extents[region] = [min_lon-0.5, max_lon+0.5, min_lat-0.5, max_lat+0.5]
    else:
        extents[region] = None  # No states found for this region

# Display the extents
for region, extent in extents.items():
    if extent:
        print(f"{region} extent: {extent}")
    else:
        print(f"{region} has no valid states.")

#this function was my first attempt at masking the regions - it worked but was very slow
# def mask_array(array, states_gdf):
#     masked_array = np.full(array.shape, np.nan)
#     points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon_grid.flatten(), lat_grid.flatten()))
#     for _, state in states_gdf.iterrows():
#         state_geom = state['geometry']
#         # Find points within the current state geometry
#         points_within_state = points[points.geometry.within(state_geom)]        
#         # Update the masked_array
#         for idx in points_within_state.index:
#             masked_array[np.unravel_index(idx, array.shape)] = array.flatten()[idx]
#     return masked_array

#built these two functions to apply the masking to the data and it works much faster
def make_region_mask(states_gdf, lat_grid, lon_grid):
    points = gpd.GeoSeries(gpd.points_from_xy(lon_grid.ravel(), lat_grid.ravel()))
    mask = np.zeros(lat_grid.size, dtype=bool)
    for _, state in states_gdf.iterrows():
        geom = state.geometry
        mask |= points.within(geom).values
    return mask.reshape(lat_grid.shape)

def apply_region_mask(array, region_mask):
    return np.where(region_mask, array, np.nan)
# =============================================================================
# 
# =============================================================================
nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
nclimgrid_prcptot = np.array(nclimgrid_data.variables['prcptot'])
#masking missing values
mask = np.isclose(nclimgrid_prcptot,1e+20)
nclimgrid_prcptot[mask] = np.nan
lat = np.array(nclimgrid_data.variables['lat'])
lon = np.array(nclimgrid_data.variables['lon'])

#30 year average
nclimgrid_prcptot_30yr_avg = np.nanmean(nclimgrid_prcptot,axis=0)
#this was hardcoded based on min and max of prcptot - can be adjusted based on variable
levels = np.arange(0,2001,10)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,nclimgrid_prcptot_30yr_avg,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("jet_r"),extend='both')
plt.savefig(image_dir+'nclimgrid_prcptot_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()


livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/ClimdexPhard_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
lat = np.array(livneh_data.variables['lat'])
lon = np.array(livneh_data.variables['lon'])
livneh_prcptot = np.array(livneh_data.variables['prcptot'])
mask = np.isclose(livneh_prcptot,1e+20)
livneh_prcptot[mask] = np.nan
#additional mask because Livneh's extent is larger than ncligmrid - so making consistent for plotting
mask = np.isnan(nclimgrid_prcptot)
livneh_prcptot[mask] = np.nan

star_ens_data = nc.Dataset("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/STAR/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_STAR_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc")
star_ens_prcptot = np.array(star_ens_data.variables['prcptotAvg_annmean'][0,:,:])
mask = np.isclose(star_ens_prcptot,1e+20)
star_ens_prcptot[mask] = np.nan


fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,star_ens_prcptot,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("jet_r"),extend='both')
plt.savefig(image_dir+'star_ens_prcptot_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

#simple bias calculation between STAR and nclimgrid
bias = star_ens_prcptot/nclimgrid_prcptot_30yr_avg
bias_levels = np.arange(0.78,1.205,0.005)

base_cmap = plt.get_cmap("bwr_r")
# Create the custom colormap by combining parts
colors = [
    (0.0, base_cmap(0.0)),       # Start with the first color of reversed 'bwr' (blue)
    (0.4, base_cmap(0.4)),       # First 40% is the first 40% of 'bwr_r'
    (0.4, 'white'),              # Start of the white section
    (0.6, 'white'),              # End of the white section (20%)
    (0.6, base_cmap(0.6)),       # Start of the last 40% of 'bwr_r'
    (1.0, base_cmap(1.0)),       # Last color of 'bwr_r' (red)
]
# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_bwr_r", colors)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,bias,bias_levels,transform=ccrs.PlateCarree(),cmap=custom_cmap,extend='both')
plt.savefig(image_dir+'nclimgrid_star_prcptot_historical_bias_blank_v2.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

# =============================================================================
# Time Series
# =============================================================================
region_names = ['Northeast','Southeast','Midwest','North Great Plains','South Great Plains','Southwest','Northwest']
regions_names_alt=['northeast','southeast','midwest','north_plains','south_plains','southwest','northwest']

#this sets up the region masking
region_masks = {}
for i, region in enumerate(regions_names_alt):
    region_masks[region] = make_region_mask(globals()[region+'_us_states'], lat_grid, lon_grid)

#applies the regional masking
for i in range(7):
    key = f'nclimgrid_{regions_names_alt[i]}'
    data = nclimgrid_prcptot
    mask = region_masks[regions_names_alt[i]]
    globals()[key] = apply_region_mask(data, mask)

#time series averaged over the NCA region
for i in range(7):
    globals()['nclimgrid_'+str(regions_names_alt[i])+'_timeseries'] = np.nanmean(globals()['nclimgrid_'+str(regions_names_alt[i])],axis=(1,2))
    
for i in range(7):
    key = f'livneh_{regions_names_alt[i]}'
    data = livneh_prcptot
    mask = region_masks[regions_names_alt[i]]
    globals()[key] = apply_region_mask(data, mask)

for i in range(7):
    globals()['livneh_'+str(regions_names_alt[i])+'_timeseries'] = np.nanmean(globals()['livneh_'+str(regions_names_alt[i])],axis=(1,2))

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']

#previous STAR data was the already calculated ensemble average, now we are getting the individual models data
star_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/ClimdexPhard_STAR-*CONUS16thD2021_1985-2014_ann.nc")
star_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['star_'+str(models[i])+'_prcptot'] = np.array(star_data.variables['prcptot'])
    mask = np.isclose(globals()['star_'+str(models[i])+'_prcptot'],1e+20)
    globals()['star_'+str(models[i])+'_prcptot'][mask] = np.nan
    mask = np.isnan(nclimgrid_prcptot)
    globals()['star_'+str(models[i])+'_prcptot'][mask]= np.nan
    
for i in range(7):
    for j in range(16):
        key = f'star_{models[j]}{regions_names_alt[i]}'
        data = globals()[f'star_{models[j]}_prcptot']
        mask = region_masks[regions_names_alt[i]]
        globals()[key] = apply_region_mask(data, mask)

for i in range(7):
    for j in range(16):
        globals()['star_'+str(models[j])+str(regions_names_alt[i])+'_timeseries'] = np.nanmean(globals()['star_'+str(models[j])+str(regions_names_alt[i])],axis=(1,2))

#probably an easier way to do this, tbh
ensemble = np.empty([7,16,30])
for i in range(7):
    for j in range(16):
        ensemble[i,j,:] = globals()['star_'+str(models[j])+str(regions_names_alt[i])+'_timeseries']

ensemble_avg = np.empty([7,30])
for i in range(7):
    ensemble_avg[i,:] = np.nanmean(ensemble[i,:,:],axis=0)

#plot the time series 
for h in range(7):
    fig = plt.figure(figsize=(14,6))
    for i in range(16):
        if i in {0}:
            plt.plot(np.arange(1985,2015,1),globals()['star_'+str(models[i])+str(regions_names_alt[h])+'_timeseries'],color='grey',ls='--',alpha=0.75,label='STAR',lw=1)
        else:
            plt.plot(np.arange(1985,2015,1),globals()['star_'+str(models[i])+str(regions_names_alt[h])+'_timeseries'],color='grey',ls='--',alpha=0.75,lw=1)
    plt.plot(np.arange(1985,2015,1),ensemble_avg[h,:],lw=3,label='Ens. Avg.',ls='--',color='k')
    plt.plot(np.arange(1985,2015,1),globals()['livneh_'+str(regions_names_alt[h])+'_timeseries'],lw=3,color='tab:orange',label='Livneh')
    plt.plot(np.arange(1985,2015,1),globals()['nclimgrid_'+str(regions_names_alt[h])+'_timeseries'],lw=3,label='nClimGrid')
    plt.legend()
    plt.xlim(1984.5,2014.5)
    plt.ylabel("Avg. Annual Total Precip. (mm)")
    plt.xticks(size=12)
    plt.savefig(image_dir+'prcptot_time_series_obs_STAR_'+str(regions_names_alt[h])+'.png',dpi=300,bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close()

# =============================================================================
# saving an image of just the colorbars
# =============================================================================
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
bias_levels = np.arange(0.78,1.205,0.005)
norm = Normalize(vmin=bias_levels.min(), vmax=bias_levels.max())

base_cmap = plt.get_cmap("bwr_r")

# Create the custom colormap by combining parts
colors = [
    (0.0, base_cmap(0.0)),       # Start with the first color of reversed 'bwr' (blue)
    (0.4, base_cmap(0.4)),       # First 40% is the first 40% of 'bwr_r'
    (0.4, 'white'),              # Start of the white section
    (0.6, 'white'),              # End of the white section (20%)
    (0.6, base_cmap(0.6)),       # Start of the last 40% of 'bwr_r'
    (1.0, base_cmap(1.0)),       # Last color of 'bwr_r' (red)
]

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_bwr_r", colors)
fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(bottom=0.5)
cb = ColorbarBase(ax, cmap=custom_cmap, norm=norm, orientation='horizontal', extend='both')
cb.set_label(r'$\frac{Model}{Observation}$')
plt.savefig(image_dir+"bias_map_colorbar.png",dpi=300,bbox_inches='tight')
plt.show()
plt.close()

levels = np.arange(10,2001,10)
norm = Normalize(vmin=levels.min(), vmax=levels.max())

fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(bottom=0.5)
cb = ColorbarBase(ax, cmap=get_cmap('jet_r'), norm=norm, orientation='horizontal', extend='both')
cb.set_label('Precipitation (mm)')
plt.savefig(image_dir+"prcptot_map_colorbar.png",dpi=300,bbox_inches='tight')
plt.show()
plt.close()

# =============================================================================
# Now LOCA2 and Livneh
# =============================================================================
livneh_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/Livneh/1985-2014/ClimdexPhard_livneh_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
lat = np.array(livneh_data.variables['lat'])
lon = np.array(livneh_data.variables['lon'])
livneh_prcptot = np.array(livneh_data.variables['prcptot'])
mask = np.isclose(livneh_prcptot,1e+20)
livneh_prcptot[mask] = np.nan
mask = np.isnan(nclimgrid_prcptot)
livneh_prcptot[mask] = np.nan

livneh_prcptot_30yr_avg = np.nanmean(livneh_prcptot,axis=0)
levels = np.arange(0,2001,10)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,livneh_prcptot_30yr_avg,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("jet_r"),extend='both')
plt.savefig(image_dir+'livneh_prcptot_historical_mean_blank_v2.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()

loca_ens_data = nc.Dataset("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Ensembles/LOCA2/1985-2014/prcptot_EnsembleStats_ClimdexPhard-ann_LOCA2_NCA16_historical_CONUS16thD2021_1985-2014_30yr_Avg_ann.nc")
loca_ens_prcptot = np.array(loca_ens_data.variables['prcptotAvg_annmean'][0,:,:])
mask = np.isclose(loca_ens_prcptot,1e+20)
loca_ens_prcptot[mask] = np.nan
mask = np.isnan(nclimgrid_prcptot[0,:,:])
loca_ens_prcptot[mask] = np.nan

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,loca_ens_prcptot,levels,transform=ccrs.PlateCarree(),cmap=get_cmap("jet_r"),extend='both')
plt.savefig(image_dir+'loca_ens_prcptot_historical_mean_blank.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()


bias = loca_ens_prcptot/livneh_prcptot_30yr_avg
bias_levels = np.arange(0.78,1.205,0.005)

base_cmap = plt.get_cmap("bwr_r")

# Create the custom colormap by combining parts
colors = [
    (0.0, base_cmap(0.0)),       # Start with the first color of reversed 'bwr' (blue)
    (0.4, base_cmap(0.4)),       # First 40% is the first 40% of 'bwr_r'
    (0.4, 'white'),              # Start of the white section
    (0.6, 'white'),              # End of the white section (20%)
    (0.6, base_cmap(0.6)),       # Start of the last 40% of 'bwr_r'
    (1.0, base_cmap(1.0)),       # Last color of 'bwr_r' (red)
]

# Create the colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_bwr_r", colors)

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,1,1,projection=ccrs.Mercator())
ax.set_extent([-125,-65,24,50],crs=ccrs.PlateCarree())
c = plt.contourf(lon,lat,bias,bias_levels,transform=ccrs.PlateCarree(),cmap=custom_cmap,extend='both')
plt.savefig(image_dir+'livneh_loca_prcptot_historical_bias_blank_v2.png',dpi=300,bbox_inches='tight', pad_inches=0,transparent=True)
plt.show()
plt.close()


loca_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/LOCA2/1985-2014/ClimdexPhard_LOCA2-*_historical_r1i1p1f1_gn_CONUS16thD2021_1985-2014_ann.nc")
loca_files.sort()

for i in range(16):
    loca_data = nc.Dataset(loca_files[i])
    globals()['loca_'+str(models[i])+'_prcptot'] = np.array(loca_data.variables['prcptot'])
    mask = np.isclose(globals()['loca_'+str(models[i])+'_prcptot'],1e+20)
    globals()['loca_'+str(models[i])+'_prcptot'][mask] = np.nan
    mask = np.isnan(nclimgrid_prcptot)
    globals()['loca_'+str(models[i])+'_prcptot'][mask]= np.nan
    
locat = time.time()
for i in range(7):
    for j in range(16):
        key = f'loca_{models[j]}{regions_names_alt[i]}'
        data = globals()[f'loca_{models[j]}_prcptot']
        mask = region_masks[regions_names_alt[i]]
        globals()[key] = apply_region_mask(data, mask)
print("Masking took", time.time() - locat, "seconds")


for i in range(7):
    for j in range(16):
        globals()['loca_'+str(models[j])+str(regions_names_alt[i])+'_timeseries'] = np.nanmean(globals()['loca_'+str(models[j])+str(regions_names_alt[i])],axis=(1,2))
        
ensemble = np.empty([7,16,30])
for i in range(7):
    for j in range(16):
        ensemble[i,j,:] = globals()['loca_'+str(models[j])+str(regions_names_alt[i])+'_timeseries']

ensemble_avg = np.empty([7,30])
for i in range(7):
    ensemble_avg[i,:] = np.nanmean(ensemble[i,:,:],axis=0)
    
for h in range(7):
    fig = plt.figure(figsize=(14,6))
    for i in range(16):
        if i in {0}:
            plt.plot(np.arange(1985,2015,1),globals()['loca_'+str(models[i])+str(regions_names_alt[h])+'_timeseries'],color='grey',ls='--',alpha=0.75,label='LOCA2',lw=1)
        else:
            plt.plot(np.arange(1985,2015,1),globals()['loca_'+str(models[i])+str(regions_names_alt[h])+'_timeseries'],color='grey',ls='--',alpha=0.75,lw=1)
    plt.plot(np.arange(1985,2015,1),ensemble_avg[h,:],lw=3,label='Ens. Avg.',ls='--',color='k')
    plt.plot(np.arange(1985,2015,1),globals()['livneh_'+str(regions_names_alt[h])+'_timeseries'],lw=3,color='tab:orange',label='Livneh')
    plt.plot(np.arange(1985,2015,1),globals()['nclimgrid_'+str(regions_names_alt[h])+'_timeseries'],lw=3,label='nClimGrid')
    plt.legend()
    plt.xlim(1984.5,2014.5)
    plt.ylabel("Avg. Annual Total Precip. (mm)")
    plt.xticks(size=12)
    plt.savefig(image_dir+'prcptot_time_series_obs_loca_'+str(regions_names_alt[h])+'.png',dpi=300,bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close()
