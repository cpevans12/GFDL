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
import geopandas as gpd

image_dir = '/Volumes/Elements/GFDL_project/images/'
# =============================================================================
# 
# =============================================================================
shape_file = "/Users/cpe28/Downloads/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp"
us_states = gpd.read_file(shape_file)

# Define the Northeast states
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


nclimgrid_data = nc.Dataset('/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/nClimGrid/1985-2014/ClimdexPhard_day_nClimGrid_historical_r0i0p0_CONUS16thD2021_1985-2014_ann.nc')
# nclimgrid_prcptot = np.array(nclimgrid_data.variables['prcptot'])
# mask = np.isclose(nclimgrid_prcptot,1e+20)
# nclimgrid_prcptot[mask] = np.nan
lat = np.array(nclimgrid_data.variables['lat'])
lon = np.array(nclimgrid_data.variables['lon'])


lon_adjust = lon-360.
lon_grid, lat_grid = np.meshgrid(lon_adjust, lat)
 # Example data


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


def mask_array(array, states_gdf):
    masked_array = np.full(array.shape, np.nan)
    points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon_grid.flatten(), lat_grid.flatten()))
    for _, state in states_gdf.iterrows():
        state_geom = state['geometry']
        # Find points within the current state geometry
        points_within_state = points[points.geometry.within(state_geom)]        
        # Update the masked_array
        for idx in points_within_state.index:
            masked_array[np.unravel_index(idx, array.shape)] = array.flatten()[idx]
    return masked_array


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
# 
# =============================================================================
region_names = ['Northeast','Southeast','Midwest','North Great Plains','South Great Plains','Southwest','Northwest']
regions_names_alt=['northeast','southeast','midwest','north_plains','south_plains','southwest','northwest']

for i in range(7):
    globals()['livneh_'+str(regions_names_alt[i])] = mask_array(livneh_slopes, globals()[str(regions_names_alt[i])+'_us_states'])
    globals()['nclimgrid_'+str(regions_names_alt[i])] = mask_array(nclimgrid_slopes, globals()[str(regions_names_alt[i])+'_us_states'])
    globals()['prism_'+str(regions_names_alt[i])] = mask_array(prism_slopes, globals()[str(regions_names_alt[i])+'_us_states'])
    
# =============================================================================
#     
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



for i in range(7):
    for j in range(16):
        globals()['loca2_'+str(models[j])+str(regions_names_alt[i])] = mask_array(globals()['loca2_slope_'+str(models[j])+'_prcptot'], globals()[str(regions_names_alt[i])+'_us_states'])

for h in range(7):
    for i in range(16):
        globals()['loca2_'+str(models[i])+str(regions_names_alt[h])+"_cleaned"]=globals()['loca2_'+str(models[i])+str(regions_names_alt[h])][~np.isnan(globals()['loca2_'+str(models[i])+str(regions_names_alt[h])])]


for h in range(7):
    globals()["loca2_"+str(regions_names_alt[h])+"_ensemble"] = np.empty([16,444,922])
    for i in range(16):
        globals()["loca2_"+str(regions_names_alt[h])+"_ensemble"][i,:,:] = globals()['loca2_'+str(models[i])+str(regions_names_alt[h])]
    #

for i in range(7):
    globals()["loca2_"+str(regions_names_alt[i])+"_ensemble_mean"] = np.nanmean(globals()["loca2_"+str(regions_names_alt[i])+"_ensemble"],axis=0)

for i in range(7):
    globals()["loca2_"+str(regions_names_alt[i])+"_ensemble_cleaned"] = globals()["loca2_"+str(regions_names_alt[i])+"_ensemble_mean"][~np.isnan(globals()["loca2_"+str(regions_names_alt[i])+"_ensemble_mean"])]


######
for i in range(7):
    globals()['nclimgrid_'+str(regions_names_alt[i])+'_cleaned'] = globals()['nclimgrid_'+str(regions_names_alt[i])][~np.isnan(globals()['nclimgrid_'+str(regions_names_alt[i])])]
    
    globals()['livneh_'+str(regions_names_alt[i])+'_cleaned'] = globals()['livneh_'+str(regions_names_alt[i])][~np.isnan(globals()['livneh_'+str(regions_names_alt[i])])]
    
    globals()['prism_'+str(regions_names_alt[i])+'_cleaned'] = globals()['prism_'+str(regions_names_alt[i])][~np.isnan(globals()['prism_'+str(regions_names_alt[i])])]

########
all_trends = [globals()['nclimgrid_'+str(regions_names_alt[0])+'_cleaned'],globals()['livneh_'+str(regions_names_alt[0])+'_cleaned'],globals()['prism_'+str(regions_names_alt[0])+'_cleaned']]
titles = ["nClimGrid","Livneh","PRISM"]
for i in range(16):
    all_trends.append(globals()['loca2_'+str(models[i])+str(regions_names_alt[0])+"_cleaned"])
    titles.append(str(models[i]))

for i in range(7):
    all_trends = [globals()['nclimgrid_'+str(regions_names_alt[i])+'_cleaned'],globals()['livneh_'+str(regions_names_alt[i])+'_cleaned'],globals()['prism_'+str(regions_names_alt[i])+'_cleaned'],globals()["loca2_"+str(regions_names_alt[i])+"_ensemble_cleaned"]]
    titles = ["nClimGrid","Livneh","PRISM","LOCA2 Ensemble"]
    for j in range(16):
        all_trends.append(globals()['loca2_'+str(models[j])+str(regions_names_alt[i])+"_cleaned"])
        titles.append(str(models[j]))
    
    fig = plt.figure(figsize=(8,6))
    plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
    # Customize plot
    plt.title(region_names[i]+" PRCPTOT Trends (Theil Slope)\nLOCA2 Models v. Obs.")
    plt.ylabel(r"$m_{\mathrm{TS}}$ (mm/year)")
    plt.xticks(np.arange(1, len(titles) + 1), titles, rotation=45, ha='right')  # Adjust X-axis label to represent dataset
    plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
    # plt.ylim(-1,3.5)
    plt.savefig(image_dir+str(regions_names_alt[i])+"_boxplot_prcptot_slopes_obs_loca2_all.png",dpi=600,bbox_inches='tight')
    plt.show()
    plt.close()


# =============================================================================
#     
# =============================================================================

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','EC-Earth3','FGOALS-g3','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM']

letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

star_files = glob.glob("/Volumes/Elements/GFDL_project/DOWNSCALING/SUBPROJECTS/LOCA2vsSTAR/diagnostics/Climdex/ClimdexPhard/Individual/STAR/1985-2014/ClimdexPhard_STAR-*CONUS16thD2021_1985-2014_ann.nc")
star_files.sort()

for i in range(16):
    star_data = nc.Dataset(star_files[i])
    globals()['star_'+str(models[i])+'_prcptot'] = np.array(star_data.variables['prcptot'])
    mask = np.isclose(globals()['star_'+str(models[i])+'_prcptot'],1e+20)
    globals()['star_'+str(models[i])+'_prcptot'][mask] = np.nan
    mask = np.isnan(nclimgrid_prcptot)
    globals()['star_'+str(models[i])+'_prcptot'][mask]= np.nan
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


for i in range(7):
    for j in range(16):
        globals()['star_'+str(models[j])+str(regions_names_alt[i])] = mask_array(globals()['star_slope_'+str(models[j])+'_prcptot'], globals()[str(regions_names_alt[i])+'_us_states'])

for h in range(7):
    for i in range(16):
        globals()['star_'+str(models[i])+str(regions_names_alt[h])+"_cleaned"]=globals()['star_'+str(models[i])+str(regions_names_alt[h])][~np.isnan(globals()['star_'+str(models[i])+str(regions_names_alt[h])])]


for h in range(7):
    globals()["star_"+str(regions_names_alt[h])+"_ensemble"] = np.empty([16,444,922])
    for i in range(16):
        globals()["star_"+str(regions_names_alt[h])+"_ensemble"][i,:,:] = globals()['star_'+str(models[i])+str(regions_names_alt[h])]
    #

for i in range(7):
    globals()["star_"+str(regions_names_alt[i])+"_ensemble_mean"] = np.nanmean(globals()["star_"+str(regions_names_alt[i])+"_ensemble"],axis=0)

for i in range(7):
    globals()["star_"+str(regions_names_alt[i])+"_ensemble_cleaned"] = globals()["star_"+str(regions_names_alt[i])+"_ensemble_mean"][~np.isnan(globals()["star_"+str(regions_names_alt[i])+"_ensemble_mean"])]

######
for i in range(7):
    globals()['nclimgrid_'+str(regions_names_alt[i])+'_cleaned'] = globals()['nclimgrid_'+str(regions_names_alt[i])][~np.isnan(globals()['nclimgrid_'+str(regions_names_alt[i])])]
    
    globals()['livneh_'+str(regions_names_alt[i])+'_cleaned'] = globals()['livneh_'+str(regions_names_alt[i])][~np.isnan(globals()['livneh_'+str(regions_names_alt[i])])]
    
    globals()['prism_'+str(regions_names_alt[i])+'_cleaned'] = globals()['prism_'+str(regions_names_alt[i])][~np.isnan(globals()['prism_'+str(regions_names_alt[i])])]

########
all_trends = [globals()['nclimgrid_'+str(regions_names_alt[0])+'_cleaned'],globals()['livneh_'+str(regions_names_alt[0])+'_cleaned'],globals()['prism_'+str(regions_names_alt[0])+'_cleaned']]
titles = ["nClimGrid","Livneh","PRISM"]
for i in range(16):
    all_trends.append(globals()['star_'+str(models[i])+str(regions_names_alt[0])+"_cleaned"])
    titles.append(str(models[i]))

for i in range(7):
    all_trends = [globals()['nclimgrid_'+str(regions_names_alt[i])+'_cleaned'],globals()['livneh_'+str(regions_names_alt[i])+'_cleaned'],globals()['prism_'+str(regions_names_alt[i])+'_cleaned'],globals()["star_"+str(regions_names_alt[i])+"_ensemble_cleaned"]]
    titles = ["nClimGrid","Livneh","PRISM","STAR Ensemble"]
    for j in range(16):
        all_trends.append(globals()['star_'+str(models[j])+str(regions_names_alt[i])+"_cleaned"])
        titles.append(str(models[j]))
    
    fig = plt.figure(figsize=(8,6))
    plt.boxplot(all_trends, vert=True, patch_artist=True, showfliers=True)
    # Customize plot
    plt.title(region_names[i]+" PRCPTOT Trends (Theil Slope)\nSTAR Models v. Obs.")
    plt.ylabel(r"$m_{\mathrm{TS}}$ (mm/year)")
    plt.xticks(np.arange(1, len(titles) + 1), titles, rotation=45, ha='right')  # Adjust X-axis label to represent dataset
    plt.axhline(0, color="gray", linestyle="--", linewidth=1)  # Reference line at zero trend
    # plt.ylim(-1,3.5)
    plt.savefig(image_dir+str(regions_names_alt[i])+"_boxplot_prcptot_slopes_obs_star_all.png",dpi=600,bbox_inches='tight')
    plt.show()
    plt.close()
