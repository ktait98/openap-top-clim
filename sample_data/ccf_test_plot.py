import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import AxesGrid

#plt.rc('font',**{'family':'serif','serif':['cmr10']})
#plt.rc('text', usetex=True)
#font = {'family' : 'normal',
#'size' : 13}

# Extract ccf dataset from nc file for further processing and plotting!
ccf = xr.open_dataset("env_processed_new.nc")

lats = ccf['latitude'].values
lons = ccf['longitude'].values
lons1,lats1 = np.meshgrid(lons,lats)
cc_lon = np.flipud(lons1)[::1, ::1]
cc_lat = np.flipud(lats1)[::1, ::1]
print(lats)
print(lats1)
print(cc_lat)

time = np.datetime64('2022-08-16')
pressure_level = 250
time_idx = np.where (ccf.time.values == time)[0][0]
pl_idx = np.where (ccf.level.values == pressure_level) [0][0]
aCCF_merged = np.flipud(ccf['aCCF_merged'].values[time_idx, pl_idx, :, :])[::1,
                                                                          ::1]
#aCCF_dcontrail = np.flipud(ccf['aCCF_dcontrail'].values[time_idx, pl_idx, :, :])[::1,
 #                                                                         ::1]

def main():
    projection = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection=projection))

    fig = plt.figure(figsize=(5,5))
    axgr = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(1,1),
                    axes_pad=1.0,
                    share_all=True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.2,
                    cbar_size='3%',
                    label_mode='') # Note empty label_mode

    for i, ax in enumerate(axgr):
        xticks = [-5, 10, 25]
        #xticks = [-20, -5, 10, 25, 40, 55]
        yticks = [40, 50, 60]
        ax.coastlines(zorder=1)
        ax.add_feature(cfeature.BORDERS, edgecolor='0', linewidth=1, alpha=0.7,
                       zorder=0.8)
        ax.set_xticks(xticks, crs=projection)
        ax.set_yticks(yticks, crs=projection)

        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.set_title(time)
        p = ax.contourf(cc_lon, cc_lat, aCCF_merged,
                        transform=projection,
                        cmap='YlOrRd', alpha=0.9, zorder=0.5)
        axgr.cbar_axes[i].colorbar(p)
        cax = axgr.cbar_axes[i]
        axis = cax.axis[cax.orientation]
        axis.label.set_text('aCCF-merged [K/kg(fuel)]')

    plt.show()

main()
"""
# Print to test output of aCCF_merged
print(ccf.level.values)
print(ccf.latitude.values)
print(ccf.longitude.values)
print(ccf.time.values)

# Select (by integer) first point in dataset (lon[i], lat[i] and level[i]
# when i = 0)
first_point = ccf.aCCF_merged.isel(latitude=0, longitude=0, level=0)


# Plot timeseries of 'first_point'
#first_point.plot()

# Do timeseries with multiple points
#ccf.aCCF_merged.sel(latitude=20, longitude=[10, 20], level=200, method="nearest").plot.line(x="time")

# Histograms
#ccf.aCCF_merged.plot()

# Maps
#Using geophysical units. 'robust' disregards outliers for colour map creation

# Reproject on globe
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=1, color='grey', alpha=0.5, linestyle='--')

#ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
#ax.set_xticks([-180, -120, -60, 0, 60, 120,180], crs=ccrs.PlateCarree())

#ax.pcolormesh(ccf.aCCF_merged.isel(time=0, level=1))
plt.show()
"""
