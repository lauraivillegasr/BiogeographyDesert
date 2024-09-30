#This is the code implemented for plotting the estimations from the different models based on the environmental variables defined as the most important to predict either genera richness or reproductive mode
#Input files used are specified in the manuscript (15arc_dem.zarr, Chelsa 1981-2010 pr and tas, NASADEM, ASTER_GDEM_DEM, GLim_wgs_500m and average_soil_and_sediment )

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import LightSource
from cmcrameri import cm
import geopandas as gp
from shapely.geometry import Point, Polygon
import matplotlib.colors as cls
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd


def sel_dem(src_dem, sel_lat=[-28, -18], sel_lon=[-72, -66]):
    """
    Crop selection of a digital elevation model (DEMs)

    :param src_dem: string,
                    source of DEM either 'nasadem', 'aster' or 'srtm'
    :param sel_lat: list,
                minimum and maximum latitude to crop part of DEM
    :param sel_lon: list,
                minimum and maximum longitude to crop part of DEM
    :return: xarray.Dataset
            Cropped selection of DEM
    """

    src_path = "/Users/acevedo/Documents/Projects/EarthEvoDryLimit/AS-LEM_scenarios/Data/"

    if src_dem == 'nasadem':
        dem = (xr.open_dataset(src_path + 'NASADEM/NASADEM_NC.001_NASADEM_HGT_doy2000042_aid0001.tif',
                               engine='rasterio')
               .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
               .band_data.sel(band=1).drop_vars('band')
               .rename('elevation')
               .rio.write_crs("epsg:4326", inplace=True)
               )
        dem = dem.where(dem >= 0, drop=True)

    elif src_dem == 'srtm':
        dem = (xr.open_zarr(src_path + '15arc_dem.zarr', consolidated=False)
               .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
               .rename({'elevation_dem15': 'elevation'})
               .rio.write_crs("epsg:4326", inplace=True)
               ).elevation

    return dem


def chelsa_data(data_name, sel_lat=[-28, -18], sel_lon=[-72, -66]):
    import glob
    src_path = "/Users/acevedo/Documents/Projects/EarthEvoDryLimit/AS-LEM_scenarios/Data/"
    lst_pr = sorted(glob.glob(src_path + 'chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/' + data_name + '/*.tif'))

    def prepros(ds):
        """
        preprocessing of CHELSA raster files
        :param ds: xr.Dataset
        :return:
        """
        if 'pr' in ds.encoding['source']:
            new_name = 'pr'
        elif 'tas' in ds.encoding['source']:
            new_name = 'tas'
        else:
            ValueError('CHELSA variable not found, must be either "pr" or "tas"')

        month = int(ds.encoding['source'][-22:-20])

        ds_out = (ds.squeeze('band')
                  .drop_vars(['band'])
                  .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
                  .rename({'band_data': new_name})
                  .assign_coords(time=month)
                  )
        return ds_out

    # CHELSA v2.1 dataset
    dts = xr.open_mfdataset(lst_pr, engine='rasterio', combine='nested', concat_dim=['time'],
                            data_vars=[data_name], preprocess=prepros)

    return dts


def lat_range_precip():
    """
    Extract predictor variables of the best model, i.e. Latitude and range of precipitation,
    CHELSA monthly climatological data.

    :param sel_lat: list,
                minimum and maximum latitude to crop part of climatic fields
    :param sel_lon: list,
                minimum and maximum longitude to crop part of climatic fields
    :return: numpy.arrays,
            latitude and range of precipitation
    """

    # CHELSA v2.1 Precipitation (mm) dataset
    dts_pr = chelsa_data('pr')
    # Range Precipitation
    minp = dts_pr['pr'].min(dim='time', skipna=True)
    maxp = dts_pr['pr'].max(dim='time', skipna=True)
    range_pre = maxp - minp

    longitude, latitude = np.meshgrid(dts_pr.x.values, dts_pr.y.values)

    return longitude, latitude, range_pre.to_numpy()


def rocktype_var(sel_lat=[-28, -18], sel_lon=[-72, -66], resampled_data=False):
    """
        Obtain rock type for a set of latitude and longitude values.

        :param sel_lat: list,
                minimum and maximum latitude to crop part of raster file
        :param sel_lon: list,
                minimum and maximum longitude to crop part of raster file
        :return: pandas.Series,
                 with rock type
                 at the provided latitude and longitude
        """
    src_path = "/Users/acevedo/Documents/Projects/EarthEvoDryLimit/AS-LEM_scenarios/Data/"
    if resampled_data:
        rocktype_file = src_path + 'RockTypes/GLim_wgs_500m_resampled_1km.tif'
    else:
        rocktype_file = src_path + 'RockTypes/GLim_wgs_500m.tif'
    rock_type_raster = (xr.open_dataarray(rocktype_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))

    return rock_type_raster


def soilthickness(sel_lat=[-28, -18], sel_lon=[-72, -66]):
    """
    Obtain soil thickness for a set of latitude and longitude values.

        :param lats: xr.DataArray,
                     latitude values to extract data
        :param lons: xr.DataArray,
                     longitude values to extract data
        :param sel_lat: list,
                minimum and maximum latitude to crop part of raster file
        :param sel_lon: list,
                minimum and maximum longitude to crop part of raster file
        :return: np.array,
                 with soil thickness
                 at the provided latitude and longitude
    """
    src_path = "/Users/acevedo/Documents/Projects/EarthEvoDryLimit/AS-LEM_scenarios/Data/"
    soilthick_file = src_path + 'SoilThickness/data/average_soil_and_sedimentary-deposit_thickness.tif'
    soilthick_raster = (xr.open_dataarray(soilthick_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))

    return soilthick_raster


def x_pred():
    precip = chelsa_data('pr')
    map = precip['pr'].sum(dim='time', skipna=True)
    temp = chelsa_data('tas')
    mint = temp['tas'].min(dim='time', skipna=True)
    maxt = temp['tas'].max(dim='time', skipna=True)
    range_tmp = maxt - mint
    soilthick = soilthickness()
    lons, lats = np.meshgrid(precip.x.values, precip.y.values)

    dem = sel_dem('srtm')
    dem = dem.where(dem > 0)
    dem_resamp = dem.coarsen(x=2).mean().coarsen(y=2).mean()
    rocktype = rocktype_var(resampled_data=True)

    X_pred = pd.DataFrame({'Latitude': lats.flatten(),
                           'dem_ele': dem_resamp.values.flatten(),
                           'MAP_chelsa': map.values.flatten(),
                           'range_tmp_chelsa': range_tmp.values.flatten(),
                           'rock_type_int': rocktype.values.flatten(),
                           'soil_thick_log': np.log(soilthick.values.flatten() + 1)

                           })

    X_pred.to_parquet('x_pred.parquet', engine='pyarrow')


def make_figure01(foutname='map_env_vars.png'):
    dem = sel_dem('nasadem')
    dem = dem.where(dem > 0)
    dem_hillshade = sel_dem('srtm')
    rocktype = rocktype_var()
    soilthick = soilthickness()

    precip = chelsa_data('pr')
    map = precip['pr'].sum(dim='time', skipna=True)

    temp = chelsa_data('tas')
    mint = temp['tas'].min(dim='time', skipna=True)
    maxt = temp['tas'].max(dim='time', skipna=True)
    range_tmp = maxt - mint

    src_path = "/Users/acevedo/Documents/Projects/EarthEvoDryLimit/AS-LEM_scenarios/Data/"
    world = gp.read_file(src_path + 'ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')

    bounds_rck = np.arange(14)

    locs = (pd.read_csv('./Location_samplespots.csv', delimiter=';')
            .drop(columns=['Unnamed: 0'])
            )
    ls = LightSource(azdeg=315, altdeg=45)
    dxy = 1e5 / 100
    ve = 10
    z = dem_hillshade.to_numpy()

    fig, axs = plt.subplots(2, 3, figsize=(10, 10))

    world.plot(ax=axs[0, 0], color='black')
    axs[0, 0].set_ylim(-57, 13)
    axs[0, 0].set_xlim(-82, -34)
    axs[0, 0].grid(True, linestyle='--')

    pcm1 = axs[0, 1].pcolormesh(map.x, map.y, map / 1000, cmap=cm.lapaz_r, vmin=0, vmax=2.00)
    pcm2 = axs[0, 2].pcolormesh(range_tmp.x, range_tmp.y, range_tmp, cmap=cm.lajolla_r, vmin=0, vmax=25)
    pcm3 = axs[1, 0].pcolormesh(dem.x, dem.y, dem / 1000, cmap=cm.bamako, vmin=0, vmax=7.00)
    pcm4 = axs[1, 1].pcolormesh(rocktype.x, rocktype.y, rocktype, cmap=cm.glasgowS,
                                norm=cls.BoundaryNorm(bounds_rck, 14))
    pcm5 = axs[1, 2].pcolormesh(soilthick.x, soilthick.y, soilthick, cmap=cm.imola)

    pcms = [pcm1, pcm2, pcm3, pcm4, pcm5]

    for i, ax in enumerate(axs.flatten()):
        if i > 0:
            ax.pcolormesh(dem_hillshade.x, dem_hillshade.y, ls.hillshade(z, vert_exag=ve, dx=dxy, dy=dxy),
                          cmap='gray', alpha=0.5)
            ax.scatter(x=locs.Longitud, y=locs.Latitude, s=50, edgecolor='white', facecolor='black', alpha=0.75)
            ax.set_xlim(-72, -66)
            ax.set_ylim(-28, -18)

            cbbox = inset_axes(ax, width='20%', height='45%', loc='upper left')
            [cbbox.spines[k].set_visible(False) for k in cbbox.spines]
            cbbox.tick_params(axis='both', left=False, top=False,
                              right=False, bottom=False, labelleft=False,
                              labeltop=False, labelright=False, labelbottom=False)
            cbbox.set_facecolor([1, 1, 1, 0.75])
            cax = inset_axes(cbbox, '16%', '90%', loc=6)
            fig.colorbar(pcms[i - 1], cax=cax, orientation='vertical')

        if 0 < i < 2:
            ax.set_xticklabels([])
            ax.yaxis.tick_right()
            ax.set_yticklabels([])
        if 2 < i < 5:
            ax.set_yticklabels([])
            ax.yaxis.tick_right()
        if i == 2:
            ax.set_xticklabels([])
            ax.yaxis.tick_right()
        if i == 5:
            ax.yaxis.tick_right()

    cross_sec_edges = gp.GeoDataFrame([['box', Point(-72.0, -18.0)],
                                       ['box', Point(-72.0, -28.0)],
                                       ['box', Point(-66.0, -28.0)],
                                       ['box', Point(-66.0, -18.0)]],
                                      columns=['shape_id', 'geometry'],
                                      geometry='geometry',
                                      crs='EPSG:4326'
                                      )
    box1 = cross_sec_edges.groupby('shape_id')['geometry'].apply(lambda x: Polygon(x.tolist())).reset_index()
    box1.plot(facecolor='none', alpha=0.75, edgecolor='red', ax=axs[0, 0], linewidth=3, aspect=None)

    axs[0, 1].set_title('MAP', weight='bold', fontsize=16)
    axs[0, 2].set_title('Temperature\nRange', weight='bold', fontsize=16)
    axs[1, 0].set_title('Elevation', weight='bold', fontsize=16)
    axs[1, 1].set_title('Rock Type', weight='bold', fontsize=16)
    axs[1, 2].set_title('Soil\nThickness', weight='bold', fontsize=16)

    fig.text(0.05, 0.38, "Latitude [degrees]", weight='bold', fontsize=16, rotation=90)
    axs[1, 1].set_xlabel("Longitude [degrees]", weight='bold', fontsize=16)

    fig.subplots_adjust(wspace=0.15, hspace=0.2)
    fig.savefig(foutname, dpi=600)


def make_figure02(foutname='map_predictions.png'):
    '''
    Map of nematode biodiversity predicted by a model that considers
    range in precipitation and latitude.
    :param foutname: string,
                     name of output file
    :return:
    '''

    dem = sel_dem('srtm')
    dem_coarse = dem.coarsen(x=2).max().coarsen(y=2).max()
    lon = dem_coarse.x
    lat = dem_coarse.y

    rich = (np.genfromtxt("./revised_models/pred_gls_richness.csv", skip_header=1,
                          delimiter=',', usecols=1)).reshape(1200, 720)
    rich = np.ma.array(rich, mask=dem_coarse.to_masked_array().mask)

    rich_rf = (np.genfromtxt("./revised_models/pred_rf_richness.csv", skip_header=1,
                             delimiter=',', usecols=1)).reshape(1200, 720)

    repmod = (np.genfromtxt("./revised_models/pred_glm_reproductivemode.csv", skip_header=1,
                            delimiter=',', usecols=1)).reshape(1200, 720)

    ls = LightSource(azdeg=315, altdeg=45)
    dxy = 1e5 / 100
    ve = 10
    z = dem.to_numpy()

    fig, axs = plt.subplots(2, 2, sharey='all', sharex='all', figsize=(8, 10))
    pcm1 = axs[0, 0].pcolormesh(lon, lat, rich, cmap=cm.nuuk,
                     vmin=0, vmax=np.ceil(np.nanmax(rich)))
    pcm2 = axs[0, 1].pcolormesh(lon, lat, rich_rf, cmap=cm.nuuk,
                                vmin=0, vmax=np.ceil(np.nanmax(rich_rf)))
    pcm3 = axs[1, 0].pcolormesh(lon, lat, repmod, cmap=cm.lipari,
                                vmin=0, vmax=1)

    pcms = [pcm1, pcm2, pcm3]

    for pcm, ax in zip(pcms, axs.flatten()):
        ax.pcolormesh(dem.x, dem.y, ls.hillshade(z, vert_exag=ve, dx=dxy, dy=dxy), cmap='gray', alpha=0.5)
        ax.set_xlim(-72, -66)
        ax.set_ylim(-28, -18)
        ax.grid(True, linestyle='--', color='black')

        cbbox = inset_axes(ax, width='20%', height='45%', loc='upper left')
        [cbbox.spines[k].set_visible(False) for k in cbbox.spines]
        cbbox.tick_params(axis='both', left=False, top=False,
                          right=False, bottom=False, labelleft=False,
                          labeltop=False, labelright=False, labelbottom=False)
        cbbox.set_facecolor([1, 1, 1, 0.75])
        cax = inset_axes(cbbox, '16%', '90%', loc=6)
        fig.colorbar(pcm, cax=cax, orientation='vertical')

    axs[0, 0].set_title('R=Latitude+MAP', weight='bold', fontsize=14)
    axs[0, 1].set_title('R=Elevation+Latitude\n+MAP+Soil+Rock', weight='bold', fontsize=14)
    axs[1, 0].set_title('P=Elevation', weight='bold', fontsize=14)

    fig.text(0.925, 0.67, "Richness [R]", weight='bold', fontsize=16, rotation=270)
    fig.text(0.925, 0.08, "Sexual reproduction probability [P]", weight='bold', fontsize=16, rotation=270)

    fig.text(0.05, 0.38, "Latitude [degrees]", weight='bold', fontsize=16, rotation=90)
    fig.text(0.38, 0.05, "Longitude [degrees]", weight='bold', fontsize=16)
    fig.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.savefig(foutname, dpi=600)


if __name__ == "__main__":
    make_figure01()
    make_figure02()
