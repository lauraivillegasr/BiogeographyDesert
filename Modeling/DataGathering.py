#This is the code implemented for gathering environmental and geomorphological variables for each of the smapling spots sampled and analyzed
#Input files used are specified in the manuscript (15arc_dem.zarr, Chelsa 1981-2010 pr and tas, NASADEM, ASTER_GDEM_DEM, GLim_wgs_500m and average_soil_and_sediment )

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import richdem as rd
import rioxarray


def sel_dem(src_dem, sel_lat=[-32, -18], sel_lon=[-72, -69]):
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

    src_path = "/data/"

    if src_dem == 'nasadem':
        dem = (xr.open_dataset(src_path + 'NASADEM/NASADEM_NC.001_NASADEM_HGT_doy2000042_aid0001.tif',
                               engine='rasterio')
               .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
               .band_data.sel(band=1).drop_vars('band')
               .rename('elevation')
               .rio.write_crs("epsg:4326", inplace=True)
               )
        dem = dem.where(dem >= 0, drop=True)

    elif src_dem == 'aster':
        dem = ((xr.open_dataset(src_path + 'ASTER_GDEM_DEM/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif',
                               engine='rasterio')
               .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
               .band_data.sel(band=1).drop_vars('band')
               .rename('elevation')
               .rio.write_crs("epsg:4326", inplace=True)
               ))
        dem = dem.where(dem >= 0, drop=True)

    elif src_dem == 'srtm':
        dem = (xr.open_zarr(src_path + '15arc_dem.zarr', consolidated=False)
               .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1]))
               .rename({'elevation_dem15': 'elevation'})
               .rio.write_crs("epsg:4326", inplace=True)
               ).elevation

    return dem


def dem_attributes(lats, lons, src_dem='srtm', res_deg=0.01, sel_lat=[-32, -18], sel_lon=[-72, -69]):
    """
    Visualization of various dem attributes
    :param src_dem: string,
                    source of DEM either 'nasadem', 'aster' or 'srtm'
    :param res_deg: float
                    rescaling factor in degrees to coarsen dem
    :return: splots
    """

    dem = sel_dem(src_dem)
    rda = rd.rdarray(dem.to_numpy(), no_data=-9999)
    slp = rd.TerrainAttribute(rda, 'slope_riserun')
    slp = xr.DataArray(slp, coords=dem.coords)
    curv = rd.TerrainAttribute(rda, 'curvature')
    curv = xr.DataArray(curv, coords=dem.coords)

    src_path = "/data/"

    soilthick_file = src_path + 'average_soil_and_sedimentary-deposit_thickness.tif'
    soilthick_raster = (xr.open_dataarray(soilthick_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))

    rocktype_file = src_path + 'RockTypes/GLim_wgs_500m.tif'
    rock_type_raster = (xr.open_dataarray(rocktype_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))

    fig, axs = plt.subplots(2, 3, sharey='row', sharex='col',  figsize=(12, 9))
    dem.plot(ax=axs[0, 0])
    slp.plot(ax=axs[0, 1])
    curv.plot(ax=axs[0, 2])
    axs[0, 0].scatter(x=lons, y=lats, c='Black')
    axs[0, 1].scatter(x=lons, y=lats, c='Black')
    axs[0, 2].scatter(x=lons, y=lats, c='Black')
    axs[0, 0].set_title('DEM')
    axs[0, 1].set_title('Slope')
    axs[0, 2].set_title('Curvature')

    x_coarse = int(round(res_deg / (abs(dem.x.min() - dem.x.max()).values.round() / dem.x.size)))
    y_coarse = int(round(res_deg / (abs(dem.y.min() - dem.y.max()).values.round() / dem.y.size)))
    #cdem = dem.coarsen(dim={'x': x_coarse, 'y': y_coarse}, boundary='pad').mean()
    #crda = rd.rdarray(cdem.to_numpy(), no_data=-9999)
    #cslp = rd.TerrainAttribute(crda, 'slope_riserun')
    #cslp = xr.DataArray(cslp, coords=cdem.coords)
    #ccurv = rd.TerrainAttribute(crda, 'curvature')
    #ccurv = xr.DataArray(ccurv, coords=cdem.coords)

    dem.coarsen(dim={'x': x_coarse, 'y': y_coarse}, boundary='pad').std().plot(ax=axs[1, 0])
    soilthick_raster.plot(ax=axs[1, 1])
    rock_type_raster.plot(ax=axs[1, 2])

    axs[1, 0].scatter(x=lons, y=lats, c='Black')
    axs[1, 1].scatter(x=lons, y=lats, c='Black')
    axs[1, 2].scatter(x=lons, y=lats, c='Black')

    #cslp.plot(ax=axs[1, 1])
    #ccurv.plot(ax=axs[1, 2])
    axs[1, 0].set_title('Topographic Complexity')
    axs[1, 1].set_title('Soil Thickness')
    axs[1, 2].set_title('Rock type')


def geomorph_vars(lats, lons, src_dem='nasadem', res_deg=0.01):
    """
    Obtain geomorphological variables for a set of latitude and
    longitude values. The computed geomorphological variables are
    elevation, curvature and topographic complexity.

    :param lats: xr.DataArray,
                 latitude values to extract data
    :param lons: xr.DataArray,
                 longitude values to extract data
    :param src_dem: string,
                source of DEM either 'nasadem', 'aster' or 'srtm'
    :param res_deg:
    :return: np.array,
             with elevation, curvature and topographic complexity values
             at the provided latitude and longitude
    """

    dem = sel_dem(src_dem)
    rda = rd.rdarray(dem.to_numpy(), no_data=-9999)
    cur = rd.TerrainAttribute(rda, 'curvature')
    cur = xr.DataArray(cur, coords=dem.coords)
    x_coarse = int(round(res_deg / (abs(dem.x.min() - dem.x.max()).values.round() / dem.x.size)))
    y_coarse = int(round(res_deg / (abs(dem.y.min() - dem.y.max()).values.round() / dem.y.size)))
    top_com = dem.coarsen(dim={'x': x_coarse, 'y': y_coarse}, boundary='pad').std()

    dem_ele = dem.sel(x=lons, y=lats, method='nearest').values
    dem_cur = cur.sel(x=lons, y=lats, method='nearest').values
    top_com = top_com.sel(x=lons, y=lats, method='nearest').values

    return dem_ele, dem_cur, top_com

def clima_vars_chelsa(lats, lons, sel_lat=[-32, -18], sel_lon=[-72, -69]):
    """
    Obtain climatical variables from CHELSA monthly climatologies and for a set of locations.
    The computed climatical variables are mean annual precipitation, range of precipitation,
    mean annual temperature, range of temperature.

    :param lats: xr.DataArray,
                 latitude values to extract data
    :param lons: xr.DataArray,
                 longitude values to extract data
    :param sel_lat: list,
                minimum and maximum latitude to crop part of climatic fields
    :param sel_lon: list,
                minimum and maximum longitude to crop part of climatic fields
    :return: numpy.arrays,
             with mean annual precipitation, range of precipitation,
             mean annual temperature, range of temperature
             at the provided latitude and longitude
    """

    import glob
    src_path = "/data/"
    lst_pr = sorted(glob.glob(src_path + 'pr/*.tif'))
    lst_tas = sorted(glob.glob(src_path + 'tas/*.tif'))

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

    # CHELSA v2.1 Precipitation (mm) dataset
    dts_pr = xr.open_mfdataset(lst_pr, combine='nested', concat_dim=['time'], data_vars=['pr'], preprocess=prepros)

    # Mean Annual Precipitation
    map = dts_pr['pr'].sum(dim='time', skipna=True)
    map_sel = map.sel(y=lats, x=lons, method='nearest')

    # Range Precipitation
    minp = dts_pr['pr'].min(dim='time', skipna=True)
    maxp = dts_pr['pr'].max(dim='time', skipna=True)
    range_pre = maxp - minp
    range_pre_sel = range_pre.sel(y=lats, x=lons, method='nearest')

    # CHELSA v2.1 Temperature (degree C) dataset
    dts_tas = xr.open_mfdataset(lst_tas, combine='nested', concat_dim=['time'], data_vars=['tas'], preprocess=prepros)

    # Mean Annual Temperature
    mat = dts_tas['tas'].mean(dim='time', skipna=True)
    mat_sel = mat.sel(y=lats, x=lons, method='nearest')

    # Range Precipitation
    mint = dts_tas['tas'].min(dim='time', skipna=True)
    maxt = dts_tas['tas'].max(dim='time', skipna=True)
    range_tmp = maxt - mint
    range_tmp_sel = range_tmp.sel(y=lats, x=lons, method='nearest')

    return map_sel, range_pre_sel, mat_sel, range_tmp_sel


def rocktype_var(lats, lons, sel_lat=[-32, -18], sel_lon=[-72, -69]):
    """
        Obtain rock type for a set of latitude and longitude values.

        :param lats: xr.DataArray,
                     latitude values to extract data
        :param lons: xr.DataArray,
                     longitude values to extract data
        :param sel_lat: list,
                minimum and maximum latitude to crop part of raster file
        :param sel_lon: list,
                minimum and maximum longitude to crop part of raster file
        :return: pandas.Series,
                 with rock type
                 at the provided latitude and longitude
        """
    src_path = "/data/"
    rocktype_file = src_path + 'GLim_wgs_500m.tif'
    rock_type_raster = (xr.open_dataarray(rocktype_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))
    rock_type_raster_sel = rock_type_raster.sel(y=lats, x=lons, method='nearest')

    def get_rock_type(val):
        rock_type_class = {
            'UncSed': [1],
            'SilRck': [2],
            'PyrRck': [3],
            'MxSRck': [4],
            'CarRck': [5],
            'EvaRck': [6],
            'AciRck': [7],
            'InVRck': [8],
            'BaVRck': [9],
            'AcPRck': [10],
            'InPRck': [11],
            'BaPRck': [12],
            'MetRck': [13],
            'WatBod': [14],
            'IceGla': [15],
            'nan': [16]
        }
        for k, v in rock_type_class.items():
            if val in v:
                return k
        return np.nan

    rock_type = rock_type_raster_sel.to_dataset(name='rock_type_int').astype(int).rock_type_int.to_pandas().apply(get_rock_type)

    rock_type = rock_type.astype('category')
    rock_type_int= rock_type_raster_sel.to_dataset(name='rock_type_int').astype(int).rock_type_int.to_pandas()

    return rock_type, rock_type_int


def soilthickness_vars(lats, lons, sel_lat=[-32, -18], sel_lon=[-72, -69]):
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
    src_path = "data/"
    soilthick_file = src_path + 'average_soil_and_sedimentary-deposit_thickness.tif'
    soilthick_raster = (xr.open_dataarray(soilthick_file, engine="rasterio")
                        .squeeze('band').drop_vars(['band'])
                        .sel(y=slice(sel_lat[1], sel_lat[0]), x=slice(sel_lon[0], sel_lon[1])))
    soilthick_raster_sel = soilthick_raster.sel(y=lats, x=lons, method='nearest')

    return soilthick_raster_sel


def main():
    """
    Main procedure to execute
    """
    import time

    start_time = time.time()

    print('computing geomorphological variables')
    dtf = pd.read_csv("/data/large_dataset.csv", delimiter='\t', decimal=',')
    lons = dtf.to_xarray().load().Longitud
    lats = dtf.to_xarray().load().Latitude
    

    dem_ele, dem_cur, top_com = geomorph_vars(lats, lons)
    dtf['dem_ele'] = dem_ele
    dtf['dem_cur'] = dem_cur
    dtf['top_com'] = top_com

    map_chelsa, range_pre_chelsa, mat_chelsa, range_tmp_chelsa = clima_vars_chelsa(lats, lons)
    dtf['MAP_chelsa'] = map_chelsa
    dtf['range_pre_chelsa'] = range_pre_chelsa
    dtf['MAT_chelsa'] = mat_chelsa
    dtf['range_tmp_chelsa'] = range_tmp_chelsa

    print('computing rock type variables')
    rock_type, rock_type_int = rocktype_var(lats, lons)
    dtf['rock_type'] = rock_type
    dtf['rock_type_int'] = rock_type_int

    print('computing soil thickness')
    soil_thick = soilthickness_vars(lats, lons)
    dtf['soil_thick'] = soil_thick

    dtf.to_csv("../results/ecovar_data_GLM.csv")

    print("Done! Time of execution " + str(np.round((time.time() - start_time) / 60, 2)) + ' minutes')


if __name__ == "__main__":
    main()
