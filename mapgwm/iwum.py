"""
Code for preprocessing output from the Irrigation Water Use Model (IWUM)
see Wilson (2020) for a description of IWUM
"""
from pathlib import Path
import os
from shapely.geometry import Point
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from gisutils import project
from gisutils.raster import get_values_at_points
from mfsetup.units import convert_length_units, convert_volume_units
from mapgwm.plot import format_xtick_labels
from mapgwm.utils import cull_data_to_active_area


def get_node(i, j, ncol):
    """Assign a unique, consecutive number for each i, j location"""
    return i * ncol + j


vol_suffix = {'feet': 'ft3',
              'meters': 'm3',
              1: 'ft3',
              2: 'm3'
              }


def plot_iwum_output(ncfile,
                     flux_variable='value', output_path='.'):
    """Make a plot of iwum output in mgal/day
    for comparison with subsequent datasets.
    """
    outpath = Path(output_path)
    ncfile = Path(ncfile)

    ds = xr.open_dataset(ncfile)
    time_variable = [k for k in ds.coords.keys() if k.lower() not in {'x', 'y'}][0]
    xydims = tuple([i for i, len in enumerate(ds[flux_variable].shape)
                    if len != ds[time_variable].shape[0]])
    ts = ds[flux_variable][:, :, :].sum(axis=xydims).to_pandas()
    if ts.index.dtype == np.object:
        ts.index = pd.to_datetime(ts.index)

    ndays = ts.index.days_in_month.tolist()
    df = pd.DataFrame(ts, columns=['m3'])
    df['m3d'] = df['m3'] / ndays # convert volumes to daily rate

    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax = df['m3d'].plot.bar(ax=ax)
    ax.set_ylabel('Cubic meters per day')
    ymin, ymax = ax.get_ylim()
    ax2 = ax.twinx()
    to_mg = convert_volume_units('m3', 'mgal')
    ax2.set_ylim(ymin * to_mg, ymax * to_mg)
    ax2.set_ylabel('Million gallons per day')

    # can't use .mean(),
    # because periods with 0 pumping may not be included
    mean_mgd = df['m3'].sum() * to_mg / np.sum(ndays)
    ax2.axhline(mean_mgd, c='r')
    ax2.text(0.75, 0.9, 'Mean: {:,.0f} mgal/day'.format(mean_mgd),
             transform=ax.transAxes)

    # format the tick labels
    maxlabels = 30
    xticklabels = df.index.strftime('%Y-%m-%d').tolist()
    stride = int(np.floor(len(xticklabels) / maxlabels))
    stride = stride if stride > 1 else 1
    formatted_labels = []
    for label in xticklabels[::stride]:
        formatted_labels += [label] + [''] * (stride - 1)
    junk = ax.set_xticklabels(formatted_labels)

    # record the file name and last modified date
    ftime = pd.Timestamp(os.path.getmtime(ncfile), unit='s')
    ax2.text(0.02, 0.98, '{}\n{}'.format(ncfile,
                                         ftime.strftime('%Y-%m-%d')),
             va='top', fontsize=8,
             transform=ax.transAxes)

    # annotate the bars with the values
    for i, p in enumerate(ax.patches):
        value = '{:,.0f}'.format(to_mg * p.get_height())
        ax.annotate(value, (p.get_x() * 1.01, p.get_height() * 1.01),
                    ha='center', fontsize=8)

    outfile = outpath / (ncfile.stem + '_{}.pdf'.format(ftime.strftime('%Y-%m-%d')))
    plt.savefig(outfile)
    print('wrote {}'.format(outfile))
    plt.close()


def preprocess_iwum_pumping(ncfile,
                            start_date=None,
                            end_date=None,
                            active_area=None,
                            active_area_id_column=None,
                            active_area_feature_id=None,
                            estimated_production_zone_top=None,
                            estimated_production_zone_botm=None,
                            flux_variable='value',
                            nc_crs=5070,
                            dest_crs=5070,
                            nc_length_units='meters',
                            estimated_production_surface_units='meters',
                            model_length_units='meters',
                            outfile=None):
    """Get pumping from the Irrigation Water Use Model (IWUM; Wilson, 2020) output and
    assign open interval information, using raster surfaces of the
    top and bottom of an estimated production zone.

    Parameters
    ----------
    ncfile : file path
        NetCDF output from Irrigation Water Use Model
    start_date : str
        Cull data before this date.
    end_date : str
        Cull data after this date.
    active_area : str
        Shapefile with polygon to cull observations to. Automatically reprojected
        to dest_crs if the shapefile includes a .prj file.
        by default, None.
    active_area_id_column : str, optional
        Column in active_area with feature ids.
        By default, None, in which case all features are used.
    active_area_feature_id : str, optional
        ID of feature to use for active area
        By default, None, in which case all features are used.
    estimated_production_zone_top : file path
        Raster surface for assigning screen tops
    estimated_production_zone_botm : file path
        Raster surface for assigning screen bottoms
    flux_variable : str
        Varible in ncfile for pumping fluxes. Fluxes are assumed to
        represent total volumes for each time period.
    nc_crs : obj
        Coordinate Reference System (CRS) of ncfile.
        A Python int, dict, str, or pyproj.crs.CRS instance
        passed to the pyproj.crs.from_user_input
        See http://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input.
        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class
    nc_length_units : str, {'meters', 'ft', etc.}
        Length units of pumped volumes in ncfile
    estimated_production_surface_units : str, {'meters', 'ft', etc.}
        Length units of elevations in estimated production surface rasters.
    model_length_units : str, {'meters', 'ft', etc.}
        Length units of model.
    outfile : csv file for output table

    Returns
    -------
    df : DataFrame
        Table of pumping rates in m3/day, location
        and open interval information.

        Columns:

        ============== ================================================
        site_no        index position of pumping rate in ncfile grid
        x              x-coordinate in `dest_crs`
        y              y-coordinate in `dest_crs`
        start_datetime start date of pumping period
        end_datetime   end date of pumping period
        screen_top     screen top elevation, in `model_length_units`
        screen_botm    screen bottom elevation, in `model_length_units`
        q              pumping rate, in model units
        geometry       shapely Point object representing location
        ============== ================================================

    Notes
    -----
    * Time units are assumed to be days.
    * Fluxes are assumed to represent total volumes for each time period
      indicated by the differences between successive values along the time axis of ncfile.
    """
    ds = xr.open_dataset(ncfile)
    time_variable = [k for k in ds.coords.keys() if k.lower() not in {'x', 'y'}][0]
    ds_x, ds_y = np.meshgrid(ds['x'], ds['y'])

    # original values are in m3, in each 1 mi2 cell
    # can leave in m3 if reassigning to 1km grid as point values
    length_conversion = convert_volume_units(nc_length_units,
                                             model_length_units) ** 3
    unit_suffix = vol_suffix[model_length_units] + 'd'
    flux_col = 'q'  # 'flux_{}'.format(unit_suffix)  # output field name for fluxes

    # get top/botm elevations
    est_screen_top = None
    est_screen_botm = None
    if estimated_production_zone_top is not None and \
            estimated_production_zone_botm is not None:
        surf_unit_conversion = convert_length_units(estimated_production_surface_units,
                                                    model_length_units)
        est_screen_top = get_values_at_points(estimated_production_zone_top, ds_x, ds_y,
                                                points_crs=nc_crs)
        est_screen_top *= surf_unit_conversion
        est_screen_botm = get_values_at_points(estimated_production_zone_botm, ds_x, ds_y,
                                                 points_crs=nc_crs)
        est_screen_botm *= surf_unit_conversion

        # in any places where screen top is less than the screen botm,
        # set both at the mean
        loc = est_screen_top < est_screen_botm
        means = np.mean([est_screen_top, est_screen_botm], axis=0)
        est_screen_top[loc] = means[loc]
        est_screen_botm[loc] = means[loc]
        print(f'Reset screen top and bottom to mean elevation at {loc.ravel().sum()} '
              f'locations where screen top was < screen bottom')

    dfs = []
    times = pd.DatetimeIndex(ds[time_variable].loc[start_date:end_date].values)
    for n, period_start_date in enumerate(times):

        # for each time entry, get the data
        kwargs = {time_variable: period_start_date}
        arr = ds[flux_variable].sel(**kwargs).values

        # make sure pumping sign is  negative
        # based on assumption that values are mostly abstraction
        if arr.sum() > 0:
            arr *= -1

        # set up a dataframe
        data = {'site_no': np.arange(ds_x.size),
                'x': ds_x.ravel(),
                'y': ds_y.ravel(),
                 }
        if est_screen_top is not None and est_screen_botm is not None:
            data.update({'screen_top': est_screen_top.ravel(),
                         'screen_botm': est_screen_botm.ravel()
                         }
                        )
        df = pd.DataFrame(data)
        df['start_datetime'] = period_start_date

        # get the end_date, handling last entry
        if n + 1 < len(times):
            period_end_date = times[n + 1]
        else:
            # set end date for last period on previous period length
            last_start = dfs[-1]['start_datetime'].values[0]
            ndays = (pd.Timestamp(period_start_date) -
                     pd.Timestamp(last_start)).days
            period_end_date = period_start_date + pd.Timedelta(ndays, unit='d')

        # convert the time units
        ndays = (pd.Timestamp(period_end_date) -
                 pd.Timestamp(period_start_date)).days
        assert ndays > 0, "period_end_date {} is before period_start_date {}"\
            .format(period_end_date, period_start_date)
        time_conversion = 1 / ndays  # original quantities are volumes for the time period

        # time indexing in pandas is through last value
        period_end_date = pd.Timestamp(period_end_date) - pd.Timedelta(1, unit='d')
        df['end_datetime'] = period_end_date
        df[flux_col] = arr.ravel() * length_conversion * time_conversion

        # only includes fluxes > 0
        df = df.loc[df[flux_col] < 0]

        dfs.append(df)
    df = pd.concat(dfs)

    # site number column (that would be unique from other integers from other data sources)
    df['site_no'] = [f'iwum_{node}' for node in df.site_no]

    # project the data to a destination crs, if provided
    # make a separate metadata dataframe with 1 row per location
    # to avoid redundant operations
    metadata = df.groupby('site_no').first().reset_index()[['site_no', 'x', 'y']]
    metadata.index = metadata['site_no']
    x_pr, y_pr = project((metadata.x.values, metadata.y.values), nc_crs, dest_crs)
    metadata['x'], metadata['y'] = x_pr, y_pr
    metadata['geometry'] = [Point(x, y) for x, y in zip(x_pr, y_pr)]

    # cull the data to the model area, if provided
    if active_area is not None:
        df, metadata = cull_data_to_active_area(df, active_area,
                                      active_area_id_column,
                                      active_area_feature_id,
                                      data_crs=dest_crs, metadata=metadata)

    # update data with x,y values projected in metadata
    x = dict(zip(metadata.site_no, metadata.x))
    y = dict(zip(metadata.site_no, metadata.y))
    df['x'] = [x[sn] for sn in df.site_no]
    df['y'] = [y[sn] for sn in df.site_no]
    if outfile is not None:
        outfile = Path(outfile)
        df.to_csv(outfile, index=False, float_format='%g')
        print('wrote {}'.format(outfile))

        # Make a plot of iwum output in mgal/day
        out_pdf_path = outfile.parent / 'plots'
        out_pdf_path.mkdir(exist_ok=True)
        plot_iwum_output(ncfile, flux_variable=flux_variable, outpath=out_pdf_path)

    return df


def plot_iwum_output(ncfile, flux_variable='value', outpath='.'):
    """Make a plot of iwum output in mgal/day
    for comparison with subsequent datasets.
    """

    ds = xr.open_dataset(ncfile)
    time_variable = [k for k in ds.coords.keys() if k.lower() not in {'x', 'y'}][0]
    xydims = tuple([i for i, len in enumerate(ds[flux_variable].shape)
                    if len != ds[time_variable].shape[0]])
    ts = ds[flux_variable][:, :, :].sum(axis=xydims).to_pandas()
    if ts.index.dtype == np.object:
        ts.index = pd.to_datetime(ts.index)

    ndays = pd.to_timedelta(np.diff(ts.index)).days.tolist()
    ndays.append(ndays[-1])  # pad the last time period
    df = pd.DataFrame(ts, columns=['m3'])
    df['m3d'] = df['m3'] / ndays # convert volumes to daily rate

    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax = df['m3d'].plot.bar(ax=ax)
    ax.set_ylabel('Cubic meters per day')
    ymin, ymax = ax.get_ylim()
    ax2 = ax.twinx()
    to_mg = convert_volume_units('m3', 'mgal')
    ax2.set_ylim(ymin * to_mg, ymax * to_mg)
    ax2.set_ylabel('Million gallons per day')

    # can't use .mean(),
    # because periods with 0 pumping may not be included
    mean_mgd = df['m3'].sum() * to_mg / np.sum(ndays)
    ax2.axhline(mean_mgd, c='r')
    ax2.text(0.75, 0.9, 'Mean: {:,.0f} mgal/day'.format(mean_mgd),
             transform=ax.transAxes)

    # format the tick labels
    format_xtick_labels(df, ax, maxlabels=30, date_format='%Y-%m-%d')
    #maxlabels = 30
    #xticklabels = df.index.strftime('%Y-%m-%d').tolist()
    #stride = max(int(np.floor(len(xticklabels) / maxlabels)), 1)
    #formatted_labels = []
    #for label in xticklabels[::stride]:
    #    formatted_labels += [label] + [''] * (stride - 1)
    #formatted_labels = formatted_labels[:len(xticklabels)]
    #junk = ax.set_xticklabels(formatted_labels)

    # record the file name and last modified date
    ftime = pd.Timestamp(os.path.getmtime(ncfile), unit='s')
    ax2.text(0.02, 0.98, '{}\n{}'.format(ncfile,
                                         ftime.strftime('%Y-%m-%d')),
             va='top', fontsize=8,
             transform=ax.transAxes)

    # annotate the bars with the values
    for i, p in enumerate(ax.patches):
        value = '{:,.0f}'.format(to_mg * p.get_height())
        ax.annotate(value, (p.get_x() * 1.01, p.get_height() * 1.01),
                    ha='center', fontsize=8)
    ncfile = Path(ncfile)
    ftime = pd.Timestamp(ncfile.stat().st_mtime, unit='s')
    outfile = Path(outpath, f'{ncfile.name}_{ftime:%Y-%m-%d}.pdf')
    plt.savefig(outfile)
    print('wrote {}'.format(outfile))
    plt.close()