import time
import numpy as np
import xarray as xr

import os
from mfsetup import load_modelgrid


def get_monthly_means(ncfile, outfile, filter=None, 
                      check_results=True, variable='net_infiltration'):
    get_monthly_values(ncfile, outfile, filter=filter, 
                      check_results=check_results, stat='mean', variable=variable)


def get_monthly_values(ncfile, outfile, filter=None, 
                       check_results=True, stat='mean', variable='net_infiltration'):
    """Get monthly mean values for a netcdf file
    with a time variable 'time' and write them to another netcdf file.
    Use dask if it is installed to limit memory required and
    distribute to cores.

    Parameters
    ----------
    ncfile : str
        Path to input netcdf file. Assumed to have variable 'net_infiltration'
        that is compressed in the output file.
    outfile : str
        Path of output netcdf file.
    filter : tuple of floats
        Bounding box of area to subset from ncfile (left, bottom, right, top).
        By default, none (entire area is retained in output file).
    check_results : bool, optional
        Check that the monthly values written to the output file match those in memory.
    variable : str
        Variable of interest in NetCDF file (e.g. 'net_infiltration')
    stat : str
        How the monthly values should be aggregated. Can use any of the pandas/xarray
        methods (e.g. 'mean' for ds.resample(time='MS').mean())
    """

    t0 = time.time()
    try:
        import dask
        print('using dask...')
        chunks = {'time': 30}
    except:
        chunks = None
    print('opening {}...'.format(ncfile))
    with xr.open_dataset(ncfile, chunks=chunks) as ds:

        if filter is not None:
            print('culling to bounding box: {}'.format(filter))
            l, b, r, t = filter
            ds = ds.sel(x=slice(l, r), y=slice(t, b))
        print('Aggregating to monthly values...')
        # average the data by month;
        # resampled data has start of each month as the timestamp
        monthly = getattr(ds.resample(time='MS'), stat)()
        if 'crs' in ds.variables:
            monthly['crs'] = ds['crs']
        monthly.attrs = ds.attrs
        monthly[variable].attrs = ds[variable].attrs
        monthly.to_netcdf(outfile, format='netcdf4', engine='netcdf4',  # engine='h5netcdf',
                          encoding={variable: {'zlib': True, 'complevel': 4,
                                                         'dtype': 'float32',  # 'scale_factor': 0.01,
                                                         '_FillValue': -9999,
                                                         # 'least_significant_digit': 8
                                                         }},
                          # encoding={'net_infiltration': {'dtype': 'int16', 'scale_factor': 0.1, '_FillValue': -9999}},
                          # encoding={'net_infiltration': {'compression': 'gzip',
                          #                               'compression_opts': 9}}
                          )
        print('wrote {}'.format(outfile))
        print("took {:.2f}s\n".format(time.time() - t0))
        if check_results:
            print("checking compressed file")
            with xr.open_dataset(outfile) as written:
                for i in range(0, len(monthly.time)):
                    np.testing.assert_equal(monthly[variable][i, :].values,
                                            written[variable][i, :].values)
    print("finished in {:.2f}s\n".format(time.time() - t0))