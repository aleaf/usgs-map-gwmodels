"""
Functions for making plots that compare model input to source data
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mfsetup.units import convert_volume_units, convert_time_units


def plot_wateruse(wel_files, perioddata, add_data=None,
                  wel_flux_col='q',
                  model_volume_units='$m^3$', model_time_units='day',
                  plot_volume_units='mgal', plot_time_units='day',
                  annual_means=False,
                  start_date=None, end_date=None,
                  outfile=None):
    """

    Parameters
    ----------
    wel_files :
        A header line with column names is assumed. For example:
        #k,i,j,q,boundname

    perioddata :
    add_data :
    model_volume_units :
    model_time_units :
    plot_volume_units :
    plot_time_units :

    Returns
    -------

    """

    # read the stress period information
    if not isinstance(perioddata, pd.DataFrame):
        perioddata = pd.read_csv(perioddata)
    else:
        perioddata = perioddata.copy()
    perioddata.index = perioddata['per']

    dfs = []
    for i, f in wel_files.items():
        df = pd.read_csv(f, delim_whitespace=True)
        df.columns = [c.strip('#') for c in df.columns]
        df['per'] = i
        df['start_datetime'] = perioddata.loc[i, 'start_datetime']
        df['end_datetime'] = perioddata.loc[i, 'end_datetime']
        dfs.append(df)
    df = pd.concat(dfs)

    # sum the model pumping by stress period
    period_sums = df.groupby('per').first()
    period_sums[wel_flux_col] = df.groupby('per')[wel_flux_col].sum()
    # fill nan values (from any periods without wel files) with 0s
    period_sums = period_sums.reindex(range(period_sums.index.max()))
    period_sums['start_datetime'] = perioddata['start_datetime']
    period_sums['end_datetime'] = perioddata['end_datetime']
    period_sums[wel_flux_col].fillna(0, inplace=True)
    period_sums.index = pd.to_datetime(period_sums['start_datetime'])
    period_sums['WEL package input'] = period_sums['q']
    period_sums = period_sums[['WEL package input', 'start_datetime', 'end_datetime']]

    # convert units
    model_vol_conv = convert_volume_units(model_volume_units, plot_volume_units)
    model_time_conv = convert_time_units(model_time_units, plot_time_units)
    model_conv = model_vol_conv * model_time_conv

    # plot any additional comparison data
    if add_data is not None:
        for label, items in add_data.items():
            # read the stress period information
            if not isinstance(items['data'], pd.DataFrame):
                items['data'] = pd.read_csv(items['data'])
            req_cols = {'q', 'start_datetime'}
            assert not req_cols.difference(items['data'].columns), \
                f"add_data: {label} data must have columns: {req_cols}"

            items['data']['start_datetime'] = pd.to_datetime(items['data']['start_datetime'])
            aux_period_sums = items['data'].groupby('start_datetime').first()
            aux_period_sums[label] = items['data'].groupby('start_datetime')['q'].sum()
            # fill nan values (from any periods without wel files) with 0s
            #aux_period_sums[label].fillna(0, inplace=True)
            aux_period_sums['start_datetime'] = aux_period_sums.index

            period_sums = period_sums.join(aux_period_sums[[label]], how='outer')
            j=2

    # forward fill nan WEL values values
    # (where other times may have been inserted)
    period_sums['WEL package input'] = period_sums['WEL package input'].ffill()
    for col in add_data.keys():
        period_sums[col] = period_sums[col].interpolate(limit_area='inside')
    #for col in 'WEL package input', '1900-2007_avg_pumping_from_meras22_m3.csv':
    #    period_sums[col] = period_sums[col].ffill()
    #period_sums = period_sums.resample('M').mean() #.ffill()

    # make a plot
    fig, ax = plt.subplots(figsize=(11, 8.5))
    if annual_means:
        period_sums = period_sums.groupby(period_sums.index.year).mean()
    ax = period_sums.loc[slice(start_date, end_date)].plot(ax=ax)
    units_text = f'{model_volume_units}/{model_time_units}'
    ax.set_ylabel(f'Pumpage, in {units_text}')
    ax.set_xlabel('')

    # second axis with another volume unit
    def second_axis_conversion(x):
        return x * model_conv

    def second_axis_conversion_r(x):
        return x * 1 / model_conv
    ax2 = ax.secondary_yaxis('right', functions=(second_axis_conversion,
                                                 second_axis_conversion_r))
    ax2.set_ylabel(f'Pumpage, in {plot_volume_units}/{plot_time_units}')
    #format_xtick_labels(period_sums, ax, maxlabels=30, date_format='%Y-%m-%d')
    h, l = ax.get_legend_handles_labels()
    means = (period_sums.mean(axis=0) * model_conv).to_dict()
    plot_units_text = f'{plot_volume_units}/{plot_time_units}'
    labels_with_means = []
    for label in l:
        new_label = label
        if label in means:
            new_label += f' (mean: {means[label]:,g} {plot_units_text})'
        labels_with_means.append(new_label)
    ax.legend(h, labels_with_means)

    if outfile is not None:
        outfile = Path(outfile)
        outfile_name = f"{outfile.stem}"
        if start_date is None:
            start_date = period_sums.index[0]
            if isinstance(period_sums.index[0], pd.Timestamp):
                start_date = start_date.strftime('%Y-%m-%d')
        outfile_name += f"_{start_date}"
        if end_date is None:
            end_date = period_sums.index[-1]
            if isinstance(period_sums.index[0], pd.Timestamp):
                end_date = end_date.strftime('%Y-%m-%d')
            outfile_name += f"_to_{end_date}"
        outfile_name = outfile_name + ".pdf"
        outfile = outfile.parent / outfile_name
        outfile.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile)
        plt.close()
        print(f'wrote {outfile}')
    else:
        return ax

    
def format_xtick_labels(df, ax, maxlabels=30, date_format='%Y-%m-%d'):
    """Clean up the xtick labels on a time axis.
    Cap the number of labels to maxlabels, and format
    dates to date_format.
    """
    xticklabels = df.index.strftime(date_format).tolist()
    stride = max(int(np.floor(len(xticklabels) / maxlabels)), 1)
    formatted_labels = []
    for label in xticklabels[::stride]:
        formatted_labels += [label] + [''] * (stride - 1)
    formatted_labels = formatted_labels[:len(xticklabels)]
    ax.set_xticklabels(formatted_labels)