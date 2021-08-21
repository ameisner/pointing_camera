import astropy.io.fits as fits
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.dates as md

fname = '/global/homes/a/ameisner/pointing/etc/sky_monitor-20210609_thru_20210710.fits'

def plot_sky_mon(night, mjdrange, xticklabels=True, title_extra='',
                 markersize=10, do_xlabel=True):

    tab = fits.getdata(fname)

    tab = tab[tab['night'] == night]

    datetimes = []
    if len(tab) == 0:
        mags = [-100, -100] # off the plot

        for val in mjdrange:
            tm = Time(val, format='mjd')
            datetimes.append(tm.to_datetime())

        plt.scatter(datetimes, mags)
        ymin = -7
        ymax = 1.5
    else:
        for row in tab:
            tm = Time(row['MJD'], format='mjd')
            datetimes.append(tm.to_datetime())

        mags = -2.5*np.log10(tab['value'])
        plt.scatter(datetimes, mags,
                    edgecolor='none', s=markersize, c='k')

        ymax = min(1.5, np.max(mags)+0.2)
        ymin = max(-7, np.min(mags)-0.2)

    plt.ylim((ymin, ymax))

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    
    title = 'sky monitor sky brightness' + title_extra
    plt.title(title)
    plt.ylabel(r'$-2.5\times$' + 'log' + r'$_{10}$' + '(value)')
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

    if do_xlabel:
        plt.xlabel(r'time (UT)', fontsize=12)

    date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
    date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
    plt.xlim((date_min, date_max))
