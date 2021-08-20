import astropy.io.fits as fits
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.dates as md

# utilities for accessing / plotting DESI GFA data

def read_gfa():

    fname = '/global/cfs/cdirs/desi/users/ameisner/GFA/conditions/offline_all_guide_ccds_SV3-thru_20210710.fits'

    # get ex=3
    tab = fits.getdata(fname, ext=3)

    return tab
    

def plot_gfa_fwhm(gfa, markersize=2, mjdrange=None, xticklabels=True,
                  do_xlabel=True):

    assert(len(np.unique(gfa['NIGHT'])) == 1)

    datetimes = []
    for row in gfa:
        tm = Time(row['MJD'], format='mjd')
        datetimes.append(tm.to_datetime())

    plt.scatter(datetimes, gfa['FWHM_ASEC'],
                edgecolor='none', s=markersize, c='k')


    if mjdrange is not None:
        date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
        date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
        plt.xlim((date_min, date_max))

    plt.ylim((0, 3.5))

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    title = 'DESI GFA' # + title_extra


    plt.ylabel('r band FWHM (asec)')
    plt.title(title)
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

    if do_xlabel:
        plt.xlabel(r'time (UT)', fontsize=12)
