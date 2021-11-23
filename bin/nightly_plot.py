#!/usr/bin/env python

import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.dates as md
import numpy as np
import argparse
import glob
import time
from astropy.table import Table, vstack
import os
from multiprocessing import Pool

default_data_dir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/reduced/v0006' # NERSC

files_processed = []
skies_table = None
zps_table = None

def quadrant_colors():
    # meant to be a small utility
    # quadrant 0 means all 4 quadrants combined

    colors = {0: 'k', 1: 'm', 2: 'r', 3: 'g', 4: 'c'}

    return colors

def sky_subplot(tab, xticklabels=True, mjdrange=None, markersize=20,
                title_extra='', skip_q0=False, do_xlabel=True):

    # build a list of datetime objects

    # factor this gathering of datetime list out into a helper
    # to be used both by the sky and zeropoints subplots
    datetimes = []
    for t in tab:
        tm = Time(t['mjd_obs'], format='mjd')
        datetimes.append(tm.to_datetime())

    #plt.scatter(tab['MJD_OBS'], tab['MEAN_ADU'], s=10, edgecolor='none')

    # forgot to add column for MEAN_ADU_PER_S during reductions !!!!!!!

    colors = quadrant_colors()

    for q in [4, 3, 2, 1, 0]:
        if (q == 0) and (skip_q0):
            continue
        colname = 'median_adu_'
        if q != 0:
            colname += 'quad' + str(q)
        else:
            # for quadrant zero, use the science (DESI) FOV footprint
            # rather than the entire pointing camera FOV
            colname += 'sci'

        colname += '_per_s'
        plt.scatter(datetimes, -2.5*np.log10(tab[colname]),
                    edgecolor='none', s=markersize, c=colors[q])

    plt.ylabel(r'$-2.5\times$' + 'log' + r'$_{10}$' + '(ADU/pix/s)',
               fontsize=12)

    if do_xlabel:
        plt.xlabel(r'time (UT)', fontsize=12)

    if mjdrange is not None:
        date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
        date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
        plt.xlim((date_min, date_max))

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #ax.legend(['sky brightness'])
    title = 'pointing camera sky brightness' + title_extra
    plt.title(title)
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

def zp_subplot(tab, xticklabels=False, mjdrange=None, markersize=20,
               title_extra='', skip_q0=False, do_xlabel=False):

    colors = quadrant_colors()

    for q in [4, 3, 2, 1, 0]:
        if (q == 0) and (skip_q0):
            continue
        _tab = tab[tab['quadrant'] == q]

        # for quadrant zero, use the science (DESI) FOV footprint
        # rather than the entire pointing camera FOV
        if (q == 0):
            _tab = _tab[_tab['science_fov_only'] == 1]

        datetimes = []
        for t in _tab:
            tm = Time(t['mjd_obs'], format='mjd')
            datetimes.append(tm.to_datetime())

        plt.scatter(datetimes,
                    _tab['zp_adu_per_s'],
                    edgecolor='none', s=markersize, c=colors[q])

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #ax.legend(['zeropoint'], loc='lower right')

    title = 'pointing camera zeropoint'
    title += title_extra
    plt.title(title)

    plt.ylabel('zeropoint (G' + r"$'$" + ' ; 1 ADU/s)')

    if mjdrange is not None:
        date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
        date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
        plt.xlim((date_min, date_max))
    
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

    if do_xlabel:
        plt.xlabel(r'time (UT)', fontsize=12)

def _twopanel(skies_table, zps_table, clobber=True, save=True,
              markersize=20, title_extra='', skip_q0=False):

    mjdrange = [min(np.min(zps_table['mjd_obs']),
                    np.min(skies_table['mjd_obs'])),
                max(np.max(zps_table['mjd_obs']),
                    np.max(skies_table['mjd_obs']))]

    plt.subplot(2, 1, 1)
    zp_subplot(zps_table, mjdrange=mjdrange, markersize=markersize,
               title_extra=title_extra, skip_q0=skip_q0)
    
    plt.subplot(2, 1, 2)
    sky_subplot(skies_table, mjdrange=mjdrange, markersize=markersize,
                title_extra=title_extra, skip_q0=skip_q0)

    # how to force the sky and zp panels to have the same
    # range of x values ... could imagine it getting confusing
    # if for whatever reason their x axes get out of alignment

    if save:
        outname = 'night.png'

    # mainly for debugging purposes
        if not clobber:
            t = time.time()
            outname = outname.replace('.png', str(round(t)) + '.png')

        plt.savefig(outname, dpi=200, bbox_inches='tight')

def _read_one_table(fname, ext):
    print('READING: ' + fname)
    assert(os.path.exists(fname))
    t, h = fits.getdata(fname, ext=ext, header=True)

    t = Table(t)

    if 'AIRMASS' in h:
        if h['AIRMASS'] is not None:
            t['airmass'] = h['AIRMASS']
        else:
            t['airmass'] = -1
    else:
        t['airmass'] = -1

    return t
    
def _read_concat_tables(flist, ext=1, nmp=None):

    assert(len(flist) > 0)

    if nmp is None:
        tables = []
        for i, f in enumerate(flist):
            t = _read_one_table(f, ext)
            if t[0]['mjd_obs'] != 0:
                tables.append(t)
    else:
        print('Reading tables in parallel...')
        p = Pool(nmp)
        args = [(f, ext) for f in flist]
        tables = p.starmap(_read_one_table, args)
        p.close()
                
    print('Attempting to append tables...')
    if len(tables) == 1:
        result = tables[0]
    else:
        result = vstack(tables)
    print('Finished appending tables...')

    return result
