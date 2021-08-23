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

default_data_dir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/reduced/v0005' # NERSC

files_processed = []
skies_table = None
zps_table = None

def quadrant_colors():
    # meant to be a small utility
    # quadrant 0 means all 4 quadrants combined

    colors = {0: 'b', 1: 'm', 2: 'r', 3: 'g', 4: 'c'}

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
        colname = 'mean_adu_'
        if q != 0:
            colname += 'quad' + str(q) + '_'
        colname += 'per_s'
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

def _proc_new_files(data_dir=default_data_dir, outdir='.', clobber=True):

    print('Checking for new reduction outputs...')

    assert(os.path.exists(data_dir))
    
    global files_processed
    global skies_table
    global zps_table

    flist_sky = glob.glob(data_dir + '/????????/*sky.fits')
    flist_zp = glob.glob(data_dir + '/????????/*zeropoints.fits')

    if (len(flist_sky) == 0) or (len(flist_zp) == 0):
        print('Nothing to plot yet...')
        return

    # figure out new list of sky.fits and zeropoints.fits files

    flist_sky_new = set(flist_sky) - set(files_processed)
    flist_zeropoints_new = set(flist_zp) - set(files_processed)

    if len(flist_sky_new) > 0:
        new_skies = _read_concat_tables(flist_sky_new)
        if skies_table is not None:
            skies_table = vstack([skies_table, new_skies])
        else:
            skies_table = new_skies
        
    if len(flist_zeropoints_new) > 0:
        new_zps = _read_concat_tables(flist_zeropoints_new)
        new_zps = new_zps[(new_zps['aper_ind'] == 1) &
                          (new_zps['quadrant'] == 0)]

        if zps_table is not None:
            zps_table = vstack([zps_table, new_zps])
        else:
            zps_table = new_zps
    
    # update skies_table and zps_table by appending new rows
    # to old rows
    
    files_processed = flist_sky + flist_zp # what if a file gets deleted though?

    if (len(flist_sky_new) > 0) or (len(flist_zeropoints_new) > 0):
        plt.cla()
        _twopanel(skies_table, zps_table, clobber=clobber, mjdrange=mjdrange)
    else:
        print('No new reduction output files found...')
    
def _watch(wait_seconds=5, data_dir=default_data_dir, outdir='.', clobber=True):

    while True:
        print('Waiting', wait_seconds, ' seconds')
        time.sleep(wait_seconds)
        _proc_new_files(data_dir=data_dir, outdir=outdir, clobber=clobber)

if __name__ == "__main__":
    descr = 'create nightly strip charts of zeropoint and sky brightness'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--start_night', default=None, type=str,
                        help="observing night to start with")

    parser.add_argument('--data_dir', default=default_data_dir, type=str,
                        help="directory with reduction pipeline output files")
    
    parser.add_argument('--outdir', default='.', type=str,
                        help="directory to write output images in")

    parser.add_argument('--wait_seconds', default=5, type=int,
                        help="plot creation polling interval in seconds")

    parser.add_argument('--no_clobber', default=False, action='store_true',
                        help="new nightly plot image for every update")

    args = parser.parse_args()

    # if start_night is None, then figure out the DESI observing night
    # based on the 

    assert(args.wait_seconds >= 1)

    _watch(wait_seconds=args.wait_seconds, data_dir=args.data_dir,
           outdir=args.outdir,
           clobber=(not args.no_clobber))
    
