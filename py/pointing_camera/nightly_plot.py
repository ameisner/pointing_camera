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

default_data_dir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/reduced/v0002' # NERSC

files_processed = []
skies_table = None
zps_table = None

def sky_subplot(tab, xticklabels=True, mjdrange=None):

    # build a list of datetime objects

    # factor this gathering of datetime list out into a helper
    # to be used both by the sky and zeropoints subplots
    datetimes = []
    for t in tab:
        tm = Time(t['MJD_OBS'], format='mjd')
        datetimes.append(tm.to_datetime())

    #plt.scatter(tab['MJD_OBS'], tab['MEAN_ADU'], s=10, edgecolor='none')

    # forgot to add column for MEAN_ADU_PER_S during reductions !!!!!!!
    plt.scatter(datetimes, -2.5*np.log10(tab['MEAN_ADU']/tab['TIME_SECONDS']),
                edgecolor='none', s=20, c='b')

    plt.ylabel(r'$-2.5\times$' + 'log' + r'$_{10}$' + '(ADU/pix/s)',
               fontsize=12)
    plt.xlabel(r'time (UT)', fontsize=12)

    if mjdrange is not None:
        date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
        date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
        plt.xlim((date_min, date_max))

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #ax.legend(['sky brightness'])
    plt.title('pointing camera sky brightness')
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

def zp_subplot(tab, xticklabels=False, mjdrange=None):

    print(tab.columns)
    datetimes = []
    for t in tab:
        tm = Time(t['MJD_OBS'], format='mjd')
        datetimes.append(tm.to_datetime())

    plt.scatter(datetimes, tab['ZP_ADU_PER_S'],
                edgecolor='none', s=20, c='b')

    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #ax.legend(['zeropoint'], loc='lower right')
    plt.title('pointing camera zeropoint')

    plt.ylabel('zeropoint (G' + r"$'$" + ' ; 1 ADU/s)')

    if mjdrange is not None:
        date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
        date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
        plt.xlim((date_min, date_max))
    
    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])
    
def _twopanel(skies_table, zps_table, clobber=True, mjdrange=None):

    plt.subplot(2, 1, 1)
    zp_subplot(zps_table, mjdrange=mjdrange)
    
    plt.subplot(2, 1, 2)
    sky_subplot(skies_table, mjdrange=mjdrange)

    # how to force the sky and zp panels to have the same
    # range of x values ... could imagine it getting confusing
    # if for whatever reason their x axes get out of alignment
    outname = 'night.png'

    # mainly for debugging purposes
    if not clobber:
        t = time.time()
        outname = outname.replace('.png', str(round(t)) + '.png')

    plt.savefig(outname, dpi=200, bbox_inches='tight')

def _read_concat_tables(flist):

    assert(len(flist) > 0)

    tables = []
    for i, f in enumerate(flist):
        print('READING: ' + f, ' ; ', i+1, ' of ', len(flist))
        assert(os.path.exists(f))
        t = Table(fits.getdata(f))
        if t[0]['MJD_OBS'] != 0:
            tables.append(t)

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
        new_zps = new_zps[(new_zps['APER_IND'] == 1) &
                          (new_zps['QUADRANT'] == 0)]

        if zps_table is not None:
            zps_table = vstack([zps_table, new_zps])
        else:
            zps_table = new_zps
    
    # update skies_table and zps_table by appending new rows
    # to old rows
    
    files_processed = flist_sky + flist_zp # what if a file gets deleted though?

    if (len(flist_sky_new) > 0) or (len(flist_zeropoints_new) > 0):

        mjdrange = [min(np.min(zps_table['MJD_OBS']),
                        np.min(skies_table['MJD_OBS'])),
                    max(np.max(zps_table['MJD_OBS']),
                        np.max(skies_table['MJD_OBS']))]
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
    
