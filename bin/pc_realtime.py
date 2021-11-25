#!/usr/bin/env python

import argparse
import glob
import time
import os
import pointing_camera.pc_proc as pipeline

default_data_dir = '/global/cscratch1/sd/ameisner/pointing_camera/nino' # NERSC

files_processed = []

def _reduce_new_files(flist_wcs, outdir='.', dont_send_redis=False):

    for f_wcs in flist_wcs:
        # check that the corresponding .fits raw image exists
        # if there's a .wcs file with no corresponding raw .fits image
        # then something is badly wrong...

        f_fits = f_wcs.replace('.wcs', '.fits')

        if not os.path.exists(f_fits):
            print('WCS file with no corresponding raw image?? skipping')
            continue

        # call the reduction pipeline
        print('Reducing ' + f_fits)

        try:
            pipeline.pc_proc(f_fits, outdir=outdir, dont_write_detrended=True,
                             skip_checkplot=False, nightly_subdir=True, send_redis=(not dont_send_redis),
                             one_aper=True, nmp=8)
        except:
            print('PROCESSING FAILURE: ' + f_fits)

def _proc_new_files(data_dir=default_data_dir, outdir='.', dont_send_redis=False):

    print('Checking for new .wcs files...')

    assert(os.path.exists(data_dir))

    global files_processed

    flist_wcs = glob.glob(data_dir + '/*.wcs')

    if len(flist_wcs) == 0:
        print('No data to reduce yet...')
        return

    # dereference symlink to get actual location on disk
    flist_wcs = [(os.readlink(_f) if os.path.islink(_f) else _f) for _f in flist_wcs]

    flist_wcs_new = set(flist_wcs) - set(files_processed)

    if len(flist_wcs_new) > 0:
        flist_wcs_new = list(flist_wcs_new)
        flist_wcs_new.sort()

        _reduce_new_files(flist_wcs_new, outdir=outdir, dont_send_redis=dont_send_redis)

        files_processed = files_processed + flist_wcs_new

def _watch(wait_seconds=5, data_dir=default_data_dir, outdir='.',
           dont_send_redis=False):

    while True:
        print('Waiting', wait_seconds, ' seconds')
        time.sleep(wait_seconds)
        _proc_new_files(data_dir=data_dir, outdir=outdir, dont_send_redis=dont_send_redis)

def _do_veto(fname):
    assert(os.path.exists(fname))

    global files_processed

    f = open(fname, 'r')

    veto_list = f.readlines()

    veto_list = [v.replace('\n', '') for v in veto_list]

    files_processed = veto_list

if __name__ == "__main__":
    descr = 'process new pointing camera astrometry mode images in real time'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('--start_night', default=None, type=str,
                        help="observing night to start with")

    parser.add_argument('--data_dir', default=default_data_dir, type=str,
                        help="directory with raw images and WCS files")

    parser.add_argument('--outdir', default='.', type=str,
                        help="directory to write pipeline outputs in")

    parser.add_argument('--wait_seconds', default=5, type=int,
                        help="polling interval in seconds")

    parser.add_argument('--veto_list', default=None, type=str,
                        help="list of files not to process")

    parser.add_argument('--dont_send_redis', default=False,
                        action='store_true',
                        help="don't send results to Redis")

    args = parser.parse_args()

    assert(args.wait_seconds >= 1)

    if args.veto_list is not None:
        _do_veto(args.veto_list)

    _watch(wait_seconds=args.wait_seconds, data_dir=args.data_dir,
           outdir=args.outdir, dont_send_redis=args.dont_send_redis)
