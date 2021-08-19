#!/usr/bin/env python

import argparse
import glob
import time
import os
import pointing_camera.pc_proc as pipeline

default_data_dir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/el_nino'

def proc_1night(year, month, day):
    # year, month, day should all be integers

    year_string = str(year)
    month_string = str(month).zfill(2)
    day_string = str(day).zfill(2)

    search_dir = default_data_dir + '/' + year_string + '-' + \
                 month_string + '/' + year_string + '-' + \
                 month_string + '-' + day_string

    assert(os.path.exists(search_dir))

    flist_wcs = glob.glob(search_dir + '/*.wcs')
    # sort this file list? shuffle it in a reproducibly random way?

    print(len(flist_wcs))

    flist_fits = [f_wcs.replace('.wcs', '.fits') for f_wcs in flist_wcs]

    outdir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/proc_1night'

    for f in flist_fits:
        assert(os.path.exists(f))
        try:
            pipeline.pc_proc(f, outdir=outdir, dont_write_detrended=True,
                             skip_checkplot=False, nightly_subdir=True,
                             send_redis=False, one_aper=True)
        except:
            print('PROCESSING FAILURE: ' + f)
