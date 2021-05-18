#!/usr/bin/env python

import argparse
from datetime import datetime
import os
import time
from pointing_camera.exposure import PC_exposure
import pointing_camera.util as util
import pointing_camera.io as io
import pointing_camera.zp as zp

def pc_proc(fname_in, outdir=None, dont_write_detrended=False,
            skip_checkplot=False, nightly_subdir=False, send_redis=False):

    print('Starting pointing camera reduction pipeline at: ' +
          str(datetime.utcnow()) + ' UTC')

    t0 = time.time()

    write_outputs = (outdir is not None)
    
    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')

    assert(os.path.exists(fname_in))

    exp = PC_exposure(fname_in)

    util.detrend_pc(exp)

    sky = util.sky_summary_table(exp)

    cat = util.pc_phot(exp)

    zps = zp.calc_many_zps(cat, exp, checkplot=(not skip_checkplot))

    if write_outputs:
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if nightly_subdir:
            # should add sanity checks on obs_night here
            outdir = os.path.join(outdir, exp.obs_night)
            if not os.path.exists(outdir):
                print('Creating subdirectory for observing night = ' + \
                      exp.obs_night)
                os.mkdir(outdir)

        if not dont_write_detrended:
            io.write_image_level_outputs(exp, outdir)

        io.write_sky_summary(sky, exp, outdir)
        io.write_source_catalog(cat, exp, outdir)
        io.write_zeropoints_table(zps, exp, outdir)

        if not skip_checkplot:
            io.save_zp_checkplot(exp, outdir)

        if send_redis:
            util.send_redis(exp, zps[(zps['aper_ind'] == 1) & (zps['quadrant'] == 0)]['zp_adu_per_s'][0], sky['mean_adu_per_s'][0])

    dt = time.time() - t0
    print('pointing camera reduction pipeline took ' + '{:.2f}'.format(dt) +
          ' seconds')
    print('pointing camera reduction pipeline completed at: ' +
          str(datetime.utcnow()) + ' UTC') 

if __name__ == "__main__":
    descr = 'run the pointing camera reduction pipeline on an exposure'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('fname_in', type=str, nargs=1,
                        help="pointing camera raw image file name")

    parser.add_argument('--outdir', default=None, type=str,
                        help="directory to write outputs in")

    parser.add_argument('--dont_write_detrended', default=False,
                        action='store_true',
                        help="don't write detrended image")

    parser.add_argument('--skip_checkplot', default=False,
                        action='store_true',
                        help="don't create a checkplot")

    parser.add_argument('--nightly_subdir', default=False, action='store_true',
                        help="create output subdirectories per observing night")

    parser.add_argument('--send_redis', default=False, action='store_true',
                        help="send results to redis")
    
    args = parser.parse_args()
    
    pc_proc(args.fname_in[0], outdir=args.outdir,
            dont_write_detrended=args.dont_write_detrended,
            skip_checkplot=args.skip_checkplot,
            nightly_subdir=args.nightly_subdir, send_redis=args.send_redis)
