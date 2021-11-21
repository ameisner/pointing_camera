#!/usr/bin/env python

import argparse
from datetime import datetime
import os
import time
from pointing_camera.exposure import PC_exposure
import pointing_camera.util as util
import pointing_camera.io as io
import pointing_camera.zp as zp
import pointing_camera.common as common

def pc_proc(fname_in, outdir=None, dont_write_detrended=False,
            skip_checkplot=False, nightly_subdir=False, send_redis=False,
            one_aper=False, bg_sigclip=False, nmp=None, max_n_stars=3000,
            pm_corr=False):

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

    exp.has_dome = util.flag_dome_vignetting(exp.detrended, exp.time_seconds)

    sky = util.sky_summary_table(exp)

    cat = util.pc_phot(exp, one_aper=one_aper, bg_sigclip=bg_sigclip,
                       nmp=nmp, max_n_stars=max_n_stars,
                       pm_corr=pm_corr)

    # intentionally don't pass nmp to zps.calc_many_zp, since doing
    # so didn't appear to provide any speed-up; could revisit later
    zps = zp.calc_many_zps(cat, exp, one_aper=one_aper, checkplot=(not skip_checkplot))

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

        # if options are added to skip certain steps
        # (like sky mag estimation), this will need
        # to be adjusted accordingly
        io.write_bintables_mef(cat, zps, sky, exp, outdir)

        if not skip_checkplot:
            io.save_zp_checkplot(exp, outdir)

        if send_redis:
            print('Attempting to send results to redis...')
            _aper_ind = 0 if one_aper else 1
            util.send_redis(exp, zps[zps['aper_ind'] == _aper_ind], sky[0])

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

    parser.add_argument('--one_aper', default=False, action='store_true',
                        help="only do aperture photometry for one aperture size")

    parser.add_argument('--bg_sigclip', default=False, action='store_true',
                        help="sigma clipping for background annulus median")

    parser.add_argument('--multiproc', default=None, type=int,
                        help="number of threads for multiprocessing")

    parser.add_argument('--max_n_stars', default=3000, type=int,
                        help="limit analysis to brightest max_n_stars Gaia stars")

    parser.add_argument('--pm_corr', default=False, action='store_true',
                        help="make Gaia proper motion corrections based on MJD")

    args = parser.parse_args()

    # basic checks on requested number of multiprocessing threads
    if args.multiproc is not None:
        par = common.pc_params()
        assert(args.multiproc > 1)
        assert(args.multiproc <= par['ncpus'])

    pc_proc(args.fname_in[0], outdir=args.outdir,
            dont_write_detrended=args.dont_write_detrended,
            skip_checkplot=args.skip_checkplot,
            nightly_subdir=args.nightly_subdir, send_redis=args.send_redis,
            one_aper=args.one_aper, bg_sigclip=args.bg_sigclip,
            nmp=args.multiproc, max_n_stars=args.max_n_stars, pm_corr=args.pm_corr)
