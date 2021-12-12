#!/usr/bin/env python

"""
pointing_camera.pc_proc
=====================

Main pipeline driver for processing pointing camera images.
"""

import argparse
from datetime import datetime
import os
import time
from pointing_camera.exposure import PC_exposure
import pointing_camera.util as util
import pointing_camera.io as io
import pointing_camera.zp as zp
import pointing_camera.common as common
import satellites

def pc_proc(fname_in, outdir=None, dont_write_detrended=False,
            skip_checkplot=False, nightly_subdir=False, send_redis=False,
            one_aper=False, bg_sigclip=False, nmp=None, max_n_stars=3000,
            pm_corr=False, skip_flatfield=False, sci_inst_name='desi',
            sci_fov_checkplot=False, check_tcs_motion=False,
            max_zp_radius=None, detect_streaks=False, plot_detrended=False,
            plot_streaks=False, plot_quiver=False):
    """
    Process one pointing camera image.

    Parameters
    ----------
        fname_in : str
            Full name of raw pointing camera image file to process.
        outdir : str, optional
            Full path of output directory. If not specified, outputs are not
            written. Default is None.
        dont_write_detrended : bool, optional
            Don't write out detrended image.
        skip_checkplot : bool, optional
            Skip making photometric zeropoint checkplot.
        nightly_subdir : bool, optional
            Create an observing night subdirectory within outdir, where
            this subdirectory name has the format YYYYMMDD.
        send_redis : bool, optional
            If True, send telemetry with reduction pipeline results to the
            Mayall's engineering database. This only makes sense to enable
            when running within the appropriate computing infrastructure
            on the mountain.
        one_aper : bool, optional
            To reduce run time, do aperture photometry for only one
            aperture radius. When only one aperture radius is used, it is
            taken to be the 'standard' aperture radius of 2.5 pixels.
        bg_sigclip : bool, optional
            Use sigma clipping when computing the background level. I believe
            that skipping the sigma clipping is meant to be an optimization
            toward minimizing run time.
        nmp : int, optional
            Number of threads for multiprocessing. Default is None,
            in which case multiprocessing is not used.
        max_n_stars : int, optional
            Maximum number of Gaia stars to recentroid/photometer. The
            idea is to cap the number of stars analyzed to avoid excessively
            long run times in dense stellar fields.
        pm_corr : bool, optional
            If True, correct Gaia stars for proper motion so that their
            positions match the epoch of pointing camera observation. This
            is False by default given that the recentroiding allowance should
            be large enough to accommodate basically all stellar proper motions.
        skip_flatfield : bool, optional
            If True, skip the flatfielding step of detrending. False by default.
            Could be useful for e.g., building a sky flat or star flat.
        sci_inst_name : str, optional
            The idea is not to hardcode 'desi' or 'DESI', in case this pipeline
            is run at some point with a different telescope and/or instrument.
        sci_fov_checkplot : bool, optional
            Set True to limit photometric zeropoint checkplot star sample to
            only those stars that fall within the science instrument's FOV.
        check_tcs_motion : bool, optional
            Set True to check header metadata ZPFLAG (which indicates
            telescope motion during the exposure), and abort the pipeline
            if this flag is set. TCS = Telescope Control System.
        max_zp_radius : int, optional
            Default is None. If set, this is the maximum radius in pixels
            relative to the detector center at which a star will be used
            in computing zeropoints. Currently, the pipeline still analyzes
            stars beyond this radius. A future optimization could be
            entirely ignoring (no centroiding or photometry) stars
            beyond this radius, in order to decrease runtime.
        detect_streaks : bool, optional
            If True, run satellite streak detection. This will increase
            the total run time by at least a couple of seconds (haven't
            checked if there are situations where the satellite streak
            detection takes an excessively long time, but there might be).
            At present, running streak detection will also require the
            astride streak detection package to be installed.
        plot_detrended : bool, optional
            If True, make and save a rendering of the detrended pointing
            camera image. Would be good to eventually allow for overplotting
            of detected satellite streaks on this rendering.
        plot_streaks : bool, optional
            If True, overplot any detected streaks on the detrended image
            rendering. If plot_detrended is False, then plot_streaks has no
            effect.
        plot_quiver : bool, optional
            If True, make and save a quiver plot of the centroid
            shifts relative to the astrometry.net WCS solution.

    """

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

    if check_tcs_motion and exp.did_telescope_move():
        print('Telescope moved during exposure, abandoning further analysis')
        return

    util.detrend_pc(exp, skip_flatfield=skip_flatfield)

    sky = util.sky_summary_table(exp)

    cat = util.pc_phot(exp, one_aper=one_aper, bg_sigclip=bg_sigclip,
                       nmp=nmp, max_n_stars=max_n_stars,
                       pm_corr=pm_corr, max_zp_radius=max_zp_radius)

    # intentionally don't pass nmp to zps.calc_many_zp, since doing
    # so didn't appear to provide any speed-up; could revisit later
    zps = zp.calc_many_zps(cat, exp, one_aper=one_aper,
                           checkplot=(not skip_checkplot),
                           sci_fov_checkplot=sci_fov_checkplot,
                           max_zp_radius=max_zp_radius)

    if detect_streaks:
    # probably want to put a timeout handler around this eventually
        streaks = satellites.detect_streaks(exp)

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

        if detect_streaks:
            io.write_streaks(exp, streaks, outdir)

        if plot_detrended:
            io.plot_detrended(exp, outdir, plot_streaks=plot_streaks)

        if plot_quiver:
            io.save_quiver_plot(exp, cat, outdir)

        if send_redis:
            print('Attempting to send results to redis...')
            _aper_ind = 0 if one_aper else 1
            util.send_redis(exp, zps[zps['aper_ind'] == _aper_ind], sky[0],
                            sci_inst_name=sci_inst_name)

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

    parser.add_argument('--skip_flatfield', default=False, action='store_true',
                        help="skip flatfielding during pixel-level detrending")

    parser.add_argument('--sci_inst_name', default='desi', type=str,
                        help='name of science instrument')

    parser.add_argument('--sci_fov_checkplot', default=False,
                        action='store_true',
                        help="restrict checkplot to science instrument FOV")

    parser.add_argument('--check_tcs_motion', default=False,
                        action='store_true',
                        help="abort reductions based on telescope motion flag")

    parser.add_argument('--max_zp_radius', default=None, type=int,
                        help="maximum radius in pixels for zeropoint stars")

    parser.add_argument('--detect_streaks', default=False,
                        action='store_true',
                        help="run satellite streak detection/cataloging")

    parser.add_argument('--plot_detrended', default=False,
                        action='store_true',
                        help="make and save plot of detrended image")

    parser.add_argument('--plot_streaks', default=False,
                        action='store_true',
                        help="overplot detected streaks on detrended image")

    parser.add_argument('--plot_quiver', default=False, action='store_true',
                       help="make and save quiver plot of centroid shifts")

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
            nmp=args.multiproc, max_n_stars=args.max_n_stars,
            pm_corr=args.pm_corr, skip_flatfield=args.skip_flatfield,
            sci_inst_name=args.sci_inst_name,
            sci_fov_checkplot=args.sci_fov_checkplot,
            check_tcs_motion=args.check_tcs_motion,
            max_zp_radius=args.max_zp_radius,
            detect_streaks=args.detect_streaks,
            plot_detrended=args.plot_detrended, plot_streaks=args.plot_streaks,
            plot_quiver=args.plot_quiver)
