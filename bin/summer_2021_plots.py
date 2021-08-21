from nightly_plot import zp_subplot, sky_subplot, _read_concat_tables
import astropy.io.fits as fits
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import gfa_utils
import sky_mon

basedir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/proc_1night'

def _multipanel_1night(night, basedir=basedir, markersize=10, save=False,
                       skip_q0=False, use_gfa=False, use_skymon=False,
                       nmp=None):

    plt.cla()

    if use_gfa:
        gfa = gfa_utils.read_gfa()
        gfa = gfa[gfa['NIGHT'] == int(night)]
        assert(len(gfa) > 0)
    
    n_panels = 2 + int(use_gfa)*2 + int(use_skymon)

    width_inches = 6.4
    height_inches = 2.8*n_panels

    dir = os.path.join(basedir, night)

    assert(os.path.exists(dir))

    # note that the sky and zp tables are now both contained within the same
    # multi-extension FITS file

    flist = glob.glob(os.path.join(dir, '*-summary.fits'))

    flist.sort()

    skies = _read_concat_tables(flist, ext=3, nmp=nmp)
    zps = _read_concat_tables(flist, ext=2, nmp=nmp)

    title_extra = ' ; ' + night
    #_twopanel(skies, zps, clobber=True, save=False, markersize=markersize,
    #          title_extra=title_extra, skip_q0=skip_q0, fig=fig)

    mjdrange = [min(np.min(zps['mjd_obs']),
                    np.min(skies['mjd_obs'])),
                max(np.max(zps['mjd_obs']),
                    np.max(skies['mjd_obs']))]

    if use_gfa:
        mjdrange[0] = min(mjdrange[0], min(gfa['MJD']))
        mjdrange[1] = max(mjdrange[1], max(gfa['MJD']))

    plt.subplot(n_panels, 1, 1)
    zp_subplot(zps, mjdrange=mjdrange, markersize=markersize,
               title_extra=title_extra, skip_q0=skip_q0)

    plt.subplot(n_panels, 1, 2)
    sky_subplot(skies, mjdrange=mjdrange, markersize=markersize,
                title_extra=title_extra, skip_q0=skip_q0,
                do_xlabel=(not use_gfa), xticklabels=(not use_gfa))

    ct = 2
    if use_gfa:
        ct += 1
        plt.subplot(n_panels, 1, ct)
        gfa_utils.plot_gfa(gfa, markersize=markersize,
                           mjdrange=mjdrange, title_extra=title_extra,
                           colname='TRANSPARENCY', xticklabels=False,
                           do_xlabel=False)
        ct += 1
        plt.subplot(n_panels, 1, ct)
        gfa_utils.plot_gfa(gfa, markersize=markersize,
                           mjdrange=mjdrange, title_extra=title_extra,
                           colname='FWHM_ASEC', do_xlabel=(not use_skymon),
                           xticklabels=(not use_skymon))

    if use_skymon:
        ct += 1
        plt.subplot(n_panels, 1, ct)
        # DO NOT FACTOR RANGE OF SKY MON MJD VALUES INTO mjdrange !
        sky_mon.plot_sky_mon(night, mjdrange, xticklabels=True,
                             title_extra=title_extra, do_xlabel=True)
        
    fig = plt.gcf()

    fig.set_size_inches(width_inches, height_inches)

    plt.subplots_adjust(hspace=0.18)
    
    if save:
        outname = 'pointing_camera-' + night
        if skip_q0:
            outname += '-no_q0'
        if use_gfa:
            outname += '-gfa'
        if use_skymon:
            outname += '-skymon'
        outname += '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')
    else:
        plt.show()

def summer_2021_nightly_plots(markersize=2, skip_q0=False, use_gfa=False,
                              use_skymon=False, nmp=None):

    nights = glob.glob(os.path.join(basedir, '????????'))

    nights = [os.path.split(night)[-1] for night in nights]

    nights.sort()

    for night in nights:
        _multipanel_1night(night, basedir=basedir, markersize=markersize,
                           save=True, skip_q0=skip_q0, use_gfa=use_gfa,
                           use_skymon=use_skymon, nmp=nmp)


