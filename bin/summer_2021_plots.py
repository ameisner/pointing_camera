from nightly_plot import zp_subplot, sky_subplot, _read_concat_tables
import astropy.io.fits as fits
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import gfa_utils
import sky_mon
from astropy.time import Time
import matplotlib.dates as md

basedir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/summer_2021/v0006'

def _plot_airmass(tab, mjdrange, xticklabels=True,
                  title_extra='', do_xlabel=True, markersize=2):

    datetimes = []
    for row in tab:
        tm = Time(row['mjd_obs'], format='mjd')
        datetimes.append(tm.to_datetime())

    print(np.min(tab['airmass']), np.max(tab['airmass']))

    plt.scatter(datetimes, tab['airmass'],
                edgecolor='none', s=markersize, c='k')

    date_min = Time(mjdrange[0] - 0.01, format='mjd').to_datetime()
    date_max = Time(mjdrange[1] + 0.01, format='mjd').to_datetime()
    plt.xlim((date_min, date_max))

    plt.ylim((0.98, 3)) # could think more about this choice...

    title = 'airmass' + title_extra

    plt.title(title)

    plt.ylabel('airmass')
    
    ax = plt.gca()
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    
    if do_xlabel:
        plt.xlabel(r'time (UT)', fontsize=12)

    if not xticklabels:
        ax.axes.xaxis.set_ticklabels([])

def _multipanel_1night(night, basedir=basedir, markersize=10, save=False,
                       skip_q0=False, use_gfa=False, use_skymon=False,
                       plot_airmass=False, nmp=None):

    plt.cla()

    if use_gfa:
        gfa = gfa_utils.read_gfa()
        gfa = gfa[gfa['NIGHT'] == int(night)]
        assert(len(gfa) > 0)
    
    n_panels = 2 + int(use_gfa)*2 + int(use_skymon) + int(plot_airmass)

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

    ct = 0

    if use_gfa:
        ct += 1
        plt.subplot(n_panels, 1, ct)
        gfa_utils.plot_gfa(gfa, markersize=markersize,
                           mjdrange=mjdrange, title_extra=title_extra,
                           colname='FWHM_ASEC', do_xlabel=(ct == n_panels),
                           xticklabels=(ct == n_panels))

        ct += 1
        plt.subplot(n_panels, 1, ct)
        gfa_utils.plot_gfa(gfa, markersize=markersize,
                           mjdrange=mjdrange, title_extra=title_extra,
                           colname='TRANSPARENCY', xticklabels=(ct == n_panels),
                           do_xlabel=(ct == n_panels))

    ct += 1
    plt.subplot(n_panels, 1, ct)
    zp_subplot(zps, mjdrange=mjdrange, markersize=markersize,
               title_extra=title_extra, skip_q0=skip_q0,
               do_xlabel=(ct == n_panels), xticklabels=(ct == n_panels))

    ct += 1
    plt.subplot(n_panels, 1, ct)
    sky_subplot(skies, mjdrange=mjdrange, markersize=markersize,
                title_extra=title_extra, skip_q0=skip_q0,
                do_xlabel=(ct == n_panels), xticklabels=(ct == n_panels))

    if use_skymon:
        ct += 1
        plt.subplot(n_panels, 1, ct)
        # DO NOT FACTOR RANGE OF SKY MON MJD VALUES INTO mjdrange !
        sky_mon.plot_sky_mon(night, mjdrange, xticklabels=(ct == n_panels),
                             title_extra=title_extra,
                             do_xlabel=(ct == n_panels))

    if plot_airmass:
        ct += 1
        plt.subplot(n_panels, 1, ct)
        _plot_airmass(skies, mjdrange, xticklabels=(ct == n_panels),
                      title_extra=title_extra, do_xlabel=(ct == n_panels))

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
        if plot_airmass:
            outname += '-airmass'
        outname += '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')
    else:
        plt.show()

def summer_2021_nightly_plots(markersize=2, skip_q0=False, use_gfa=False,
                              use_skymon=False, nmp=None, plot_airmass=False):

    nights = glob.glob(os.path.join(basedir, '????????'))

    nights = [os.path.split(night)[-1] for night in nights]

    nights.sort()

    for night in nights:
        _multipanel_1night(night, basedir=basedir, markersize=markersize,
                           save=True, skip_q0=skip_q0, use_gfa=use_gfa,
                           use_skymon=use_skymon, nmp=nmp,
                           plot_airmass=plot_airmass)


