"""
pointing_camera.zp
==================

Calculations and plots related to determining pointing camera image zeropoints.
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pointing_camera.common as common
import copy
from astropy.table import Table, vstack
from multiprocessing import Pool
import pointing_camera.util as util

def calc_zp(_cat, aper_ind, time_seconds, fname_im, quadrant=0,
            one_aper=False, checkplot=True, science_fov_only=False,
            sci_fov_checkplot=False):
    """
    Calculate a zeropoint value for a particular aperture size.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Source catalog with refined centroids and aperture photometry.
        aper_ind : int
            Index of aperture to use within array of apertures photometered
            per object.
        time_seconds : float
            Exposure time in seconds. Zeropoints computed correspond to
            the magnitude of a source with flux of 1 ADU/second, so
            the exposure time is needed...
        fname_im : str
            File name of the raw pointing camera image. This gets packaged
            into the zeropoints table and also is used for constructing the
            plot title.
        quadrant : int, optional
            Should be one of 0, 1, 2, 3, 4 where 0 means the entire image
            (all four quadrants combined).
        one_aper : bool, optional
            If True, aperture photometry only exists for the 'standard' aperture,
            rather than the whole series of apertures.
        checkplot : bool, optional
            Create zeropoint scatterplot checkplot. Set False to skip checkplot.
        science_fov_only : bool, optional
            If True, restrict zeropoint computation to only those stars that
            fall within the science instrument field of view footprint.
        sci_fov_checkplot : bool, optional
             If True, restrict the sample of stars plotted in the checkplot to
             only the science field of view footprint.

    Returns
    -------
        result : astropy.table.table.Table
            Single-row table of zeropoints measurements and associated
            metadata/metrics.

    """

    print('Computing zeropoint for quadrant : ', quadrant, \
          ' , aper ', aper_ind, ' , science FOV only = ', science_fov_only)

    assert(time_seconds > 0)
    assert(quadrant in [0, 1, 2, 3, 4]) # note the 0 option here...

    if science_fov_only:
        # for now, I don't want to compute per-quadrant zeropoints
        # when restricting to the science instrument field of view
        # (not sure if there'd be enough stars to do that well per-quadrant?)
        # may revisit this later...
        assert(quadrant == 0)

    par = common.pc_params()

    n_aper = len(par['aper_phot_objrad'])

    aper_ind = int(aper_ind)
    assert(aper_ind in np.arange(n_aper))

    cat = copy.deepcopy(_cat)

    good = np.logical_not(cat['centroid_pixel_saturated']) & \
           (cat['centroid_raw_pixel_val'] < 15300) & \
           np.logical_not(cat['centroid_shift_flag']) & \
           np.logical_not(cat['wrong_source_centroid']) & \
           np.logical_not(cat['radius_too_large']) & \
           np.isfinite(cat['PHOT_BP_MEAN_MAG']) & \
           np.isfinite(cat['PHOT_RP_MEAN_MAG']) & \
           np.isfinite(cat['PHOT_G_MEAN_MAG'])

    if science_fov_only:
        good = np.logical_and(good, cat['in_science_fov'])

    if quadrant != 0:
        good = good & (cat['quadrant'] == quadrant)

    if np.sum(good) == 0:
        return None

    cat = cat[good]
    n = len(cat)

    m_inst = cat['m_inst'][:, aper_ind]

    diff = cat['G_PRIME'] - m_inst
    nf = np.sum(np.isfinite(diff))

    zp = np.nanmedian(diff)

    resid = diff - zp

    # now calculate robust sigma about the median zeropoint offset

    # at first glance argsort appears to put NaN's at end of sorted array
    sind = np.argsort(resid) # should look into what exactly happens with NaN's

    ind_l = max(int(round(0.16*nf)), 0)
    ind_u = min(int(round(0.84*nf)), n-1)

    resid_l = resid[sind[ind_l]]
    resid_u = resid[sind[ind_u]]

    sig_robust = (np.abs(resid_l) + np.abs(resid_u))/2.0

    result = Table()

    result['quadrant'] = [quadrant]
    result['aper_ind'] = [aper_ind]
    result['zp_adu_per_s'] = [zp]
    result['n_sources_for_zp'] = [n]
    result['time_seconds'] = [time_seconds]
    result['bp_rp_median'] = [np.nanmedian(cat['BP_RP'])]
    result['gaia_g_median'] = [np.nanmedian(cat['PHOT_G_MEAN_MAG'])]
    result['robust_sigma_mag'] = [sig_robust]
    result['fname_raw'] = [fname_im]
    result['science_fov_only'] = [int(science_fov_only)]

    # checkplot
    best_aper_ind = 1 if not one_aper else 0

    if checkplot and (quadrant == 0) and (aper_ind == best_aper_ind) and \
       (science_fov_only == sci_fov_checkplot):
        plt.cla()
        plt.figure(1)
        xtitle = 'G + 0.25' + r'$\times$' + '(BP-RP) [' + \
                 str(len(cat)) + ' stars]'
        ytitle = '-2.5' + r'$\times$' + 'log' + r'$_{10}$' + '(ADU/sec)'
        title = fname_im.split('/')[-1]
        title = title.replace('.fits', '')
        title += '; aper' + str(best_aper_ind) + '; '
        title += 'sci FOV' if sci_fov_checkplot else 'full FOV'

        plt.scatter(cat['G_PRIME'], m_inst, s=20, edgecolor='none',
                    facecolor='k')

        xmin = np.nanmin(cat['G_PRIME'])
        xmax = np.nanmax(cat['G_PRIME'])
        ymin = np.nanmin(m_inst)
        ymax = np.nanmax(m_inst)

        xsamp = np.array([xmin, xmax])
        ysamp = xsamp - zp

        plt.plot(xsamp, ysamp, linewidth=2, c='r')

        ax = plt.gca()

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        xtext = xlim[0] + (xlim[1] - xlim[0])*0.125
        ytext = ylim[0] + (ylim[1] - ylim[0])*0.875

        # this shouldn't crash for NaN zeropoint value...
        plt.text(xtext, ytext, 'ZP = ' + '{:.2f}'.format(zp), color='r')

        plt.title(title)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

    return result

def calc_many_zps(cat, exp, one_aper=False, checkplot=True, nmp=None,
                  sci_fov_checkplot=False, max_zp_radius=None):
    """
    Calculate multiple photometric zeropoint measurements.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Source catalog with refined centroids and aperture photometry.
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.
        one_aper : bool, optional
            If True, aperture photometry only exists for the 'standard' aperture,
            rather than the whole series of apertures.
        checkplot : bool, optional
            Create zeropoint scatterplot checkplot. Set False to skip checkplot.
        nmp : int, optional
            Number of processes to use for multiprocessing. Should be an integer
            greater than 1 but less than or equal to the number of CPU's. Default
            of None means no multiprocessing.
        sci_fov_checkplot : bool, optional
             If True, restrict the sample of stars plotted in the checkplot to
             only the science field of view footprint.
        max_zp_radius : float, optional
            Maximum radius in pixels from center of image for inclusion of
            stars in zeropoint determination. The idea is that beyond some
            radius from the image center, the PSF may degrade significantly
            such that we don't want to take stars far from the image center
            into account when computing the zeropoints. Default of None
            means that no such maximum image center distance cut is made.

    Returns
    -------
        results : astropy.table.table.Table
            Table of zeropoints measurements and associated metadata/metrics.

    """

    print('Attempting to calculate zeropoints')

    par = common.pc_params()

    aper_radii = par['aper_phot_objrad'] if not one_aper else [par['aper_phot_objrad_best']]

    args = []
    for q in [0, 1, 2, 3, 4]:
        for aper_ind in range(len(aper_radii)):
            args.append((cat, aper_ind, exp.time_seconds, exp.fname_im, q, one_aper, checkplot, False, sci_fov_checkplot))
            if q == 0:
                args.append((cat, aper_ind, exp.time_seconds, exp.fname_im, q, one_aper, checkplot, True, sci_fov_checkplot))

    if nmp is not None:
        p = Pool(min(nmp, len(args)))
        results = p.starmap(calc_zp, args)
        p.close()
        p.join()
    else:
        results = [calc_zp(*_arg) for _arg in args]

    results = vstack(results)
    results['mjd_obs'] = exp.header['MJD-OBS']
    results['obs_night'] = exp.obs_night

    if 'DEV_FNUM' in exp.header:
        results['DEV_FNUM'] = exp.header['DEV_FNUM']

    util.add_field_center_cols(results, exp.header)

    if exp.has_dome is not None:
        results['has_dome'] = int(exp.has_dome)

    sind = np.argsort(1000*results['quadrant'] + results['aper_ind'] + \
                      0.5*results['science_fov_only'])
    results = results[sind]

    # note that this could be confusing for rows that *further* restrict
    # to stars within the science FOV when computing the zeropoint
    results['max_zp_radius'] = 99999 if max_zp_radius is None else max_zp_radius

    return results
