"""
pointing_camera.util
====================

A collection of pointing camera related utility functions.
"""

import astropy.io.fits as fits
from astropy import wcs
import os
import glob
import pointing_camera.common as common
import pointing_camera.io as io
import pointing_camera.analysis.djs_maskinterp as djs_maskinterp
from astropy.table import Table, hstack, vstack
import numpy as np
from pointing_camera.gaia import read_gaia_cat
import time
from pointing_camera.analysis.djs_photcen import djs_photcen
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
import photutils
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from multiprocessing import Pool
from functools import lru_cache
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

def get_wcs_filename(fname_im, verify=True):
    """
    Construct the astrometry.net WCS file name corresponding to a raw image.

    Parameters
    ----------
        fname_im : str
            Full path of raw image file name for which to construct the
            corresponding WCS file name.
        verify : bool (optional)
            If True, check whether the WCS file name constructed exists.

    Returns
    -------
        fname_wcs : str
            WCS file name constructed including full path.

    Notes
    -----
        Currently this is pretty trivial but perhaps it could
        become more complicated in the future.

    """

    fname_wcs = fname_im.replace('.fits', '.wcs')

    if verify:
        assert(os.path.exists(fname_wcs))

    return fname_wcs

def load_wcs(fname_wcs):
    """
    Load astrometry.net WCS file.

    Parameters
    ----------
        fname_wcs : str
            Full file name of the WCS file to be read in.

    Returns
    -------
        w : astropy.wcs.wcs.WCS
            Astropy WCS object.
        header : astropy.io.fits.header.Header
            WCS file as a header object.

    """

    assert(os.path.exists(fname_wcs))

    hdul = fits.open(fname_wcs)

    header = hdul[0].header

    # somehow this both averts a warning and saves ~0.25 seconds
    # when creating the wcs.WCS object
    header['NAXIS'] = 2

    w = wcs.WCS(header)

    return w, header

def load_exposure_image(fname):
    """
    Read in raw pointing camera image and its header.

    Parameters
    ----------
        fname : str
            Full file path to the raw image.

    Returns
    -------
        im : numpy.ndarray
            Raw image pixel data.
        h  : astropy.io.fits.header.Header
            Raw image FITS header object.

    Notes
    -----
        This function gets the pixel data and header for the actual pointing
        camera image -- the WCS solution is loaded separately via the
        load_wcs function.

    """

    assert(os.path.exists(fname))

    print('Attempting to read exposure : ' + fname)

    im, h = fits.getdata(fname, header=True)

    print('Successfully read raw pointing camera image file')

    return im, h

def check_image_dimensions(image):
    """
    Verify that raw image dimensions are correct.

    Parameters
    ----------
        image : numpy.ndarray
            Raw pointing camera image.

    Notes
    -----
        Of particular relevance is not getting fooled by downbinned data...

    """

    print('Checking raw pointing camera image dimensions...')

    par = common.pc_params()

    sh = image.shape

    assert(sh[0] == par['ny'])
    assert(sh[1] == par['nx'])

    print('Raw pointing camera image has correct dimensions')

def get_exptime(h_im, milliseconds=False):
    """
    Retrieve exposure time from raw pointing camera image header.

    Parameters
    ----------
        h_im : astropy.io.fits.header.Header
            Header of a raw pointing camera image.
        milliseconds : bool
            If True, return exposure time in milliseconds (default
            units for return value are seconds).

    Returns
    -------
        float
            Exposure time in either seconds (when milliseconds=False) or
            ms (when milliseconds=True).

    """

    if 'EXPOSURE' in h_im:
        exptime_ms = h_im['EXPOSURE']
    elif 'EXPTIME' in h_im:
        exptime_ms = h_im['EXPTIME']
    else:
        print('Could not find an exposure time in the raw image header!')
        assert(False)

    exptime_ms = float(exptime_ms)

    # think that for El Nino this ais ctually impossible since neither
    # readout mode offers 0 ms exposure times...
    assert(exptime_ms != 0)

    print('Exposure time is ' + str(round(exptime_ms)) + ' ms')
    if milliseconds:
        return exptime_ms
    else:
        return exptime_ms/1000.0

def _check_bitpix(h_im):
    """
    Check that the raw FITS header indicates the correct data type.

    Parameters
    ----------
        h_im : astropy.io.fits.header.Header
            Header of a raw pointing camera image

    Notes
    -----
        Concern is to flag/avoid cases where the camera may be reading
        out in some non-standard mode that won't result in sensible
        pipeline outputs.

    """

    par = common.pc_params()

    print('Checking raw image BITPIX value')

    assert(h_im['BITPIX'] == par['bitpix'])

def quad_pix_limits(quad):
    """
    Get bounding pixel indices for each raw image quadrant.

    Parameters
    ----------
        quad : int
            Integer quadrant number, using the standard convention
            (when the image is rendered as in DS9 and IDL/ATV. Should
            be one of [1, 2, 3, 4]. Not intended to work for multi-element
            inputs.

    Returns
    -------
        result : dict
            Dictionary containing minimum and maximum pixel indices
            along each axis.

    """

    assert(quad in [1, 2, 3, 4])

    par = common.pc_params()

    half_x = par['nx'] // 2
    half_y = par['ny'] // 2

    if (quad == 1) or (quad == 4):
        xmin = half_x
        xmax = par['nx']
    else:
        xmin = 0
        xmax = half_x

    if (quad == 1) or (quad == 2):
        ymin = half_y
        ymax = par['ny']
    else:
        ymin = 0
        ymax = half_y

    result = {'xmin' : xmin, 'xmax' : xmax, 'ymin' : ymin, 'ymax' : ymax}

    return result

def get_quadrant(im, q):
    """
    Extract a quadrant from a full-sized pointing camera image.

    Parameters
    ----------
        im : numpy.ndarray
            2D numpy array. Should have the same dimensions as a raw
            pointing camera image.
        q : int
            Integer quadrant number, using the standard convention
            (when the image is rendered as in DS9 and IDL/ATV. Should
            be one of [1, 2, 3, 4]. Not intended to work for multi-element
            inputs.

    Returns
    -------
        quad : numpy.ndarray
            Extracted image quadrant.

    """

    assert(q in [1, 2, 3, 4])

    check_image_dimensions(im)

    par = common.pc_params()

    p = quad_pix_limits(q)

    quad = im[p['ymin']:p['ymax'], p['xmin']:p['xmax']]

    return quad

def quadrant_from_xy(x, y):
    """
    Convert from image pixel coordinates to quadrant number.

    Parameters
    ----------
        x : numpy.ndarray
            x pixel coordinate values of interest (0-indexed).
            Should have the same dimensions as y.

        y : numpy.ndarray
            y pixel coordinate values of interest (0-indexed).
            Should have the same dimensions as x.

    Returns
    -------
        quadrant : numpy.ndarray
            Integer image quadrant values corresponding to the input
            (x, y) pairs.

    """

    par = common.pc_params()

    half_x = par['nx']/2 - 0.5
    half_y = par['ny']/2 - 0.5

    quadrant = np.zeros(len(x), dtype=int)

    quadrant[(x <= half_x) & (y <= half_y)] = 3

    quadrant[(x <= half_x) & (y > half_y)] = 2

    quadrant[(x > half_x) & (y <= half_y)] = 4

    quadrant[(x > half_x) & (y > half_y)] = 1

    return quadrant

def min_edge_dist_pix(x, y):
    """
    Compute minimum distance to any image edge.

    Parameters
    ----------
        x : numpy.ndarray
            x pixel coordinate values.
            Should have the same dimensions as y.
        y : numpy.ndarray
            y pixel coordinate values of interest (0-indexed).
            Should have the same dimensions as x.

    Returns
    -------
        min_edge_dist : numpy.ndarray
            Minimum distance to any image edge. 'Edge'
            is taken to be the outer edge of a bounding pixel,
            for instance on the left edge x would be -0.5.

    """

    min_edge_dist = 20000

    par = common.pc_params()

    min_edge_dist = np.minimum(x + 0.5, y + 0.5)
    min_edge_dist = np.minimum(min_edge_dist, par['nx'] - 0.5 - x)
    min_edge_dist = np.minimum(min_edge_dist, par['ny'] - 0.5 - y)

    return min_edge_dist

def subtract_dark_current(im, time_seconds):
    """
    Subtract per-quadrant offsets representing accumulated dark current.

    Parameters
    ----------
        im : numpy.ndarray
            Image from which to subtract the dark current. Should
            be float data type rather than integer.

        time_seconds : float
            Exposure time in seconds. Gets multiplied by the dark current
            rate to determine total dark current to subtract.

    Returns
    -------
        result : numpy.ndarray
            Dark current subtracted version of input image.

    """

    print('Subtracting dark current')

    assert(time_seconds < 30)

    par = common.pc_params()

    result = im.astype('float32')

    for q in [1, 2, 3, 4]:
        dark_adu = par['dark_adu_per_s_quad'][q-1]*time_seconds
        print('quadrant : ', q, ', dark counts/pix : ', '{:.3f}'.format(dark_adu))

        p = quad_pix_limits(q)

        result[p['ymin']:p['ymax'], p['xmin']:p['xmax']] -= dark_adu

    return result

def subtract_master_dark(im, time_seconds):
    """
    Subtract image-level master dark as part of pixel-level detrending.

    Parameters
    ----------
        im : numpy.ndarray
            Image from which to subtract the dark current. Should
            be float data type rather than integer.

        time_seconds : float
            Exposure time in seconds. Gets multiplied by the dark current
            rate to determine total dark current to subtract.

    Returns
    -------
        result : numpy.ndarray
            Dark current subtracted version of input image.

    """

    dark = io.load_master_dark()

    print('Subtracting image-level master dark')

    im -= dark*time_seconds

    return im

def subtract_quad_offs(im):
    """
    Subtract a scalar offset per image quadrant.

    Parameters
    ----------
        im : numpy.ndarray
            2D image array that should have the standard dimensions
            of a raw pointing camera image.

    Returns
    -------
        im : numpy.ndarray
            Version of the input image with per-quadrant
            scalar offsets subtracted. Output is floating point
            even if input was integer data type.

    """

    par = common.pc_params()

    # just in case this hasn't already been taken care of...
    if im.dtype.name != 'float32':
        im = im.astype('float32')

    check_image_dimensions(im)

    for q in [1, 2, 3, 4]:
        p = quad_pix_limits(q)

        im[p['ymin']:p['ymax'], p['xmin']:p['xmax']] -= par['bias_med'][q-1]

    return im

def subtract_master_bias(im):
    """
    Subtract master bias as part of pixel level detrending.

    Parameters
    ----------
        im : numpy.ndarray
            2D array holding pointing camera image.

    Returns
    -------
        im : numpy.ndarray
            Pointing camera image after bias subtraction.

    Notes
    -----
        This function subtracts off an image-level master bias template,
        rather than merely subtracting off scalar per-quadrant offsets.

    """

    bias = io.load_master_bias()

    print('Subtracting image-level master bias')

    im -= bias

    return im

def apply_flatfield(im):
    """
    Apply flatfield correction.

    Parameters
    ----------
        im : numpy.ndarray
            2D array holding pointing camera image. Should already
            have bias and dark current removed.

    Returns
    -------
        im : numpy.ndarray
            Pointing camera image after flat field correction.

    Notes
    -----
        Flatfield should have a median value of 1. Flatfield built
        with ~1,600 serendipitous 'sky flat' exposures from
        summer 2021 taken with bright sky and no clouds.

    """

    flat = io.load_master_flat()

    print('Applying flat field correction')

    im /= flat

    return im

def badpix_interp(im):
    """
    Interpolate over static bad pixels.

    Parameters
    ----------
        im : numpy.ndarray
            2D array holding pointing camera image. Should have the
            usual dimensions of a raw pointing camera image. Should
            be floating point data type.

    Returns
    -------
        result : numpy.ndarray
            2D array holding modified version of input image, where
            pixels in the static bad pixel mask have been interpolated over.

    """

    mask = io.load_static_badpix()

    result = djs_maskinterp.average_bilinear(im, (mask != 0))

    return result

def detrend_pc(exp, skip_flatfield=False, ml_dome_flag=False):
    """
    Driver for detrending a raw pointing camera image.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.
        skip_flatfield : bool, optional
            Set True to skip the flat field application step.
        ml_dome_flag : bool, optional
            If True, use machine learning classifier to flag dome vignetting
            (otherwise a heuristic method is employed).

    """

    print('Attempting to detrend the raw pointing camera image')

    assert(not exp.is_detrended)

    # try to avoid accidentally using units of ms instead of seconds
    assert(exp.time_seconds < 30)

    # final thing is to set detrended == true

    im = exp.raw_image.astype('float32')

    im = subtract_master_bias(im)

    im = subtract_master_dark(im, exp.time_seconds)

    if not skip_flatfield:
        im = apply_flatfield(im)

    im = badpix_interp(im)

    exp.is_detrended = True

    exp.detrended = im

    exp.update_dome_flag(use_ml=ml_dome_flag)

    print('Finished detrending the raw pointing camera image')

def sky_metrics(im, mask=None):
    """
    Compute sky background pixel value statistics.

    Parameters
    ----------
        im : numpy.ndarray
            2D numpy array, intended to be either a full-frame image or one
            quadrant of a full-frame image.
        mask : numpy.ndarray
            A 2D boolean image array with the same dimensions as im.

    Returns
    -------
        result : dict
            Various sky background value statistics such as robust standard
            deviation, sigma clipped mean, and median.

    """

    sz = im.shape
    assert(len(sz) == 2)

    if mask is None:
        im_sorted = np.ravel(im)
    else:
        im_sorted = im[mask]

    sind = np.argsort(im_sorted)

    im_sorted = im_sorted[sind]

    n = len(im_sorted)

    ind_med = int(round(0.5*n))
    ind_u = int(round(0.84*n))
    ind_l = int(round(0.16*n))

    ind_med = min(max(ind_med, 0), n-1)
    ind_u = min(max(ind_u, 0), n-1)
    ind_l = min(max(ind_l, 0), n-1)

    med = im_sorted[ind_med]

    val_u = im_sorted[ind_u]
    val_l = im_sorted[ind_l]

    sig = (val_u - val_l)/2.0

    nsig_thresh = 5.0

    val_l = med - nsig_thresh*sig
    val_u = med + nsig_thresh*sig

    good = (im_sorted >= val_l) & (im_sorted <= val_u)
    ngood = int(np.sum(good))

    assert(ngood > 0)

    result = {'median': med,
              'sky_clipped_mean': np.mean(im_sorted[good]),
              'sky_sig': sig,
              'n_pixels_used': ngood}

    return result

def sky_summary_table(exp):
    """
    Driver for sky analysis and metrics.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.

    Returns
    -------
        t : astropy.table.table.Table
            Table of sky metrics/metadata. Should have one row with many
            columns.

    """

    print('Computing sky background statistics')

    par = common.pc_params()

    t = Table()
    t['time_seconds'] = [exp.time_seconds]
    t['fname_raw'] = [exp.fname_im]

    metrics = sky_metrics(exp.detrended)
    t['median_adu'] = [metrics['median']]
    t['median_adu_per_s'] = [metrics['median']/exp.time_seconds]
    t['mean_adu'] = [metrics['sky_clipped_mean']]
    t['mean_adu_per_s'] = [metrics['sky_clipped_mean']/exp.time_seconds]

    sci_mask = circular_mask(par['science_radius_pix'])
    metrics_sci = sky_metrics(exp.detrended, mask=sci_mask)
    t['median_adu_sci'] = [metrics_sci['median']]
    t['median_adu_sci_per_s'] = [metrics_sci['median']/exp.time_seconds]
    t['mean_adu_sci'] = [metrics_sci['sky_clipped_mean']]
    t['mean_adu_sci_per_s'] = [metrics_sci['sky_clipped_mean']/exp.time_seconds]

    qmeds = []
    for q in [1, 2, 3, 4]:
        qmetrics = sky_metrics(get_quadrant(exp.detrended, q))
        qmed = qmetrics['median']
        qmean = qmetrics['sky_clipped_mean']
        qmeds.append(qmed)

        t['median_adu_quad' + str(q)] = [qmed]

        t['median_adu_quad' + str(q) + '_per_s'] = [qmed/exp.time_seconds]

        t['mean_adu_quad' + str(q)] = [qmean]

        t['mean_adu_quad' + str(q) + '_per_s'] = [qmean/exp.time_seconds]

    diff = np.array(qmeds) - np.median(qmeds)

    rms_adu = np.sqrt(np.mean(np.power(diff, 2)))

    t['inter_quad_rms_adu'] = [rms_adu]
    t['inter_quad_rms_adu_per_s'] = [rms_adu/exp.time_seconds]
    t['mjd_obs'] = exp.header['MJD-OBS']
    t['obs_night'] = exp.obs_night

    add_field_center_cols(t, exp.header)

    if 'DEV_FNUM' in exp.header:
        t['DEV_FNUM'] = exp.header['DEV_FNUM']

    if exp.has_dome is not None:
        t['has_dome'] = int(exp.has_dome)

    return t

def _validate_ctype(wcs):
    """
    Check that ctype values in WCS are as expected.

    Parameters
    ----------
        wcs : astropy.wcs.wcs.WCS
            An astropy WCS object.

    """

    valid_ctypes = ['RA---TAN', 'DEC--TAN', 'RA---TAN-SIP', 'DEC--TAN-SIP']
    for ctype in wcs.wcs.ctype:
        assert(ctype in valid_ctypes)

def max_gaia_mag(time_seconds):
    """
    Maximum Gaia magnitude to retain as a function of exposure time.

    Parameters
    ----------
        time_seconds : float
            Exposure time in units of seconds. Should work for both array
            and scalar inputs.

    Returns
    -------
        val : float
            The maximum Gaia magnitude to retain given the exposure time. If
            time_seconds is an array, then the return value should also be an
            array with the same dimensions.

    """

    fac = 1.21

    val = 15.5 + 2.5*np.log10(time_seconds/26.0)*fac

    val = np.minimum(np.maximum(val, 10), 15.5)

    return val

def xy_subsamp_grid():
    """
    Subsampled grid of pixel coordinate locations spanning the full image extent.

    Returns
    -------
        xgrid : numpy.ndarray
            2D array of x pixel coordinate values. Same shape as ygrid.
        ygrid : numpy.ndarray
            2D array of y pixel coordinate values. Same shape as xgrid.

    Notes
    -----
        The values hardcoded here are going to give me problems when trying
        to run on data with different dimensions e.g., la Nina. Should
        generalize this function in the future.

    """

    x = np.arange(104)*32
    y = np.arange(104)*24

    xgrid, ygrid = np.meshgrid(x, y)

    return xgrid, ygrid

def pc_gaia_cat(exp, mag_thresh=None, edge_pad_pix=0, nmp=None,
                max_n_stars=3000, pm_corr=False):
    """
    Assemble Gaia catalog relevant to a particular pointing camera exposure.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.

        mag_thresh : float, optional
            Maximum Gaia G mag for which to retain Gaia stars. Default value
            of None means that no Gaia minimum brightness threshold will be applied.

        edge_pad_pix : float, optional
            Minimum distance from image boundaries to enforce when deciding whether
            to retain each Gaia star. Default of 0 simply means that stars within the
            boundaries are retained and those outside of the image boundaries are not
            retained.

        nmp : int, optional
            Number of processes to use for multiprocessing. Should be an integer
            greater than 1 but less than the number of CPU's. Default of None means
            that read-in of the Gaia chunk files will happen in serial.

        max_n_stars : int, optional
            Maximum number of Gaia stars to retain. The idea is to avoid running
            into a situation where the pipeline takes forever trying to photometer
            a huge number of stars in dense fields. If the total number of stars
            exceeds max_n_stars, then the brightest max_n_stars are retained.

        pm_corr : bool, optional
            If set True, apply proper motion corrections to Gaia star positions,
            so that their positions match the epoch of the pointing camera
            exposure. By default this is False, on the assumption that only a
            tiny fraction of stars will have high enough proper motions for this
            to matter (El Nino, for instance, has ~8.6" pixels), especially given
            that this pipeline does empirical centroid refinement for all stars.

    Returns
    -------
        cat : astropy.table.table.Table
            Gaia catalog for stars relevant to the pointing camera exposure.

    Notes
    -----
        Would be good to add a check that the number of requested multiprocessing
        processes does not exceed the number of CPU's.

        Upgrade Gaia catalogs to eDR3/DR[3+] at some point? This has to do
        with what files are on disk though, and wouldn't affect the code in
        this function at all.

    """

    print('Reading Gaia DR2 catalogs...')

    xgrid, ygrid = xy_subsamp_grid()

    wcs = exp.wcs

    # last arg is 0 rather than 1 because that's what agrees with IDL
    ra, dec = wcs.all_pix2world(xgrid, ygrid, 0)

    cat = read_gaia_cat(ra, dec, nmp=nmp)

    if mag_thresh is not None:
        keep = (cat['PHOT_G_MEAN_MAG'] <= mag_thresh)
        assert(np.sum(keep) > 0)
        cat = cat[keep]

    if pm_corr:
        gaia_pm_corr(cat, exp.header['MJD-OBS'])

    x_gaia_guess, y_gaia_guess = wcs.all_world2pix(cat['RA'], cat['DEC'], 0)

    par = common.pc_params()

    # should switch to calling min_edge_dist utility here?
    keep  = (x_gaia_guess > edge_pad_pix) & (y_gaia_guess > edge_pad_pix) & (x_gaia_guess < (par['nx'] - 1 - edge_pad_pix)) & (y_gaia_guess < (par['ny'] - 1 - edge_pad_pix))

    assert(np.sum(keep) > 0)

    cat = cat[keep]

    x_gaia_guess = x_gaia_guess[keep]
    y_gaia_guess = y_gaia_guess[keep]

    cat = Table(cat)
    cat['x_gaia_guess'] = x_gaia_guess
    cat['y_gaia_guess'] = y_gaia_guess

    if len(cat) > max_n_stars:
        # retain brightest max_n_stars
        # it'd be better to do this cut based on
        # 'G_PRIME', the color-corrected pointing
        # camera Gaia-based mag
        # in the future can evaluate trying to spread the
        # selected max_n_stars evenly across quadrants
        # (could imagine a pathological case with e.g., a
        # globular cluster in the FOV)
        # might also be worth considering throwing out stars that
        # are excessively bright, though those are probably always
        # a very small fraction of the retained max_n_stars (?)
        print('Restricting to the brightest ' + str(max_n_stars) + \
              ' of ' + str(len(cat)) + ' stars')
        sind = np.argsort(cat['PHOT_G_MEAN_MAG'])
        cat = cat[sind[0:max_n_stars]]

    return cat

def pc_recentroid(im, cat):
    """
    Perform star recentroiding.

    Parameters
    ----------
        im : numpy.ndarray
            Pointing camera image. Should be the detrended image.
        cat : astropy.table.table.Table
            Table of stars to recentroid. Columns 'x_gaia_guess',
            'y_gaia_guess' provide the starting locations adopted
            for each star.

    Returns
    -------
        cat : astropy.table.table.Table
            A modified version of the input catalog with more columns.

    """

    par = common.pc_params()

    im = im.astype(float)

    x_center = par['nx']/2 + 0.5
    y_center = par['ny']/2 + 0.5

    dist_pix = np.sqrt(np.power(cat['x_gaia_guess'] - x_center, 2) + \
                       np.power(cat['y_gaia_guess'] - y_center, 2))

    dist_max = np.sqrt(x_center**2 + y_center**2)

    slope = 2.0/dist_max

    cmaxshift = np.minimum(np.maximum(dist_pix*slope + 2.5, 2.5), 4.5)

    assert(np.sum(cmaxshift < 2.5) == 0)
    assert(np.sum(cmaxshift > 4.5) == 0)

    xcen = np.zeros(len(cat), dtype=float)
    ycen = np.zeros(len(cat), dtype=float)

    qmaxshift = np.zeros(len(cat), dtype=int)

    for i in range(len(cat)):
        _xcen, _ycen, q = djs_photcen(cat['x_gaia_guess'][i],
                                      cat['y_gaia_guess'][i], im, cbox=8,
                                      cmaxiter=10, cmaxshift=cmaxshift[i])
        xcen[i] = _xcen
        ycen[i] = _ycen
        qmaxshift[i] = q

    result = Table()

    result['xcentroid'] = xcen
    result['ycentroid'] = ycen
    result['x_shift'] = xcen - cat['x_gaia_guess']
    result['y_shift'] = ycen - cat['y_gaia_guess']
    result['cmaxshift'] = cmaxshift
    result['qmaxshift'] = qmaxshift

    result['centroid_shift_flag'] = (np.abs(result['x_shift']) > cmaxshift) | (np.abs(result['y_shift']) > cmaxshift) | (qmaxshift != 0)

    assert(len(result) == len(cat))
    cat = hstack([cat, result])

    return cat

def flag_wrong_centroids(cat, full_cat):
    """
    Flag cases of a centroid wandering off to an entirely different star.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Table with columns including xcentroid, ycentroid, SOURCE_ID. Can be a
            subset of the rows for the entire pointing camera exposure's
            catalog, with the idea being that partial lists can be run
            in parallel to reduce run time.
        full_cat : astropy.table.table.Table
            Table with columns including x_gaia_guess, y_gaia_guess, SOURCE_ID.
            Needs to be the full star catalog for this pointing camera exposure.

    Returns
    -------
        cat : astropy.table.table.Table
            Input catalog cat but with an added column called wrong_source_centroid,
            which has a value of 1 (centroid has wandered too far) or 0.

    """

    wrong_source_centroid = np.zeros(len(cat), dtype=bool)

    for i in range(len(cat)):
        _dist = np.sqrt(np.power(full_cat['x_gaia_guess'] - cat['xcentroid'][i], 2) + np.power(full_cat['y_gaia_guess'] - cat['ycentroid'][i], 2))
        indmin = np.argmin(_dist)
        wrong_source_centroid[i] = (full_cat[indmin]['SOURCE_ID'] != cat[i]['SOURCE_ID'])

    cat['wrong_source_centroid'] = wrong_source_centroid.astype(int)

    return cat

def _get_area_from_ap(ap):
    """
    Get the area of a photutils aperture in pixels.

    Parameters
    ----------
        ap : photutils.aperture.circle.CircularAperture
            photutils aperture object. Should be a single object rather
            than a list of them.

    Returns
    -------
        area : float
            Area of aperture in pixels.

    Notes
    -----
        This is to try and work around the photutils API change
        between versions 0.6 and 0.7.

    """

    if (photutils.__version__.find('0.7') != -1) or (photutils.__version__.find('1.0') != -1) or (photutils.__version__.find('1.2') != -1):
        area = ap.area # 0.7
    else:
        area = ap.area() # 0.6

    return area

def pc_aper_phot(im, cat, one_aper=False, bg_sigclip=False):
    """
    Driver for aperture photometry.

    Parameters
    ----------
        im : numpy.ndarray
            2D array representing the detrended image.
        cat : astropy.table.table.Table
            Source catalog, with each row representing a Gaia star.
        one_aper : bool, optional
            If True, only compute aperture photometry in the 'standard' aperture,
            rather than a series of apertures. This is meant to be a way
            of reducing run time.
        bg_sigclip : bool, optional
            If True, use sigma clipping when computing background level median
            value in sky annulus. False by default for speed, as the
            sigma_clipped_stats call is slow. May be worth investigating
            other sigma clipped median utilities to see if those are
            somehow faster.

    Returns
    -------
        cat : astropy.table.table.Table
            Version of the input catalog with added columns pertaining to
            the aperture photometry results/outputs.

    """

    im = im.astype(float)

    par = common.pc_params()

    positions = list(zip(cat['xcentroid'], cat['ycentroid']))

    radii = par['aper_phot_objrad'] if not one_aper else [par['aper_phot_objrad_best']]
    ann_radii = par['annulus_radii'] # should have 2 elements - inner and outer

    apertures = [CircularAperture(positions, r=r) for r in radii]
    annulus_apertures = CircularAnnulus(positions, r_in=ann_radii[0],
                                        r_out=ann_radii[1])
    annulus_masks = annulus_apertures.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(im)
        annulus_data_1d = annulus_data[mask.data > 0]
        if bg_sigclip:
            # this sigma_clipped_stats call is actually the slow part !!
            _, median_sigclip, std_bg = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        else:
            bkg_median.append(np.median(annulus_data_1d))

    bkg_median = np.array(bkg_median)
    phot = aperture_photometry(im, apertures)

    for i, aperture in enumerate(apertures):
        aper_bkg_tot = bkg_median*_get_area_from_ap(aperture)
        cat['aper_sum_bkgsub_' + str(i)] = phot['aperture_sum_' + str(i)] - aper_bkg_tot

        cat['aper_bkg_' + str(i)] = aper_bkg_tot

    cat['sky_annulus_area_pix'] = _get_area_from_ap(annulus_apertures)
    cat['sky_annulus_median'] = bkg_median

    flux_adu = np.zeros((len(cat), len(radii)), dtype=float)

    for i in range(len(radii)):
        flux_adu[:, i] = cat['aper_sum_bkgsub_' + str(i)]

    cat['flux_adu'] = flux_adu

    return cat

def get_g_prime(G, BP_RP):
    """
    Color-corrected Gaia magnitudes to match the pointing camera bandpass.

    Parameters
    ----------
        G : numpy.ndarray
            Gaia G magnitude(s). Should work for either array or
            scalar input. Input dimensions need to be compatible with those
            of BP_RP input.
        BP_RP :
            Gaia BP-RP color(s). Should work for either array or
            scalar input. Input dimensions need to be compatible with those
            of G input.

    Notes
    -----
        Right now this is pretty much trivial but in the future it
        could become more complex. For instance, the color correction
        could get updated to have a second order term or to be
        piecewise linear. Other potential areas for future work are
        handling NaN BP_RP values and/or bounding BP_RP to some range.

    """

    par = common.pc_params()

    g_prime = G + par['bp_rp_coeff']*BP_RP

    return g_prime

def source_raw_pixel_metrics(cat, raw):
    """
    Tabulate raw pixel values at central locations of sources.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Source catalog.
        raw : numpy.ndarray
            Raw pointing camera image.

    """

    par = common.pc_params()

    ixcen = np.round(cat['xcentroid']).astype(int)
    iycen = np.round(cat['ycentroid']).astype(int)

    ixcen = np.minimum(np.maximum(ixcen, 0), par['nx']-1)
    iycen = np.minimum(np.maximum(iycen, 0), par['ny']-1)

    centroid_pixel_vals = raw[iycen, ixcen] # ordering of indices !

    centroid_pixel_saturated = (centroid_pixel_vals == par['raw_satur_val']).astype(int)

    # modify the input catalog
    cat['centroid_raw_pixel_val'] = centroid_pixel_vals
    cat['centroid_pixel_saturated'] = centroid_pixel_saturated

def gaia_pm_corr(cat, pc_mjd):
    """
    Compute proper motion corrected Gaia coordinates.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Source catalog.

        pc_mjd : float
            MJD of pointing camera exposure epoch.

    Notes
    -----
        Leave Gaia (RA, DEC) unchanged in cases lacking full
        five-parameter astrometric solution.
        When Gaia (PMRA, PMDEC) are available, correct
        Gaia (RA, DEC) to pc_mjd epoch.
        Don't do anything about parallax given how large the
        pointing camera pixels are. Input 'cat' gets modified.

    """

    print('Correcting Gaia positions for proper motion when possible')

    # epoch for Gaia DR2 (RA, DEC) values
    mjd_gaia = 57205.625

     # should probably add some sort of check(s) on pc_mjd
     # to catch e.g., cases of it somehow being 0 or a JD somehow

    ra_corr = cat['RA'] + ((pc_mjd - mjd_gaia)/365.25)*cat['PMRA']/(np.cos(cat['DEC']/(180.0/np.pi))*3600.0*1000.0)
    dec_corr = cat['DEC'] + ((pc_mjd - mjd_gaia)/365.25)*cat['PMDEC']/(3600.0*1000.0)

    full_solution = np.isfinite(cat['PMRA'])

    # RA, DEC columns get overwritten
    cat['RA'][full_solution] = ra_corr[full_solution]
    cat['DEC'][full_solution] = dec_corr[full_solution]

def recentroid_and_photometer(im, cat, one_aper=False, bg_sigclip=False):
    """
    Driver for running recentroiding followed by aperture photometry.

    Parameters
    ----------
        im : numpy.ndarray
            2D array with the detrended pointing camera image.
        cat : astropy.table.table.Table
            Source catalog.
        one_aper : bool, optional
            If True, only compute aperture photometry in the 'standard' aperture,
            rather than a series of apertures. This is meant to be a way
            of reducing run time.
        bg_sigclip : bool, optional
            If True, use sigma clipping when computing background level median
            value in sky annulus. False by default for speed, as the
            sigma_clipped_stats call is slow.

    Returns
    -------
        cat : astropy.table.table.Table
            Modified version of input table cat.

    """

    cat = pc_recentroid(im, cat)

    cat = pc_aper_phot(im, cat, one_aper=one_aper, bg_sigclip=bg_sigclip)

    return cat

def pc_phot(exp, one_aper=False, bg_sigclip=False, nmp=None, max_n_stars=3000,
            pm_corr=False, max_zp_radius=None):
    """
    Main photometry driver.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.
        one_aper : bool, optional
            If True, only compute aperture photometry in the 'standard' aperture,
            rather than a series of apertures. This is meant to be a way
            of reducing run time.
        bg_sigclip : bool, optional
            If True, use sigma clipping when computing background level median
            value in sky annulus. False by default for speed, as the
            sigma_clipped_stats call is slow.
        nmp : int, optional
            Number of processes to use for multiprocessing. Should be an integer
            greater than 1 but less than or equal to the number of CPU's. Default
            of None means no multiprocessing.
        max_n_stars : int, optional
            Maximum number of Gaia stars to retain. The idea is to avoid running
            into a situation where the pipeline takes forever trying to photometer
            a huge number of stars in dense fields. If the total number of stars
            exceeds max_n_stars, then the brightest max_n_stars are retained.
        pm_corr : bool, optional
            If set True, apply proper motion corrections to Gaia star positions,
            so that their positions match the epoch of the pointing camera
            exposure. By default this is False, on the assumption that only a
            tiny fraction of stars will have high enough proper motions for this
            to matter (El Nino, for instance, has ~8.6" pixels), especially given
            that this pipeline does empirical centroid refinement for all stars.
        max_zp_radius : float, optional
            Maximum radius in pixels from center of image for inclusion of
            stars in zeropoint determination. The idea is that beyond some
            radius from the image center, the PSF may degrade significantly
            such that we don't want to take stars far from the image center
            into account when computing the zeropoints. Default of None
            means that no such maximum image center distance cut is made.

    Returns
    -------
        cat : astropy.table.table.Table
            Source catalog with refined centroids and aperture photometry.

    """

    mag_thresh = max_gaia_mag(exp.time_seconds)

    cat = pc_gaia_cat(exp, mag_thresh=mag_thresh, nmp=nmp,
                      max_n_stars=max_n_stars, pm_corr=pm_corr)

    print('Recentroiding and aperture photometering...')
    t0 = time.time()
    if nmp is None:
        cat = recentroid_and_photometer(exp.detrended, cat, one_aper=one_aper,
                                        bg_sigclip=bg_sigclip)
    else:
        p =  Pool(nmp)
        parts = split_table(cat, nmp)
        args = [(exp.detrended, _cat, one_aper, bg_sigclip) for _cat in parts]
        cats = p.starmap(recentroid_and_photometer, args)
        cat = vstack(cats)
        p.close()
        p.join()

    dt = time.time()-t0
    print('recentroiding and photometry took ', '{:.3f}'.format(dt), ' seconds')

    print('Attempting to flag wrong centroids...')
    t0 = time.time()
    cat = flag_wrong_centroids_kdtree(cat)
    dt = time.time()-t0
    print('flagging wrong centroids took ', '{:.3f}'.format(dt), ' seconds')

    # add columns for quadrant, min_edge_dist_pix, BP-RP, m_inst, g_prime
    cat['BP_RP'] = cat['PHOT_BP_MEAN_MAG'] - cat['PHOT_RP_MEAN_MAG']

    cat['G_PRIME'] = get_g_prime(cat['PHOT_G_MEAN_MAG'], cat['BP_RP'])

    cat['quadrant'] = quadrant_from_xy(cat['xcentroid'], cat['ycentroid'])

    cat['min_edge_dist_pix'] = min_edge_dist_pix(cat['xcentroid'],
                                                 cat['ycentroid'])

    cat['flux_adu_per_s'] = cat['flux_adu']/exp.time_seconds

    cat['m_inst'] = -2.5*np.log10(cat['flux_adu_per_s'])

    par = common.pc_params()

    cat['radius_pix'] = \
        np.hypot(cat['xcentroid'] - central_pixel_coord(par['nx']),
                 cat['ycentroid'] - central_pixel_coord(par['ny']))

    cat['in_science_fov'] = \
        (cat['radius_pix'] < par['science_radius_pix']).astype('int16')

    if max_zp_radius is None:
        cat['radius_too_large'] = False
    else:
        cat['radius_too_large'] = \
        (cat['radius_pix'] > max_zp_radius)
    cat['radius_too_large'] = cat['radius_too_large'].astype('int16')

    if exp.has_dome is not None:
        cat['has_dome'] = int(exp.has_dome)

    source_raw_pixel_metrics(cat, exp.raw_image)

    return cat

def get_obs_night(date_string_local, time_string_local):
    """
    Determine observing night from local timestamp.

    Parameters
    ----------
        date_string_local : str
            Should be something like '2020/11/08'.
        time_string_local : str
            Should be something like '04:44:49'.

    Returns
    -------
        str
        Observing night string formatted like 'YYYYMMDD'

    Notes
    -----
        'local' for KPNO means MST.

    """

    # strip spaces from date_string_local and time_string_local
    date_string_local = date_string_local.replace(' ', '')
    time_string_local = time_string_local.replace(' ', '')

    assert(len(date_string_local) == 10)
    assert(len(time_string_local) == 8)

    hours = int(time_string_local[0:2])

    # stipulate that observing night rolls over at noon local time
    if hours >= 12:
        return date_string_local.replace('/', '')
    else:
        # figure out what was the previous calendar date
        fiducial_time_string = date_string_local.replace('/', '-') + 'T01:00:00'

        t = Time(fiducial_time_string, scale='utc')

        mjd_yesterday = t.mjd - 1.0

        t_yesterday = Time(mjd_yesterday, scale='utc', format='mjd')

        date_string_yesterday = np.datetime_as_string(t_yesterday.datetime64)[0:10].replace('-', '')

        return date_string_yesterday

def send_redis(exp, zp_info, sky_info, sci_inst_name='desi'):
    """
    Send pointing camera photometry telemetry to redis.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object.
        zp_info : astropy.table.table.Table
            A table with five rows of zeropoints information (quadrants 0-4
            inclusive with zero representing whole image).
        sky_info : astropy.table.table.Table
            A table with one row of sky background level information.
        sci_inst_name : str, optional
            Science instrument name. Meant to allow this to be configurable.
            Propagates into some of the field names that get sent to Redis.

    Returns
    -------
        This is fairly specific to the Mayall telescope.
        Should only be run within the appropriate computing
        environment where the relevant Redis connection can
        be established.

    """

    import redis

    timestamp = exp.header['DATE'].replace(' ', '').replace('/', '-') + '/' + \
                exp.header['TIME'].replace(' ', '') + '/MST/'

    # should think about potential failure modes here,
    # like MJD-OBS being empty in the header
    mjd_obs = float(exp.header['MJD-OBS'])

    zp_adu_per_s = [zp_info[(zp_info['quadrant'] == q) & (zp_info['science_fov_only'] == 0)]['zp_adu_per_s'][0] for q in range(5)]

    n_sources_for_zp = [int(zp_info[(zp_info['quadrant'] == q) & (zp_info['science_fov_only'] == 0)]['n_sources_for_zp'][0]) for q in range(5)]

    sky_adu_per_s = sky_info['median_adu_per_s']

    print('Redis timestamp = ' + timestamp)

    host = os.environ['REDISHOST']
    port = int(os.environ['REDISPORT'])
    db = int(os.environ['REDISDBNUM'])
    key = os.environ['REDISKEY']

    print('Attempting to connect with Redis...')
    r = redis.Redis(host=host, port=port, db=db)

    # should be careful about e.g., NaN values in quantities
    # sent to Redis; not sure what happens in such cases
    data = {'timestamp': timestamp,
            'zp_adu_per_s': zp_adu_per_s[0],
            'zp_adu_per_s_q1': zp_adu_per_s[1],
            'zp_adu_per_s_q2': zp_adu_per_s[2],
            'zp_adu_per_s_q3': zp_adu_per_s[3],
            'zp_adu_per_s_q4': zp_adu_per_s[4],
            'sky_adu_per_s': sky_adu_per_s,
            'mjd_obs': mjd_obs,
            'n_stars_for_zp': n_sources_for_zp[0],
            'n_stars_for_zp_q1': n_sources_for_zp[1],
            'n_stars_for_zp_q2': n_sources_for_zp[2],
            'n_stars_for_zp_q3': n_sources_for_zp[3],
            'n_stars_for_zp_q4': n_sources_for_zp[4],
            'sky_adu_per_s_q1': sky_info['median_adu_quad1_per_s'],
            'sky_adu_per_s_q2': sky_info['median_adu_quad2_per_s'],
            'sky_adu_per_s_q3': sky_info['median_adu_quad3_per_s'],
            'sky_adu_per_s_q4': sky_info['median_adu_quad4_per_s'],
            'sky_adu_per_s_' + sci_inst_name : sky_info['median_adu_sci_per_s'],
            'zp_adu_per_s_' + sci_inst_name : zp_info[zp_info['science_fov_only'] == 1]['zp_adu_per_s'][0],
            'n_stars_for_zp_' + sci_inst_name : int(zp_info[zp_info['science_fov_only'] == 1]['n_sources_for_zp'][0]),
            'flag' : int(sky_info['has_dome'])}

    print(data)

    r.hset(key, mapping=data)
    print('Redis data sent...')

def split_table(tab, n_parts):
    """
    Split a table into equal or near-equal sized chunks.

    Parameters
    ----------
        tab : astropy.table.table.Table
            Table to split into n_parts equal or nearly equal parts.
        n_parts : int
            Number of equal or nearly equal parts into which to split the
            input table. n_parts should be less than or equal to the
            number of rows in tab...

    Returns
    -------
        parts : list
            List with n_parts elements, with each element being a table.

    Notes
    -----
        Basically meant to be a wrapper for numpy.array_split that applies to
        an astropy Table rather than a numpy array. Quite possibly/likely
        this was a reinvention of something that could be done easily
        with existing tools.

    """

    nrows = len(tab)

    ind = np.arange(nrows)

    # list of arrays of indices into tab
    inds_per_part = np.array_split(ind, n_parts)

    # list of tables
    parts = [tab[_ind] for _ind in inds_per_part]

    # consistency check
    _nrows = 0
    for part in parts:
        _nrows += len(part)

    assert(nrows == _nrows)

    return parts

def flag_dome_vignetting(detrended, exptime_seconds):
    """
    Boolean classification of whether dome is vignetting any of the FOV.

    Parameters
    ----------
        detrended : numpy.ndarray
            Detrended image in units of ADU/pix (not ADU/pix/s).
        exptime_seconds : float
            Exposure time in seconds.

    Returns
    -------
        has_dome : bool
            Boolean label (True = has any amount of dome vignetting,
            False = no dome vignetting)

    Notes
    -----
        Counts number of pixels in detrended image (no downbinning)
        with < 0.5 ADU/pix/second (this value is tuned for El Nino).
        Nominal dark sky corresponding to r ~ 21 AB generates
        3.2 ADU/pix/second. This approach seems to work well
        for long El Nino exposure times of 20 seconds (the standard/routine
        value), but is likely to not work as well for shorter exposure times,
        where any variations of the bias levels will matter much more.

    """

    print('Attempting to flag dome vignetting...')

    par = common.pc_params()

    if exptime_seconds < 20:
        print('WARNING: dome flagging works better for longer exposures ' + \
              '(exposure time is only ' + \
              '{:.2f}'.format(exptime_seconds) + ' seconds)')

    frac = np.sum(detrended/exptime_seconds < par['dome_thresh_adu_per_s'])/detrended.size

    # note the 0.01 threshold here -- possible this could use more tuning
    has_dome = frac >= 0.01

    return has_dome

def central_pixel_coord(sidelen):
    """
    Return central pixel coordinate along an axis.

    Parameters
    ----------
        sidelen : int
            Integer side length in pixels.

    Returns
    -------
        coord : float
            Pixel coordinate of center. Will be a half-integer for
            even side length and integer for odd side length.

    Notes
    -----
        Convention is that center of first pixel corresponds to
        coordinate value of 0.

    """

    coord = float(sidelen)/2 - 0.5

    return coord

@lru_cache(maxsize=1)
def circular_mask(radius_pix):
    """
    Create a boolean mask representing a circle centered at the image center.

    Parameters
    ----------
        radius_pix : float
            Radius value to use for the circular mask.

    Returns
    -------
        mask : numpy.ndarray
            Boolean mask as a 2D numpy array. True means that a pixel
            location is within radius_pix pixels of the center.

    Notes
    -----
        Need to work on speeding this up, as there are definitely
        possible optimizations that remain to be made.

    """

    par = common.pc_params()

    sh = (par['ny'], par['nx'])

    Y, X = np.ogrid[:par['ny'], :par['nx']]

    x_center = central_pixel_coord(par['nx'])
    y_center = central_pixel_coord(par['ny'])

    dist = np.hypot(X - x_center, Y - y_center)

    dist = dist.reshape(sh)

    mask = (dist <= radius_pix)

    return mask

def add_field_center_cols(tab, header):
    """
    Add columns giving the RA, Dec coordinates of the field being observed.

    Parameters
    ----------
        tab : astropy.table.table.Table
            Table to add columns to.
        header : astropy.io.fits.header.Header
            FITS header that should include cards RADEG, DECDEG,
            REAL_RA, REAL_DEC

    Notes
    -----
        Input table is modified via the addition of columns (though I don't
        currently check whether the columns to be added already exist within the
        input table).

    """

    card2col = {'RADEG' : 'target_ra_deg',
                'DECDEG' : 'target_dec_deg',
                'REAL_RA' : 'real_ra_deg',
                'REAL_DEC' : 'real_dec_deg'}

    for card in card2col.keys():
        if card in header:
            tab[card2col[card]] = header[card]

def quiver_plot(_cat):
    """
    Make a quiver plot of recentroiding shifts.

    Parameters
    ----------
        _cat : astropy.table.table.Table
            Pointing camera source catalog including refined centroids.

    Returns
    -------
        status : bool
            Whether a quiver plot was made.

    """

    good = np.logical_not(_cat['centroid_pixel_saturated']) & \
           (_cat['centroid_raw_pixel_val'] < 15300) & \
           np.logical_not(_cat['centroid_shift_flag']) & \
           np.logical_not(_cat['wrong_source_centroid'])

    cat = _cat[good]

    if len(cat) == 0:
        print('no sources with good centroids to make quiver plot?')
        return False

    dx = cat['xcentroid'] - cat['x_gaia_guess']
    dy = cat['ycentroid'] - cat['y_gaia_guess']

    dist = np.hypot(dx, dy)

    plt.cla()
    plt.quiver(cat['xcentroid'], cat['ycentroid'], dx, dy,
               scale=10, scale_units='inches')

    plt.xticks([])
    plt.yticks([])

    par = common.pc_params()

    plt.xlim((0, par['nx']-1))
    plt.ylim((0, par['ny']-1))

    return True

def rebin(a, shape):
    """
    Rebin an array by via simple averaging.

    Parameters
    ----------
        a : numpy.ndarray
            2D numpy array to rebin.
        shape : tuple
            2-element tuple listing the target dimensions of the rebinned
            output array. Each element of shape should be an integer
            factor of its corresponding element in the dimensions of a.

    Notes
    -----
        Meant to be a Python port of the IDL routine REBIN.

        https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
    """

    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def ml_feature_vector(detrended):
    """
    Convert detrended image to a feature vector for machine learning classifier.

    Parameters
    ----------
        detrended : numpy.ndarray
            2D array holding the detrended image.

    Returns
    -------
        X : numpy.ndarray
            Feature vector corresponding to the detrended image;
            a flattened (1D) version of the downbinned and normalized
            detrended image.

    Notes
    -----
        The detrended image gets normalized based upon its mean, so it
        doesn't matter whether the detrended image has units of ADU versus
        ADU/second.

    """

    sh = detrended.shape

    binfac = 103 # extract this special number somewhere else...

    assert((sh[0] % binfac == 0) and (sh[1] % binfac == 0))

    ny_small = sh[0] // binfac
    nx_small = sh[1] // binfac

    reb = rebin(detrended, (ny_small, nx_small))

    reb /= np.mean(reb)

    X = np.ravel(reb)

    X = X.reshape(1, -1) # recommended by sklearn error/warning message...

    return X

def ml_dome_flag(detrended):
    """
    Machine learning boolean prediction of dome vignetting.

    Parameters
    ----------
        detrended : numpy.ndarray
            Detrended image as a 2D numpy array.

    Returns
    -------
        has_dome : bool
            Boolean label (True = has any amount of dome vignetting,
            False = no dome vignetting)

    Notes
    -----
        Loads pre-trained classifier, converts detrended image into feature
        vector, and computes/returns the classifier's boolean prediction.

    """

    print('Flagging dome vignetting with machine learning approach...')
    t0 = time.time()

    X = ml_feature_vector(detrended)

    svc = io.load_dome_clf()
    has_dome = bool(svc.predict(X)[0])

    dt = time.time()-t0
    print('Done with feature extraction and classification...took ' + \
          '{:.2f}'.format(dt) + ' seconds')

    return has_dome

@lru_cache()
def pointing_camera_index(night, require_standard_exptime=True,
                          raw_data_basedir=None):
    """
    Get a tabulation of pointing camera exposures and their MJD values.

    Parameters
    ----------
        night : str
            Eight element observing night string, YYYYMMDD format.
        require_standard_exptime : bool, optional
            Whether or not to require the standard exposure time, which
            for El Nino has been 20 seconds.
        raw_data_basedir : str, optional
            Full path to raw pointing camera data base directory. If no
            value is specified, it will be determined from PC_RAW_DATA_DIR.

    Returns
    -------
        t : astropy.table.table.Table
            Table of pointing camera image file names and their metadata
            for the relevant observing night.

    Notes
    -----
        Should add a column for ZPFLAG as well, being careful about the
        fact that ZPFLAG was only added in winter of 2021, so won't always
        be present.

        Why does this gather (HA, Dec) in sexagesimal format??

        Add multiprocessing option for speed-up?

        Caches result to avoid wasting time on repeated I/O.

    """

    par = common.pc_params()

    year = night[0:4]
    month = night[4:6]
    day = night[6:8]

    # pointing camera raw data convention
    monthdir = year + '-' + month
    nightdir = year + '-' + month + '-' + day

    if raw_data_basedir is None:
        raw_data_basedir = \
            os.environ[par['raw_data_env_var']]

    pattern = os.path.join(raw_data_basedir, monthdir, nightdir, '*.fits')

    flist = glob.glob(pattern)

    if len(flist) == 0:
        print('NO POINTING CAMERA EXPOSURES FOR NIGHT ' + night)
        return None

    flist.sort()

    tai_utc_offs = 37.0/(24.0*3600.0) # in days

    results = []
    for i, f in enumerate(flist):
        print('working on ', i+1, ' of ', len(flist), ' : ' + f)
        try:
            h = fits.getheader(f)
        except:
            print('corrupt raw pointing camera image??' + f)
            continue

        if 'MJD-OBS' in h:
            # would it be better to replace missing ZPFLAG with 0
            # rather than None?
            zpflag = h['ZPFLAG'] if 'ZPFLAG' in h else None
            result = (f, h['MJD-OBS'] - tai_utc_offs, h['HA'], h['DEC'], \
                      h['EXPTIME'], zpflag)
            results.append(result)
        else:
            print(f + ' does not have MJD-OBS??')

    t = Table()
    t['FNAME'] = [r[0] for r in results]
    t['MJD'] = [r[1] for r in results]
    t['HA'] = [r[2] for r in results]
    t['DEC'] = [r[3] for r in results]
    t['EXPTIME_SECONDS'] = [r[4]/1000.0 for r in results]
    t['ZPFLAG'] = [r[5] for r in results]

    if require_standard_exptime:
        t = t[t['EXPTIME_SECONDS'] == par['standard_exptime_seconds']]

    return t

def dome_slit_edge_azel(domeaz):
    """
    Return list of (az, el) coordinates running along Mayall dome slit edges.


    Parameters
    ----------
        domeaz : float
            Dome azimuth in degrees. Should be a scalar in the
            interval [0, 360).

    Returns
    -------
        az : numpy.ndarray
            List of azimuth coordinates in degrees along the Mayall dome slit
            edges. Should be within the interval [0, 360).
        el : numpy.ndarray
            List of elevation values in degrees along the Mayall dome slit
            edges.

    Notes
    -----
        Based on D. Joyce's Domeplot.ipynb notebook.

    """

    x = np.arange(0, 90, 0.01) # degrees

    rad2deg = 180.0/np.pi

    # Semi-width of dome slit (dome crane rails) = 3962 mm
    phi = rad2deg*np.arctan(3962/(15975*np.cos(np.radians(x))))

    y1 = domeaz - phi
    y2 = domeaz + phi

    az = np.concatenate((y1, y2))
    el = np.concatenate((x, x))

    az[az >= 360] -= 360.0
    az[az < 0] += 360

    return az, el

def flag_wrong_centroids_kdtree(cat):
    """
    Flag cases of a centroid wandering off to an entirely different star.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Table with columns including xcentroid, ycentroid, MY_BSC_ID.

    Returns
    -------
        cat : astropy.table.table.Table
            Input catalog cat but with an added column called wrong_source_centroid,
            which has a value of 1 (centroid has wandered too far) or 0.

    Notes
    -----
        The idea is that using a KDTree should make this fast...

    """

    n = len(cat)

    tree = KDTree(np.c_[cat['x_gaia_guess'], cat['y_gaia_guess']])

    dists, inds = tree.query(np.array((cat['xcentroid'],cat['ycentroid'])).T, k=1)

    inds = inds.reshape(n)

    cat['wrong_source_centroid'] = (inds != np.arange(n)).astype(int)

    return cat
