import astropy.io.fits as fits
from astropy import wcs
import os
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

def get_wcs_filename(fname_im, verify=True):
    # right now this is pretty trivial but in perhaps it could
    # become more complicated in the future

    fname_wcs = fname_im.replace('.fits', '.wcs')

    if verify:
        assert(os.path.exists(fname_wcs))

    return fname_wcs

def load_wcs(fname_wcs):
    assert(os.path.exists(fname_wcs))

    hdul = fits.open(fname_wcs)

    header = hdul[0].header

    # somehow this both averts a warning and saves ~0.25 seconds
    # when creating the wcs.WCS object
    header['NAXIS'] = 2

    w = wcs.WCS(header)

    return w, header

def load_exposure_image(fname):
    # this gets the pixel data and header for the actual pointing camera
    # image -- the WCS solution is loaded separately via load_wcs

    assert(os.path.exists(fname))

    print('Attempting to read exposure : ' + fname)

    im, h = fits.getdata(fname, header=True)

    print('Successfully read raw pointing camera image file')

    return im, h

def check_image_dimensions(image):
    # check that image dimensions make sense
    # of particular relevance is not getting fooled by downbinned data

    print('Checking raw pointing camera image dimensions...')

    par = common.pc_params()

    sh = image.shape

    assert(sh[0] == par['ny'])
    assert(sh[1] == par['nx'])

    print('Raw pointing camera image has correct dimensions')

def get_exptime(h_im, milliseconds=False):
    # h_im should be a header of a raw pointing camera image

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
    par = common.pc_params()

    print('Checking raw image BITPIX value')

    assert(h_im['BITPIX'] == par['bitpix'])

def quad_pix_limits(quad):
    # q is an integer quadrant number, one of [1, 2, 3, 4]

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
    assert(q in [1, 2, 3, 4])

    check_image_dimensions(im)

    par = common.pc_params()

    p = quad_pix_limits(q)

    quad = im[p['ymin']:p['ymax'], p['xmin']:p['xmax']]

    return quad

def quadrant_from_xy(x, y):
    # works for array-valued x, y

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
    # minimum distance to any image edge
    # works for array-valued x, y

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
    # interpolate over static badpixels

    mask = io.load_static_badpix()

    result = djs_maskinterp.average_bilinear(im, (mask != 0))

    return result

def detrend_pc(exp):
    # exp is a PC_exposure object

    print('Attempting to detrend the raw pointing camera image')

    assert(not exp.is_detrended)

    # try to avoid accidentally using units of ms instead of seconds
    assert(exp.time_seconds < 30)

    # final thing is to set detrended == true

    im = exp.raw_image.astype('float32')

    im = subtract_master_bias(im)

    im = subtract_master_dark(im, exp.time_seconds)

    im = apply_flatfield(im)
    
    im = badpix_interp(im)

    exp.is_detrended = True

    exp.detrended = im

    exp.update_dome_flag()

    print('Finished detrending the raw pointing camera image')

def sky_metrics(im):
    # im is a 2D numpy array, intended to be either a full-frame
    # image or one quadrant of a full-frame image

    sz = im.shape
    assert(len(sz) == 2)

    im_sorted = np.ravel(im)

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

    print('Computing sky background statistics')

    metrics = sky_metrics(exp.detrended)
    med = metrics['median']
    mean = metrics['sky_clipped_mean']

    t = Table()

    t['median_adu'] = [med]
    t['median_adu_per_s'] = [med/exp.time_seconds]
    t['mean_adu'] = [mean]
    t['mean_adu_per_s'] = [mean/exp.time_seconds]
    t['time_seconds'] = [exp.time_seconds]
    t['fname_raw'] = [exp.fname_im]

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

    if exp.has_dome is not None:
        t['has_dome'] = exp.has_dome.astype('int16')

    return t

def _validate_ctype(wcs):
    # wcs should be an astropy WCS object

    valid_ctypes = ['RA---TAN', 'DEC--TAN', 'RA---TAN-SIP', 'DEC--TAN-SIP']
    for ctype in wcs.wcs.ctype:
        assert(ctype in valid_ctypes)

def max_gaia_mag(time_seconds):
    # trying to make this work for both array and scalar inputs

    fac = 1.21

    val = 15.5 + 2.5*np.log10(time_seconds/26.0)*fac

    val = np.minimum(np.maximum(val, 10), 15.5)

    return val

def xy_subsamp_grid():
    # the values hardcoded here are going to give me problems
    # when trying to run on data with different dimensions e.g., la Nina

    x = np.arange(104)*32
    y = np.arange(104)*24

    xgrid, ygrid = np.meshgrid(x, y)

    return xgrid, ygrid


def pc_gaia_cat(exp, mag_thresh=None, edge_pad_pix=0, nmp=None,
                max_n_stars=3000, pm_corr=False):

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
        print('Restricting to the brightest ' + str(max_n_stars) + \
              ' of ' + str(len(cat)) + ' stars')
        sind = np.argsort(cat['PHOT_G_MEAN_MAG'])
        cat = cat[sind[0:max_n_stars]]

    return cat

def pc_recentroid(im, cat):
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
    wrong_source_centroid = np.zeros(len(cat), dtype=bool)

    for i in range(len(cat)):
        _dist = np.sqrt(np.power(full_cat['x_gaia_guess'] - cat['xcentroid'][i], 2) + np.power(full_cat['y_gaia_guess'] - cat['ycentroid'][i], 2))
        indmin = np.argmin(_dist)
        wrong_source_centroid[i] = (full_cat[indmin]['SOURCE_ID'] != cat[i]['SOURCE_ID'])

    cat['wrong_source_centroid'] = wrong_source_centroid.astype(int)

    return cat

def _get_area_from_ap(ap):
    # this is to try and work around the photutils API change
    # between versions 0.6 and 0.7
    if (photutils.__version__.find('0.7') != -1) or (photutils.__version__.find('1.0') != -1) or (photutils.__version__.find('1.2') != -1):
        area = ap.area # 0.7
    else:
        area = ap.area() # 0.6

    return area

def pc_aper_phot(im, cat, one_aper=False, bg_sigclip=False):

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
    # right now this is pretty much trivial but in the future it
    # could become more complex

    par = common.pc_params()

    g_prime = G + par['bp_rp_coeff']*BP_RP

    return g_prime

def source_raw_pixel_metrics(cat, raw):

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
    # leave Gaia (RA, DEC) unchanged in cases lacking full
    # five-parameter astrometric solution
    # when Gaia (PMRA, PMDEC) are available, correct
    # Gaia (RA, DEC) to pc_mjd epoch
    # don't do anything about parallax given how large the
    # the pointing camera pixels are...

    # 'cat' gets modified

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

    # wrapper that sequentially does recentroiding
    # and aperture photometry by calling pc_recentroid
    # and then pc_aper_phot
    # expect that im would typically be the full detrended pointing
    #    camera image

    cat = pc_recentroid(im, cat)

    cat = pc_aper_phot(im, cat, one_aper=one_aper, bg_sigclip=bg_sigclip)

    return cat

def pc_phot(exp, one_aper=False, bg_sigclip=False, nmp=None, max_n_stars=3000,
            pm_corr=False):
    # main photometry driver; exp is a PC_exposure object

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
    if nmp is None:
        cat = flag_wrong_centroids(cat, cat)
    else:
        p =  Pool(nmp)
        parts = split_table(cat, nmp)
        args = [(_cat, cat) for _cat in parts]
        cats = p.starmap(flag_wrong_centroids, args)
        cat = vstack(cats)
        p.close()
        p.join()

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

    if exp.has_dome is not None:
        cat['has_dome'] = exp.has_dome.astype('int16')

    source_raw_pixel_metrics(cat, exp.raw_image)

    return cat

def get_obs_night(date_string_local, time_string_local):
    # 'local' for KPNO means MST
    # date_string_local should be something like 2020/11/08
    # time_string_local should be something like 04:44:49

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

def send_redis(exp, zp_info, sky_info):

    # zp_info is a table with five rows:
    #     quadrants 0-4 inclusive with zero representing whole image
    # sky_info is one table row

    import redis

    timestamp = exp.header['DATE'].replace(' ', '').replace('/', '-') + '/' + \
                exp.header['TIME'].replace(' ', '') + '/MST/'

    # should think about potential failure modes here,
    # like MJD-OBS being empty in the header
    mjd_obs = float(exp.header['MJD-OBS'])

    zp_adu_per_s = [zp_info[zp_info['quadrant'] == q]['zp_adu_per_s'][0] for q in range(5)]

    n_sources_for_zp = [int(zp_info[zp_info['quadrant'] == q]['n_sources_for_zp'][0]) for q in range(5)]

    sky_adu_per_s = sky_info['mean_adu_per_s']

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
            'sky_adu_per_s_q1': sky_info['mean_adu_quad1_per_s'],
            'sky_adu_per_s_q2': sky_info['mean_adu_quad2_per_s'],
            'sky_adu_per_s_q3': sky_info['mean_adu_quad3_per_s'],
            'sky_adu_per_s_q4': sky_info['mean_adu_quad4_per_s']}

    print(data)

    r.hset(key, mapping=data)
    print('Redis data sent...')

def split_table(tab, n_parts):
    # n_parts is number of equal or nearly equal
    # parts into which to split tab

    # basically meant to be a wrapper for
    # numpy.array_split that applies to
    # an astropy Table rather than a numpy array

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
