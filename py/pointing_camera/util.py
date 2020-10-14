import astropy.io.fits as fits
from astropy import wcs
import os
import common
import pointing_camera.io as io
import pointing_camera.analysis.djs_maskinterp as djs_maskinterp
from astropy.table import Table, hstack
import numpy as np
from pointing_camera.gaia import read_gaia_cat
import time
from pointing_camera.analysis.djs_photcen import djs_photcen
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
import photutils
from astropy.stats import sigma_clipped_stats

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

    w = wcs.WCS(hdul[0].header)

    return w, hdul[0].header

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
    # right now not meant to be vectorized...

    par = common.pc_params()

    half_x = par['nx']/2 - 0.5
    half_y = par['ny']/2 - 0.5

    if (x <= half_x) & (y <= half_y):
        return 3
    if (x <= half_x) & (y > half_y):
        return 2
    if (x > half_x) & (y <= half_y):
        return 4
    if (x > half_x) & (y > half_y):
        return 1

def _loop_quadrant_from_xy(x, y):

    assert(len(x) == len(y))
    assert(np.sum(np.logical_not(np.isfinite(x))) == 0)
    assert(np.sum(np.logical_not(np.isfinite(y))) == 0)
    
    quadrants = [quadrant_from_xy(*c) for c in zip(x, y)]

    return np.array(quadrants)

def min_edge_dist_pix(x, y):
    # minimum distance to any image edge
    # for now inputs are meant to be scalar, not array/list

    min_edge_dist = 20000

    par = common.pc_params()

    min_edge_dist = min(min_edge_dist, x + 0.5)
    min_edge_dist = min(min_edge_dist, y + 0.5)
    min_edge_dist = min(min_edge_dist, par['nx'] - 0.5 - x)
    min_edge_dist = min(min_edge_dist, par['ny'] - 0.5 - y)

    return min_edge_dist

def _loop_min_edge_dist_pix(x, y):

    assert(len(x) == len(y))

    _dists = [min_edge_dist_pix(*c) for c in zip(x, y)]

    return _dists

def subtract_dark_current(im, time_seconds):

    print('Subtracting dark current')

    assert(time_seconds < 30)

    par = common.pc_params()

    result = im.astype(float)

    for q in [1, 2, 3, 4]:
        dark_adu = par['dark_adu_per_s_quad'][q-1]*time_seconds
        print('quadrant : ', q, ', dark counts/pix : ', dark_adu)

        p = quad_pix_limits(q)

        result[p['ymin']:p['ymax'], p['xmin']:p['xmax']] -= dark_adu

    return result

def subtract_quad_offs(im):
    par = common.pc_params()

    # just in case this hasn't already been taken care of...
    im = im.astype(float)

    check_image_dimensions(im)

    for q in [1, 2, 3, 4]:
        p = quad_pix_limits(q)

        im[p['ymin']:p['ymax'], p['xmin']:p['xmax']] -= par['bias_med'][q-1]

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

    im = exp.raw_image.astype(float)

    im = subtract_quad_offs(im)

    im = subtract_dark_current(im, exp.time_seconds)

    im = badpix_interp(im)

    exp.is_detrended = True

    exp.detrended = im

    print('Finished detrending the raw pointing camera image')

def sky_summary_table(exp):

    print('Computing sky background statistics')

    med = np.median(exp.detrended)

    t = Table()

    t['median_adu'] = [med]
    t['median_adu_per_s'] = [med/exp.time_seconds]
    t['time_seconds'] = [exp.time_seconds]
    t['fname_raw'] = [exp.fname_im]

    qmeds = []
    for q in [1, 2, 3, 4]:
        qmed = np.median(get_quadrant(exp.detrended, q))
        qmeds.append(qmed)
        
        t['median_adu_quad' + str(q)] = [qmed]

        t['median_adu_quad' + str(q) + '_per_s'] = [qmed/exp.time_seconds]

    diff = np.array(qmeds) - np.median(qmeds)

    rms_adu = np.sqrt(np.mean(np.power(diff, 2)))

    t['inter_quad_rms_adu'] = [rms_adu]
    t['inter_quad_rms_adu_per_s'] = [rms_adu/exp.time_seconds]

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

def pc_gaia_cat(wcs, mag_thresh=None, edge_pad_pix=0):
    # wcs should be an astropy WCS object

    print('Reading Gaia DR2 catalogs...')

    xgrid, ygrid = xy_subsamp_grid()

    # last arg is 0 rather than 1 because that's what agrees with IDL
    ra, dec = wcs.all_pix2world(xgrid, ygrid, 0)

    cat = read_gaia_cat(ra, dec)

    if mag_thresh is not None:
        keep = (cat['PHOT_G_MEAN_MAG'] <= mag_thresh)
        assert(np.sum(keep) > 0)
        cat = cat[keep]

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

    print('Recentroiding...')
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

    wrong_source_centroid = np.zeros(len(result), dtype=bool)

    print('Attempting to flag wrong centroids...')

    t0 = time.time()
    for i in range(len(result)):
        _dist = np.sqrt(np.power(cat['x_gaia_guess'] - result['xcentroid'][i], 2) + np.power(cat['y_gaia_guess'] - result['ycentroid'][i], 2))
        indmin = np.argmin(_dist)
        wrong_source_centroid[i] = (indmin != i)

    dt = time.time()-t0

    result['wrong_source_centroid'] = wrong_source_centroid.astype(int)

    return result

def _get_area_from_ap(ap):
    # this is to try and work around the photutils API change
    # between versions 0.6 and 0.7
    if (photutils.__version__.find('0.7') != -1) or (photutils.__version__.find('1.0') != -1):
        area = ap.area # 0.7
    else:
        area = ap.area() # 0.6

    return area

def pc_aper_phot(im, cat):

    print('Attempting to do aperture photometry')

    im = im.astype(float)

    par = common.pc_params()

    positions = list(zip(cat['xcentroid'], cat['ycentroid']))

    radii = par['aper_phot_objrad']
    ann_radii = par['annulus_radii'] # should have 2 elements - inner and outer

    apertures = [CircularAperture(positions, r=r) for r in radii]
    annulus_apertures = CircularAnnulus(positions, r_in=ann_radii[0],
                                        r_out=ann_radii[1])
    annulus_masks = annulus_apertures.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(im)
        annulus_data_1d = annulus_data[mask.data > 0]
        # this sigma_clipped_stats call is actually the slow part !!
        _, median_sigclip, std_bg = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

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

def pc_phot(exp):
    # main photometry driver; exp is a PC_exposure object

    mag_thresh = max_gaia_mag(exp.time_seconds)

    cat = pc_gaia_cat(exp.wcs, mag_thresh=mag_thresh)

    centroids = pc_recentroid(exp.detrended, cat)

    assert(len(centroids) == len(cat))

    cat = hstack([cat, centroids])

    pc_aper_phot(exp.detrended, cat)

    # add columns for quadrant, min_edge_dist_pix, BP-RP, m_inst, g_prime
    cat['BP_RP'] = cat['PHOT_BP_MEAN_MAG'] - cat['PHOT_RP_MEAN_MAG']

    cat['G_PRIME'] = get_g_prime(cat['PHOT_G_MEAN_MAG'], cat['BP_RP'])

    cat['quadrant'] = _loop_quadrant_from_xy(cat['xcentroid'], cat['ycentroid'])

    cat['min_edge_dist_pix'] = _loop_min_edge_dist_pix(cat['xcentroid'],
                                                       cat['ycentroid'])

    cat['flux_adu_per_s'] = cat['flux_adu']/exp.time_seconds

    cat['m_inst'] = -2.5*np.log10(cat['flux_adu_per_s'])

    source_raw_pixel_metrics(cat, exp.raw_image)

    return cat
