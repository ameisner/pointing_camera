import astropy.io.fits as fits
from astropy import wcs
import os
import common

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
