import astropy.io.fits as fits
from astropy import wcs
import os

def load_wcs(fname_wcs):
    assert(os.path.exists(fname_wcs))

    w = wcs.WCS(hdul[0].header)

    return w
