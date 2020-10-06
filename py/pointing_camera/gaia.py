import pointing_camera.common as common
import astropy.io.fits as fits
import healpy
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

# this is intended to mirror how the DESI imaging surveys access
# Gaia, namely through the HEALPix-elized full-sky catalog at:
#     /global/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom
#
# each catalog contains one Nside = 32 HEALPix pixel worth of Gaia souces
# the HEALPix indices are determined using RA/Dec as longitude/latitude
# HEALPix indexing is ring-ordered

nside = 32

def gaia_chunknames(ipix):
    # could add checks to make sure that all ipix values are 
    # sane HEALPix pixel indices
    # RIGHT NOW THIS ASSUMES IPIX IS AN ARRAY !! 
    # should eventually make this also work for scalar ipix

    par = common.pc_params()

    gaia_dir = os.environ[par['gaia_env_var']]

    flist = [os.path.join(gaia_dir, 'chunk-' + str(i).zfill(5) + 
                                    '.fits') for i in ipix]
    return flist

def read_gaia_cat(ra, dec):
    # should add checks to make sure that ra and dec have compatible dimensions
    # should also check that this works for both scalar and array ra/dec

    ipix_all = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)

    ipix_u = np.unique(ipix_all)

    flist = gaia_chunknames(ipix_u)

    tablist = []
    for f in flist:
        print('READING : ', f)
        tab = fits.getdata(f)
        tablist.append(tab)

    return np.hstack(tuple(tablist))
