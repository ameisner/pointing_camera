import pointing_camera.common as common
import astropy.io.fits as fits
import healpy
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import time
from multiprocessing import Pool
from functools import lru_cache

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

    flist.sort()

    return flist

# in the densest fields the cached catalog
# for one pointing camera FOV could be up to ~0.5 GB, so
# don't want to keep many such catalogs cached
@lru_cache(maxsize=1)
def read_gaia_chunknames(flist, nmp=None):

    tablist = []

    if nmp is not None:
        p = Pool(nmp)
        tablist = p.map(fits.getdata, flist)
        p.close()
        p.join()
    else:
        for f in flist:
            print('READING : ', f)
            tab = fits.getdata(f)
            tablist.append(tab)

    return tablist

def read_gaia_cat(ra, dec, nmp=None):
    # should add checks to make sure that ra and dec have compatible dimensions
    # should also check that this works for both scalar and array ra/dec

    ipix_all = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)

    ipix_u = np.unique(ipix_all)

    flist = gaia_chunknames(ipix_u)

    t0 = time.time()

    # tuple so that caching decorator works...
    tablist = read_gaia_chunknames(tuple(flist), nmp=nmp)

    dt = time.time()-t0

    print('took ' + '{:.3f}'.format(dt) + ' seconds to read ' + \
          str(len(flist)) + ' Gaia files')

    return np.hstack(tuple(tablist))
