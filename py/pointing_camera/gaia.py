"""
pointing_camera.util
==================

Utilities for reading in the Gaia catalog.

    Notes
    -----
        These utilities are intended to mirror how the DESI imaging surveys
        access Gaia, namely through a HEALPix-elized full-sky catalog. Each
        catalog file contains one Nside = 32 HEALPix pixel worth of Gaia
        sources. The HEALPix indices are determined using RA/Dec as
        longitude/latitude. HEALPix indexing is ring-ordered.

"""

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

def gaia_chunknames(ipix):
    """
    Get the full file names of Gaia chunk files from their HEALPix indices.

        Parameters
        ----------
            ipix : np.ndarray
                List or numpy array of integer HEALPix pixel indices. Currently
                does not work for scalar ipix.

        Returns
        -------
            flist : list
                Sorted list of full HEALPix chunk file names.

        Notes
        -----
            Could add checks to make sure that all ipix values are sensible
            HEALPix pixel indices. Assumes that ipix input is *not* a scalar.
            Should eventually make this also work for scalar ipix input.

    """

    par = common.pc_params()

    gaia_dir = os.environ[par['gaia_env_var']]

    flist = [os.path.join(gaia_dir, 'chunk-' + str(i).zfill(5) + 
                                    '.fits') for i in ipix]

    flist.sort()

    return flist

@lru_cache(maxsize=1)
def read_gaia_chunknames(flist, nmp=None):
    """
    Read in a list of Gaia chunk files.

    Parameters
    ----------
        flist : tuple
            Full file names to read in. Note that the variable name
            flist is misleading, in that the caching decorator won't
            tolerate list input, but will tolerate tuple input.
        nmp : int (optional)
            Number of multiprocessing processes. Should be > 1 if set.
            Default is None, in which case the Gaia catalogs are read in
            serially.

    Returns
    -------
        tablist : list
            List of astropy Table objects, with one element per element of
            the input file list.

    Notes
    -----
        In the densest fields the cached catalog for one pointing camera FOV
        (many individual chunks combined) could be up to ~0.5 GB, so don't
        want to keep many such per-FOV catalogs cached.

        Caching is meant to address the common real-time use case
        where many poniting camera exposures will be taken while tracking
        at fixed sky location. In that case it's inefficient to re-read the
        same Gaia chunk files from disk for every single pointing camera
        exposure.

    """

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

    # should probably just do the table stacking here...
    return tablist

def read_gaia_cat(ra, dec, nmp=None):
    """
    Read in set of Gaia catalogs corresponding to a list of RA, Dec coords.

    Parameters
    ----------
        ra : numpy.ndarray
            RA values of coordinate pairs that need to be encompassed by
            the Gaia files read in. Units should be degrees.
        dec : numpy.ndarray
            Dec values of coordinate pairs that need to be encompassed by
            the Gaia files read in. Units should be degrees.
        nmp : int (optional)
            Number of multiprocessing processes. Should be > 1 if set.
            Default is None, in which case the Gaia catalogs are read in
            serially.

    Returns
    -------
        astropy.table.table.Table
            The Gaia catalog for the requested set of (ra, dec) coodinates.

    Notes
    -----
        Should add checks to make sure that ra and dec have compatible
        dimensions. Should also check that this works for both scalar and array
        ra/dec. Think that ra/dec inputs could equally well be e.g., lists
        rather than numpy arrays.

    """

    par = common.pc_params()

    ipix_all = healpy.pixelfunc.ang2pix(par['gaia_nside'], ra, dec, nest=False,
                                        lonlat=True)

    ipix_u = np.unique(ipix_all)

    flist = gaia_chunknames(ipix_u)

    t0 = time.time()

    # tuple so that caching decorator works...
    tablist = read_gaia_chunknames(tuple(flist), nmp=nmp)

    dt = time.time()-t0

    print('took ' + '{:.3f}'.format(dt) + ' seconds to read ' + \
          str(len(flist)) + ' Gaia files')

    return np.hstack(tuple(tablist))
