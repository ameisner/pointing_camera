from nightly_plot import _twopanel, _read_concat_tables
import astropy.io.fits as fits
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

basedir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/proc_1night'

def _twopanel_1night(night, basedir=basedir, markersize=10):

    dir = os.path.join(basedir, night)

    assert(os.path.exists(dir))

    # note that the sky and zp tables are now both contained within the same
    # multi-extension FITS file

    flist = glob.glob(os.path.join(dir, '*-summary.fits'))

    skies = _read_concat_tables(flist, ext=3)
    zps = _read_concat_tables(flist, ext=2)

    _twopanel(skies, zps, clobber=True, save=False, markersize=markersize)

    plt.show()
    


