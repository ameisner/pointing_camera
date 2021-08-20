from nightly_plot import _twopanel, _read_concat_tables
import astropy.io.fits as fits
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

basedir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/proc_1night'

def _twopanel_1night(night, basedir=basedir, markersize=10, save=False,
                     skip_q0=False):

    plt.cla()

    dir = os.path.join(basedir, night)

    assert(os.path.exists(dir))

    # note that the sky and zp tables are now both contained within the same
    # multi-extension FITS file

    flist = glob.glob(os.path.join(dir, '*-summary.fits'))

    flist.sort()

    skies = _read_concat_tables(flist, ext=3)
    zps = _read_concat_tables(flist, ext=2)

    title_extra = ' ; ' + night
    _twopanel(skies, zps, clobber=True, save=False, markersize=markersize,
              title_extra=title_extra, skip_q0=skip_q0)

    if save:
        outname = 'pointing_camera-' + night
        if skip_q0:
            outname += '-no_q0'
        outname += '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')
    else:
        plt.show()

def summer_2021_nightly_plots(markersize=2, skip_q0=False):

    nights = glob.glob(os.path.join(basedir, '????????'))

    nights = [os.path.split(night)[-1] for night in nights]

    nights.sort()
    
    for night in nights:
        _twopanel_1night(night, basedir=basedir, markersize=markersize,
                         save=True, skip_q0=skip_q0)


