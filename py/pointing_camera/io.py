"""
pointing_camera.io
================

I/O functions for pointing camera reduction pipeline.
"""

import astropy.io.fits as fits
import os
import pointing_camera.common as common
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from functools import lru_cache
import pickle
from scipy.stats import scoreatpercentile
import numpy as np

def write_image_level_outputs(exp, outdir):
    """
    Write image-level output of pointing camera reduction pipeline.

    Parameters
    ----------
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        outdir : str
            Full path of output directory.

    Notes
    -----
        For now the only image-level output that gets written is the
        the detrended image. Could imagine other image-level reduction
        pipeline outputs in the future...

    """

    print('Attempting to write image level outputs')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-detrended.fits')

    outname = os.path.join(outdir, outname)

    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    hdu = fits.PrimaryHDU(exp.detrended.astype('float32'), header=exp.header)
    hdu.writeto(outname_tmp)
    os.rename(outname_tmp, outname)


def write_bintables_mef(cat, zps, sky, exp, outdir):
    """
    Write binary table pipeline outputs to a multi-extension FITS file.

    Parameters
    ----------
        cat : astropy.table.table.Table
            Table containing the source catalog with centroid, photometry
            measurements and Gaia cross-match columns.
        zps : astropy.table.table.Table
            Table containing the summary of photometric zeropoint measurements.
        sky : astropy.table.table.Table
            Table containing the sky brightness measurements.
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        outdir : str
            Full path of output directory.

    """

    print('Attempting to write binary tables as multi-extension FITS')

    # for now assume that cat, zps, sky tables all exist

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-summary.fits')

    outname = os.path.join(outdir, outname)

    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    hdul = fits.HDUList(hdus=[fits.PrimaryHDU(header=exp.header),
                              fits.BinTableHDU(data=cat, header=exp.header),
                              fits.BinTableHDU(data=zps, header=exp.header),
                              fits.BinTableHDU(data=sky, header=exp.header)])

    hdul[1].header['EXTNAME'] = 'CATALOG'
    hdul[2].header['EXTNAME'] = 'ZEROPOINTS'
    hdul[3].header['EXTNAME'] = 'SKY'

    hdul.writeto(outname_tmp)

    os.rename(outname_tmp, outname)

@lru_cache(maxsize=1)
def load_static_badpix():
    """
    Read in the static bad pixel mask.

    Returns
    -------
        mask : numpy.ndarray
            Static bad pixel mask image. Should have the same dimensions
            as a standard raw image would. Should be an integer data type.

    """
    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['static_mask_filename'])

    assert(os.path.exists(fname))

    mask = fits.getdata(fname)

    return mask

@lru_cache(maxsize=1)
def load_master_bias():
    """
    Read in the master bias.

    Returns
    -------
        bias : numpy.ndarray
            Master bias image. Should have the same dimensions
            as a standard raw image would. Should be floating point
            data type (not integer).

    """

    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['master_bias_filename'])

    assert(os.path.exists(fname))

    print('READING MASTER BIAS: ' + fname)
    bias = fits.getdata(fname)

    return bias

@lru_cache(maxsize=1)
def load_master_dark():
    """
    Read in the master dark.

    Returns
    -------
        dark : numpy.ndarray
            Master dark image. Should have the same dimensions
            as a standard raw image would. Should be floating point
            data type (not integer).

    """

    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['master_dark_filename'])

    assert(os.path.exists(fname))

    print('READING MASTER DARK: ' + fname)
    dark = fits.getdata(fname)

    return dark

@lru_cache(maxsize=1)
def load_master_flat():
    """
    Read in the master flat.

    Returns
    -------
        flat : numpy.ndarray
            Master flat image. Should have the same dimensions
            as a standard raw image would. Should be floating point
            data type (not integer). Should have an overall median
            value of 1.

    """

    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['master_flat_filename'])

    assert(os.path.exists(fname))

    print('READING MASTER FLAT: ' + fname)
    flat = fits.getdata(fname)

    return flat

@lru_cache(maxsize=1)
def load_dome_clf():
    """
    Load dome vignetting machine learning classifier.

    Returns
    -------
        clf : sklearn.svm._classes.SVC
            Machine learning classifier for dome vignetting. Exact type
            of classifier may evolve over time based on further
            optimization/experimentation.

    """

    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['fname_clf'])

    assert(os.path.exists(fname))

    print('READING DOME VIGNETTING ML CLASSIFIER : ' + fname)
    clf = pickle.load(open(fname,"rb"))

    return clf

def save_zp_checkplot(exp, outdir):
    """
    Save a file with the zeropoint check plot.

    Parameters
    ----------
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        outdir : str
            Full path of output directory.

    """

    print('Attempting to save zeropoint check plot')

    if not plt.fignum_exists(1):
        return

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-zp.png')
    outname_tmp = 'tmp.' + outname

    outname = os.path.join(outdir, outname)
    outname_tmp = os.path.join(outdir, outname_tmp)

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    plt.savefig(outname_tmp, dpi=200, bbox_inches='tight')
    os.rename(outname_tmp, outname)

def write_streaks(exp, streaks, outdir):
    """
    Write a file with a summary of detected streaks.

    Parameters
    ----------
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        streaks : list
            List of streaks, each of which is a dictionary with data
            defining one detected streak.
        outdir : str
            Full path of output directory.

    Notes
    -----
        If input streaks variable is an empty list (no streaks detected) then
        no output is written.

        streaks list contains numpy arrays as dictionary values, which
        makes it not possible to dump to JSON, hence the pickle output file
        type.

    """

    if not len(streaks):
        print('Streaks file not written because no streaks were cataloged.')
        return

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-streaks.pkl')

    outname = os.path.join(outdir, outname)

    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    print('Writing satellite streaks...')

    pickle.dump(streaks, open(outname_tmp, "wb" ) )

    os.rename(outname_tmp, outname)

def plot_detrended(exp, outdir, plot_streaks=False):
    """
    Save a rendering of the detrended pointing camera image.

    Parameters
    ----------
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        outdir : str
            Full path of output directory.
        plot_streaks : bool, optional
            If True, overplot satellite streak outlines (if any were
            detected).

    Notes
    -----
        Downbinning?
        Should parcel out the plotting code to its own function...

    """

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-detrended.png')
    outname_tmp = 'tmp.' + outname

    outname = os.path.join(outdir, outname)
    outname_tmp = os.path.join(outdir, outname_tmp)

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))


    limits = scoreatpercentile(np.ravel(exp.detrended), [1, 99])

    plt.cla()
    plt.imshow(exp.detrended, vmin=limits[0], vmax=limits[1], origin='lower',
               interpolation='nearest', cmap='gray')

    plt.xticks([])
    plt.yticks([])

    title = exp.fname_im.split('/')[-1]
    title = title.replace('.fits', '')

    plt.title(title)

    if plot_streaks:
        # overplot satellite streaks
        if exp.streaks is not None:
            for streak in exp.streaks:
                plt.plot(streak['x'], streak['y'], linewidth=0.25)

    plt.savefig(outname_tmp, dpi=200, bbox_inches='tight')
    os.rename(outname_tmp, outname)

def save_quiver_plot(exp, cat, outdir):
    """
    Make and write out a quiver plot of the centroid refinement shifts.

    Parameters
    ----------
        exp : exposure.PC_exposure
            Pointing camera exposure object.
        cat : astropy.table.table.Table
            Table containing the source catalog with centroid, photometry
            measurements and Gaia cross-match columns.
        outdir : str
            Full path of output directory.

    """

    from pointing_camera.util import quiver_plot

    status = quiver_plot(cat)

    if not status:
        return

    title = exp.fname_im.split('/')[-1]
    title = title.replace('.fits', '')

    title += ' ; ' + str(exp.wcs.wcs.ctype)

    plt.title(title)

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-quiver.png')
    outname_tmp = 'tmp.' + outname

    outname = os.path.join(outdir, outname)
    outname_tmp = os.path.join(outdir, outname_tmp)

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    plt.savefig(outname_tmp, dpi=200, bbox_inches='tight')
    os.rename(outname_tmp, outname)
