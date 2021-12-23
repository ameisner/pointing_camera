"""
pointing_camera.desi
======================

DESI-specific utilities, including those related to the DESI telemetry database.
"""

import pointing_camera.util as util
import pointing_camera.common as common
from scipy.stats import scoreatpercentile
import imageio
import numpy as np
from pointing_camera.exposure import PC_exposure
import os

def desi_exposures_1night(night):
    """
    Gather list of DESI exposures for a given observing night.

    Parameters
    ----------
        night : str
            Eight element observing night string, YYYYMMDD format.

    Returns
    -------
        data : list
            List of dictionary-like objects, one per row of output returned
            by the SQL query.

    Notes
    -----
        What happens if there are no rows of output from the SQL query?

        This routine is not currently intended to downselect to exposures that
        were full-fledged DESI sequences. Might change this in the future.

    """

    assert(isinstance(night, str))
    assert(len(night) == 8)

    import DOSlib.exposure as exp
    import psycopg2 as psycopg
    import psycopg2.extras

    # downselect to full-fledged DESI sequences?
    sql = "SELECT id, mjd_obs, night, exptime, reqra, reqdec, skyra, skydec, targtra, targtdec FROM exposure WHERE (night = " + night + ") AND (flavor = 'science') AND (sequence = 'DESI')"

    conn =  psycopg.connect(exp.dsn)

    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cursor.execute(sql)

    data = cursor.fetchall()

    cursor.close()

    conn.close()

    return data

def desi_exp_movie(_pc_index, expid, mjdmin, mjdmax, outdir='.'):
    """
    Make a movie of pointing camera images during one DESI exposure.

    Parameters
    ----------
        _pc_index : astropy.table.table.Table
            Index table of pointing camera exposures spanning at least the
            (mjdmin, mjdmax) time interval. Expect that _pc_index will
            typically cover a full observing night.
        expid : int
            Eight digit integer representing the observing night, YYYYMMDD
            format.
        mjdmin : float
            Minimum MJD of the DESI exposure.
        mjdmax : float
            Maximum MJD of the DESI exposure.
        outdir : str, optional
            Full path of output directory.

    Notes
    -----
        Writes out an animated GIF movie of pointing camera images acquired
        during the DESI exposure.

    """

    # figure out which subset of pointing camera images
    # falls in the correct MJDRANGE

    # eventually get more sophisticated by taking into account ZPFLAG, to
    # avoid frames taken during a pointing offset or partially while slewing

    # could be careful about pointing camera timestamps being beginning
    # versus middle versus end of pointing camera exposure
    good = (_pc_index['MJD'] > mjdmin) & (_pc_index['MJD'] < mjdmax)

    if np.sum(good) == 0:
        return None

    pc_index = _pc_index[good]

    ims = []
    for row in pc_index:
        im = one_pc_rendering(row['FNAME'])
        ims.append(im)

    # what happens in the case of exactly one pointing camera image?

    outname = str(expid).zfill(8) + '-pointing_camera.gif'
    outname = os.path.join(outdir, outname)

    # use imageio to create the GIF, specifying an FPS value
    imageio.mimsave(outname, ims, fps=10)
    

def all_movies_1night(night, outdir='.'):
    """
    Generate all pointing camera animations for one observing night.

    Parameters
    ----------
        night : str
            Eight element observing night string, YYYYMMDD format.
        outdir : str, optional
            Full path of output directory. Defaults to current directory.

    """

    exp_desi = desi_exposures_1night(night)
    exp_pc = util.pointing_camera_index(night)

    seconds_per_day = 86400.0

    for exposure in exp_desi:
        if exposure['exptime'] is None:
            continue

        mjdmin = exposure['mjd_obs']
        mjdmax = mjdmin + exposure['exptime']/seconds_per_day
        desi_exp_movie(exp_pc, exposure['id'], mjdmin, mjdmax, outdir=outdir)

def one_pc_rendering(fname, dome_flag_ml=False):
    """
    Make a low-resolution rendering of one pointing camera image.

    Parameters
    ----------
        fname : str
            Raw pointing camera image file name.
        dome_flag_ml : bool, optional
            If True, use machine learning approach when flagging dome
            vignetting.

    Returns
    -------
        im : numpy.ndarray
            Downbinned rendering of the detrended pointing camera image that
            can be used as one frame in an animation.

    Notes
    -----
        Eventually upgrade to use REQRA, REQDEC and pointing camera WCS to most
        accurately overplot the DESI FOV.

    """

    exp = PC_exposure(fname)

    util.detrend_pc(exp)

    exp.update_dome_flag(use_ml=dome_flag_ml)

    has_dome = exp.has_dome

    # downsample via averaging
    binfac = 8

    sh = exp.detrended.shape
    im = util.rebin(exp.detrended, (int(sh[0]/binfac), int(sh[1]/binfac)))
    sh = im.shape
    
    # figure out the stretch (should this be held constant across frames??)
    limits = scoreatpercentile(np.ravel(im), [1, 99])

    im[im < limits[0]] = limits[0]
    im[im > limits[1]] = limits[1]

    # overplot approximate DESI FOV to guide the eye
    par = common.pc_params()
    
    ybox = np.arange(sh[0]*sh[1], dtype=int) // sh[1]
    xbox = np.arange(sh[0]*sh[1], dtype=int) % sh[1]

    xbox = xbox.astype('float')
    ybox = ybox.astype('float')

    x_center = (sh[1] // 2) - 0.5*((sh[1] % 2) == 0)
    y_center = (sh[0] // 2) - 0.5*((sh[0] % 2) == 0)

    xbox -= x_center
    ybox -= y_center

    xbox = xbox.reshape(sh)
    ybox = ybox.reshape(sh)

    # do this with numpy.hypot instead
    dist = np.sqrt(np.power(xbox, 2) + np.power(ybox, 2))

    mask = np.abs(dist - par['science_radius_pix']/binfac) < 0.5

    im[mask] = limits[1]

    if has_dome:
        # denote dome vignetting flag as a second white circle surrounding
        # the DESI FOV circle
        _mask = np.abs(dist - 1.05*par['science_radius_pix']/binfac) < 0.5
        im[_mask] = limits[1]

    im -= limits[0]

    return im
