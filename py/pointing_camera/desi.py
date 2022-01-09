"""
pointing_camera.desi
====================

DESI-specific utilities, including those related to the DESI telemetry database.
"""

import pointing_camera.util as util
import pointing_camera.common as common
from scipy.stats import scoreatpercentile
import imageio
import numpy as np
from pointing_camera.exposure import PC_exposure
import os
from multiprocessing import Pool
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob

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

    """

    assert(isinstance(night, str))
    assert(len(night) == 8)

    import DOSlib.exposure as exp
    import psycopg2 as psycopg
    import psycopg2.extras

    sql = "SELECT id, mjd_obs, night, exptime, reqra, reqdec, skyra, skydec, targtra, targtdec FROM exposure WHERE (night = " + night + ") AND (flavor = 'science') AND ((sequence = 'DESI') OR (sequence = '_Split'))"

    conn =  psycopg.connect(exp.dsn)

    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cursor.execute(sql)

    data = cursor.fetchall()

    cursor.close()

    conn.close()

    if len(data) == 0:
        print('NO DESI EXPOSURES FOUND FOR NIGHT ' + night)

    return data

def desi_exp_movie(_pc_index, expid, mjdmin, mjdmax, outdir='.',
                   reqra=None, reqdec=None):
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
        reqra : float, optional
            Actual desired field center RA for DESI sequences (degrees).
        reqdec : float, optional
            Actual desired field center Dec for DESI sequences (degrees).

    Notes
    -----
        Writes out an animated GIF movie of pointing camera images acquired
        during the DESI exposure.

    """

    print('WORKING ON EXPID = ' + str(expid))

    if _pc_index is None:
        return None

    # figure out which subset of pointing camera images
    # falls in the correct MJDRANGE

    # eventually get more sophisticated by taking into account ZPFLAG, to
    # avoid frames taken during a pointing offset or partially while slewing

    # could be careful about pointing camera timestamps being beginning
    # versus middle versus end of pointing camera exposure
    good = (_pc_index['MJD'] > mjdmin) & (_pc_index['MJD'] < mjdmax)

    if np.sum(good) == 0:
        print('NO POINTING CAMERA IMAGES CORRESPONDING TO DESI EXPID ' + str(expid))
        return None

    pc_index = _pc_index[good]

    ims = []
    for row in pc_index:
        im = one_pc_rendering(row['FNAME'], reqra=reqra, reqdec=reqdec)
        ims.append(im)

    # what happens in the case of exactly one pointing camera image?

    outname = str(expid).zfill(8) + '-pointing_camera.gif'
    outname = os.path.join(outdir, outname)

    # use imageio to create the GIF, specifying an FPS value
    imageio.mimsave(outname, ims, fps=10)

def all_movies_1night(night, outdir='.', nmp=None):
    """
    Generate all pointing camera animations for one observing night.

    Parameters
    ----------
        night : str
            Eight element observing night string, YYYYMMDD format.
        outdir : str, optional
            Full path of output directory. Defaults to current directory.
        nmp : int, optional
            Number of processes for multiprocessing. Values > 1 but less than
            the total number of CPUs on the machine make sense. Default
            of None means that the movies will be generated in serial.

    """

    print('WORKING ON NIGHT : ' + night)

    exp_desi = desi_exposures_1night(night)
    exp_pc = util.pointing_camera_index(night)

    seconds_per_day = 86400.0

    args = []
    for exposure in exp_desi:
        if (exposure['exptime'] is None) or (exposure['mjd_obs'] is None):
            continue

        mjdmin = exposure['mjd_obs']
        mjdmax = mjdmin + exposure['exptime']/seconds_per_day
        args.append((exp_pc, exposure['id'], mjdmin, mjdmax, outdir,
                     exposure['reqra'], exposure['reqdec']))

    if (nmp is None) or (nmp == 1):
        for arg in args:
            print('WORKING ON DESI EXPID = ' + str(arg[1]))
            desi_exp_movie(*arg)
    else:
        p =  Pool(nmp)
        p.starmap(desi_exp_movie, args)
        p.close()
        p.join()

def one_pc_rendering(fname, dome_flag_ml=False, reqra=None, reqdec=None):
    """
    Make a low-resolution rendering of one pointing camera image.

    Parameters
    ----------
        fname : str
            Raw pointing camera image file name.
        dome_flag_ml : bool, optional
            If True, use machine learning approach when flagging dome
            vignetting.
        reqra : float, optional
            Actual desired field center RA for DESI sequences (degrees).
        reqdec : float, optional
            Actual desired field center Dec for DESI sequences (degrees).

    Returns
    -------
        im : numpy.ndarray
            Downbinned rendering of the detrended pointing camera image that
            can be used as one frame in an animation.

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

    if (reqra is None) or (reqdec is None):
        x_center = (sh[1] // 2) - 0.5*((sh[1] % 2) == 0)
        y_center = (sh[0] // 2) - 0.5*((sh[0] % 2) == 0)
    else:
        x_center, y_center = exp.wcs.all_world2pix(reqra, reqdec, 0)
        x_center = float(x_center)/binfac
        y_center = float(y_center)/binfac

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

def radec_to_petal_loc(ra_deg, dec_deg, reqra, reqdec):
    """
    Compute DESI petal location values for sky locations given a field center.

    Parameters
    ----------
        ra_deg : numpy.ndarray
            Set of RA values (in degrees) for which to compute the corresponding
            DESI petal_loc values. Should have the same size as dec_deg.
        dec_deg : numpy.ndarray
            Set of Dec values (in degrees) for which to compute the
            corresponding DESI petal_loc values. Should have the same
            size as ra_deg.
        reqra : float
            DESI field center RA in degrees. Should be a scalar value.
        reqdec : float
            DESI field center Dec in degrees. Should be a scalar value.

    Returns
    -------
        petal_loc : numpy.ndarray
            Array of integer petal location values (0-9, inclusive) with
            same number of elements as input ra_deg and dec_deg coordinate
            lists.

    Notes
    -----
        -1 is used as a dummy/placeholder value to indicate when a sky location
        falls outside of the DESI FOV.

    """

    desi_fov_radius_deg = 1.6 # factor out this special number?

    sc_center = SkyCoord(reqra*u.deg, reqdec*u.deg, frame='icrs')
    sc = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame='icrs')

    pos_angle = sc_center.position_angle(sc).to(u.deg)

    petal_loc = np.floor((np.mod(pos_angle + 180.0*u.deg + 18.0*u.deg, 360.0*u.deg))/36.0*u.deg).astype(int)

    petal_loc = np.array(petal_loc)

    ang_sep = sc_center.separation(sc)

    petal_loc[ang_sep.deg > desi_fov_radius_deg] = -1

    return petal_loc

def movies_nightly_webpage(_dir):
    """
    Generate a nightly summary webpage displaying a grid of movie thumbnails.

    Parameters
    ----------
        _dir : str
            This should be one full directory path representing both the input
            and output directory (same directory for both).

    Notes
    -----
        Remember to handle the case of no movies for a night without crashing.
        Writes a file called summary.html in the _dir directory.

    """

    if not os.path.exists(_dir):
        return

    flist = glob.glob(os.path.join(_dir, '*.gif'))

    if len(flist) == 0:
        return

    flist.sort()

    outname = os.path.join(_dir, 'summary.html')

    lines = []
    lines.append('<HTML>')
    lines.append('<HEAD>')
    lines.append('</HEAD>')
    lines.append('')

    for i,f in enumerate(flist):
       if (i % 2) == 0:
           lines.append('<tr>')
       url = os.path.basename(f)
       lines.append('<td align="center"><a href="' + url + \
                    '"><img src="' + url + '"></a></td>')
       if ((i % 2) != 0) or (i == (len(flist)-1)):
           lines.append('</tr>')

    lines.append('<table border="0" width="1020">')
    lines.append('</table>')
    lines.append('')
    lines.append('</body>')
    lines.append('')
    lines.append('</HTML>')

    f = open(outname, 'w')

    f.writelines('\n'.join(lines))

    f.close()
