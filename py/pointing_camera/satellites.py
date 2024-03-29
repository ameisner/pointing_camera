"""
pointing_camera.satellites
==========================

Satellite streak detection.
"""

from astride import Streak
import time
import os
import timeout_decorator

def streak_radec(streak, wcs):
    """
    Compute streak boundary sky coordinates.

    Parameters
    ----------
        streak : dict
            Dictionary representing the parameters of a single streak. Must
            have keys 'x' and 'y' containing the pixel coordinate boundaries
            of the streak.
        wcs : astropy.wcs.wcs.WCS
            Astropy WCS object.

    Notes
    -----
        Input dictionary gets modified via the addition of keys called
        'ra' and 'dec' holding the (ra, dec) coordinates of the streak in
        units of degrees.

    """

    ra, dec = wcs.all_pix2world(streak['x'], streak['y'], 0)

    streak['ra'] = ra
    streak['dec'] = dec

class TookTooLongError(Exception):
    """Exception for computation that has taken too long to run."""
    pass

def detect_streaks(exp):
    """
    Run streak detection.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object with detrended image available.

    Returns
    -------
        streaks : list
            List of streaks, each of which is a dictionary with data
            defining one detected streak.

    Notes
    -----
        A lot more work could be done tuning the astride parameters
        that presumably exist.

        astride discards contours that aren't closed, which is a problem
        for long pointing camera exposures where the streaks often extend
        off of one or more image boundaries.

        Streak detection should operate on the detrended (rather than raw)
        pointing camera image.

    """

    print('Running streak detection...')

    t0 = time.time()

    exp._write_tmp_detrended(null_edge=True)

    fname = exp._tmp_detrended_filename()

    streak = Streak(fname)

    streak.detect()

    exp._del_tmp_detrended()

    streaks = streak.streaks

    for streak in streaks:
        streak_radec(streak, exp.wcs)

    dt = time.time() - t0

    print('Found ' + str(len(streaks)) + ' streaks')
    print('Streak detection took ' + '{:.2f}'.format(dt) + ' seconds ; ' +
          os.path.basename(exp.fname_im))

    exp.streaks = streaks

    return streaks

# allow up to 10 seconds to finish, could revisit this value
@timeout_decorator.timeout(10, timeout_exception=TookTooLongError)
def detect_streaks_time_limit(exp):
    """
    Wrapper for streak detection with an execution time limit.

    Parameters
    ----------
        exp : pointing_camera.exposure.PC_exposure
            Pointing camera exposure object with detrended image available.

    Returns
    -------
        streaks : list
            List of streaks, each of which is a dictionary with data
            defining one detected streak.

    """

    streaks = detect_streaks(exp)

    return streaks
