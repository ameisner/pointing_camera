import movies_1night

def pc_file_list(night, expid):
    """
    Get pointing camera file list corresponding to a particular DESI EXPID.

    Parameters
    ----------
        night : str
            DESI observing night as an 8-character string (YYYYMMDD).
        expid : int
            DESI EXPID (doesn't need to be a full DESI sequence, can be e.g.,
            just a DESI guide cube with no corresponding spectroscopy)

    Returns
    -------
        flist : list
            List of raw pointing camera FITS image file names that 
            correspond to the requested DESI EXPID based on timestamp.

    """

    pc = movies_1night.pointing_camera_index(night)
    cubes = movies_1night.guide_cube_mjd_ranges(night)
    cubes['EXPID'] = cubes['EXPID'].astype(int)

    keep = cubes[cubes['EXPID'] == expid]

    if not len(keep):
        print('EXPID ' + str(expid) + ' does not exist on night ' + \
              str(night) + '?')
        return None


    mjdrange = [keep['MJDMIN'], keep['MJDMAX']]

    print(mjdrange)

    pc = movies_1night.pc_index_for_expid(pc, mjdrange)

    if pc is None:
        print('no pointing camera exposures for exposure ' + str(expid))
        return None

    return pc['FNAME']
