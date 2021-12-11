"""
pointing_camera.exposure
======================

A class representing a pointing camera exposure.
"""

import pointing_camera.util as util
import pointing_camera.common as common
import os
import astropy.io.fits as fits
import copy

class PC_exposure:
    """Object encapsulating the contents of a single pointing camera exposure"""

    def __init__(self, fname_im):
        """
        Create a pointing camera exposure object.

        Parameters
        ----------
            fname_in : str
                File name of the raw pointing camera FITS image.

        """

        im, h = util.load_exposure_image(fname_im)

        util.check_image_dimensions(im)

        util._check_bitpix(h)

        self.obs_night = util.get_obs_night(h['DATE'], h['TIME'])

        # image file name
        self.fname_im = fname_im

        # pixel data
        self.raw_image = im

        # image header
        self.header = h

        # wcs filename
        self.fname_wcs = util.get_wcs_filename(fname_im)

        print('WCS filename: ' + self.fname_wcs)

        wcs, wcs_header = util.load_wcs(self.fname_wcs)

        util._validate_ctype(wcs)

        # wcs_header
        self.wcs_header = wcs_header

        # wcs object
        self.wcs = wcs

        # exposure time in seconds
        self.time_seconds = util.get_exptime(self.header)

        self.is_detrended = False

        self.has_dome = None

        self.streaks = None

    def update_dome_flag(self):
        """
        Compute and store boolean flag for dome vignetting.

        Notes
        -----
            Detrended image needs to have been computed/stored prior to
            running this.

        """

        assert(self.is_detrended)

        self.has_dome = util.flag_dome_vignetting(self.detrended,
                                                  self.time_seconds)

        self.header['DOMEFLAG'] = (self.has_dome, 'potential dome vignetting')

    def did_telescope_move(self):
        """
        Check raw image header metadata for telescope motion flag.

        Returns
        -------
            bool
                Returns True if ZPFLAG is 1, False if ZPFLAG is 0, None if
                ZPFLAG is not present.

        Notes
        -----
            First observing night of El Nino data with ZPFLAG in raw image
            headers should be 20211201.

        """

        if 'ZPFLAG' in self.header:
            return bool(self.header['ZPFLAG'])
        else:
            return None

    def _tmp_detrended_filename(self):
        """
        Get name of temporary detrended file.

        Returns
        -------
            outname : str
                Full temporary file name.

        """

        outname = os.path.split(self.fname_im)[-1].replace('.fits',
            '-detrended.fits.tmp')

        outname = os.path.join('/tmp', outname)

        return outname

    def _write_tmp_detrended(self, null_edge=False):
        """
        Write detrended image and header to /tmp as a FITS image file.

            null_edge : bool, optional
                Null out edge pixel values (the one edge-most row/column
                along each boundary) so that astride will 'close' the
                contours of satellite streaks that extend off the image edges.
                Could imagine trying other choices, like some value less
                than the minimum pixel value across the entire detrended image.
                Or maybe NaN values? Could also try nulling N > 1 edge-most
                rows/columns, rather than just N = 1.

        Notes
        -----
            This is a workaround to accommodate the astride Streak
            constructor.

        """

        outname = self._tmp_detrended_filename()

        if null_edge:
            par = common.pc_params()
            im = copy.deepcopy(self.detrended)
            im[:, 0] = 0
            im[:, par['nx']-1] = 0
            im[0, :] = 0
            im[par['ny']-1, :] = 0
            hdu = fits.PrimaryHDU(im, self.header)
        else:
            hdu = fits.PrimaryHDU(self.detrended, self.header)

        hdu.writeto(outname)

    def _del_tmp_detrended(self):
        """
        Delete temporary detrended image FITS file.

        """

        fname = self._tmp_detrended_filename()

        if os.path.exists(fname):
            os.remove(fname)
