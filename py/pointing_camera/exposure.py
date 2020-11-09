import pointing_camera.util as util

class PC_exposure:
    """Object encapsulating the contents of a single pointing camera exposure"""

    def __init__(self, fname_im):

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
        
        wcs, wcs_header = util.load_wcs(self.fname_wcs)

        util._validate_ctype(wcs)

        # wcs_header
        self.wcs_header = wcs_header
        
        # wcs object
        self.wcs = wcs

        # exposure time in seconds
        self.time_seconds = util.get_exptime(self.header)
        
        self.is_detrended = False
