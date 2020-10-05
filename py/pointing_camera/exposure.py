import pointing_camera.util as util

class PC_exposure:
    """Object encapsulating the contents of a single pointing camera exposure"""

    def __init__(self, fname_im):

        im, h = util.load_exposure_image(fname_im)

        # image file name
        self.fname_im = fname_im

        # pixel data
        self.image = im
        
        # image header
        self.header = h

        # wcs filename
        self.fname_wcs = util.get_wcs_filename(fname_im)
        
        wcs, wcs_header = util.load_wcs(self.fname_wcs)

        # wcs_header
        self.wcs_header = wcs_header
        
        # wcs object
        self.wcs = wcs
        
        self.detrended = False
