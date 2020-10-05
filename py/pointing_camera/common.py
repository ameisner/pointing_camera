def pc_params():
    # eventually generalize this to return correct parameters
    # for either el Nino or la Nina -- right now it's just for El Nino

    par = {'nx' : 3296,
           'ny' : 2472,
           'bitpix' : 16,
           'bias_med' : [87, 152, 118, 110],
           'bias_clipped_mean' : [87.5545, 151.814, 118.435, 109.746],
           'readnoise' : [12.5401, 13.0523, 12.7788, 12.6456],
           'dark_adu_per_s' : 0.089665301,
           'dark_adu_per_s_quad': [0.096094832, 0.099679574, 0.083682373,
                                   0.083619997],
           'meta_env_var' : 'POINTING_CAMERA_META',
           'static_mask_filename' : 'pc_badpix_mask.fits.gz'}

    return par
