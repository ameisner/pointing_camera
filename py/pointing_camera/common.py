import numpy as np

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
           'static_mask_filename' : 'pc_badpix_mask.fits.gz',
           'master_bias_filename' : 'pc_master_bias.fits.gz',
           'master_dark_filename' : 'pc_master_dark.fits.gz',
           'master_flat_filename' : 'pc_master_flat.fits.gz',
           'gaia_env_var' : 'PC_GAIA_DIR',
           'aper_phot_objrad' : 2.0 + 0.5*np.arange(7),
           'aper_phot_objrad_best' :  2.5,
           'annulus_radii' : [12.0, 20.0],
           'bp_rp_coeff' : 0.25,
           'raw_satur_val': 16383,
           'ncpus': 8,
           'science_radius_pix': 665.49217,
           'zp_adu_per_s': 17.58,
           'dome_thresh_adu_per_s': 0.5}

    return par
