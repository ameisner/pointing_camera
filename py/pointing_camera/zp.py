import numpy as np
import matplotlib.pyplot as plt
import pointing_camera.common as common
import copy
from astropy.table import Table

def calc_zp(_cat, aper_ind, time_seconds, fname_im, quadrant=0):
    # quadrant = 0 means whole image (all quadrants combined)

    assert(time_seconds > 0)
    assert(quadrant in [0, 1, 2, 3, 4]) # note the 0 option here...
    
    par = common.pc_params()

    n_aper = len(par['aper_phot_objrad'])

    aper_ind = int(aper_ind)
    assert(aper_ind in np.arange(n_aper))

    cat = copy.deepcopy(_cat)

    good = np.logical_not(cat['centroid_pixel_saturated']) & \
           (cat['centroid_raw_pixel_val'] < 15300) & \
           np.logical_not(cat['centroid_shift_flag']) & \
           np.logical_not(cat['wrong_source_centroid']) & \
           np.isfinite(cat['PHOT_BP_MEAN_MAG']) & \
           np.isfinite(cat['PHOT_RP_MEAN_MAG']) & \
           np.isfinite(cat['PHOT_G_MEAN_MAG'])

    if quadrant is not 0:
        good = good & (cat['quadrant'] == quadrant)

    if np.sum(good) == 0:
        return None

    cat = cat[good]
    n = len(cat)

    m_inst = cat['m_inst'][:, aper_ind]

    diff = cat['g_prime'] - m_inst
    nf = np.sum(np.isfinite(diff))
    
    zp = np.nanmedian(diff)

    resid = diff - zp
    
    # now calculate robust sigma about the median zeropoint offset

    # at first glance argsort appears to put NaN's at end of sorted array
    sind = np.argsort(resid) # should look into what exactly happens with NaN's

    ind_l = max(round(0.16*nf), 0)
    ind_u = min(round(0.84*nf), n-1)

    resid_l = resid[sind[ind_l]]
    resid_u = resid[sind[ind_u]]

    sig_robust = (np.abs(resid_l) + np.abs(resid_u))/2.0

    result = Table()

    result['quadrant'] = [quadrant]
    result['zp_adu_per_s'] = [zp]
    result['n_sources_for_zp'] = [n]
    result['time_seconds'] = [time_seconds]
    result['aper_ind'] = [aper_ind]
    result['bp_rp_median'] = [np.nanmedian(cat['BP_RP'])]
    result['gaia_g_median'] = [np.nanmedian(cat['PHOT_G_MEAN_MAG'])]
    result['robust_sigma_mag'] = [sig_robust]
    result['fname_raw'] = [fname_im]

    return result
    
    
