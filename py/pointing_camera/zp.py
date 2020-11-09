import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pointing_camera.common as common
import copy
from astropy.table import Table, vstack

def calc_zp(_cat, aper_ind, time_seconds, fname_im, quadrant=0,
            checkplot=True):
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

    if quadrant != 0:
        good = good & (cat['quadrant'] == quadrant)

    if np.sum(good) == 0:
        return None

    cat = cat[good]
    n = len(cat)

    m_inst = cat['m_inst'][:, aper_ind]

    diff = cat['G_PRIME'] - m_inst
    nf = np.sum(np.isfinite(diff))
    
    zp = np.nanmedian(diff)

    resid = diff - zp
    
    # now calculate robust sigma about the median zeropoint offset

    # at first glance argsort appears to put NaN's at end of sorted array
    sind = np.argsort(resid) # should look into what exactly happens with NaN's

    ind_l = max(int(round(0.16*nf)), 0)
    ind_u = min(int(round(0.84*nf)), n-1)

    resid_l = resid[sind[ind_l]]
    resid_u = resid[sind[ind_u]]

    sig_robust = (np.abs(resid_l) + np.abs(resid_u))/2.0

    result = Table()

    result['quadrant'] = [quadrant]
    result['aper_ind'] = [aper_ind]
    result['zp_adu_per_s'] = [zp]
    result['n_sources_for_zp'] = [n]
    result['time_seconds'] = [time_seconds]
    result['bp_rp_median'] = [np.nanmedian(cat['BP_RP'])]
    result['gaia_g_median'] = [np.nanmedian(cat['PHOT_G_MEAN_MAG'])]
    result['robust_sigma_mag'] = [sig_robust]
    result['fname_raw'] = [fname_im]

    # checkplot (eventually make this optional)
    if checkplot and (quadrant == 0) and (aper_ind == 1):
        plt.cla()
        plt.figure(1)
        xtitle = 'G + 0.25*(BP-RP)'
        ytitle = '-2.5' + r'$\times$' + 'log' + r'$_{10}$' + '(ADU/sec)'
        title = fname_im.split('/')[-1]
        title = title.replace('.fits', '')
        title += '; aper1; all quads'

        plt.scatter(cat['G_PRIME'], m_inst, s=20, edgecolor='none',
                    facecolor='k')

        xmin = np.nanmin(cat['G_PRIME'])
        xmax = np.nanmax(cat['G_PRIME'])
        ymin = np.nanmin(m_inst)
        ymax = np.nanmax(m_inst)

        xsamp = np.array([xmin, xmax])
        ysamp = xsamp - zp

        plt.plot(xsamp, ysamp, linewidth=2, c='r')

        ax = plt.gca()

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        xtext = xlim[0] + (xlim[1] - xlim[0])*0.125
        ytext = ylim[0] + (ylim[1] - ylim[0])*0.875

        # this shouldn't crash for NaN zeropoint value...
        plt.text(xtext, ytext, 'ZP = ' + '{:.2f}'.format(zp), color='r')

        plt.title(title)
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

    return result
    
def calc_many_zps(cat, exp, checkplot=True):

    print('Attempting to calculate zeropoints')

    par = common.pc_params()

    results = []
    for q in [0, 1, 2, 3, 4]:
        for aper_ind in range(len(par['aper_phot_objrad'])):
            print('Computing zeropoint for quadrant : ', q, ' , aper ',
                  aper_ind)
            result = calc_zp(cat, aper_ind, exp.time_seconds, exp.fname_im,
                             quadrant=q, checkplot=checkplot)
            result['mjd_obs'] = exp.header['MJD-OBS']
            result['obs_night'] = exp.obs_night
            results.append(result)

    results = vstack(results)

    return results
