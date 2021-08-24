#!/usr/bin/env python

import astropy.io.fits as fits
import numpy as np
import os
import glob
from astropy.table import Table
import imageio
from pointing_camera.exposure import PC_exposure
import pointing_camera.util as util
from scipy.stats import scoreatpercentile

desidir = '/global/cfs/cdirs/desi/spectro/data'
lostfound = '/global/cfs/cdirs/desi/spectro/staging/lost+found'
pcdir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/el_nino'

def rebin(a, shape):
    # https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def guide_cube_flist(night):
    '''Get list of guide cubes'''

    flist = []
    for basedir in [desidir, lostfound]:
        pattern = os.path.join(basedir, night, '?'*8, 'guide-????????.fits.fz')
        flist += glob.glob(pattern)

    # any reason to sort the file list?

    # try to be careful about possible duplicates in both spectro/data and
    # staging/lost+found?

    return flist

def guide_cube_mjd_ranges(night):
    '''get a tabulation of guide cube file names and corresponding MJD ranges'''
    # check both spectro/data and staging/lost+found

    flist = guide_cube_flist(night)

    if len(flist) == 0:
        return None

    extnames = ['GUIDE' + str(gfa_loc) + 'T' for gfa_loc in [0, 2, 3, 5, 7, 8]]

    results = []
    for i, f in enumerate(flist):
        print('working on ', i+1, ' of ', len(flist), ' : ' + f)
        try:
            hdul = fits.open(f)
        except:
            print('corrupt raw guide cube??' + f)
            continue

        for extname in extnames:
            if extname in hdul:
                tab = hdul[extname].data
                expid = (os.path.split(f)[-1]).replace('guide-', '').replace('.fits.fz', '')
                # seems like I should discard the acq image when
                # determining MJDMIN
                result = (f, np.min(tab['MJD-OBS']), np.max(tab['MJD-OBS']),
                          extname, expid)
                results.append(result)
                break

    t = Table()
    t['FNAME'] = [r[0] for r in results]
    t['MJDMIN'] = [r[1] for r in results]
    t['MJDMAX'] = [r[2] for r in results]
    t['EXTNAME'] = [r[3] for r in results]
    t['EXPID'] = [r[4] for r in results]

    return t

def pointing_camera_index(night):
    '''get a tabulation of pointing camera exposures and their MJD values'''

    year = night[0:4]
    month = night[4:6]
    day = night[6:8]

    # pointing camera raw data convention
    monthdir = year + '-' + month
    nightdir = year + '-' + month + '-' + day
    
    pattern = os.path.join(pcdir, monthdir, nightdir, '*.fits')
    
    flist = glob.glob(pattern)

    if len(flist) == 0:
        return None
    
    flist.sort()

    results = []
    for i, f in enumerate(flist):
        print('working on ', i+1, ' of ', len(flist), ' : ' + f)
        try:
            h = fits.getheader(f)
        except:
            print('corrupt raw pointing camera image??' + f)
            continue

        if 'MJD-OBS' in h:
            result = (f, h['MJD-OBS'])
            results.append(result)
        else:
            print(f + ' does not have MJD-OBS??')

    t = Table()
    t['FNAME'] = [r[0] for r in results]
    t['MJD'] = [r[1] for r in results]

    return t

def one_pc_rendering(fname):
    '''make a low-resolution rendering of one pointing camera image'''

    # read in
    exp = PC_exposure(fname)

    # detrend
    util.detrend_pc(exp)

    # downsample via averaging
    binfac = 8

    sh = exp.detrended.shape
    im = rebin(exp.detrended, (int(sh[0]/binfac), int(sh[1]/binfac)))
    
    # figure out the stretch (should this be held constant across frames??)
    limits = scoreatpercentile(np.ravel(im), [1, 99])

    im[im < limits[0]] = limits[0]
    im[im > limits[1]] = limits[1]

    im -= limits[0]
    
    # overplot approximate DESI FOV to guide the eye

    return im

def desi_exp_movie(mjdrange, expid, pc_index, fps=5):
    '''make an animation of pointing camera images during one DESI exposure'''

    outdir = '.' # generalize this later...

    print('Working on DESI exposure ', expid)

    keep = (pc_index['MJD'] >= mjdrange[0]) & (pc_index['MJD'] < mjdrange[1])

    if np.sum(keep) == 0:
        print('no pointing camera exposures for exposure ' + expid)
        return

    pc_index = pc_index[keep]
    pc_index.sort('MJD')

    ims = []
    for row in pc_index:
        im = one_pc_rendering(row['FNAME'])
        ims.append(im)

    # what happens iin the case of exactly one pointing camera image?

    outname = expid + '-pointing_camera.gif'
    outname = os.path.join(outdir, outname)

    # use imageio to create the GIF, specifying an FPS value
    imageio.mimsave(outname, ims, fps=10)

def movies_1night(night):
    '''generate all the per DESI exposure animations for one night'''

    pc = pointing_camera_index(night)
    cubes = guide_cube_mjd_ranges(night)
    
    # if no guide cubes or no pointing camera images on requested night
    # then give up

    if (pc is None) or (cubes is None):
        print('no pointing camera animations to generate')
        return

    for cube in cubes:
        # at some point could try to take into account the
        # *ending* MJD (rather than starting MJD) of the final guide frame
        mjdrange = (cube['MJDMIN'], cube['MJDMAX'])
        desi_exp_movie(mjdrange, cube['EXPID'], pc)
